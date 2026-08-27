#!/usr/bin/env python3
"""Losslessly repack an unsharded Zarr-v2/v3 VCZ as sharded Zarr v3."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import itertools
import json
import math
import shutil
import threading
import time
from pathlib import Path

import numpy as np
from numcodecs import blosc as numcodecs_blosc
import zarr
from zarr.codecs import BloscCodec
from zarr.core.buffer import default_buffer_prototype


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument(
        "--variant-shard-chunks",
        type=int,
        default=25,
        help="Number of logical variant chunks packed into each physical shard.",
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--validate", action="store_true")
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--compression-level", type=int, default=7)
    parser.add_argument(
        "--codec-threads",
        type=int,
        default=1,
        help="Blosc threads per shard task; keep at 1 with outer parallelism.",
    )
    parser.add_argument(
        "--compressed-passthrough",
        action="store_true",
        help=(
            "Pack compatible v2/v3 Blosc chunk payloads directly into v3 "
            "shards. Incompatible arrays automatically use decoded copying."
        ),
    )
    return parser.parse_args()


def digest(array: zarr.Array, block_rows: int) -> str:
    checksum = hashlib.sha256()
    if array.ndim == 0:
        checksum.update(np.asarray(array[...]).tobytes())
        return checksum.hexdigest()
    for start in range(0, array.shape[0], block_rows):
        stop = min(start + block_rows, array.shape[0])
        values = np.asarray(array[start:stop])
        if values.dtype.kind in {"O", "T", "U", "S"}:
            for value in values.flat:
                encoded = str(value).encode("utf-8")
                checksum.update(len(encoded).to_bytes(8, "little"))
                checksum.update(encoded)
        else:
            checksum.update(values.tobytes(order="C"))
    return checksum.hexdigest()


def v2_passthrough_config(
    zarray_metadata: dict[str, object], compression_level: int
) -> str | None:
    compressor = zarray_metadata.get("compressor")
    if not isinstance(compressor, dict):
        return None
    filters = zarray_metadata.get("filters")
    compatible_filters = filters is None or filters == [{"id": "vlen-utf8"}]
    compatible = (
        zarray_metadata.get("order") == "C"
        and compatible_filters
        and compressor.get("id") == "blosc"
        and compressor.get("cname") == "zstd"
        and compressor.get("clevel") == compression_level
        and compressor.get("shuffle") in {0, 1, 2}
        and compressor.get("blocksize") == 0
    )
    if not compatible:
        return None
    return {0: "noshuffle", 1: "shuffle", 2: "bitshuffle"}[compressor["shuffle"]]


def v3_passthrough_config(
    zarray_metadata: dict[str, object], compression_level: int
) -> str | None:
    codecs = zarray_metadata.get("codecs")
    if not isinstance(codecs, list) or len(codecs) != 2:
        return None
    array_codec, compressor = codecs
    if not isinstance(array_codec, dict) or not isinstance(compressor, dict):
        return None
    if array_codec.get("name") not in {"bytes", "vlen-utf8"}:
        return None
    configuration = compressor.get("configuration")
    if compressor.get("name") != "blosc" or not isinstance(configuration, dict):
        return None
    shuffle = configuration.get("shuffle")
    compatible = (
        configuration.get("cname") == "zstd"
        and configuration.get("clevel") == compression_level
        and shuffle in {"noshuffle", "shuffle", "bitshuffle"}
        and configuration.get("blocksize") == 0
    )
    return shuffle if compatible else None


def main() -> None:
    args = arguments()
    if args.variant_shard_chunks < 1:
        raise ValueError("--variant-shard-chunks must be positive")
    if args.workers < 1:
        raise ValueError("--workers must be positive")
    if not 0 <= args.compression_level <= 9:
        raise ValueError("--compression-level must be between 0 and 9")
    if args.codec_threads < 1:
        raise ValueError("--codec-threads must be positive")
    if args.destination.exists():
        if not args.overwrite:
            raise FileExistsError(args.destination)
        shutil.rmtree(args.destination)

    started = time.monotonic()
    root_v3_metadata = args.source / "zarr.json"
    source_format = 3 if root_v3_metadata.exists() else 2
    source = zarr.open_group(args.source, mode="r", zarr_format=source_format)
    destination = zarr.open_group(
        args.destination,
        mode="w",
        zarr_format=3,
        attributes=dict(source.attrs),
    )
    variant_count = source["variant_position"].shape[0]
    # Zarr shard writes are already parallelized below. Without this limit each
    # task can create its own Blosc thread team and heavily oversubscribe a host.
    numcodecs_blosc.set_nthreads(args.codec_threads)
    report: dict[str, object] = {
        "source": str(args.source),
        "destination": str(args.destination),
        "source_zarr_format": source_format,
        "variant_count": variant_count,
        "variant_shard_chunks": args.variant_shard_chunks,
        "workers": args.workers,
        "compression_level": args.compression_level,
        "codec_threads": args.codec_threads,
        "compressed_passthrough": args.compressed_passthrough,
        "arrays": {},
    }

    tasks: list[tuple[str, zarr.Array, zarr.Array, slice]] = []
    raw_tasks: list[
        tuple[
            str,
            Path,
            Path,
            object,
            tuple[int, ...],
            tuple[int, ...],
            tuple[int, ...],
            tuple[int, ...],
        ]
    ] = []
    for name, source_array in sorted(source.arrays()):
        source_array_dir = args.source / name
        metadata_name = "zarr.json" if source_format == 3 else ".zarray"
        zarray_metadata = json.loads((source_array_dir / metadata_name).read_text())
        if source_format == 3:
            passthrough_shuffle = v3_passthrough_config(
                zarray_metadata, args.compression_level
            )
        else:
            passthrough_shuffle = v2_passthrough_config(
                zarray_metadata, args.compression_level
            )
        shuffle_name = passthrough_shuffle or "shuffle"
        array_compressor = BloscCodec(
            cname="zstd",
            clevel=args.compression_level,
            shuffle=shuffle_name,
        )
        logical_chunks = tuple(source_array.chunks)
        shard_shape = list(logical_chunks)
        is_variant_array = source_array.ndim > 0 and source_array.shape[0] == variant_count
        if is_variant_array:
            shard_shape[0] = min(
                source_array.shape[0],
                logical_chunks[0] * args.variant_shard_chunks,
            )
        shard_shape = tuple(
            min(axis_size, shard_axis)
            for axis_size, shard_axis in zip(source_array.shape, shard_shape, strict=True)
        )
        output = destination.create_array(
            name,
            shape=source_array.shape,
            dtype=source_array.dtype,
            chunks=logical_chunks,
            shards=shard_shape,
            fill_value=source_array.fill_value,
            attributes=dict(source_array.attrs),
            dimension_names=getattr(source_array.metadata, "dimension_names", None),
            compressors=[array_compressor],
        )

        use_passthrough = (
            args.compressed_passthrough
            and source_array.ndim > 0
            and passthrough_shuffle is not None
        )

        if source_array.ndim == 0:
            output[...] = source_array[...]
        elif use_passthrough:
            chunks_per_shard = tuple(
                math.ceil(shard / chunk)
                for shard, chunk in zip(shard_shape, logical_chunks, strict=True)
            )
            shards_per_axis = tuple(
                math.ceil(size / shard)
                for size, shard in zip(source_array.shape, shard_shape, strict=True)
            )
            sharding_codec = output.metadata.codecs[0]
            for shard_coords in itertools.product(
                *(range(count) for count in shards_per_axis)
            ):
                raw_tasks.append(
                    (
                        name,
                        source_array_dir,
                        args.destination / name,
                        sharding_codec,
                        shard_coords,
                        chunks_per_shard,
                        source_array.shape,
                        logical_chunks,
                        source_format,
                    )
                )
        else:
            rows_per_write = shard_shape[0]
            for start in range(0, source_array.shape[0], rows_per_write):
                stop = min(start + rows_per_write, source_array.shape[0])
                tasks.append((name, source_array, output, slice(start, stop)))

        entry: dict[str, object] = {
            "shape": source_array.shape,
            "dtype": str(source_array.dtype),
            "logical_chunks": logical_chunks,
            "shards": shard_shape,
            "expected_shards": math.prod(
                math.ceil(size / shard)
                for size, shard in zip(source_array.shape, shard_shape, strict=True)
            ),
            "copy_tasks": math.ceil(source_array.shape[0] / shard_shape[0])
            if source_array.ndim
            else 1,
            "copy_mode": "compressed_passthrough" if use_passthrough else "decoded",
        }
        report["arrays"][name] = entry

    timing_lock = threading.Lock()
    task_seconds: dict[str, float] = {name: 0.0 for name in report["arrays"]}
    task_max_seconds: dict[str, float] = {name: 0.0 for name in report["arrays"]}

    def copy_shard(task: tuple[str, zarr.Array, zarr.Array, slice]) -> None:
        name, source_array, output, row_slice = task
        task_started = time.monotonic()
        output[row_slice] = source_array[row_slice]
        elapsed = time.monotonic() - task_started
        with timing_lock:
            task_seconds[name] += elapsed
            task_max_seconds[name] = max(task_max_seconds[name], elapsed)

    prototype = default_buffer_prototype()

    def pack_raw_shard(
        task: tuple[
            str,
            Path,
            Path,
            object,
            tuple[int, ...],
            tuple[int, ...],
            tuple[int, ...],
            tuple[int, ...],
            int,
        ]
    ) -> None:
        (
            name,
            source_array_dir,
            destination_array_dir,
            sharding_codec,
            shard_coords,
            chunks_per_shard,
            array_shape,
            logical_chunks,
            task_source_format,
        ) = task
        task_started = time.monotonic()
        shard_chunks: dict[tuple[int, ...], object] = {}
        for relative_coords in itertools.product(
            *(range(count) for count in chunks_per_shard)
        ):
            global_coords = tuple(
                shard_coord * per_shard + relative_coord
                for shard_coord, per_shard, relative_coord in zip(
                    shard_coords,
                    chunks_per_shard,
                    relative_coords,
                    strict=True,
                )
            )
            if any(
                coordinate * chunk >= size
                for coordinate, chunk, size in zip(
                    global_coords, logical_chunks, array_shape, strict=True
                )
            ):
                continue
            if task_source_format == 3:
                chunk_path = source_array_dir / "c"
                for coordinate in global_coords:
                    chunk_path /= str(coordinate)
            else:
                chunk_path = source_array_dir / ".".join(map(str, global_coords))
            if chunk_path.exists():
                shard_chunks[relative_coords] = prototype.buffer.from_bytes(
                    chunk_path.read_bytes()
                )

        encoded = sharding_codec._encode_shard_dict_sync(
            shard_chunks,
            chunks_per_shard=chunks_per_shard,
            buffer_prototype=prototype,
        )
        if encoded is not None:
            shard_path = destination_array_dir / "c"
            for coordinate in shard_coords:
                shard_path /= str(coordinate)
            shard_path.parent.mkdir(parents=True, exist_ok=True)
            shard_path.write_bytes(encoded.to_bytes())
        elapsed = time.monotonic() - task_started
        with timing_lock:
            task_seconds[name] += elapsed
            task_max_seconds[name] = max(task_max_seconds[name], elapsed)

    copy_started = time.monotonic()
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        decoded_future = executor.map(copy_shard, tasks)
        passthrough_future = executor.map(pack_raw_shard, raw_tasks)
        list(decoded_future)
        list(passthrough_future)
    report["copy_elapsed_seconds"] = time.monotonic() - copy_started

    for name, entry in report["arrays"].items():
        entry["copy_task_seconds"] = task_seconds[name]
        entry["copy_max_task_seconds"] = task_max_seconds[name]

    if args.validate:
        for name, source_array in sorted(source.arrays()):
            output = destination[name]
            block_rows = output.shards[0] if source_array.ndim else 1
            entry = report["arrays"][name]
            entry["source_sha256"] = digest(source_array, block_rows)
            entry["destination_sha256"] = digest(output, block_rows)
            entry["valid"] = entry["source_sha256"] == entry["destination_sha256"]
            if not entry["valid"]:
                raise RuntimeError(f"validation failed for {name}")

    zarr.consolidate_metadata(args.destination)
    report["elapsed_seconds"] = time.monotonic() - started
    report["valid"] = all(
        entry.get("valid", True) for entry in report["arrays"].values()
    )
    report_path = args.destination.with_suffix(args.destination.suffix + ".validation.json")
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
