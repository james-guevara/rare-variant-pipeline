#!/usr/bin/env python3
"""Infer XX-like/XY-like genotype patterns from chrX and chrY VCZ masks.

This is a cohort QC inference, not a declaration of biological sex. It samples
GRCh38 non-PAR blocks, deliberately avoiding the pseudoautosomal regions.
"""

import argparse
import csv
import json
from pathlib import Path

import numpy as np
import zarr


def sampled_block(root, start: int, end: int, maximum: int) -> tuple[np.ndarray, np.ndarray]:
    positions = root["variant_position"][:]
    lo = int(np.searchsorted(positions, start))
    hi = int(np.searchsorted(positions, end + 1))
    if hi <= lo:
        raise ValueError(f"no variants in requested interval {start}-{end}")
    center = (lo + hi) // 2
    block_start = max(lo, center - maximum // 2)
    block_end = min(hi, block_start + maximum)
    genotype = root["call_genotype"][block_start:block_end, :, :]
    mask = root["call_genotype_mask"][block_start:block_end, :, :]
    called = (~mask) & (genotype >= 0)
    return genotype, called


def sampled_depth(root, start: int, end: int, maximum: int) -> np.ndarray:
    positions = root["variant_position"][:]
    lo = int(np.searchsorted(positions, start))
    hi = int(np.searchsorted(positions, end + 1))
    if hi <= lo:
        raise ValueError(f"no variants in requested interval {start}-{end}")
    center = (lo + hi) // 2
    block_start = max(lo, center - maximum // 2)
    block_end = min(hi, block_start + maximum)
    depth = root["call_DP"][block_start:block_end, :].astype(float)
    depth[depth < 0] = np.nan
    return np.nanmedian(depth, axis=0)


DEFAULT_THRESHOLDS = {
    "x_haploid_rate_min": 0.8,
    "x_diploid_rate_min": 0.8,
    "y_present_call_rate_min": 0.5,
    "y_absent_call_rate_max": 0.2,
}


def classify(
    x_haploid_rate: float,
    x_diploid_rate: float,
    y_call_rate: float,
    thresholds: dict | None = None,
) -> str:
    values = DEFAULT_THRESHOLDS if thresholds is None else thresholds
    if (
        x_haploid_rate >= values["x_haploid_rate_min"]
        and y_call_rate >= values["y_present_call_rate_min"]
    ):
        return "XY-like"
    if (
        x_diploid_rate >= values["x_diploid_rate_min"]
        and y_call_rate < values["y_absent_call_rate_max"]
    ):
        return "XX-like"
    return "ambiguous"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--x-zarr", required=True, type=Path)
    parser.add_argument("--y-zarr", required=True, type=Path)
    parser.add_argument(
        "--autosome-zarr",
        type=Path,
        help="optional autosomal store used to normalize chrX/chrY median DP",
    )
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--config",
        required=True,
        type=Path,
        help="checksum-pinned cohort calibration JSON",
    )
    parser.add_argument("--variants", default=2000, type=int)
    args = parser.parse_args()
    if args.variants <= 0:
        parser.error("--variants must be positive")
    config = json.loads(args.config.read_text())
    if config.get("schema_version") != 1:
        raise ValueError("sex-chromosome config must use schema_version 1")
    regions = config["regions"]
    thresholds = config["thresholds"]

    x = zarr.open_group(str(args.x_zarr), mode="r")
    y = zarr.open_group(str(args.y_zarr), mode="r")
    x_samples = [str(value) for value in x["sample_id"][:]]
    y_samples = [str(value) for value in y["sample_id"][:]]
    if x_samples != y_samples:
        raise ValueError("chrX and chrY sample order differs")

    # Central non-PAR windows avoid both GRCh38 PARs and the X-transposed
    # region. Multiple thousands of sites make the rates stable to local noise.
    x_region = regions["x_nonpar"]
    y_region = regions["y_discriminating"]
    x_gt, x_called = sampled_block(
        x, x_region["start"], x_region["end"], args.variants
    )
    # The 11-13 Mb region has extensive XX-like off-target calls in G2MH.
    # The 3-5 Mb block cleanly separates Y-supported from Y-absent samples.
    _, y_called = sampled_block(
        y, y_region["start"], y_region["end"], args.variants
    )
    x_ploidy = x_called.sum(axis=2)
    y_ploidy = y_called.sum(axis=2)
    x_haploid_rate = (x_ploidy == 1).mean(axis=0)
    x_diploid_rate = (x_ploidy == 2).mean(axis=0)
    y_call_rate = (y_ploidy > 0).mean(axis=0)
    x_het_rate = (
        (x_ploidy == 2)
        & x_called[:, :, 0]
        & x_called[:, :, 1]
        & (x_gt[:, :, 0] != x_gt[:, :, 1])
    ).mean(axis=0)

    x_depth = y_depth = autosome_depth = None
    x_depth_ratio = y_depth_ratio = None
    if args.autosome_zarr:
        autosome = zarr.open_group(str(args.autosome_zarr), mode="r")
        autosome_samples = [str(value) for value in autosome["sample_id"][:]]
        if autosome_samples != x_samples:
            raise ValueError("autosome and sex-chromosome sample order differs")
        baseline = regions["autosome_baseline"]
        x_depth = sampled_depth(
            x, x_region["start"], x_region["end"], args.variants
        )
        y_depth = sampled_depth(
            y, y_region["start"], y_region["end"], args.variants
        )
        autosome_depth = sampled_depth(
            autosome, baseline["start"], baseline["end"], args.variants
        )
        x_depth_ratio = np.divide(
            x_depth, autosome_depth,
            out=np.full_like(x_depth, np.nan), where=autosome_depth > 0,
        )
        y_depth_ratio = np.divide(
            y_depth, autosome_depth,
            out=np.full_like(y_depth, np.nan), where=autosome_depth > 0,
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    counts = {"XX-like": 0, "XY-like": 0, "ambiguous": 0}
    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        header = [
            "sample_id", "inferred_karyotype", "x_haploid_rate",
            "x_diploid_rate", "x_het_rate", "y_call_rate",
        ]
        if args.autosome_zarr:
            header.extend([
                "x_median_dp", "y_median_dp", "autosome_median_dp",
                "x_autosome_dp_ratio", "y_autosome_dp_ratio",
                "x_depth_concordant",
            ])
        writer.writerow(header)
        for index, sample_id in enumerate(x_samples):
            label = classify(
                float(x_haploid_rate[index]),
                float(x_diploid_rate[index]),
                float(y_call_rate[index]),
                thresholds,
            )
            counts[label] += 1
            row = [
                sample_id,
                label,
                f"{x_haploid_rate[index]:.6f}",
                f"{x_diploid_rate[index]:.6f}",
                f"{x_het_rate[index]:.6f}",
                f"{y_call_rate[index]:.6f}",
            ]
            if args.autosome_zarr:
                ratio = float(x_depth_ratio[index])
                concordant = (
                    label == "XX-like"
                    and thresholds["xx_x_autosome_dp_ratio_min"] <= ratio
                    <= thresholds["xx_x_autosome_dp_ratio_max"]
                ) or (
                    label == "XY-like"
                    and thresholds["xy_x_autosome_dp_ratio_min"] <= ratio
                    <= thresholds["xy_x_autosome_dp_ratio_max"]
                )
                row.extend([
                    f"{x_depth[index]:.6f}",
                    f"{y_depth[index]:.6f}",
                    f"{autosome_depth[index]:.6f}",
                    f"{x_depth_ratio[index]:.6f}",
                    f"{y_depth_ratio[index]:.6f}",
                    "1" if concordant else "0",
                ])
            writer.writerow(row)
    print(
        "samples={} XX-like={} XY-like={} ambiguous={}".format(
            len(x_samples), counts["XX-like"], counts["XY-like"], counts["ambiguous"]
        )
    )
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
