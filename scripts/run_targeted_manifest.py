#!/usr/bin/env python3
"""Validate a portable targeted-run manifest and execute its chromosome worker."""

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


ENVIRONMENT_KEYS = {
    "zarr_store": "ZARR_STORE",
    "target_bed": "TARGET_BED",
    "annotation_root": "ANNOTATION_ROOT",
    "loftee_root": "LOFTEE_ROOT",
    "genebayes": "GENEBAYES",
    "missense_candidates": "MISSENSE_CANDIDATES",
    "postprocess_config": "POSTPROCESS_CONFIG",
    "sample_sex_qc": "SAMPLE_SEX_QC",
    "sex_chromosome_regions": "SEX_CHROMOSOME_REGIONS",
}
REQUIRED_RESOURCES = {
    "zarr_store", "annotation_root", "loftee_root", "genebayes"
}


def read_json(path: Path) -> dict:
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"cannot read JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise ValueError(f"expected a JSON object: {path}")
    return value


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_string(document: dict, key: str, source: str) -> str:
    value = document.get(key)
    if not isinstance(value, str) or not value:
        raise ValueError(f"{source}.{key} must be a nonempty string")
    return value


def resolve(manifest: dict, bindings: dict) -> tuple[dict[str, str], list[str]]:
    if manifest.get("schema_version") != 1 or bindings.get("schema_version") != 1:
        raise ValueError("manifest and bindings must use schema_version 1")
    chromosome = require_string(manifest, "chromosome", "manifest")
    contig = require_string(manifest, "contig", "manifest")
    cohort = require_string(manifest, "cohort", "manifest")
    contig_length = manifest.get("contig_length")
    if not isinstance(contig_length, int) or contig_length <= 0:
        raise ValueError("manifest.contig_length must be a positive integer")
    specifications = manifest.get("resources")
    paths = bindings.get("resources")
    if not isinstance(specifications, dict) or not isinstance(paths, dict):
        raise ValueError("manifest.resources and bindings.resources must be objects")
    missing = sorted(REQUIRED_RESOURCES - specifications.keys())
    if missing:
        raise ValueError(f"manifest is missing required resources: {','.join(missing)}")

    environment = {
        "RUN_ROOT": require_string(bindings, "run_root", "bindings"),
        "COHORT": cohort,
        "CHROMOSOME": chromosome,
        "CONTIG": contig,
        "CONTIG_LENGTH": str(contig_length),
    }
    observations = []
    for logical_name, specification in specifications.items():
        if logical_name not in ENVIRONMENT_KEYS:
            raise ValueError(f"unsupported resource identity: {logical_name}")
        if not isinstance(specification, dict):
            raise ValueError(f"manifest.resources.{logical_name} must be an object")
        binding_name = specification.get("binding", logical_name)
        raw_path = paths.get(binding_name)
        if not isinstance(raw_path, str) or not raw_path:
            raise ValueError(f"bindings.resources.{binding_name} is required")
        root = Path(raw_path)
        kind = specification.get("kind", "file")
        if kind == "directory":
            if not root.is_dir():
                raise ValueError(f"resource directory is absent: {logical_name}={root}")
        elif kind == "file":
            if not root.is_file():
                raise ValueError(f"resource file is absent: {logical_name}={root}")
        else:
            raise ValueError(f"unsupported resource kind for {logical_name}: {kind}")
        check_path = root / specification["sentinel"] if specification.get("sentinel") else root
        if not check_path.is_file():
            raise ValueError(f"resource sentinel is absent: {logical_name}={check_path}")
        expected = specification.get("sha256")
        if expected:
            observed = sha256(check_path)
            if observed != expected:
                raise ValueError(
                    f"resource checksum differs: {logical_name}: expected={expected} observed={observed}"
                )
            observations.append(f"checksum.{logical_name}={observed}")
        environment[ENVIRONMENT_KEYS[logical_name]] = str(root)

    thresholds = manifest.get("thresholds", {})
    if not isinstance(thresholds, dict):
        raise ValueError("manifest.thresholds must be an object")
    for source_key, env_key in (
        ("population_af_max", "POPULATION_AF_MAX"),
        ("cohort_af_max", "COHORT_AF_MAX"),
    ):
        if source_key in thresholds:
            value = thresholds[source_key]
            if not isinstance(value, (int, float)) or not 0 <= value <= 1:
                raise ValueError(f"manifest.thresholds.{source_key} must be between 0 and 1")
            environment[env_key] = str(value)
    optional_outputs = manifest.get("optional_outputs", {})
    if not isinstance(optional_outputs, dict):
        raise ValueError("manifest.optional_outputs must be an object")
    synonymous = optional_outputs.get("synonymous_tiered_controls", True)
    if not isinstance(synonymous, bool):
        raise ValueError("manifest.optional_outputs.synonymous_tiered_controls must be boolean")
    environment["SYNONYMOUS_TIERED_CONTROLS"] = "1" if synonymous else "0"
    expected_qc = manifest.get("qc")
    if expected_qc is not None:
        if "POSTPROCESS_CONFIG" not in environment:
            raise ValueError("manifest.qc requires the postprocess_config resource")
        postprocess = read_json(Path(environment["POSTPROCESS_CONFIG"]))
        if postprocess.get("qc") != expected_qc:
            raise ValueError(
                "postprocess QC differs from the scientific manifest: "
                f"expected={expected_qc} observed={postprocess.get('qc')}"
            )
    if chromosome in {"chrX", "chrY"}:
        required_sex = {"SAMPLE_SEX_QC", "SEX_CHROMOSOME_REGIONS"}
        missing_sex = sorted(required_sex - environment.keys())
        if missing_sex:
            raise ValueError(
                "sex-chromosome manifest requires resources: "
                + ",".join(missing_sex)
            )
    return environment, observations


def preflight_run_root(run_root: Path) -> None:
    run_root.mkdir(parents=True, exist_ok=True)
    if not os.access(run_root, os.W_OK):
        raise ValueError(f"run root is not writable: {run_root}")
    with tempfile.NamedTemporaryFile(prefix=".preflight-", dir=run_root, delete=True) as handle:
        handle.write(b"ok\n")
        handle.flush()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--bindings", required=True, type=Path)
    parser.add_argument("--run-root", type=Path, help="override the binding's unique run root")
    parser.add_argument(
        "--lof-only", action="store_true",
        help="run only the LoF branch even when missense resources are declared",
    )
    default_root = Path(os.environ.get("RVP_REPO", Path(__file__).resolve().parents[1]))
    parser.add_argument(
        "--worker", type=Path,
        default=default_root / "scripts/run_targeted_chromosome.sh",
    )
    parser.add_argument("--preflight-only", action="store_true")
    parser.add_argument("--skip-runtime-checks", action="store_true")
    parser.add_argument(
        "--regression", type=Path,
        help="validate completed outputs against this pinned regression document",
    )
    args = parser.parse_args()

    manifest = read_json(args.manifest)
    bindings = read_json(args.bindings)
    if args.lof_only:
        manifest.get("resources", {}).pop("missense_candidates", None)
    regression = args.regression
    if regression is None and not args.lof_only and bindings.get("regression"):
        regression = Path(require_string(bindings, "regression", "bindings"))
    if args.run_root:
        bindings["run_root"] = str(args.run_root)
    environment, observations = resolve(manifest, bindings)
    if "TARGET_BED" not in environment:
        environment["ALL_OBSERVED"] = "1"
    if args.lof_only:
        environment.pop("MISSENSE_CANDIDATES", None)
    preflight_run_root(Path(environment["RUN_ROOT"]))
    if not args.skip_runtime_checks:
        for command in ("bcftools", "fastvep"):
            if shutil.which(command) is None:
                raise ValueError(f"required runtime command is absent: {command}")
    print("preflight=PASS")
    print(f"cohort={environment['COHORT']}")
    print(f"chromosome={environment['CHROMOSOME']}")
    for observation in observations:
        print(observation)
    sys.stdout.flush()
    if args.preflight_only:
        return
    worker = args.worker.resolve()
    if not worker.is_file():
        raise ValueError(f"worker is absent: {worker}")
    if regression is None:
        os.execvpe(str(worker), [str(worker)], os.environ | environment)
    if not regression.is_file():
        raise ValueError(f"regression document is absent: {regression}")
    subprocess.run([str(worker)], env=os.environ | environment, check=True)
    validator = default_root / "scripts/validate_targeted_regression.py"
    subprocess.run(
        [
            sys.executable,
            str(validator),
            "--regression", str(regression),
            "--bindings", str(args.bindings),
            "--run-root", environment["RUN_ROOT"],
        ],
        check=True,
    )


if __name__ == "__main__":
    main()
