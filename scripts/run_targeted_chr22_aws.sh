#!/usr/bin/env bash
# Backward-compatible entry point. New deployments use run_targeted_manifest.py.
set -euo pipefail
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
export COHORT=${COHORT:-g2mh}
export CHROMOSOME=${CHROMOSOME:-chr22}
export CONTIG=${CONTIG:-22}
export CONTIG_LENGTH=${CONTIG_LENGTH:-50818468}
export RVP_REPO=${RVP_REPO:-/home/ubuntu/work/rare-variant-pipeline-targeted-aws}
export RUN_ROOT=${RUN_ROOT:-/fsx/rare-variant-pilot/targeted-workflows/g2mh/chr22-lof-v1}
export ZARR_STORE=${ZARR_STORE:-/fsx/rare-variant-pilot/g2mh-vcz-v3/v1/chr22.sharded-v3.zarr}
export TARGET_BED=${TARGET_BED:-/home/ubuntu/work/rare-variant-missense-scope/targets/lof-plus-missense-candidates.chr22.bed}
export ANNOTATION_ROOT=${ANNOTATION_ROOT:-/home/ubuntu/work/fastvep-runtime/ensembl-115}
export LOFTEE_ROOT=${LOFTEE_ROOT:-/fsx/loftee-parity/resources}
export GENEBAYES=${GENEBAYES:-/home/ubuntu/work/standalone-loftee-test/port/GeneBayes.Supplementary_Table_1.tsv}
export FASTVEP_ROOT=${FASTVEP_ROOT:-/home/ubuntu/work/fastvep-runtime}
exec "$script_dir/run_targeted_chromosome.sh" "$@"
