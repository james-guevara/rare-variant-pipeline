# Container for the Python stages of the pipeline: PREPARE_VARIANTS,
# MERGE_ANNOTATIONS, CONVERT_PARQUET and the POSTPROCESS (PP_*) stages.
#
# WHY THIS EXISTS IN THE REPO
# Those stages previously ran against a micromamba env at an absolute path in one
# user's home directory, so nobody else could run the pipeline. duckdb, polars and
# polars-bio are conda-forge/PyPI packages with no bioconda equivalent, so there is
# no prebuilt biocontainer to pull. Expanse cannot build images locally either (no
# /etc/subuid mapping, so --fakeroot is refused), which is why this is built by CI
# and pulled as a published image instead.
#
# VERSIONS ARE PINNED to what the original micromamba env `rvp_env` ran, not to
# latest. polars and polars_bio both churn APIs, and prepare_variants.py is written
# against polars_bio 0.28 specifically (that release replaced the use_zero_based
# kwarg with metadata-based coordinate detection). Upgrading should be a deliberate,
# reviewed change — that is the point of pinning them here rather than in an env.
FROM python:3.12-slim

LABEL org.opencontainers.image.source="https://github.com/james-guevara/rare-variant-pipeline"
LABEL org.opencontainers.image.description="Python stages of the rare-variant pipeline (duckdb, polars, polars-bio)"

# procps is required by Nextflow, not by any of our code: the task wrapper shells
# out to `ps` to collect per-task metrics, and aborts the task if it is missing.
# python:*-slim omits it, which fails every process with only
# "Command 'ps' required by nextflow to collect task metrics cannot be found"
# in .command.err and an empty .command.out — the script never runs at all.
RUN apt-get update \
 && apt-get install -y --no-install-recommends procps \
 && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
        duckdb==1.5.2 \
        polars==1.39.3 \
        polars-bio==0.28.0 \
        pyarrow==22.0.0

# Fail the build, not a 2-hour pipeline run, if a pin resolves to something that
# lacks the API these scripts actually call.
RUN command -v ps >/dev/null || { echo "procps missing — Nextflow needs ps"; exit 1; }

RUN python - <<'PY'
import duckdb, polars, polars_bio, pyarrow
print("duckdb", duckdb.__version__, "| polars", polars.__version__,
      "| polars_bio", polars_bio.__version__, "| pyarrow", pyarrow.__version__)
for attr in ("count_overlaps", "set_option", "POLARS_BIO_COORDINATE_SYSTEM_CHECK"):
    assert hasattr(polars_bio, attr), f"polars_bio.{attr} missing"
print("API surface OK")
PY

CMD ["python"]
