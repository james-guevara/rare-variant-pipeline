FROM public.ecr.aws/docker/library/rust:1.89-bookworm AS fastvep-builder

ARG FASTVEP_REPOSITORY=https://github.com/james-guevara/fastVEP.git
ARG FASTVEP_COMMIT=3bdb862b3153a90ef8cca0b07f02f357a15a3eb0
RUN apt-get update \
    && apt-get install -y --no-install-recommends ca-certificates git \
    && rm -rf /var/lib/apt/lists/* \
    && git clone "${FASTVEP_REPOSITORY}" /src/fastvep \
    && git -C /src/fastvep checkout --detach "${FASTVEP_COMMIT}" \
    && cargo build --locked --release --manifest-path /src/fastvep/Cargo.toml \
    && /src/fastvep/target/release/fastvep --version

FROM ghcr.io/astral-sh/uv:0.12.5 AS uv

FROM public.ecr.aws/docker/library/python:3.12-slim-bookworm

LABEL org.opencontainers.image.source="https://github.com/james-guevara/rare-variant-pipeline"
LABEL org.opencontainers.image.description="Targeted Zarr, FastVEP, picker, standalone LOFTEE, and carrier workflow"

RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        bcftools \
        ca-certificates \
        procps \
        tabix \
    && rm -rf /var/lib/apt/lists/*

COPY --from=uv /uv /uvx /usr/local/bin/
COPY --from=fastvep-builder /src/fastvep/target/release/fastvep /opt/fastvep/fastvep

WORKDIR /opt/rvp
COPY pyproject.toml uv.lock ./
COPY src ./src
RUN uv sync --frozen --no-dev \
    && .venv/bin/python -c "import duckdb, pyarrow, pyBigWig, pysam, zarr" \
    && /opt/fastvep/fastvep --version \
    && bcftools --version | head -n 1

COPY scripts ./scripts

ENV PATH="/opt/rvp/.venv/bin:/opt/fastvep:${PATH}" \
    RVP_REPO="/opt/rvp" \
    RVP_PYTHON="/opt/rvp/.venv/bin/python" \
    FASTVEP_ROOT="/opt/fastvep" \
    ANNOTATION_ROOT="/resources/ensembl-115" \
    LOFTEE_ROOT="/resources/loftee" \
    GENEBAYES="/resources/GeneBayes.Supplementary_Table_1.tsv" \
    ZARR_STORE="/inputs/chr22.sharded-v3.zarr" \
    TARGET_BED="/inputs/targets.chr22.bed" \
    RUN_ROOT="/output"

ENTRYPOINT []
CMD ["/opt/rvp/scripts/run_targeted_chr22_aws.sh"]
