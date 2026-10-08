FROM python:3.12.12-slim-bookworm@sha256:593bd06efe90efa80dc4eee3948be7c0fde4134606dd40d8dd8dbcade98e669c
ARG BENCHMARK_REVISION=unknown
LABEL org.opencontainers.image.source="https://github.com/james-guevara/rare-variant-pipeline"
LABEL org.opencontainers.image.title="Standalone carrier backend benchmark"
LABEL org.opencontainers.image.revision="${BENCHMARK_REVISION}"
ENV PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
WORKDIR /opt/rvp
COPY benchmarks/carrier/requirements.lock /tmp/requirements.lock
RUN python -m pip install --no-cache-dir --only-binary=:all: --require-hashes -r /tmp/requirements.lock
COPY scripts/extract_exact_carriers.py scripts/extract_filtered_carriers.py /opt/rvp/scripts/
COPY benchmarks/carrier /opt/rvp/benchmarks/carrier/
ENTRYPOINT ["python", "/opt/rvp/benchmarks/carrier/benchmark.py"]
