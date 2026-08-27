FROM ghcr.io/astral-sh/uv:0.12.5 AS uv

FROM ubuntu:24.04

RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        bcftools \
        ca-certificates \
        python3 \
        tabix \
    && rm -rf /var/lib/apt/lists/*

COPY --from=uv /uv /uvx /usr/local/bin/

RUN uv venv /opt/vcz \
    && uv pip install --python /opt/vcz/bin/python \
        "bio2zarr[vcf]==0.2.1" \
        numpy==2.5.2 \
        zarr==3.3.0 \
    && /opt/vcz/bin/python -c "import cyvcf2, numpy, zarr; import bio2zarr" \
    && /opt/vcz/bin/vcf2zarr --version \
    && bcftools --version | head -n 1

COPY scripts/repack_vcz_zarr_v3.py /usr/local/bin/repack_vcz_zarr_v3.py
COPY docker/vcz-batch/run_chromosome.sh /usr/local/bin/run_chromosome.sh

ENV PATH="/opt/vcz/bin:${PATH}"
ENTRYPOINT []
CMD ["/usr/local/bin/run_chromosome.sh"]
