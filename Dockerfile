# mageckinator runtime: MAGeCK, Python (numpy/pandas/scipy), aria2.
#
#   docker build -t ghcr.io/ckangas2/mageckinator:0.1.0 .
#
# Versions are pinned in environment.yml, shared with `-profile conda`.
# The pipeline's own scripts are NOT copied in: Nextflow puts bin/ on PATH at
# run time, and a baked-in copy would sit earlier on PATH and silently run
# stale code after every edit.
FROM mambaorg/micromamba:1.5.10-bookworm-slim

LABEL org.opencontainers.image.source="https://github.com/ckangas2/mageckinator" \
      org.opencontainers.image.description="MAGeCK reanalysis of a pooled CRISPR screen, with checks"

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
 && micromamba clean --all --yes \
 && rm /tmp/environment.yml

ENV PATH="/opt/conda/bin:${PATH}"

# Fail the build, not the first pipeline run, if anything is missing.
RUN mageck --version 2>&1 | grep -q 0.5.9 \
 && mageck count --help >/dev/null 2>&1 \
 && mageck test --help  >/dev/null 2>&1 \
 && python -c "import numpy, pandas, scipy" \
 && command -v ps zcat aria2c >/dev/null
