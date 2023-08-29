# Builds in ~ 5 min and is ~ 3 GB on a linux laptop
FROM mambaorg/micromamba:0.24.0

# Label image following opencontainers image-spec annotations recommendation:
# https://github.com/opencontainers/image-spec/blob/main/annotations.md
LABEL org.opencontainers.image.description="Container for InSAR time series analysis with MintPy"
LABEL org.opencontainers.image.authors="Forrest Williams <forrestfwilliams@icloud.com>, Joseph H Kennedy <me@jhkennedy.org>, Andre Theron <andretheronsa@gmail.com>"
LABEL org.opencontainers.image.url="https://github.com/insarlab/MintPy"
LABEL org.opencontainers.image.source="https://github.com/insarlab/MintPy"
LABEL org.opencontainers.image.documentation="https://mintpy.readthedocs.io/en/latest/"
LABEL org.opencontainers.image.licenses="GPL-3.0-or-later"

# Dynamic labels to define at build time via `docker build --label`
# LABEL org.opencontainers.image.created=""
# LABEL org.opencontainers.image.version=""
# LABEL org.opencontainers.image.revision=""

USER root

ARG DEBIAN_FRONTEND=noninteractive
ENV PYTHONDONTWRITEBYTECODE=true

RUN apt-get update && \
    apt-get install -y --no-install-recommends git vim wget && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

USER mambauser
WORKDIR /home/mambauser

ENV PATH=/opt/conda/bin:${PATH}

ARG MINTPY_HOME=/home/mambauser/tools/MintPy
COPY --chown=mambauser:mambauser . ${MINTPY_HOME}/

ARG PYTHON_VERSION="3.9"
RUN micromamba install -y -n base -c conda-forge python=${PYTHON_VERSION}  \
      jupyterlab ipympl gdal">=3" isce2 -f ${MINTPY_HOME}/requirements.txt && \
    python -m pip install --no-cache-dir ${MINTPY_HOME} && \
    micromamba clean --all --yes

# Jupyter setup
COPY --chown=mambauser:mambauser scripts/jupyter_notebook_config.py /home/mambauser/.jupyter/
EXPOSE 8888
