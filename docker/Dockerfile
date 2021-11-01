# Use mambaforge as the base image
# Builds in ~9.5 min and is ~ 2.5 GB on a windows laptop
FROM mambaorg/micromamba:0.15.3

# Label image
LABEL \
    "Description"="Container for open source time series InSAR processing with Mintpy" \
    "Github Source"="https://github.com/insarlab/MintPy/" \
    "Installation"="https://github.com/insarlab/MintPy/blob/main/docs/installation.md" \
    "Dockerfile Author"="Forrest Williams" \
    "Email"="forrestfwilliams@icloud.com"

# Set user name (the micromamba image comes with the micromamba user)
ARG user=micromamba

# Clone MintPy and PyAPS repositories
USER root
RUN apt-get update && \
    apt-get install -y git vim && \
    rm /var/lib/apt/lists/* && \
    mkdir /home/micromamba/tools && \
    git clone https://github.com/insarlab/MintPy.git /home/micromamba/tools/MintPy && \
    git clone https://github.com/yunjunz/PyAPS.git /home/micromamba/tools/PyAPS

# Prepare environment for MintPy
ARG MINTPY_HOME=/home/$user/tools/MintPy
ARG PYAPS_HOME=/home/$user/tools/PyAPS
ARG PYTHON3DIR=/opt/conda

ENV MINTPY_HOME=${MINTPY_HOME}
ENV PYAPS_HOME=${PYAPS_HOME}
ENV PATH=${MINTPY_HOME}/mintpy:${PYTHON3DIR}/bin:${PATH}
ENV PYTHONPATH=${MINTPY_HOME}:${PYAPS_HOME}

# # # Optionally add Jupyter Lab to environment file
# # RUN echo "  - jupyterlab\n  - ipympl" >> ${MINTPY_HOME}/docs/environment.yml

# ADD mintpy.yml /tmp
RUN micromamba install -y -n base -f ${MINTPY_HOME}/docs/environment.yml python=3.6 && \
    micromamba clean --all --yes

# Create user and set work directory
USER $user
WORKDIR /home/$user

# # Have the container start with a Jupyter Lab instance
# CMD ["jupyter", "lab", "--port=8888", "--no-browser", "--ip=0.0.0.0","--NotebookApp.token=mintpy"]