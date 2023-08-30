## The MintPy Docker container ##

[Docker](https://docs.docker.com/get-started/) allows you to run MintPy in a dedicated container, which is essentially an efficient virtual machine. Check [here](https://docs.docker.com/install/) for the installation instruction.

### 1. Pulling the mintpy Docker image ###

We publish the mintpy Docker images in the [GitHub Container Registry](https://docs.github.com/en/packages/working-with-a-github-packages-registry/working-with-the-container-registry) at [ghcr.io/insarlab/mintpy](https://github.com/insarlab/MintPy/pkgs/container/mintpy).

The latest stable released version can be pulled to your local machine via the `latest` tag as:

```shell
docker pull ghcr.io/insarlab/mintpy:latest
```

The latest development version (the current HEAD of the main branch) can be pulled via the `develop` tag as:

```shell
docker pull ghcr.io/insarlab/mintpy:develop
```

Note that both `latest` and `develop` are *rolling* tags, meaning they change as MintPy evolves. Thus, in a production system, one may want to use a specific version for reproducebility. This is available (since version 1.3.3) via the version tag as:

```shell
docker pull ghcr.io/insarlab/mintpy:v1.3.3
```

### 2. Running the mintpy Docker container ###

Run the following to start an interactive shell session in the container with a host path to the data directory using [volumes](https://docs.docker.com/storage/volumes/):

```shell
docker run -it -v </path/to/data/dir>:/home/mambauser/data ghcr.io/insarlab/mintpy:latest
# use "docker run --name" option to name the container, e.g. "--name mintpy"
# then enter the running container as "docker exec -it mintpy"

# now inside the container
cd data/FernandinaSenDT128/mintpy
smallbaselineApp.py FernandinaSenDT128.txt
```

Or run mintpy executables directly as:

```shell
docker run -it -v </path/to/data/dir>:/home/mambauser/data ghcr.io/insarlab/mintpy:latest smallbaselineApp.py --help
docker run -it -v </path/to/data/dir>:/home/mambauser/data ghcr.io/insarlab/mintpy:latest smallbaselineApp.py /home/mambauser/data/FernandinaSenDT128/mintpy/FernandinaSenDT128.txt
```

Or run the following to launch the Jupyter Lab server, then copy and paste the printed `http://localhost:8888/lab?token=` url in a browser.

```shell
# to launch a Jupyter Notebook frontend, replace "lab" with "notebook" in the command below
docker run -p 8888:8888 -it ghcr.io/insarlab/mintpy:latest jupyter lab
```

Or launch the Jupyter server with custom startup options as:

```shell
docker run -p 8888:8888 -it ghcr.io/insarlab/mintpy:latest jupyter {lab,notebook} [JUPYTER_OPTIONS]
# to see all the custom startup option:
docker run -p 8888:8888 -it ghcr.io/insarlab/mintpy:latest jupyter {lab,notebook} --help-all
```

### Notes ###

+ The container image is built using the [mambaorg/micromamba](https://hub.docker.com/r/mambaorg/micromamba) as a base. To manage conda environments inside the container use the `micromamba` command. For more information on micromamba, see: https://github.com/mamba-org/mamba#micromamba

+ Docker tightly maps user/group ids (uid/gid) inside and outside the container. By default, a `mambauser` with `uid=1000` and `gid=1000` will run inside the container and write files as that user. If you mount in a volume, files written to that volume will be owned by the *user on your local machine* with `uid=1000` and `gid=1000`. On linux and mac these are the default uid/gid values, but on a shared or managed system, these may not be *your* uid/gid values. You can override the users running inside the container with the `--user` argument to `docker run`, see: https://docs.docker.com/engine/reference/run/#user
