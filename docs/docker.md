## The MintPy Docker container

[Docker](https://docs.docker.com/get-started/) allows you to run MintPy in a dedicated container (essentially an efficient virtual machine). [Here](https://docs.docker.com/install/) is the instruction to install docker.

### Pulling Mintpy container images

MintPy publishes Docker container images to [ghcr.io/insarlab/MintPy](https://github.com/insarlab/MintPy/pkgs/container/mintpy) in the [GitHub Container Registry](https://docs.github.com/en/packages/working-with-a-github-packages-registry/working-with-the-container-registry). To pull a MintPy container to your local machine with the latest stable MintPy release, run:

```shell
docker pull docker pull ghcr.io/insarlab/mintpy:latest
```

You can also pull a specific MintPy version (>=1.3.3) using the version tag:
```shell
docker pull docker pull ghcr.io/insarlab/mintpy:v1.3.3
```

Or the latest development version, which is the current HEAD of the main branch, with the `develop` tag:
```shell
docker pull docker pull ghcr.io/insarlab/mintpy:develop
```

*Note:* `latest` and `develop` are *rolling* tags, meaning they change as MintPy is developed. It's a good idea to use a specific version of MintPy in production systems as `docker pull` will always pull the most current version of those tags.

### Running the Mintpy container

To start an interactive shell session in the container, run:

```shell
docker run -it docker pull ghcr.io/insarlab/mintpy:latest
```

You can also run MintPy executables directly like:
```shell
docker run -it ghcr.io/insarlab/mintpy:latest smallbaselineApp.py --help
```

To map data on the host (local) machine to the container use [volumes](https://docs.docker.com/storage/volumes/):

```shell
docker run -it -v </path/to/data/dir>:/home/mambauser/work ghcr.io/insarlab/mintpy:latest
# Now inside the container
cd work
```

Which would allow you to process the mapped data like:

```
docker run -it -v </path/to/data/dir>:/home/mambauser/work ghcr.io/insarlab/mintpy:latest smallbaselineApp.py /home/mambauser/work/smallbaselineApp.cfg
```

### Launching Jupyter Server

The MintPy container includes a Jupyter Server with both the Jupyter Lab frontend and the classic Jupyter Notebook frontend. To launch Jupyter Lab, run:
```shell
 docker run -p 8888:8888 -it ghcr.io/insarlab/mintpy:latest jupyter lab
```
You can connect to the Jupyter Lab server at the `http://127.0.0.1:8888...` link printed to stdout in your terminal.

If you'd prefer to specify custom startup options to Jupyter, or prefer to use the Jupyter Notebook frontend instead of Jupyter Lab, you can specify them as part of the docker run command like:

```shell
docker run -p 8888:8888 -it ghcr.io/insarlab/mintpy:latest jupyter {lab,notebook} [JUPYTER_OPTIONS]
```
To see all the custom startup options, run:
```shell
docker run -p 8888:8888 -it ghcr.io/insarlab/mintpy:latest jupyter {lab,notebook} --help-all
```

### Notes ###

+ The container image is built using the [mambaorg/micromamba](https://hub.docker.com/r/mambaorg/micromamba) as a base. To manage conda environments inside the container use the `micromamba` command. For more information on micromamba, see: https://github.com/mamba-org/mamba#micromamba

+ Docker tightly maps user/group ids (uid/gid) inside and outside the container. By default, a `mambauser` with `uid=1000` and `gid=1000` will run inside the container and write files as that user. If you mount in a volume, files written to that volume will be owned by the *user on your local machine* with `uid=1000` and `gid=1000`. On linux and mac these are the default uid/gid values, but on a shared or managed system, these may not be *your* uid/gid values. You can override the users running inside the container with the `--user` argument to `docker run`, see: https://docs.docker.com/engine/reference/run/#user
