## 1. Install the released version ##

#### a. via conda ####

MintPy is available on the [conda-forge](https://anaconda.org/conda-forge/mintpy) channel. The latest released version can be installed via `conda` as:

```bash
conda install -c conda-forge mintpy
```

#### b. via docker ####

Docker allows one to run MintPy in a dedicated container, which is essentially an efficient virtual machine, and to be independent of platform OS. First, install [docker](https://docs.docker.com/install/) if you have not already done so. Then run the following to pull the latest stable released constainer image version from [MintPy GitHub Container Registry](https://github.com/insarlab/MintPy/pkgs/container/mintpy) to your local machine:

```bash
docker pull ghcr.io/insarlab/mintpy:latest
```

Check [here](docker.md) for more details on Docker container image usage, e.g. pulling development version and running in shell or Jupyter server.

#### c. via apt (Linux Debian) ####

MintPy is available in the main archive of the [Debian](https://tracker.debian.org/pkg/mintpy) GNU/Linux OS. It can be installed by using your favourite package manager or running the following command:

```bash
apt install mintpy
```

The same procedure, in priciple, can be used in [Ubuntu](https://ubuntu.com) and all the other [Debian derivatives](https://wiki.debian.org/Derivatives/Census). Check [here](https://salsa.debian.org/debian-gis-team/mintpy/-/blob/master/debian/README.Debian) for more detailed usage.


#### d. [Post-Installation Setup](#3-post-installation-setup) ####

## 2. Install the development version ##

Note: The installation note below is tested on Linux and macOS, and is still experimental on Windows (may have bugs).

MintPy is written in Python 3 and relies on several Python modules, check the [requirements.txt](https://github.com/insarlab/MintPy/blob/main/requirements.txt) file for details. We recommend using [conda](https://docs.conda.io/en/latest/miniconda.html) or [macports](https://www.macports.org/install.php) to install the python environment and the prerequisite packages, because of the convenient management and default [performance setting with numpy/scipy](http://markus-beuckelmann.de/blog/boosting-numpy-blas.html) and [pyresample](https://pyresample.readthedocs.io/en/latest/installation.html#using-pykdtree).

Quick links:

+ [Install on Linux](#21-install-on-linux)
+ [Install on macOS](#22-install-on-macos)
+ [Install on Windows](#23-install-on-windows)

### 2.1 Install on Linux ###

#### a. Download source code ####

Run the following in your terminal to download the development version of MintPy:

```bash
cd ~/tools
git clone https://github.com/insarlab/MintPy.git
```

#### b. Install dependencies via conda ####

Install [miniconda](https://docs.conda.io/en/latest/miniconda.html) if you have not already done so. You may need to close and restart the shell for changes to take effect.

```bash
# use wget or curl to download in command line or click from the web browser
# for macOS, use Miniconda3-latest-MacOSX-x86_64.sh instead.
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p ~/tools/miniconda3
~/tools/miniconda3/bin/conda init bash
```

Install the dependencies into a custom existing environment [recommended] by running:

```bash
# To create a new custom environment, e.g. named "insar", run "conda create --name insar; conda activate insar"
# To speedup, try "conda install mamba", then use "mamba install" to replace "conda install" below

# Add "gdal'>=3'" below to install extra dependencies if you use ARIA, FRInGE, HyP3 or GMTSAR
# Add "isce2"     below to install extra dependencies if you use ISCE-2
conda install -c conda-forge --file ~/tools/MintPy/requirements.txt
```

Or install the dependencies to a new environment named "mintpy" by running:

```bash
conda env create -f ~/tools/MintPy/docs/environment.yml
conda activate mintpy
```

#### c. Install MintPy ####

Install MintPy into the current environment with pip by running:

```bash
cd ~/tools
python -m pip install MintPy
```

Or install MintPy with pip in development mode [recommended] as below. The development mode allows one to install the package without copying files to your interpreter directory (e.g. the `site-packages` directory), thus, one could "edit" the source code and have changes take effect immediately without having to rebuild and reinstall.

```bash
cd ~/tools
python -m pip install -e MintPy
```

Or simply set up the environment variables as below in your source file, e.g. `~/.bash_profile` for _bash_ users or `~/.cshrc` for _csh/tcsh_ users:

```bash
if [ -z ${PYTHONPATH+x} ]; then export PYTHONPATH=""; fi
export MINTPY_HOME=~/tools/MintPy
export PATH=${PATH}:${MINTPY_HOME}/mintpy/cli
export PYTHONPATH=${PYTHONPATH}:${MINTPY_HOME}
```

#### d. [Post-Installation Setup](#3-post-installation-setup) ####

### 2.2 Install on macOS ###

Install Xcode with command line tools, if you have not already done so.

+ Install `Xcode` from App store

+ Install `command line tools` within XCode and agree to the terms of license.

  ```bash
  xcode-select --install -s /Applications/Xcode.app/Contents/Developer/
  sudo xcodebuild -license
  ```

+ Install [XQuartz](https://www.xquartz.org), then restart the terminal.

#### a. Install MintPy via conda ####

Same as the [instruction for Linux](#21-install-on-linux).

#### b. Install MintPy via MacPorts ####

Same as the [instruction for Linux](#21-install-on-linux), except for the dependencies' installation, which is as below.

Install [macports](https://www.macports.org/install.php) if you have not done so. Add the following at the bottom of your `~/.bash_profile` file:

```bash
# MacPorts Installer addition on 2017-09-02_at_01:27:12: adding an appropriate PATH variable for use with MacPorts.
export PATH=/opt/local/bin:/opt/local/sbin:${PATH}
export MANPATH=/opt/local/share/man:${MANPATH}
# Finished adapting your PATH environment variable for use with MacPorts.
```

Update the port tree with the following command. If your network prevent the use of rsync or svn via http of port tree, try [Portfile Sync via a Snapshot Tarball](https://trac.macports.org/wiki/howto/PortTreeTarball).

```
sudo port selfupdate
```

Install the dependencies by running:

```bash
# install dependencies with macports
# use "port -N install" to use the safe default for prompt questions
sudo port install $(cat MintPy/docs/ports.txt)

# install dependencies not available on macports: pysolid, pykml, pykdtree, pyresample, cdsapi
sudo -H /opt/local/bin/pip install git+https://github.com/insarlab/PySolid.git
sudo -H /opt/local/bin/pip install git+https://github.com/insarlab/PyAPS.git
sudo -H /opt/local/bin/pip install git+https://github.com/tylere/pykml.git
sudo -H /opt/local/bin/pip install git+https://github.com/storpipfugl/pykdtree.git
sudo -H /opt/local/bin/pip install git+https://github.com/pytroll/pyresample.git
sudo -H /opt/local/bin/pip install git+https://github.com/ecmwf/cdsapi.git
```

### 2.3 Install on Windows ###

Same as the [instruction for Linux](#21-install-on-linux), except for the "c. Install MintPy" section, only the `pip install` approaches are recommended, as the "setup environment variable" approach is not tested.

## 3. Post-Installation Setup

#### a. ERA5 for tropospheric correction ####

Set up an account for ERA5 to download weather re-analysis datasets for tropospheric delay correction as described in [insarlab/PyAPS](https://github.com/insarlab/pyaps#2-account-setup-for-era5).

`WEATHER_DIR`: Optionally, if you defined an environment variable named `WEATHER_DIR` to contain the path to a directory, MintPy will download the GAM files into the indicated directory. Also, MintPy will look for the GAM files in the directory before downloading a new one to prevent downloading multiple copies if you work with different dataset that cover the same date/time.

#### b. Dask for parallel processing ####

We recommend setting the `temporary-directory` in your [Dask configuration file](https://docs.dask.org/en/stable/configuration.html), e.g. `~/.config/dask/dask.yaml`, by adding the following line, to avoid potential [workspace lock issue](https://github.com/insarlab/MintPy/issues/725). Check more details on parallel processing with Dask [here](./dask.md).

```yaml
temporary-directory: /tmp  # Directory for local disk like /tmp, /scratch, or /local

# If you are sharing the same machine with others, use the following instead to avoid permission issues with others.
# temporary-directory: /tmp/{replace_this_with_your_user_name}
```
