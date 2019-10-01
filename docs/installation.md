## Install MintPy

Tested on macOS and Linux, not sure about Windows.

### Notes for Mac users ###

Install Xcode with command line tools, if you have not already done so.

+ Install `Xcode` from App store

+ Install `command line tools` within XCode and agree to the terms of license.

```
xcode-select --install -s /Applications/Xcode.app/Contents/Developer/
sudo xcodebuild -license
```

+ Install [XQuartz](https://www.xquartz.org), then restart the terminal.

### Notes for Docker users ###

Docker allows one to run MintPy in a dedicated container (essentially an efficient virtual machine) and to be independent of platform OS. After installing [docker](https://docs.docker.com/install/), run the following to pull the [MintPy container from DockerHub](https://hub.docker.com/r/andretheronsa/mintpy) to your local machine, check more details at [here](docker.md).

```
docker pull andretheronsa/mintpy:latest
```

### 1. Download and setup MintPy ###

To use the package, you need to setup the environment a) by adding _${MINTPY_HOME}_ to your _$PYTHONPATH_ to make mintpy importable in Python and b) by adding _${MINTPY_HOME}/mintpy_ to your _$PATH_ to make application scripts executable in command line, as shown below.

Add to your **_~/.bash_profile_** file for _bash_ user. For _tcsh_ user, check the example [here](https://github.com/yunjunz/macOS_Setup/blob/master/.tcshrc). Source the file for the first time. It will be sourced automatically next time when you login.

```bash
if [ -z ${PYTHONPATH+x} ]; then export PYTHONPATH=""; fi

##--------- MintPy ------------------##
export MINTPY_HOME=~/python/MintPy
export PYTHONPATH=${PYTHONPATH}:${MINTPY_HOME}
export PATH=${PATH}:${MINTPY_HOME}/mintpy

##--------- PyAPS ------------------##
export PYAPS_HOME=~/python/PyAPS
export PYTHONPATH=${PYTHONPATH}:${PYAPS_HOME}
```

Run the following in your terminal to download the development version of MintPy and PyAPS:

```
# download MintPy and PyAPS
git clone https://github.com/insarlab/MintPy.git $MINTPY_HOME
git clone https://github.com/yunjunz/pyaps3.git $PYAPS_HOME/pyaps3
```

### 2. Install dependencies ###

MintPy is written in Python3 and relies on several Python modules, check the [requirements.txt](https://github.com/insarlab/MintPy/blob/master/docs/requirements.txt) file for details. We recommend using [conda](https://conda.io/miniconda.html) or [macports](https://www.macports.org/install.php) to install the python environment and the prerequisite packages, because of the convenient managenment and default [performance setting with numpy/scipy](http://markus-beuckelmann.de/blog/boosting-numpy-blas.html) and [pyresample](https://pyresample.readthedocs.io/en/latest/installation.html#using-pykdtree).

#### Installing via conda ####

Add to your **_~/.bash_profile_** file:

```bash
export PYTHON3DIR=~/python/miniconda3
export PATH=${PATH}:${PYTHON3DIR}/bin
```

Run the following in your terminal to install miniconda:

```
# download and install miniconda
# use wget or curl to download in command line or from anaconda's web brower
# curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -o Miniconda3-latest-MacOSX-x86_64.sh
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
chmod +x Miniconda3-latest-MacOSX-x86_64.sh
./Miniconda3-latest-MacOSX-x86_64.sh -b -p $PYTHON3DIR
```

Run the following in your terminal to install the dependencies to the default environment _base_:
```
# install dependencies with conda
$PYTHON3DIR/bin/conda config --add channels conda-forge
$PYTHON3DIR/bin/conda install --yes --file $MINTPY_HOME/docs/conda.txt

# install dependencies not compatiable from conda: basemap, pykml
# run "conda uninstall basemap" if basemap was installed with conda
$PYTHON3DIR/bin/pip install git+https://github.com/matplotlib/basemap.git#egg=mpl_toolkits
$PYTHON3DIR/bin/pip install git+https://github.com/tylere/pykml.git
```

Or run the following in your terminal to install the dependencies to a new environment _mintpy_:

```
$PYTHON3DIR/bin/conda env create -f $MINTPY_HOME/docs/conda_env.yml
$PYTHON3DIR/bin/conda activate mintpy
```

#### Installing via MacPorts ####

Install [macports](https://www.macports.org/install.php) if you have not done so. Add the following at the bottom of your **_~/.bash_profile_** file:

```bash
# MacPorts Installer addition on 2017-09-02_at_01:27:12: adding an appropriate PATH variable for use with MacPorts.
export PATH=/opt/local/bin:/opt/local/sbin:${PATH}
export MANPATH=/opt/local/share/man:${MANPATH}
# Finished adapting your PATH environment variable for use with MacPorts.

#For py36-pyhdf in macports
export INCLUDE_DIRS=/opt/local/include
export LIBRARY_DIRS=/opt/local/lib
```

Update the port tree with the following command. If your network prevent the use of rsync or svn via http of port tree, try [Portfile Sync via a Snapshot Tarball](https://trac.macports.org/wiki/howto/PortTreeTarball).

```
sudo port selfupdate
```

Run the following in your terminal in _bash_:

```bash
# install dependencies with macports
# use "port -N install" to use the safe default for prompt questions
sudo port install $(cat $MINTPY_HOME/docs/ports.txt)

# install dependencies not available on macports: pykml, pykdtree, pyresample, cdsapi, pyhdf
sudo /opt/local/bin/pip install git+https://github.com/tylere/pykml.git
sudo /opt/local/bin/pip install git+https://github.com/storpipfugl/pykdtree.git
sudo /opt/local/bin/pip install git+https://github.com/pytroll/pyresample.git
sudo /opt/local/bin/pip install git+https://github.com/ecmwf/cdsapi.git
sudo /opt/local/bin/pip install git+https://github.com/fhs/pyhdf.git

# install the basemap-dev version [remove after basemap-v1.2 release]
# run "sudo port uninstall py36-matplotlib-basemap" if basemap was installed with port
sudo /opt/local/bin/pip install git+https://github.com/matplotlib/basemap.git#egg=mpl_toolkits
```

### Notes on [PyAPS](https://github.com/yunjunz/pyaps3) ###

+ We use PyAPS for tropospheric delay correction calculated from Global Atmospheric Models (GAMs) such as ERA-5, ERA-Interim, HRES-ECMWF, MERRA and NARR.

+ Check [Caltech Earthdef](http://earthdef.caltech.edu) for accounts setup information for various GAM datasets.

+ If you defined an environment variable named `WEATHER_DIR` to contain the path to a
directory, MintPy applications will download the GAM files into the indicated directory. Also MintPy
application will look for the GAM files in the directory before downloading a new one to prevent downloading
multiple copies if you work with different dataset that cover the same date/time.

### Notes on parallel processing ###

We use [Dask](https://www.dask.org) for parallel processing on High Performance Compute (HPC) cluster, it can be setup as below:

```
mkdir -p ~/.config/dask
cp $MINTPY_HOME/mintpy/defaults/dask_mintpy.yaml ~/.config/dask/dask_mintpy.yaml
```

Edit `~/.config/dask/dask_mintpy.yaml` file according to your HPC settings. Currently, only `LSFCluster` job scheduler is tested, `PBSCluster` should also work after minor adjustment in `ifgram_inversion.py`.

### Notes on vim ###

[Here](https://github.com/yunjunz/macOS_Setup/blob/master/vim.md) is some useful setup of Vim editor for general use and Python.
