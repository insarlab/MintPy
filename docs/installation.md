## Install MintPy

### Install the released version via conda

Mintpy is available on the [conda-forge](https://anaconda.org/conda-forge/mintpy) channel. The latest released version can be installed via `conda` as:

```bash
conda install -c conda-forge mintpy
```

### Install the released version via docker

Docker allows one to run MintPy in a dedicated container (essentially an efficient virtual machine) and to be independent of platform OS. After installing [docker](https://docs.docker.com/install/), run the following to pull the [MintPy container from DockerHub](https://hub.docker.com/r/forrestwilliams/mintpy) to your local machine, check more details at [here](docker.md).

```bash
docker pull forrestwilliams/mintpy:1.3.1
```

### Install the development version

The installation note below is tested on Linux and macOS, and is still experimental on Windows (may has bugs).

MintPy is written in Python 3 and relies on several Python modules, check the [requirements.txt](https://github.com/insarlab/MintPy/blob/main/docs/requirements.txt) file for details. We recommend using [conda](https://docs.conda.io/en/latest/miniconda.html) or [macports](https://www.macports.org/install.php) to install the python environment and the prerequisite packages, because of the convenient management and default [performance setting with numpy/scipy](http://markus-beuckelmann.de/blog/boosting-numpy-blas.html) and [pyresample](https://pyresample.readthedocs.io/en/latest/installation.html#using-pykdtree). You can control the number of threads used by setting the _environment variables_, e.g. `OMP_NUM_THREADS`.



### Notes for Mac users ###

Install Xcode with command line tools, if you have not already done so.

+ Install `Xcode` from App store

+ Install `command line tools` within XCode and agree to the terms of license.

  ```bash
  xcode-select --install -s /Applications/Xcode.app/Contents/Developer/
  sudo xcodebuild -license
  ```

+ Install [XQuartz](https://www.xquartz.org), then restart the terminal.


### 1. Download MintPy ###

Run the following in your terminal to download the development version of MintPy:

```bash
git clone https://github.com/insarlab/MintPy.git
```

### 2. Install dependencies and MintPy ###


#### a. via conda ####

Install [miniconda](https://docs.conda.io/en/latest/miniconda.html) if you have not already done so. You may need to close and restart the shell for changes to take effect.

```bash
# download and install miniconda
# use wget or curl to download in command line or click from the web browser
# curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -o Miniconda3-latest-MacOSX-x86_64.sh
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
bash Miniconda3-latest-MacOSX-x86_64.sh
```

You may set up a new conda environment (recommended) for MintPy by running:

```bash
conda env create -f MintPy/docs/environment.yml
conda activate mintpy
```

Or you can install the dependencies into an existing environment by running:

```bash
# Add "gdal'>=3'" below to install extra dependencies if you use ARIA, FRInGE, HyP3 or GMTSAR
# Add "isce2"     below to install extra dependencies if you use ISCE-2
conda install -c conda-forge --file ~/tools/MintPy/docs/requirements.txt
```

Then install MintPy into the environment by running:
```bash
python -m pip install MintPy
```

For development of MintPy, you may want to use an "editable" install so changes MintPy will be immediately reflected in the environment by running:
```bash
python -m pip install -e MintPy
```

#### b. via MacPorts ####

Install [macports](https://www.macports.org/install.php) if you have not done so. Add the following at the bottom of your `~/.bash_profile` file:

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

Run the following in your terminal in `bash` to install the dependencies:

```bash
# install dependencies with macports
# use "port -N install" to use the safe default for prompt questions
sudo port install $(cat MintPy/docs/ports.txt)

# install dependencies not available on macports: pysolid, pykml, pykdtree, pyresample, cdsapi, pyhdf
sudo -H /opt/local/bin/pip install git+https://github.com/insarlab/PySolid.git
sudo -H /opt/local/bin/pip install git+https://github.com/tylere/pykml.git
sudo -H /opt/local/bin/pip install git+https://github.com/storpipfugl/pykdtree.git
sudo -H /opt/local/bin/pip install git+https://github.com/pytroll/pyresample.git
sudo -H /opt/local/bin/pip install git+https://github.com/ecmwf/cdsapi.git
sudo -H /opt/local/bin/pip install git+https://github.com/fhs/pyhdf.git
```

Then install MintPy into the environment by running:
```bash
sudo -H /opt/local/bin/pip install MintPy
```

For development of MintPy, you may want to use an "editable" install so changes MintPy will be immediately reflected in the environment by running:
```bash
sudo -H /opt/local/bin/pip install -e MintPy
```

### Notes on [PySolid](https://github.com/insarlab/PySolid) ###

We use PySolid for solid Earth tides correction, which is available on conda-forge for Linux, macOS, and Windows. If the conda-forge version install from above does not work, run the following to compile from source:

```bash
# install Fortran compiler via conda
conda install -c conda-forge fortran-compiler
python -m pip install git+https://github.com/insarlab/PySolid
```

### Notes on [PyAPS](https://github.com/insarlab/PyAPS) ###

+ We use PyAPS (Jolivet et al., 2011; 2014) for tropospheric delay correction calculated from Global Atmospheric Models (GAMs) such as ERA-5, ERA-Interim, HRES-ECMWF, MERRA and NARR.

+ Check [Earthdef/PyAPS](http://earthdef.caltech.edu/projects/pyaps/wiki/Main#) for accounts setup information for ERA-Interim and MERRA.

+ Check [GitHub/PyAPS](https://github.com/insarlab/PyAPS) for account setup for ERA-5. **Make sure that you:**

  -   accept the data license in the Terms of use on ECMWF website and 
  -   run `examples/TestECMWF.ipynb` to test the data downloading and running.

+ If you defined an environment variable named `WEATHER_DIR` to contain the path to a directory, MintPy applications will download the GAM files into the indicated directory. Also, MintPy application will look for the GAM files in the directory before downloading a new one to prevent downloading multiple copies if you work with different dataset that cover the same date/time.
