## Install MintPy

Tested on macOS and Linux, not sure about Windows.

### Notes for Mac users ###

Install Xcode with command line tools, if you have not already done so.

+   Install `Xcode` from App store

+   Install `command line tools` within XCode and agree to the terms of license.

```
xcode-select --install -s /Applications/Xcode.app/Contents/Developer/
sudo xcodebuild -license
```

+   Install [XQuartz](https://www.xquartz.org), then restart the terminal.

### Notes for Docker users ###

Docker allows one to run MintPy in a dedicated container (essentially an efficient virtual machine) and to be independent of platform OS. After installing [docker](https://docs.docker.com/install/), run the following to pull the [MintPy container from DockerHub](https://hub.docker.com/r/andretheronsa/mintpy) to your local machine, check more details at [here](docker.md).

```
docker pull andretheronsa/mintpy:latest
```

### 1. Download and Setup ###

Run the following in your terminal to download the development version of MintPy and PyAPS:

```bash
cd ~/tools
git clone https://github.com/insarlab/MintPy.git
git clone https://github.com/yunjunz/PyAPS.git
```

Set the following environment variables in your source file (e.g. **_~/.bash_profile_** for _bash_ users or **_~/.cshrc_** for _csh/tcsh_ users).

```bash
if [ -z ${PYTHONPATH+x} ]; then export PYTHONPATH=""; fi

##--------- MintPy ------------------##
export MINTPY_HOME=~/tools/MintPy
export PATH=${PATH}:${MINTPY_HOME}/mintpy
export PYTHONPATH=${PYTHONPATH}:${MINTPY_HOME}:~/tools/PyAPS
```

### 2. Install dependencies ###

MintPy is written in Python3 and relies on several Python modules, check the [requirements.txt](https://github.com/insarlab/MintPy/blob/main/docs/requirements.txt) file for details. We recommend using [conda](https://docs.conda.io/en/latest/miniconda.html) or [macports](https://www.macports.org/install.php) to install the python environment and the prerequisite packages, because of the convenient managenment and default [performance setting with numpy/scipy](http://markus-beuckelmann.de/blog/boosting-numpy-blas.html) and [pyresample](https://pyresample.readthedocs.io/en/latest/installation.html#using-pykdtree). You can control the number of threads used by setting the _environment variables_, e.g. `OMP_NUM_THREADS`.

#### a. via conda ####

Install [miniconda](https://docs.conda.io/en/latest/miniconda.html) if you have not already done so. You may need to close and restart the shell for changes to take effect.

```bash
# download and install miniconda
# use wget or curl to download in command line or click from the web brower
# curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -o Miniconda3-latest-MacOSX-x86_64.sh
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
chmod +x Miniconda3-latest-MacOSX-x86_64.sh
./Miniconda3-latest-MacOSX-x86_64.sh -b -p ~/tools/miniconda3
~/tools/miniconda3/bin/conda init bash
```

Run the following in your terminal to install the dependencies to your conda environment (recommended). The default is _**base**_; a new custom environment is recommended.

```
# To create a new conda environment, e.g. named "insar", run "conda create --name insar; conda activate insar"

# Add "gdal>=3" below to install extra dependencies if you use ARIA, FRInGE, HyP3 or GMTSAR
# Add "isce2"   below to install extra dependencies if you use ISCE-2
conda install --yes -c conda-forge --file ~/tools/MintPy/docs/requirements.txt

$CONDA_PREFIX/bin/pip install git+https://github.com/insarlab/PySolid.git
$CONDA_PREFIX/bin/pip install git+https://github.com/tylere/pykml.git
```

Or run the following in your terminal to install the dependencies to a new environment _**mintpy**_:

```
conda env create -f $MINTPY_HOME/docs/environment.yml
conda activate mintpy
```

#### b. via MacPorts ####

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

Run the following in your terminal in _bash_ to install the dependencies:

```bash
# install dependencies with macports
# use "port -N install" to use the safe default for prompt questions
sudo port install $(cat $MINTPY_HOME/docs/ports.txt)

# install dependencies not available on macports: pysolid, pykml, pykdtree, pyresample, cdsapi, pyhdf
sudo -H /opt/local/bin/pip install git+https://github.com/insarlab/PySolid.git
sudo -H /opt/local/bin/pip install git+https://github.com/tylere/pykml.git
sudo -H /opt/local/bin/pip install git+https://github.com/storpipfugl/pykdtree.git
sudo -H /opt/local/bin/pip install git+https://github.com/pytroll/pyresample.git
sudo -H /opt/local/bin/pip install git+https://github.com/ecmwf/cdsapi.git
sudo -H /opt/local/bin/pip install git+https://github.com/fhs/pyhdf.git
```

### Notes on [PySolid](https://github.com/insarlab/PySolid) ###

We use PySolid for solid Earth tides correction. If the pre-compiled version install from above does not work, run the following to compile from source:

```bash
# install Fortran compiler via conda
conda install -c conda-forge fortran-compiler

# compile Fortran code into a Python interface using f2py
cd ~/tools/PySolid/pysolid
f2py -c -m solid solid.for
```

### Notes on [PyAPS](https://github.com/yunjunz/PyAPS) ###

+   We use PyAPS (Jolivet et al., 2011; 2014) for tropospheric delay correction calculated from Global Atmospheric Models (GAMs) such as ERA-5, ERA-Interim, HRES-ECMWF, MERRA and NARR.

+   Check [Earthdef/PyAPS](http://earthdef.caltech.edu/projects/pyaps/wiki/Main#) for accounts setup information for ERA-Interim and MERRA.

+   Check [GitHub/PyAPS](https://github.com/yunjunz/PyAPS) for account setup for ERA-5. **Make sure that you:**

    -   accept the data license in the Terms of use on ECMWF website and 
    -   run `examples/TestECMWF.ipynb` to test the data downloading and running.

+   If you defined an environment variable named `WEATHER_DIR` to contain the path to a
directory, MintPy applications will download the GAM files into the indicated directory. Also MintPy
application will look for the GAM files in the directory before downloading a new one to prevent downloading
multiple copies if you work with different dataset that cover the same date/time.
