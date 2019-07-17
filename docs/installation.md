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

### 1. Setup Paths ###

To use the package, you need to setup the environment:

+ To make mintpy importable in python, by adding the path to MintPy directory to your _$PYTHONPATH_    
+ To make utility scripts available in command line, by adding _${MINTPY_HOME}/mintpy_ and _${MINTPY_HOME}/sh_ to your _$path_.

Add to your **_~/.bash_profile_** file for _bash_ user:

```bash
############################  Python  ###############################
if [ -z ${PYTHONPATH+x} ]; then export PYTHONPATH=""; fi

##--------- MintPy ------------------##
export MINTPY_HOME=~/python/MintPy
export PYTHONPATH=${PYTHONPATH}:${MINTPY_HOME}
export PATH=${PATH}:${MINTPY_HOME}/mintpy:${MINTPY_HOME}/sh

##--------- PyAPS ------------------##
export PYAPS_HOME=~/python/PyAPS
export PYTHONPATH=${PYAPS_HOME}:${PYTHONPATH}
```

Source the file for the first time. It will be sourced automatically next time when you login. [Here](https://github.com/yunjunz/macOS_Setup/blob/master/cshrc.md) is an example _.cshrc_ file for _csh/tcsh_ user.

### 2. Install Python dependecies ###

MintPy is written in Python3 and relies on several Python modules, check the [requirements.txt](https://github.com/insarlab/MintPy/blob/master/docs/requirements.txt) file for details. We recommend using [conda](https://conda.io/miniconda.html) or [macports](https://www.macports.org/install.php) to install the python environment and the prerequisite packages, because of the convenient managenment and default [performance setting with numpy/scipy](http://markus-beuckelmann.de/blog/boosting-numpy-blas.html) and [pyresample](https://pyresample.readthedocs.io/en/latest/installation.html#using-pykdtree).

#### Installing via conda ####

Add to your **_~/.bash_profile_** file:

```bash
export PYTHON3DIR=~/python/miniconda3
export PATH=${PATH}:${PYTHON3DIR}/bin
export PROJ_LIB=`python3 -c "import pyproj; print(pyproj.pyproj_datadir)"`   #Temporary fix for basemap import error
```

Run the following in your terminal:

```
cd ~/python

# download MintPy
git clone https://github.com/insarlab/MintPy.git

# download PyAPS
mkdir -p PyAPS; cd PyAPS
git clone https://github.com/yunjunz/pyaps3.git

# install miniconda
wget https://repo.continuum.io/miniconda/Miniconda3-4.5.4-MacOSX-x86_64.sh
chmod +x Miniconda3-4.5.4-MacOSX-x86_64.sh
./Miniconda3-4.5.4-MacOSX-x86_64.sh -b -p $PYTHON3DIR

# install dependencies with conda
$PYTHON3DIR/bin/conda config --add channels conda-forge
$PYTHON3DIR/bin/conda install --yes --file $MINTPY_HOME/docs/conda.txt

# install pykml
$PYTHON3DIR/bin/pip install git+https://github.com/yunjunz/pykml.git
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
cd ~/python

# download MintPy
git clone https://github.com/insarlab/MintPy.git

# download PyAPS
mkdir -p PyAPS; cd PyAPS
git clone https://github.com/yunjunz/pyaps3.git

# install dependencies with macports
# use "port -N install" to use the safe default for prompt questions
sudo port install $(cat $MINTPY_HOME/docs/ports.txt)

# install pykml
git clone https://github.com/yunjunz/pykml.git; cd pykml
python3 setup.py build
python3 setup.py install

# intall CDS API for PyAPS with ERA-5
wget https://github.com/ecmwf/cdsapi/archive/v0.1.4.tar.gz -O cdsapi-0.1.4.tar.gz
tar -xvf cdsapi-0.1.4.tar.gz; cd cdsapi-0.1.4
python3 setup.py build
python3 setup.py install

# install PyHDF for PyAPS with MERRA
git clone https://github.com/fhs/pyhdf.git; cd pyhdf
python3 setup.py build
python3 setup.py install
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

### Notes on pyKML ###

The pykml installed through conda/pip supports python2 only, the python2/3 compatible version is available [here](https://github.com/yunjunz/pykml.git) and installed throught the commands above by default.

### Notes on vim ###

[Here](https://github.com/yunjunz/macOS_Setup/blob/master/vim.md) is some useful setup of Vim editor for general use and Python.
