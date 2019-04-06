## Install PySAR

Tested on macOS and Linux, not sure about Windows.

### Notes for Mac users

Install Xcode with command line tools, if you have not already done so.

1. Install Xcode from App store

2. Install command line tools from within XCode and agree to the terms of license.

```
xcode-select --install -s /Applications/Xcode.app/Contents/Developer/ 
sudo xcodebuild -license 
```

3. Install [XQuartz](https://www.xquartz.org), then restart the terminal.

### 1. Setup Paths    

To use the package, you need to setup the environment:

1. To make pysar importable in python, by adding the path to PySAR directory to your _$PYTHONPATH_    
2. To make utility scripts available in command line, by adding _${PYSAR_HOME}/pysar_ and _${PYSAR_HOME}/sh_ to your _$path_.   

Add to your **_~/.bash_profile_** file for bash user:

```bash
############################  Python  ###############################
if [ -z ${PYTHONPATH+x} ]; then export PYTHONPATH=""; fi

##--------- PySAR ------------------## 
export PYSAR_HOME=~/python/PySAR
export PYTHONPATH=${PYTHONPATH}:${PYSAR_HOME}   
export PATH=${PATH}:${PYSAR_HOME}/pysar:${PYSAR_HOME}/sh   
```

Source the file for the first time. It will be sourced automatically next time when you login. [Here](https://github.com/yunjunz/macOS_Setup/blob/master/cshrc.md) is an example _.cshrc_ file for csh/tcsh user.


### 2. Install Python dependecies
PySAR is written in Python3 (3.5+) and it relies on several Python modules, check the [requirements.txt](./requirements.txt) file for details. We recommend using [conda](https://conda.io/miniconda.html) or [macports](https://www.macports.org/install.php) to install the python environment and the prerequisite packages, because of the convenient managenment and default [performance setting with numpy/scipy](http://markus-beuckelmann.de/blog/boosting-numpy-blas.html) and [pyresample](https://pyresample.readthedocs.io/en/latest/installation.html#using-pykdtree).

#### Installing via conda

Add to your **_~/.bash_profile_** file:

```bash
export PYTHON3DIR=~/python/miniconda3
export PATH=${PATH}:${PYTHON3DIR}/bin
export PROJ_LIB=${PYTHON3DIR}/share/proj   #Temporary fix for basemap import error
```

Run the following in your terminal:

```
cd ~/python
wget https://repo.continuum.io/miniconda/Miniconda3-4.5.4-MacOSX-x86_64.sh
chmod +x Miniconda3-4.5.4-MacOSX-x86_64.sh
./Miniconda3-4.5.4-MacOSX-x86_64.sh -b -p $PYTHON3DIR

$PYTHON3DIR/bin/conda config --add channels conda-forge
$PYTHON3DIR/bin/conda install --yes --file $PYSAR_HOME/docs/conda.txt

git clone https://github.com/yunjunz/pykml.git; cd pykml
$PYTHON3DIR/bin/python3 setup.py build
$PYTHON3DIR/bin/python3 setup.py install
```

#### Installing via MacPorts

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
sudo port install $(cat $PYSAR_HOME/docs/ports.txt)

git clone https://github.com/fhs/pyhdf.git; cd pyhdf
python3 setup.py build
python3 setup.py install

git clone https://github.com/yunjunz/pykml.git; cd pykml
python3 setup.py build
python3 setup.py install
```

#### Notes on pyKML

The pykml installed through conda/pip supports python2 only, the python2/3 compatible version is available [here](https://github.com/yunjunz/pykml.git) and installed throught the command line above by default.

#### Notes on [PyAPS](http://earthdef.caltech.edu/projects/pyaps/wiki/Main)

+ We use PyAPS for tropospheric delay correction calculated from Global Atmospheric Models (GAMs) such as ERA-Interim, MERRA and NARR. Check [Caltech Earthdef](http://earthdef.caltech.edu) for the code download and account setup information for various GAM datasets.    

+ PyAPS supports python2 only, thus you will need to install a python2 environment to be able to use it. py3-PySAR is compatible with py2-PyAPS. 

+ If you defined an environment variable named `WEATHER_DIR` to contain the path to a 
directory, PySAR applications will download the GAM files into the indicated directory. Also PySAR
application will look for the GAM files in the directory before downloading a new one to prevent downloading
multiple copies if you work with different dataset that cover the same date/time.

#### Notes on vim

[Here](https://github.com/yunjunz/macOS_Setup/blob/master/vim.md) is some useful setup of Vim editor for general use and Python.
