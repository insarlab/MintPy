### Install PySAR

#### 0. For Mac users     

Install Xcode with command line tools, if you have not already done so.

1. Install Xcode from App store
2. Install command line tools from within XCode and agree to the terms of license.

```tcsh   
xcode-select --install -s /Applications/Xcode.app/Contents/Developer/ 
sudo xcodebuild -license 
```   

3. Install [XQuartz](https://www.xquartz.org), then restart the terminal.

4. Install [macports](https://www.macports.org/install.php) if you want to use Macports to install Python environment (if you use Anaconda instead of Macports, skip now to the next section "Setup Paths"); then update the port tree with the following command. If your network prevent the use of rsync or svn via http of port tree, try [Portfile Sync via a Snapshot Tarball](https://trac.macports.org/wiki/howto/PortTreeTarball).
```tcsh
sudo port selfupdate
```

    

5. Restart terminal

#### 1. Setup Paths    

To use the package, you need to setup the environment. Depending on your shell, you may use commands below to setup PySAR, by adding the following to your source file. They are for:   
1. To make pysar importable in python, by adding the path to PySAR directory to your _$PYTHONPATH_    
2. To make utility scripts available in command line, by adding _${PYSAR_HOME}/pysar_ and _${PYSAR_HOME}/sh_ to your _$path_.   
   
For csh/tcsh user, add to your **_~/.cshrc_** file for example:   

    ############################  Python  ###############################
    if ( ! $?PYTHONPATH ) then
        setenv PYTHONPATH ""
    endif
    
    ##--------- conda ------------------## 
    setenv PYTHON3DIR    ~/python/miniconda3
    setenv PATH          ${PATH}:${PYTHON3DIR}/bin
    setenv PROJ_LIB      ${PYTHON3DIR}/share/proj   #Temporary fix for basemap import error
    
    ##--------- PySAR ------------------## 
    setenv PYSAR_HOME    ~/python/PySAR       #for released version, "~/python/PySAR-0.4.0"
    setenv PYTHONPATH    ${PYTHONPATH}:${PYSAR_HOME}
    setenv PATH          ${PATH}:${PYSAR_HOME}/pysar:${PYSAR_HOME}/sh
   
For bash user, add to your **_~/.bashrc_** file for example:   

    ############################  Python  ###############################
    if [ -z ${PYTHONPATH+x} ]; then export PYTHONPATH=""; fi
    
    ##--------- conda ------------------## 
    export PYTHON3DIR=~/python/miniconda3
    export PATH=${PATH}:${PYTHON3DIR}/bin
    export PROJ_LIB=${PYTHON3DIR}/share/proj   #Temporary fix for basemap import error
    
    ##--------- PySAR ------------------## 
    export PYSAR_HOME=~/python/PySAR        #for released version, "~/python/PySAR-0.4.0"
    export PYTHONPATH=${PYTHONPATH}:${PYSAR_HOME}   
    export PATH=${PATH}:${PYSAR_HOME}/pysar:${PYSAR_HOME}/sh   

Source the file for the first time. It will be sourced automatically next time when you login. An example of .cshrc is [here](https://github.com/yunjunz/macOS_Setup/blob/master/cshrc.md).
   
   
#### 2. Install Python dependecies
PySAR is written in Python3 (3.5+) and it relies on several Python modules, check the [requirements.txt](./requirements.txt) file for details. We recommend using [conda](https://conda.io/miniconda.html) or [macports](https://www.macports.org/install.php) to install the python environment and the prerequisite packages, because of the convenient managenment and default [performance setting with numpy/scipy](http://markus-beuckelmann.de/blog/boosting-numpy-blas.html) and [pyresample](https://pyresample.readthedocs.io/en/latest/installation.html#using-pykdtree).


For conda user, run the following in your terminal in _bash/tcsh_:   

    cd ~/python
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
    chmod +x Miniconda3-latest-MacOSX-x86_64.sh
    ./Miniconda3-latest-MacOSX-x86_64.sh -b -p $PYTHON3DIR
    
    $PYTHON3DIR/bin/conda config --add channels conda-forge
    $PYTHON3DIR/bin/conda install --yes --file $PYSAR_HOME/docs/conda.txt
    
    git clone https://github.com/yunjunz/pykml.git; cd pykml
    $PYTHON3DIR/bin/python3 setup.py build     
    $PYTHON3DIR/bin/python3 setup.py install    
   
For macports user, install modules in the [ports.txt](https://github.com/yunjunz/PySAR/blob/master/docs/ports.txt) file as below in _bash_:       
    
    sudo port install $(cat $PYSAR_HOME/docs/ports.txt)
    
    git clone https://github.com/yunjunz/pykml.git; cd pykml
    python3 setup.py build     
    python3 setup.py install   

Note that pykml through conda supports python2 only, we provide a python2/3 compatible version [here](https://github.com/yunjunz/pykml.git) and installed throught the command line above by default.
  
We use [PyAPS](http://earthdef.caltech.edu/projects/pyaps/wiki/Main) for tropospheric delay correction calculated from weather re-analysis datadset such as ERA-Interim, MERRA and NARR. Check [caltech's website](http://earthdef.caltech.edu/projects/pyaps/wiki/Main) for the code download and account setup.

[Here](https://github.com/yunjunz/macOS_Setup/blob/master/vim.md) is some useful setup of Vim editor for general use and Python.
