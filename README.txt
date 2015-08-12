pysar is an InSAR time-series analysis package. To use the package add the path to PySAR directory to your $PYTHONPATH and add PySAR/pysar to your $path

Depending on your shell you may use command such as following example to setup pysar:

export PYTHONPATH=/nethome/hfattahi/development/PySAR/pysar:${PYTHONPATH}
export PATH="/nethome/hfattahi/development/PySAR:$PATH"

OR:

setenv PYTHONPATH "/nethome/hfattahi/development/PySAR" 
set path = (/nethome/hfattahi/development/PySAR/pysar $path)

Run pysarApp.py to see the examples of processing options.

##########################################################

The current version of PySAR is compatible with roi_pac outputs. pysar reads
unwrapped interefrograms (at the same coordinate system: radar or geo) and 
the baseline files for each interefrogram. You need to give the path to where
the interferograms are and pysar takes care of the rest!

Run pysarApp.py to see examples of processing options. 

##########################################################

pykml installation:

website:http://pythonhosted.org/pykml/

wget https://pypi.python.org/packages/source/p/pykml/pykml-0.1.0.tar.gz
tar -xvf pykml-0.1.0.tar.gz
cd pykml-0.1.0
easy_install pykml

##########################################################
GDAL installation:
%%%%%%%%%%%
wget ftp://ftp.remotesensing.org/gdal/gdal-1.9.1.tar.gz

./configure --with-python --prefix=/nethome/hfattahi/development/utilities/gdal-1.9.1
make 
make install

%%%%%%%%%%%%
setenv GDALHOME /nethome/hfattahi/development/utilities/gdal-1.9.1
set path= ( $path $GDALHOME/bin )
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${GDALHOME}/lib

##########################################################
PySAR uses cvxopt-1.1.6 for L1 norm minimization
See http://cvxopt.org to download and installation

This package is used if user choose to use L1 norm minimization 
for inversion of interferograms or to estimate the velocity field.

link to download:
https://github.com/cvxopt/cvxopt/archive/1.1.6.tar.gz

To install:
Untar the package
cd cvxopt-1.1.6
python setup.py install
##########################################################


