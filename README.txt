pysar is an InSAR time-series analysis package. To use the package add the path to PySAR directory to your $PYTHONPATH and add PySAR/pysar to your $path

Depending on your shell you may use command such as following example to setup pysar:

export PYTHONPATH=/nethome/hfattahi/development/PySAR/pysar:${PYTHONPATH}
export PATH="/nethome/hfattahi/development/PySAR:$PATH"
export TSSARDIR=/nethome/timeseries/
OR:

setenv PYTHONPATH "/nethome/hfattahi/development/PySAR" 
set path = (/nethome/hfattahi/development/PySAR/pysar $path)
setenv TSSARDIR "/nethome/timeseries/"

Run pysarApp.py to see the examples of processing options.

##########################################################

The current version of PySAR is compatible with roi_pac outputs. pysar reads
unwrapped interefrograms (at the same coordinate system: radar or geo) and 
the baseline files for each interefrogram. You need to give the path to where
the interferograms are and pysar takes care of the rest!

Run pysarApp.py to see examples of processing options. 

##########################################################
Ho to run pysar:

When you have a stack of interferograms processed with roi_pac, make a pysar processing file (a text file) in your shell using for example vi or any other text editor:

eg: vi YourProjectName.template

and include the followng pysar processing options in your template:

#**********************

pysar.inputdata=/scratch/hfattahi/PROCESS/SanAndreasT356EnvD/DONE/IFG*/filt*0*c10.unw
pysar.CorFiles = /scratch/hfattahi/PROCESS/SanAndreasT356EnvD/DONE/IFG*/filt*0*.cor
pysar.wraped = /scratch/hfattahi/PROCESS/SanAndreasT356EnvD/DONE/IFG*/filt*0*.int
pysar.geomap = /scratch/hfattahi/PROCESS/SanAndreasT356EnvD/GEO/geomap_12/geomap_8rlks.trans
pysar.dem = /scratch/hfattahi/PROCESS/SanAndreasT356EnvD/DONE/IFG_20050102_20070809/radar_8lks.hgt
pysar.topo_error = yes # [no]
pysar.orbit_error = yes # [np]
pysar.orbit_error.method = plane  #['quadratic', 'plane', 'quardatic_range', 'quadratic_azimiuth', 'plane_range', 'plane_azimuth','baselineCor','BaseTropCor']
pysar.mask=yes
pysar.mask.threshold = 0.7

#**********************

Save your template file and run pysar as:
pysarApp.py YourProjectName.template

##########################################################

You may need to install some more packages including, pyaps, pykml, GDAL to get full advantage of PySAR. Basic time-series analysis does not need these packages though. However you need python with numpy, scipy, h5py and matplotlib installed.

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


