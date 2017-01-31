## Welcome to PySAR!

PySAR is an InSAR time-series package to produce three dimensional (space and time) ground displacement from InSAR data. To use the package add the path to PySAR directory to your $PYTHONPATH and add PySAR/pysar to your $path   

Depending on your shell you may use commands such as the following examples to setup pysar:   

Using bash:   
export PYTHONPATH=/nethome/hfattahi/development/PySAR:${PYTHONPATH}   
export PATH="/nethome/hfattahi/development/PySAR/pysar:$PATH"   
export TSSARDIR=/nethome/timeseries/   

Using csh:   
setenv PYTHONPATH "/nethome/hfattahi/development/PySAR"    
set path = (/nethome/hfattahi/development/PySAR/pysar $path)   
setenv TSSARDIR "/nethome/timeseries/"   

Run pysarApp.py to see the examples of processing options.   

##########################################################   

The current version of PySAR is compatible with roi_pac outputs. pysar reads unwrapped interefrograms (at the same coordinate system: radar or geo) and the baseline files for each interefrogram. You need to give the path to where the interferograms are and pysar takes care of the rest!   

Run pysarApp.py to see examples of processing options.    

##########################################################   
How to run pysar:   

When you have a stack of interferograms processed with roi_pac, make a pysar processing file (a text file) in your shell using for example vi or any other text editor:   

eg: vi YourProjectName.template   

and include the following pysar processing options in your template:   

########################   

pysar.inputdata=/scratch/hfattahi/PROCESS/SanAndreasT356EnvD/DONE/IFG*/filt*0*c10.unw   
pysar.CorFiles = /scratch/hfattahi/PROCESS/SanAndreasT356EnvD/DONE/IFG*/filt*0*.cor   
pysar.wrapped = /scratch/hfattahi/PROCESS/SanAndreasT356EnvD/DONE/IFG*/filt*0*.int   
pysar.geomap = /scratch/hfattahi/PROCESS/SanAndreasT356EnvD/GEO/geomap_12/geomap_8rlks.trans   
pysar.dem = /scratch/hfattahi/PROCESS/SanAndreasT356EnvD/DONE/IFG_20050102_20070809/radar_8lks.hgt   
pysar.topo_error = yes # [no]   
pysar.orbit_error = yes # [np]   
pysar.orbit_error.method = plane  #['quadratic', 'plane', 'quardatic_range', 'quadratic_azimiuth', 'plane_range',    'plane_azimuth','baselineCor','BaseTropCor']   
pysar.mask=yes   
pysar.mask.threshold = 0.7   

########################   

Save your template file and run pysar as:   
pysarApp.py YourProjectName.template   

pysar reads the unwrapped interferograms, refernces all of them to the same coherent pixel (a seed point point), calculates the phase closure and estimates the unwrapping errors (if it has been asked for), inverts the interferograms, calculates a parameter called "temporal_coherence" which can be used to evaluate the quality of inversion, removes ramps or surface from time-series epochs, corrects dem errors, corrects local oscilator drift (for Envisat only), corrects stratified tropospheric delay (using pyaps and using phase-elevation approach), ... and finally estimates the velocity.   

use view.py to view any pysar output.   
use tsviewer.py to plot the time-series for each point (relative to the refernce point and epoch!).    

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


Github Page: https://yunjunz.github.io/PySAR/     
