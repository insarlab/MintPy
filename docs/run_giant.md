## Notes to run GIAnT and compare it with MintPy

GIAnT is developed at Caltech, check the website for download information and more: [link](http://earthdef.caltech.edu/projects/giant/wiki).

     cd ~/insarlab/Galapagos/GalapagosSenDT128
     mkdir GIANT; cd GIANT
     save_ifg_list4giant.py  ../MintPy/INPUTS/ifgramStack.h5  --sensor SEN     
     cp ../ISCE/merged/interferograms/20150307_20150319/filt_fine.unw.rsc .
     
Edit _userfn.py_ file to setup the path for input unwrapped interferograms and coherence.    
Edit _prepxml.py_ file to setup the path for geometry files, cropped area for processing, reference area, and parameters for SBAS and minTS approach. Use absolute path.    
Then run prepxml.py to write XML files: data.xml, sbas.xml and mints.xml    
     
     ./prepxml.py     
     PrepImageStack.py
     ProcessStack.py
     SBASInvert.py
     NSBASInvert.py
     TimefnInvert.py
     
Then run prep_giant.py to prepare metadata needed for MintPy.     
     
     prep_giant.py  LS-PARAMS.h5 -x ../data.xml ../sbas.xml ../mints.xml
     prep_giant.py  TS-PARAMS.h5 -x ../data.xml ../sbas.xml ../mints.xml
     prep_giant.py  NSBAS-PARAMS.h5 -x ../data.xml ../sbas.xml ../mints.xml

Use `view.py` and `tsview.py` in MintPy for data visualiztion.
