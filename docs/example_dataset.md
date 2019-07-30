Run the following command to test the demo datasets:

```
cd ${PROJECT_NAME}/mintpy
smallbaselineApp.py ${PROJECT_NAME}.txt
```

#### Sentinel-1 on Fernandina with ISCE ####

Area: Fernandina volcano at Galápagos Islands, Ecuador     
Data: Sentinel-1 A/B descending track 128 during Dec 2014 - June 2018 (98 acquisitions; [Zenodo](https://zenodo.org/record/2748487))      
Size: ~750 MB     

```
wget https://zenodo.org/record/2748487/files/FernandinaSenDT128.tar.xz
tar -xvJf FernandinaSenDT128.tar.xz
cd FernandinaSenDT128/mintpy
smallbaselineApp.py ${MINTPY_HOME}/docs/examples/input_files/FernandinaSenDT128.txt     
```

<p align="left">
  <img width="650" src="https://yunjunzhang.files.wordpress.com/2019/06/fernandinasendt128_poi.jpg">
</p>

#### Envisat of the 2008 Wells earthquake with Gamma ####

Area: Wells, Nevada, USA       
Data: Envisat ASAR descending track 399 during July 2007 - September 2008 (11 acquisitions; [Zenodo](https://zenodo.org/record/2748560))      
Size: ~280 MB      

```
wget https://zenodo.org/record/2748560/files/WellsEnvD2T399.tar.xz
tar -xvJf WellsEnvD2T399.tar.xz
cd WellsEnvD2T399/mintpy
smallbaselineApp.py ${MINTPY_HOME}/docs/examples/input_files/WellsEnvD2T399.txt
```

<p align="left">
  <img width="650" src="https://yunjunzhang.files.wordpress.com/2019/06/wellsenvd2t399_co_poi.jpg">
</p>

#### ALOS-1 on Kuju with ROI_PAC ####

Area: Kuju volcano at Kyushu island, SW Japan     
Data: ALOS-1 PALSAR ascending track 422 during January 2007 - January 2011 (24 acquisitions; [Zenodo](https://zenodo.org/record/2748170))     
Size: ~240 MB

```
wget https://zenodo.org/record/2748170/files/KujuAlosAT422F650.tar.xz
tar -xvJf KujuAlosAT422F650.tar.xz
cd KujuAlosAT422F650/mintpy
smallbaselineApp.py ${MINTPY_HOME}/docs/examples/input_files/KujuAlosAT422F650.txt     
```

<p align="left">
  <img width="650" src="https://yunjunzhang.files.wordpress.com/2019/06/kujualosat422f650_vel.jpg">
</p>

<br><br>
Check the auto plotted figures in **mintpy/pic** folder, and modify the plotting parameters to adjust your plotting, and re-run the plotting script to update your plotting result:   

```
./plot_smallbaselineApp.sh
```
