Here are example interferogram stacks pre-processed using different InSAR processors.

### Sentinel-1 on Fernandina with ISCE2/topsStack ###

+ Area: Fernandina volcano at Galápagos Islands, Ecuador
+ Data: Sentinel-1 A/B descending track 128 during Dec 2014 - June 2018 (98 acquisitions; [Zenodo](https://zenodo.org/record/3952953))
+ Size: ~750 MB

```bash
wget https://zenodo.org/record/3952953/files/FernandinaSenDT128.tar.xz
tar -xvJf FernandinaSenDT128.tar.xz
cd FernandinaSenDT128/mintpy
smallbaselineApp.py ${MINTPY_HOME}/docs/templates/FernandinaSenDT128.txt
```

<p align="left">
  <img width="650" src="https://github.com/insarlab/design-docs/blob/main/docs/FernandinaSenDT128-ISCE2.jpg">
</p>

Relevant literature:

+ Yunjun, Z., H. Fattahi, and F. Amelung (2019), Small baseline InSAR time series analysis: Unwrapping error correction and noise reduction, _Computers & Geosciences, 133,_ 104331, doi:10.1016/j.cageo.2019.104331.

### Sentinel-1 on San Francisco Bay with ARIA ###

+ Area: San Francisco Bay, California, USA
+ Data: Sentinel-1 A/B descending track 42 during May 2015 - March 2020 (114 acquisitoins; [Zenodo](https://zenodo.org/record/4265413))
+ Size: ~2.7 GB

```bash
wget https://zenodo.org/record/4265413/files/SanFranSenDT42.tar.xz
tar -xvJf SanFranSenDT42.tar.xz
cd SanFranSenDT42/mintpy
smallbaselineApp.py ${MINTPY_HOME}/docs/templates/SanFranSenDT42.txt
```

<p align="left">
  <img width="650" src="https://github.com/insarlab/design-docs/blob/main/docs/SanFranSenDT42-ARIA.jpg">
</p>

Relevant literature:

+ Chaussard, E., R. Bürgmann, H. Fattahi, R. M. Nadeau, T. Taira, C. W. Johnson, and I. Johanson (2015), Potential for larger earthquakes in the East San Francisco Bay Area due to the direct connection between the Hayward and Calaveras Faults, _Geophysical Research Letters,_ 42(8), 2734-2741, doi:10.1002/2015GL063575.

### Sentinel-1 of the 2019 Ridgecrest, California earthquake sequence with ASF HyP3 ###

+ Area: Owens Valley, California, USA ([USGS event page](https://earthquake.usgs.gov/earthquakes/eventpage/ci38457511/executive))
+ Data: Sentinel-1 descending track 71 during June - August 2019 (7 acquisitions; [Zenodo](https://zenodo.org/record/11049257))
+ Size: ~240 MB

```bash
wget https://zenodo.org/record/11049257/files/RidgecrestSenDT71.tar.xz
tar -xvJf RidgecrestSenDT71.tar.xz
cd RidgecrestSenDT71
smallbaselineApp.py ${MINTPY_HOME}/docs/templates/RidgecrestSenDT71.txt
```

### Sentinel-1 on San Francisco Bay with GMTSAR ###

+ Area: San Francisco Bay, California, USA
+ Data: Sentinel-1 A/B descending track 42 during December 2014 - June 2024 (333 acquisitoins; [Zenodo](https://zenodo.org/records/12773014))
+ Size: ~2.3 GB

```bash
wget https://zenodo.org/records/12773014/files/SanFranBaySenD42.tar.xz
tar -xvJf SanFranBaySenD42.tar.xz
cd SanFranBaySenD42
smallbaselineApp.py ${MINTPY_HOME}/docs/templates/SanFranBaySenD42.txt
```

<p align="left">
  <img width="600" src="https://github.com/insarlab/design-docs/blob/main/docs/SanFranBaySenD42-GMTSAR.jpg">
</p>

Relevant literature:

+ Chaussard, E., R. Bürgmann, H. Fattahi, R. M. Nadeau, T. Taira, C. W. Johnson, and I. Johanson (2015), Potential for larger earthquakes in the East San Francisco Bay Area due to the direct connection between the Hayward and Calaveras Faults, _Geophysical Research Letters,_ 42(8), 2734-2741, doi:10.1002/2015GL063575.

### Envisat of the 2008 Wells, Nevada earthquake with Gamma ###

+ Area: Wells, Nevada, USA ([USGS event page](https://earthquake.usgs.gov/earthquakes/eventpage/nn00234425/executive))
+ Data: Envisat ASAR descending track 399 during July 2007 - September 2008 (11 acquisitions; [Zenodo](https://zenodo.org/record/3952950))
+ Size: ~280 MB

```bash
wget https://zenodo.org/record/3952950/files/WellsEnvD2T399.tar.xz
tar -xvJf WellsEnvD2T399.tar.xz
cd WellsEnvD2T399/mintpy
smallbaselineApp.py ${MINTPY_HOME}/docs/templates/WellsEnvD2T399.txt
```

<p align="left">
  <img width="650" src="https://github.com/insarlab/design-docs/blob/main/docs/WellsEnvD2T399-Gamma.jpg">
</p>

Relevant literature:

+ Nealy, J. L., H. M. Benz, G. P. Hayes, E. A. Bergman, and W. D. Barnhart (2017), The 2008 Wells, Nevada, Earthquake Sequence: Source Constraints Using Calibrated Multiple‐Event Relocation and InSAR, _Bulletin of the Seismological Society of America_, 107(3), 1107-1117, doi:10.1785/0120160298.

### Sentinel-1 on Western Cape, South Africa with SNAP ###

+ Area: West coast of Western Cape province, South Africa
+ Data: Sentinel-1 ascending track 29 during March - June 2019 (10 acquisitions; [Zenodo](https://zenodo.org/record/4318134))
+ Size: ~520 MB

```bash
wget https://zenodo.org/record/4318134/files/WCapeSenAT29.tar.xz
tar -xvJf WCapeSenAT29.tar.xz
cd WCapeSenAT29
smallbaselineApp.py ${MINTPY_HOME}/docs/templates/WCapeSenAT29.txt
```

### ALOS on Kuju with ROI_PAC ###

+ Area: Kuju volcano at Kyushu island, SW Japan
+ Data: ALOS PALSAR ascending track 422 during January 2007 - January 2011 (24 acquisitions; [Zenodo](https://zenodo.org/record/3952917))
+ Size: ~240 MB

```bash
wget https://zenodo.org/record/3952917/files/KujuAlosAT422F650.tar.xz
tar -xvJf KujuAlosAT422F650.tar.xz
cd KujuAlosAT422F650/mintpy
smallbaselineApp.py ${MINTPY_HOME}/docs/templates/KujuAlosAT422F650.txt
```

<p align="left">
  <img width="650" src="https://github.com/insarlab/design-docs/blob/main/docs/KujuAlosAT422F650-ROIPAC.jpg">
</p>

Relevant literature:

+ Nakaboh, M., H. Ono, M. Sako, Y. Sudo, T. Hashimoto, and A. W. Hurst (2003), Continuing deflation by fumaroles at Kuju Volcano, Japan, _Geophysical Research Letters_, 30(7), doi:10.1029/2002gl016047.
+ Ishitsuka, K., T. Tsuji, T. Matsuoka, J. Nishijima, and Y. Fujimitsu (2016), Heterogeneous surface displacement pattern at the Hatchobaru geothermal field inferred from SAR interferometry time-series, _International Journal of Applied Earth Observation and Geoinformation_, 44, 95-103, doi:10.1016/j.jag.2015.07.006.
