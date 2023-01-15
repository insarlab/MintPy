## Plate boundary files

### GSRM

The Global Strain Rate Map (GSRM) version 2.1 (http://geodesy.unr.edu/GSRM/) is defined in Kreemer et al. (2014). It includes 50 rigid plates. The no-net rotation (NNR) Euler poles for the IGS08 reference frame of all plates are in Table S2 of Kreemer et al. (2014). The boundary file of all plates are in GMT (lon, lat) format and is downloaded as below:

```bash
mkdir -p GSRM
wget http://geodesy.unr.edu/GSRM/GSRM_plate_outlines.gmt -O GSRM/plate_outlines.lola
```

### MORVEL

The Mid-Ocean Ridge VELocity (MORVEL) (http://www.geology.wisc.edu/~chuck/MORVEL/) is defined in DeMets et al. (2010) with the NNR plate motion model (NNR-MORVEL56) defined in Argus et al. (2011). It includes 56 plates. The NNR Euler poles are in Table 1 of Argus et al. (2011). The boundary file of all plates are in GMT (lat, lon) format and is downloaded as below:

```bash
mkdir -p MORVEL
wget http://www.geology.wisc.edu/~chuck/MORVEL/NnrMRVL_PltBndsLatLon.zip
unzip -p NnrMRVL_PltBndsLatLon.zip All_boundaries >MORVEL/plate_outlines.lalo
```

### References

+   Argus, D. F., Gordon, R. G., & DeMets, C. (2011). Geologically current motion of 56 plates relative to the no-net-rotation reference frame. _Geochemistry, Geophysics, Geosystems, 12_(11). https://doi.org/10.1029/2011GC003751
+   Bird, P. (2003). An updated digital model of plate boundaries. _Geochemistry, Geophysics, Geosystems, 4_(3). https://doi.org/10.1029/2001GC000252
+   DeMets, C., Gordon, R. G., & Argus, D. F. (2010). Geologically current plate motions. _Geophysical journal international, 181_(1), 1-80. https://doi.org/10.1111/j.1365-246X.2009.04491.x
+   Kreemer, C., Blewitt, G., & Klein, E. C. (2014). A geodetic plate motion and Global Strain Rate Model. _Geochemistry, Geophysics, Geosystems, 15_(10), 3849-3889. https://doi.org/10.1002/2014GC005407
