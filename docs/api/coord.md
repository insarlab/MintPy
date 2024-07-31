There are two coordination systems in MintPy: **radar coordinate** and **geo coordinate**. Geo coordinate is defined in WGS84 coordination for horizontal direction, and determined by the following [ROI_PAC attributes](https://github.com/insarlab/MintPy/wiki/Attributes) in latitude and longitude. The following shows examples from *KujuAlosAT422F650/mintpy/geo/geo_velocity.h5*:

```
X_FIRST    131.02409876
Y_FIRST    33.63756779
X_STEP     0.00033333
Y_STEP     -0.00033333
X_UNIT     degrees
Y_UNIT     degrees
```

X/Y_FIRST are the longitude/latitude value of the first (upper left corner) pixel’s upper left corner, as shown below:

<img src="https://insarlab.github.io/figs/docs/mintpy/coord_index.png" width="400">
