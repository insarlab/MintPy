## Custom colormaps

MintPy support the following colormaps:

+ [Matplotlib colormaps](https://matplotlib.org/stable/tutorials/colors/colormaps.html)
+ Custom colormaps: `cmy`, `dismph`, and `romanian`
+ Custom colormaps in **.cpt** (color palette tables) format. To add your own colormap, drop the corresponding .cpt file in `$MINTPY/mintpy/data/colormaps`.

We recommend to use cyclic colormap `cmy` for wrapped phase/displacement measurement.

<p align="left">
  <img width="280" src="https://insarlab.github.io/figs/docs/mintpy/cmap_cmy.png">
</p>

To use colormap `cmy` in view.py:

```bash
view.py velocity.h5 -c cmy
```

To use colormap `cmy` in python:

```python
from mintpy.colors import ColormapExt
cmap = ColormapExt('cmy').colromap
```

### Colormaps from [GMT](http://www.soest.hawaii.edu/gmt/) ###

All GMT cpt files, e.g. the 20 built-in colormaps shown below, can be recognized by setting the variable `GMT_CPT_DIR` in `$MINTPY_HOME/src/mintpy/objects/colors.py`. The default hardwired value is `/opt/local/share/gmt/cpt` for macOS users with GMT installed using [MacPorts](https://www.macports.org).

<p align="left">
  <img width="600" src="https://docs.generic-mapping-tools.org/5.4/_images/GMT_App_M_1a.png">
  <img width="600" src="https://docs.generic-mapping-tools.org/5.4/_images/GMT_App_M_1b.png">
</p>

### Colormaps from [cpt-city](http://soliton.vm.bytemark.co.uk/pub/cpt-city/views/totp-cpt.html) ###

The following colormaps is included by default:

+ BlueWhiteOrangeRed
+ DEM_print
+ differences
+ GMT_haxby
+ GMT_no_green
+ seminf-haxby
+ temp-c
+ temperature
+ wiki-2.0
+ wiki-schwarzwald-d050
+ wiki-scotland
+ More at [cpt-city](http://soliton.vm.bytemark.co.uk/pub/cpt-city/views/totp-cpt.html)

### Colormaps from [Scientific Color-Maps](http://www.fabiocrameri.ch/colourmaps.php) by Fabio Crameri ###

The following colormaps is included by default:

+ batlow (the scientific rainbow)
+ hawaii
+ oleron (surface topography)
+ roma (seismic tomography)
+ vik (diverging)
+ vikO (cyclic diverging)
+ More at [Scientific Color-Maps](http://www.fabiocrameri.ch/colourmaps.php) ([Crameri, 2018](https://doi.org/10.5194/gmd-11-2541-2018))

<p align="left">
  <img src="https://insarlab.github.io/figs/docs/mintpy/cmap_scientific_colour_maps_fabiocrameri.jpg">
</p>

### Interactive [web tool](https://jdherman.github.io/colormap/) to generate custom colormaps by Jon Herman ###

This web tool creates a custom colormap (for Matplotlib/Matlab) by dragging points on the RGB intensity curves.

+ Choose output as *plaintext* style and *RGB* format.
+ Copy and save the RGB table to a text file
+ Use this script [rgb2cpt.py](https://github.com/yuankailiu/utils/blob/main/trivia/rgb2cpt.py) to convert the RGB table to 8-column CPT file with heading, overrule background, foreground, and NaN colors.
