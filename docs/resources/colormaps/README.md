## Custom colormaps

MintPy support the following colormaps:

+ [Matplotlib colormaps](https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html)
+ Custom colormaps: `cmy` and `dismph`
+ Custom colormaps in **.cpt** (color palette tables) format. To add your own colormap, drop the corresponding .cpt file in `$MINTPY/docs/resources/colormaps`.

We recommend to use cyclic colormap `cmy` for wrapped phase/displacement measurement.

<p align="left">
  <img width="280" src="https://yunjunzhang.files.wordpress.com/2020/01/cmap_cmy-1.png">
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

For macOS users, if GMT is installed using [MacPorts](https://www.macports.org), all GMT cpt files (located at `/opt/local/share/gmt/cpt`) will be recognized and supported, i.e. the 20 built-in colormaps as shown below:

<p align="left">
  <img width="600" src="https://docs.generic-mapping-tools.org/5.4/_images/GMT_App_M_1a.png">
  <img width="600" src="https://docs.generic-mapping-tools.org/5.4/_images/GMT_App_M_1b.png">
</p>

### Colormaps from [cpt-city](http://soliton.vm.bytemark.co.uk/pub/cpt-city/views/totp-cpt.html) ###

The following colormaps is included by default:

+ BlueWhiteOrangeRed
+ DEM_print
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
  <img width="600" src="http://www.fabiocrameri.ch/resources/ScientificColourMaps_FabioCrameriCompact.png">
</p>
