## Custom colormaps

MintPy support all colormaps from [Matplotlib](https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html) and custom colormap files in **.cpt** (color palette tables) format. To add your own colormap, drop the corresponding .cpt file in `$MINTPY/docs/resources/colormaps`.

To use vik colormap in view.py for example:

```bash
view.py velocity.h5 -c vik
```

### Colormaps from [GMT](http://www.soest.hawaii.edu/gmt/) ###

For macOS users, if GMT is installed using [MacPorts](https://www.macports.org), all GMT cpt files (located at `/opt/local/share/gmt/cpt`) will be recognized and supported, i.e. the 20 built-in colormaps as shown below:

<p align="left">
  <img width="600" src="https://www-k12.atmos.washington.edu/~ovens/gmt/doc/html/GMT_Docs/img277.png">
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
+ More at [Scientific Color-Maps](http://www.fabiocrameri.ch/colourmaps.php)

<p align="left">
  <img width="600" src="http://www.fabiocrameri.ch/resources/ScientificColourMaps_FabioCrameriCompact.png">
</p>
