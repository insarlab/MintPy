#!/bin/sh

gmt set FONT_ANNOT_PRIMARY 16p
gmt set FORMAT_GEO_MAP dddF
gmt set MAP_FRAME_TYPE inside

ps='coast.ps'

#gmt pscoast -R113/147/21/47 -JM6i -Dl -B10wens -Na -W0.5 -Swhite -Ggray -P -V -K > $outName
#gmt pscoast -R113/147/21/47 -Js115/12/45/1:200000000 -Baf -Dl -N1 -W0.5 -Swhite -Ggray -P -V -K > $ps
gmt pscoast -Rd -JA-91.1/0.4/3i -Dl -Ba -A1000 -Ggray30 -P -V -K > $ps

#gmt pstext -R -J -F+f16,Helvetica-Oblique,black -V -O -K <<END>> $ps
#138.0    30.0    Pacific Ocean
#134.0    39.5    Sea of Japan
#END

#gmt pstext -R -J -F+f16,Helvetica,black -V -O -K <<END>> $ps
#118.0    42.0    China
#138.3    36.2    Japan
#END



#gmt psxy -R -J -W2,red -V -O <<END>> $ps
#127.8      36.5
#137.8      36.5
#137.8      30
#127.8      30
#127.8      36.5
#END

gmt psxy -R -J -Sc0.3 -W2,red -V -O <<END>> $ps
-91.1    0.4
END



##### Output - Raster
gmt psconvert $ps -E600 -Tt
