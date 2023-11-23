#!/bin/sh
# -*- coding: utf-8 -*-
###############################################################################
# $Id$
#
#  Project:  GDAL
#  Purpose:  (Interactive) script to identify and fix typos
#  Author:   Even Rouault <even.rouault at spatialys.com>
#
###############################################################################
#  Copyright (c) 2016, Even Rouault <even.rouault at spatialys.com>
#
#  Permission is hereby granted, free of charge, to any person obtaining a
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included
#  in all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#  DEALINGS IN THE SOFTWARE.
###############################################################################

set -eu

SCRIPT_DIR=$(dirname "$0")
case $SCRIPT_DIR in
    "/"*)
        ;;
    ".")
        SCRIPT_DIR=$(pwd)
        ;;
    *)
        SCRIPT_DIR=$(pwd)"/"$(dirname "$0")
        ;;
esac
GDAL_ROOT=$SCRIPT_DIR/..
cd "$GDAL_ROOT"

if [ ! -f fix_typos/gdal_dict.txt ]; then
    # Get our fork of codespell that adds --words-white-list and full filename support for -S option
    mkdir -p fix_typos
    (
        cd fix_typos
        if [ ! -d codespell ]; then
          git clone https://github.com/rouault/codespell
        fi
        (cd codespell && git checkout gdal_improvements)
        # Aggregate base dictionary + QGIS one + Debian Lintian one
        curl https://github.com/qgis/QGIS/blob/master/scripts/spell_check/spelling.dat | sed "s/:/->/" | sed "s/:%//" | grep -v "colour->" | grep -v "colours->" > qgis.txt
        curl https://salsa.debian.org/lintian/lintian/-/raw/master/data/spelling/corrections | grep "||" | grep -v "#" | sed "s/||/->/" > debian.txt
        cat codespell/data/dictionary.txt qgis.txt debian.txt | awk 'NF' > gdal_dict.txt
        echo "difered->deferred" >> gdal_dict.txt
        echo "differed->deferred" >> gdal_dict.txt
        grep -v 404 < gdal_dict.txt > gdal_dict.txt.tmp
        mv gdal_dict.txt.tmp gdal_dict.txt
    )
fi

EXCLUDED_FILES="*/.svn*,*/.git/*,configure,config.log,config.status,config.guess,config.sub,*/autom4te.cache/*,*.ai,*.svg"
if [ -z ${AUTHORIZED_LIST+x} ]; then export AUTHORIZED_LIST=""; fi
AUTHORIZED_LIST="$AUTHORIZED_LIST,te" # gdalwarp switch
AUTHORIZED_LIST="$AUTHORIZED_LIST,LaTeX,BibTeX"
AUTHORIZED_LIST="$AUTHORIZED_LIST,ALOS,Alos"
AUTHORIZED_LIST="$AUTHORIZED_LIST,lon,Lon,LON"
# New Mintpy ones
AUTHORIZED_LIST="$AUTHORIZED_LIST,alos,ALOS,alosStack"
AUTHORIZED_LIST="$AUTHORIZED_LIST,NED,LOD,Lod,lod"
AUTHORIZED_LIST="$AUTHORIZED_LIST,waterMask,watermask"
AUTHORIZED_LIST="$AUTHORIZED_LIST,smallbaselineApp"
AUTHORIZED_LIST="$AUTHORIZED_LIST,Nealy" # Author in reference

python fix_typos/codespell/codespell.py -w -i 3 -q 2 -S "$EXCLUDED_FILES,./autotest/*,./build*/*,./fix_typos/*" \
    --words-white-list="$AUTHORIZED_LIST" \
    -D ./fix_typos/gdal_dict.txt .
