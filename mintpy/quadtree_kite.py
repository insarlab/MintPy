#!/usr/bin/env python3
#################################################################
# Program using kite software to do quadtree downsampling       #
# Author: Lv Xiaoran                                            #
# Created: June 2020                                            #
#################################################################

import os
import argparse
import numpy as np
import logging
from kite import Scene
import matplotlib.pyplot as plt
from matplotlib import cm, colors

import mintpy
from mintpy.utils import readfile, writefile
######################################################################################
EXAMPLE = """example:
  quadtree_kite.py geo_20171117_20200205_kite.unw -e 0.007 --nan 0.9 -max 0.5922566019916117 -min 0.002 --plot -o geo_20171117_20200205_kite_quadtree
  quadtree_kite.py geo_20171117_20200205_kite.unw -e 0.007 --nan 0.9 -max 0.5922566019916117 -min 0.002 -o geo_20171117_20200205_kite_quadtree
  
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Prepare data for Kite software',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', nargs=1, type=str, help='geocoded unw or h5 files to be converted\n')

    parser.add_argument('-e', '--epsilon', dest='epsilon', type=float, nargs=1,
                        help='variance threshold')
    parser.add_argument('--nan', dest='nan', nargs=1, type=float,
                        help='Percentage of NaN values allowed per tile/leave')
    parser.add_argument('-max','--tile_size_max', dest='tile_max', nargs=1, type=float,
                        help='Maximum leave edge length in [m] or [deg]')
    parser.add_argument('-min','--tile_size_min', dest='tile_min', nargs=1, type=float,
                        help='Minimum leave edge length in [m] or [deg]')
    parser.add_argument('--plot', action='store_true', default=False, 
                        help='whether plot quadtree leaves') 
    parser.add_argument('-o','--outfile', dest = 'outfile', nargs=1, type=str,
                        help='output file name')

    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)  
    
    return inps

def quadtree_kite(inps):
    """using kite do quadtree donwsampling"""
    logging.basicConfig(level=logging.DEBUG)

    file = inps.file[0]
    sc = Scene.import_data(file)
    
    # For convenience we set an abbreviation to the quadtree
    qt = sc.quadtree
    
    # Parametrisation of the quadtree
    qt.epsilon = inps.epsilon[0]        # Variance threshold
    qt.nan_allowed = inps.nan[0]        # Percentage of NaN values allowed per tile/leave
    
    # Be careful here, if you scene is referenced in degree use decimal values!
    qt.tile_size_max = inps.tile_max[0]    # Maximum leave edge length in [m] or [deg]
    qt.tile_size_min = inps.tile_min[0]    # Minimum leave edge length in [m] or [deg]
    
    print('the reduction rsm is %f' % qt.reduction_rms)   # In units of [m] or [deg]
    
    # We save the scene in kite's format
    outname = inps.outfile[0]
    sc.save(outname)
    
    # export to csv format
    export_csv(sc, outname)
    # Or export the quadtree to CSV file
    #qt.export_csv(outname + '.csv')
    #print(sc.phi)
    #qt.export_geojson(outname + '.json')

    return

def export_csv(sc,outfile):
    """export quadtree leaves parameters"""
    llLat = sc.frame.llLat
    llLon = sc.frame.llLon
    dE = sc.frame.dE
    dN = sc.frame.dN

    qt = sc.quadtree
    
    filename = outfile + '.csv'
    print('Exporting Quadtree as to %s' % filename)
    with open(filename, mode='w') as f:
        f.write(
            '# node_id, focal_point_E, focal_point_N, lat, lon,'
            ' mean_displacement, median_displacement, absolute_weight\n')
        for lf in qt.leaves:
            lat, lon = rowcolm2latlon(lf.focal_point[1] * 1000, lf.focal_point[0] * 1000, llLat, llLon, dE, dN)
            f.write(
                '{lf.id}, {lf.focal_point[0]}, {lf.focal_point[1]}, '
                '{lat}, {lon}, {lf.mean}, {lf.median}, {lf.weight}\n'.format(lf=lf,lat=lat,lon=lon))
    f.close()

    return

def rowcolm2latlon(row, colm, llLat, llLon, dE, dN):
    """covert row/colm to lat/lon"""
    
    lat = row * dN + llLat 
    lon = colm * dE + llLon

    return lat, lon 

def quadtree_plot(inps):
    """plot the leaves"""
    outfile = inps.outfile[0]
    sc = Scene.load(outfile + '.yml')
    qt = sc.quadtree
    
    fig = plt.figure()
    ax = fig.gca()
    
    limit = np.abs(qt.leaf_medians).max()
    color_map = cm.ScalarMappable(
        norm=colors.Normalize(vmin=-limit, vmax=limit),
        cmap=cm.get_cmap('jet'))
   
    for rect, leaf in zip(qt.getMPLRectangles(), qt.leaves):
        color = color_map.to_rgba(leaf.median)
        rect.set_facecolor(color)
        ax.add_artist(rect)
    
    ax.set_xlim(qt.leaf_eastings.min(), qt.leaf_eastings.max())
    ax.set_ylim(qt.leaf_northings.min(), qt.leaf_northings.max())
    
    fig.savefig(outfile + '.png', dip=300, bbox_inches='tight')
######################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)   

    # generate dem.jpeg
    quadtree_kite(inps)

    if inps.plot:
        quadtree_plot(inps)

######################################################################################
if __name__ == '__main__':
    main()
