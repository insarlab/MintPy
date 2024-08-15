"""Class wrapped around matplotlib.colors for colormaps."""
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2019                               #
############################################################
# Recommend import:
#     from mintpy.objects.colors import ColormapExt
#     from mintpy.utils import plot as pp
#     cmap = pp.ColormapExt('cmy').colormap


import colorsys
import glob
import os
import re

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, to_rgb

import mintpy

MINTPY_CPT_DIR = os.path.join(os.path.dirname(mintpy.__file__), 'data', 'colormaps')
GMT_CPT_DIR = '/opt/local/share/gmt/cpt'  #location of GMT colormap files, default for macOS with GMT installed via MacPorts

# To manually create custom diverging colormaps, check the link below:
# diverging_map.py in https://github.com/ethankruse/kepler_orrery
# https://www.kennethmoreland.com/color-maps/


def isnumber(n):
    try:
        float(n)
    except ValueError:
        return False
    return True


################################## ColormapExt class begin #####################################
class ColormapExt(ScalarMappable):
    """Extended colormap class inherited from matplotlib.cm.ScalarMappable class

    Colormaps priority:
        1. (user_input)
        2. MINTPY_CPT_DIR (cpt-city + sci colormap)
        3. Matplotlib
        4. (GMT_CPT_DIR)

    Example:
        from mintpy.objects.colors import ColormapExt
        cmap = ColormapExt('cmy').colormap          #for cyclic phase
        cmap = ColormapExt('jet').colormap          #from matplotlib
        cmap = ColormapExt('RdBu').colormap         #from matplotlib
        cmap = ColormapExt('haxby').colormap        #from GMT
        cmap = ColormapExt('vik').colormap          #from Scientific Color-Maps
        cmap = ColormapExt('temperature').colormap  #from cpt-city

        ## derivative names
        # reverse  colormap by adding suffix "_r"
        cmap = ColormapExt('RdBu_r').colormap

        # truncate colormap by adding suffix "_truncate"
        cmap = ColormapExt('RdBu_truncate', vlist=[0.2, 0.4, 1.0]).colormap

        # repeat   colormap by adding suffix "_{}".format(int)
        cmap = ColormapExt('jet_5').colormap

        ## combined derivative names has to follow the order: _r, _truncate and _{int}:
        i.e. 'RdBu_r_truncate'
             'RdBu_truncate_3'
             'jet_r_5'
    """

    def __init__(self, cmap_name, cmap_lut=256, vlist=[0.0, 0.7, 1.0], cpt_dir=None):
        """ Initiate an ColormapExt object
        Parameters: cmap_name - str, colormap name. Default: viridis
                    cmap_lut  - int, number of increment in the color lookup table
                    vlist     - list of 3 float numbers, for truncated colormap only
                    cpt_dir   - list of str, extra directories of cpt files to be recognized
        """
        # default setup
        self.reverse_colormap = False
        self.truncate_colormap = False
        self.num_repeat = 1

        # initiate member variables
        self.cmap_name = cmap_name
        self.cmap_lut = cmap_lut
        self.vlist = vlist

        # initiate cpt_dirs for custom colormaps with:
        # 1. cpt_dir from custom input argument during initiation
        # 2. MINTPY_CPT_DIR, including some from cpt-city and scientific color maps
        # 3. GMT_CPT_DIR, if GMT is installed
        self.cpt_dirs = [i for i in [cpt_dir, MINTPY_CPT_DIR, GMT_CPT_DIR]
                         if i and os.path.isdir(i)]

        # initiate member functions
        self.get_colormap_name_list()

        # generate colormap object
        self.check_input_colormap_name()
        self.get_colormap()
        return


    def get_colormap_name_list(self):
        """list of colormap supported in string for name of colormap, from two sources:
            1) local GMT cpt files
            2) matlotlib
        """
        self.cpt_cmap_name_list = self.get_cpt_colormap(cmap_name=None)
        self.plt_cmap_name_list = sorted(m for m in plt.colormaps() if not m.endswith('_r'))
        self.cmap_name_list = self.cpt_cmap_name_list + self.plt_cmap_name_list + ['dismph','cmy']
        return self.cmap_name_list


    def check_input_colormap_name(self):
        """Check 1) input colormap name is supported or not
                 2) derivative settings, from the suffix
        """
        if self.cmap_name in self.cmap_name_list:
            return
        else:
            # check repeat number if file ends with "_{int}"
            re_num = re.search(r'_\d+$', self.cmap_name)
            if re_num is not None:
                suffix = re_num[0]
                self.num_repeat = int(suffix.split('_')[-1])
                self.cmap_name = self.cmap_name.split(suffix)[0]

            # check truncate setting
            if self.cmap_name.endswith('_truncate'):
                self.truncate_colormap = True
                self.cmap_name = self.cmap_name.split('_truncate')[0]

            # check reverse setting
            if self.cmap_name.endswith('_r'):
                self.reverse_colormap = True
                self.cmap_name = self.cmap_name.split('_r')[0]

            # check if input colormap name is supported
            if self.cmap_name not in self.cmap_name_list:
                msg = f'un-recognized input colormap name: {self.cmap_name}\n'
                msg += f'supported colormap from cpt files:\n{self.cpt_cmap_name_list}\n'
                msg += f'supported colormap from matplotlib:\n{self.plt_cmap_name_list}\n'
                raise ValueError(msg)
        return


    def get_colormap(self):
        self.colormap = self.get_single_colormap(cmap_name=self.cmap_name, cmap_lut=self.cmap_lut)

        # reverse setting
        if self.reverse_colormap:
            self.colormap = self.colormap.reversed()

        # truncate setting
        if self.truncate_colormap:
            self.cmap_lut = 2560 #higher color resolution to distinguish colors near the jump value

            n1_ratio = (self.vlist[1] - self.vlist[0]) / (self.vlist[2] - self.vlist[0])
            n1 = np.rint(self.cmap_lut * n1_ratio).astype('int')
            n2 = self.cmap_lut - n1
            colors1 = self.colormap(np.linspace(0.0, 0.3, n1))
            colors2 = self.colormap(np.linspace(0.6, 1.0, n2))
            self.colormap = LinearSegmentedColormap.from_list(
                name=self.cmap_name+'_truncate',
                colors=np.vstack((colors1, colors2)),
                N=self.cmap_lut,
            )

        # repeat setting
        if self.num_repeat > 1:
            colors = np.tile(self.colormap(np.linspace(0., 1., self.cmap_lut)), (self.num_repeat,1))
            self.colormap = LinearSegmentedColormap.from_list(
                name=self.cmap_name+f'_{self.num_repeat}',
                colors=colors,
                N=self.cmap_lut*self.num_repeat,
            )

        return self.colormap


    def get_single_colormap(self, cmap_name, cmap_lut=256):
        if cmap_name == 'dismph':
            # color list from bakerunavco/pygmtsar:
            # reference: showintf.py in https://github.com/bakerunavco/pygmtsar
            clist = ['#f579cd', '#f67fc6', '#f686bf', '#f68cb9', '#f692b3', '#f698ad',
                     '#f69ea7', '#f6a5a1', '#f6ab9a', '#f6b194', '#f6b78e', '#f6bd88',
                     '#f6c482', '#f6ca7b', '#f6d075', '#f6d66f', '#f6dc69', '#f6e363',
                     '#efe765', '#e5eb6b', '#dbf071', '#d0f477', '#c8f67d', '#c2f684',
                     '#bbf68a', '#b5f690', '#aff696', '#a9f69c', '#a3f6a3', '#9cf6a9',
                     '#96f6af', '#90f6b5', '#8af6bb', '#84f6c2', '#7df6c8', '#77f6ce',
                     '#71f6d4', '#6bf6da', '#65f6e0', '#5ef6e7', '#58f0ed', '#52e8f3',
                     '#4cdbf9', '#7bccf6', '#82c4f6', '#88bdf6', '#8eb7f6', '#94b1f6',
                     '#9aabf6', '#a1a5f6', '#a79ef6', '#ad98f6', '#b392f6', '#b98cf6',
                     '#bf86f6', '#c67ff6', '#cc79f6', '#d273f6', '#d86df6', '#de67f6',
                     '#e561f6', '#e967ec', '#ed6de2', '#f173d7']
            colormap = LinearSegmentedColormap.from_list('dismph', clist, N=cmap_lut)
            #colormap = self.cmap_map(lambda x: x/2 + 0.5, colormap)  # brighten colormap
            #colormap = self.cmap_map(lambda x: x*0.75, colormap)     # darken colormap
            colormap.set_bad('w', 0.0)

        elif cmap_name == 'cmy':
            # Default cyclic colormap from isce/mdx, provided by Piyush Agram, Jan 2020
            # generate the color list
            rgbs = np.zeros((256,3), dtype=np.uint8)

            for kk in range(85):
                rgbs[kk,0] = kk*3
                rgbs[kk,1] = 255-kk*3
                rgbs[kk,2] = 255

            rgbs[85:170,0] = rgbs[0:85,2]
            rgbs[85:170,1] = rgbs[0:85,0]
            rgbs[85:170,2] = rgbs[0:85,1]

            rgbs[170:255,0] = rgbs[0:85,1]
            rgbs[170:255,1] = rgbs[0:85,2]
            rgbs[170:255,2] = rgbs[0:85,0]

            rgbs[255,0] = 0
            rgbs[255,1] = 255
            rgbs[255,2] = 255

            rgbs = np.roll(rgbs, int(256/2-214), axis=0)  #shift green to the center
            rgbs = np.flipud(rgbs)   #flip up-down so that orange is in the later half (positive)

            # color list --> colormap object
            colormap = LinearSegmentedColormap.from_list('cmy', rgbs/255., N=cmap_lut)

        elif cmap_name in self.cpt_cmap_name_list:
            colormap = self.get_cpt_colormap(cmap_name, cmap_lut=cmap_lut)

        else:
            colormap = plt.get_cmap(cmap_name, lut=cmap_lut)

        return colormap


    def get_cpt_colormap(self, cmap_name=None, cmap_lut=256):
        """Load GMT .cpt colormap file.
        Modified from Scipy Cookbook originally written by James Boyle.
        Link: http://scipy-cookbook.readthedocs.io/items/Matplotlib_Loading_a_colormap_dynamically.html

        Parameters: cmap_name : string, colormap name, e.g. temperature
        Returns:    colormap : matplotlib.colors.LinearSegmentedColormap object
        Example:    colormap = get_cpt_colormap('temperature')
                    colormap = get_cpt_colormap('temperature_r')
                    cpt_cm_list = get_cpt_colormap(None)
        """
        # Return list of existing colormaps, if cmap_name is None.
        if not cmap_name:
            cm_list = []
            for cpt_dir in self.cpt_dirs:
                cpt_files = sorted(glob.glob(os.path.join(cpt_dir, '*.cpt')))
                cm_list += [os.path.splitext(os.path.basename(i))[0] for i in cpt_files]
            return cm_list

        # support _r for reversed colormap
        reverse_colormap = False
        if cmap_name.endswith('_r'):
            reverse_colormap = True
            cmap_name = cmap_name[0:-2]

        # search for cpt_file
        cpt_file = None
        for cpt_dir in self.cpt_dirs:
            cpt_file = os.path.join(cpt_dir, f"{cmap_name}.cpt")
            if os.path.isfile(cpt_file):
                break

        colormap = self.read_cpt_file(cpt_file, cmap_lut=cmap_lut)

        if reverse_colormap:
            colormap = colormap.reversed()
        return colormap


    @staticmethod
    def read_cpt_file(cpt_file, cmap_lut=256):
        """Read *.cpt file into colorDict
        Modified from Scipy Cookbook originally written by James Boyle.
        Link: http://scipy-cookbook.readthedocs.io/items/Matplotlib_Loading_a_colormap_dynamically.html
        """
        if not os.path.isfile(cpt_file):
            raise FileNotFoundError(f"file {cpt_file} not found")

        # read file into list of strings
        with open(cpt_file) as f:
            lines = f.readlines()

        # list of string --> x/r/g/b
        x, r, g, b = [], [], [], []
        colorModel = "RGB"
        for line in lines:
            ls = re.split(' |\t|\n|/', line)

            # skip empty lines
            if not ls:
                continue

            # remove empty element
            ls = [i for i in ls if i]

            # parse header info
            if line[0] == "#":
                if ls[-1] == "HSV":
                    colorModel = "HSV"
                    continue
                else:
                    continue

            # skip BFN info
            if ls[0] in ["B", "F", "N"]:
                continue

            # convert color name (in GMT cpt file sometimes) to rgb values
            if not isnumber(ls[1]):
                ls0 = list(ls) + [0,0]
                ls0[1:4] = [i*255. for i in to_rgb(ls[1])]
                ls0[4:] = ls[2:]
                ls = list(ls0)

            if not isnumber(ls[5]):
                ls0 = list(ls) + [0,0]
                ls0[5:8] = [i*255. for i in to_rgb(ls[5])]
                ls = list(ls0)

            # convert str to float
            ls = [float(i) for i in ls]

            # parse color vectors
            x.append(ls[0])
            r.append(ls[1])
            g.append(ls[2])
            b.append(ls[3])

            # save last row
            xtemp = ls[4]
            rtemp = ls[5]
            gtemp = ls[6]
            btemp = ls[7]
        x.append(xtemp)
        r.append(rtemp)
        g.append(gtemp)
        b.append(btemp)

        x = np.array(x, np.float32)
        r = np.array(r, np.float32)
        g = np.array(g, np.float32)
        b = np.array(b, np.float32)

        if colorModel == "HSV":
            # convert HSV to RGB
            for i in range(r.shape[0]):
                r[i], g[i], b[i] = colorsys.hsv_to_rgb(r[i]/360., g[i], b[i])
        elif colorModel == "RGB":
            r /= 255.
            g /= 255.
            b /= 255.

        # x/r/g/b --> colorDict
        red, blue, green = [], [], []
        xNorm = (x - x[0]) / (x[-1] - x[0])
        for i in range(len(x)):
            red.append((xNorm[i], r[i], r[i]))
            green.append((xNorm[i], g[i], g[i]))
            blue.append((xNorm[i], b[i], b[i]))

        # return colormap
        cmap_name = os.path.splitext(os.path.basename(cpt_file))[0]
        colorDict = {"red":tuple(red), "green":tuple(green), "blue":tuple(blue)}
        colormap = LinearSegmentedColormap(cmap_name, colorDict, N=cmap_lut)
        return colormap


    @staticmethod
    def cmap_map(function, cmap):
        """ Applies function (which should operate on vectors of shape 3: [r, g, b]), on colormap cmap.
        This routine will break any discontinuous points in a colormap.
        Link: http://scipy-cookbook.readthedocs.io/items/Matplotlib_ColormapTransformations.html
        """
        cdict = cmap._segmentdata
        step_dict = {}

        # First get the list of points where the segments start or end
        for key in ('red', 'green', 'blue'):
            step_dict[key] = list(map(lambda x: x[0], cdict[key]))
        step_list = sum(step_dict.values(), [])
        step_list = np.array(list(set(step_list)))

        # Then compute the LUT, and apply the function to the LUT
        reduced_cmap = lambda step : np.array(cmap(step)[0:3])
        old_LUT = np.array(list(map(reduced_cmap, step_list)))
        new_LUT = np.array(list(map(function, old_LUT)))

        # Now try to make a minimal segment definition of the new LUT
        cdict = {}
        for i, key in enumerate(['red','green','blue']):
            this_cdict = {}
            for j, step in enumerate(step_list):
                if step in step_dict[key]:
                    this_cdict[step] = new_LUT[j, i]
                elif new_LUT[j,i] != old_LUT[j, i]:
                    this_cdict[step] = new_LUT[j, i]
            colorvector = list(map(lambda x: x + (x[1], ), this_cdict.items()))
            colorvector.sort()
            cdict[key] = colorvector
        return LinearSegmentedColormap('colormap',cdict,1024)

################################## ColormapExt class end ########################################
