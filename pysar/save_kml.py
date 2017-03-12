#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
#
# Yunjun, Jul 2015: add 'timeseries'/'wrapped' option
# Yunjun, Oct 2015: merge all HDF5 option into one
#                   add support for ROI_PAC product
# Yunjun, Nov 2015: support different fig unit
#                   update colorbar


import os
import sys
import getopt

try:
    from pykml.factory import KML_ElementMaker as KML
except:
    sys.exit('pykml should be installed!')

import h5py
from lxml import etree
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import pysar._readfile as readfile


############################################################
def rewrap(unw):
    rewrapped = unw - np.round(unw/(2*np.pi)) * 2*np.pi
    return rewrapped

############################################################
def usage():
    print '''
***************************************************************
  generating  kml kmz files. (needs geocoded files )

  Usage:
      save_kml.py file
      save_kml.py -f file -m min -M max -d epoch_date -c color_map -i no(yes)
  
      -f : geocoded PySAR / ROI_PAC product
      -m : minimum value
      -M : maximum value
      -d : date of interferogram or time-series epoch to be converted to kml
           for interferogram, like 971220-990703; 
           for timeseries, like 060924 or 20060924
      --displacement : convert to displacement (for radian unit file type)

      -c : colormap, jet as default
      -i : inverse the colormap
      -w : re-wrapping the interferogram [default : no]
      -r : dpi (dots per inch) [default = 300]
      --fig-size     : figure size in inch, default is [8.0,12.0]
      --noreference  : do not show reference point
      --ref-size     : reference point marker size in points

      --cbar-bin-num : colorbar bin number [9 by default]
      --cbar-label   : colorbar label
                       LOS displacement velocity (default)

  Example:
      save_kml.py -f geo_velocity.h5 -m -0.05 -M 0.05
      save_kml.py -f geo_velocity.h5 -m -0.05 -M 0.05 --ref-size 2
      save_kml.py -f geo_velocity.h5 -m -0.05 -M 0.05 -i yes -c jet -r 250
      save_kml.py -f unwrapIfgram.h5 -d 971220-990703 
      save_kml.py -f timeseries.h5 -d 20060924         
      save_kml.py -f geo_filt_100820-101120-sim_HDR_4rlks_c10.unw
      save_kml.py gsi10m.dem

***************************************************************
    '''


############################################################
def main(argv):
    plt.switch_backend('Agg')

    cbar_bin_num  = 9
    cbar_label    = 'Mean LOS velocity'
    color_map     = 'jet'
    data_alpha    = 0.7
    disp_opposite = 'no'
    disp_colorbar = 'yes'
    rewrapping    = 'no'
    fig_dpi       = 500
    #fig_size      = [6.0,9.0]
    fig_unit      = 'mm/yr'
    disp_ref      = 'yes'
    ref_size      = 5
    dispDisplacement = 'no'

    if len(sys.argv)>2:
        try:   opts, args = getopt.getopt(argv,"f:m:M:d:c:w:i:r:",['noreference','fig-size',\
                                               'ref-size=','cbar-label=','displacement','cbar-bin-num='])
        except getopt.GetoptError:  usage() ; sys.exit(1)

        for opt,arg in opts:
            if   opt == '-f':        File = arg
            elif opt == '-m':        Vmin = float(arg)
            elif opt == '-M':        Vmax = float(arg)
            elif opt == '-d':        epoch_date    = arg
            elif opt == '-c':        color_map     = arg
            elif opt == '-i':        disp_opposite = arg
            elif opt == '-w':        rewrapping    = arg
            elif opt == '-r':        fig_dpi = int(arg)
            elif opt == '--cbar-bin-num' :   cbar_bin_num     = int(arg)
            elif opt == '--cbar-label'   :   cbar_label       = arg
            elif opt == '--displacement' :   dispDisplacement = 'yes'
            elif opt == '--fig-size'     :   fig_size = [float(i) for i in arg.split(',')][0:2]
            elif opt == '--ref-size'     :   ref_size = int(arg)
            elif opt == '--noreference'  :   disp_ref = 'no'

    elif len(sys.argv)==2:
        if argv[0]=='-h':               usage(); sys.exit(1)
        elif os.path.isfile(argv[0]):   File = argv[0]
        else:                           usage(); sys.exit(1)
    else:                             usage(); sys.exit(1)

    #######################################################
    ###################  Prepare Data  ####################
    ## prepare: data, North, East, South, West

    ext = os.path.splitext(File)[1].lower()
    atr = readfile.read_attribute(File)
    k = atr['FILE_TYPE']
    #print '\n*************** Output to KMZ file ****************'
    print 'Input file is '+k

    if ext == '.h5':
        try:      h5file=h5py.File(File,'r')
        except:   usage() ; sys.exit(1)
        outName=File.split('.')[0]

        if k in ('interferograms','wrapped','coherence'):
            ifgramList=h5file[k].keys()
            for i in range(len(ifgramList)):
                if epoch_date in ifgramList[i]:
                    epoch_number = i
            print ifgramList[epoch_number]
            outName = ifgramList[epoch_number]
            #outName=epoch_date

            dset = h5file[k][ifgramList[epoch_number]].get(ifgramList[epoch_number])
            data = dset[0:dset.shape[0],0:dset.shape[1]]

            if k == 'wrapped':
                print 'No wrapping for wrapped interferograms. Set rewrapping=no'
                rewrapping = 'no'
                Vmin = -np.pi    
                Vmax = np.pi

        elif 'timeseries' in k:
            epochList=h5file['timeseries'].keys()
            for i in range(len(epochList)):
                if epoch_date in epochList[i]:
                    epoch_number = i

            #### Out name
            try:    ref_date = atr['ref_date']
            except: ref_date = ut.yyyymmdd(atr['DATE'])[0]
            #ref_date=h5file['timeseries'].attrs['ref_date']
            if len(epoch_date)==8:  outName=ref_date[2:]+'-'+epoch_date[2:]
            else:                   outName=ref_date[2:]+'-'+epoch_date

            dset = h5file['timeseries'].get(epochList[epoch_number])
            data = dset[0:dset.shape[0],0:dset.shape[1]]

        ### one dataset format: velocity, mask, temporal_coherence, rmse, std, etc.
        else:
            dset = h5file[k].get(k)
            data=dset[0:dset.shape[0],0:dset.shape[1]]
            if disp_opposite in('yes','Yes','Y','y','YES'):
                data=-1*data

            try:
                xref=h5file[k].attrs['ref_x']
                yref=h5file[k].attrs['ref_y']
            except: pass

    elif ext in ['.unw','.cor','.hgt','.trans','.dem']:
        if   ext in ['.unw','.cor','.hgt','.trans']:
            a,data,atr = readfile.read_float32(File)
            outName = File
            if ext in ['.unw']:
                if dispDisplacement == 'yes':
                    print 'show displacement'
                    phase2range = -float(atr['WAVELENGTH']) / (4*np.pi)
                    data *= phase2range
                    atr['UNIT'] = 'm'
                    rewrapping == 'no'
                    fig_unit = 'mm'
                if rewrapping == 'yes':
                    data = rewrap(data,atr)
                    fig_unit = 'radian'
        elif ext == '.dem':
            data,atr = readfile.read_real_int16(File)
            outName = File
        if   ext in ['.hgt','.dem']:     fig_unit = 'm'
        elif ext in ['.cor','.trans']:   fig_unit = ' '
    else: sys.exit('Do not support '+ext+' file!')


    ########################################################
    if rewrapping=='yes':
        data=rewrap(data)
        Vmin = -np.pi    #[-pi,pi] for wrapped interferograms
        Vmax =  np.pi
    else:
        try:     Vmin
        except:  Vmin = np.nanmin(data)
        try:     Vmax
        except:  Vmax = np.nanmax(data)

    try:
        lon_step = float(atr['X_STEP'])
        lat_step = float(atr['Y_STEP'])
        lon_unit = atr['Y_UNIT']
        lat_unit = atr['X_UNIT']
        West     = float(atr['X_FIRST'])
        North    = float(atr['Y_FIRST'])
        South    = North+lat_step*(data.shape[0]-1)
        East     = West +lon_step*(data.shape[1]-1)
        geocoord = 'yes'
        print 'Geocoded'
    except:
        print '%%%%%%%%%%'
        print 'Error:\nThe input file is not geocoded\n'
        print '%%%%%%%%%%'
        usage();sys.exit(1)


    #######################################################
    ###################  Output KMZ  ######################

    ############### Make PNG file
    print 'Making png file ...'   
    length = data.shape[0]
    width  = data.shape[1]
    try:fig_size
    except:
        fig_size_0 = 6.0           ## min figure dimension: 6.0
        ratio = float(length)/float(width)
        fig_size = [fig_size_0,fig_size_0*ratio]
    print 'figure size:  %.1f, %.1f'%(fig_size[0],fig_size[1])
    ccmap = plt.get_cmap(color_map)
    fig = plt.figure(figsize=fig_size,frameon=False)
    ax = fig.add_axes([0., 0., 1., 1.])
    ax.set_axis_off()

    aspect = width/(length*1.0)
    try:     ax.imshow(data,aspect='auto',cmap=ccmap,vmax=Vmax,vmin=Vmin)
    except:  ax.imshow(data,aspect='auto',cmap=ccmap)

    if disp_ref == 'yes':
        try:
            xref = int(atr['ref_x'])
            yref = int(atr['ref_y'])
            ax.plot(xref,yref,'ks',ms=ref_size)
            print 'show reference point'
        except: print 'Cannot find reference point info!'

    ax.set_xlim([0,width])
    ax.set_ylim([length,0])

    figName = outName + '.png'
    print 'writing '+figName
    plt.savefig(figName, pad_inches=0.0, transparent=True, dpi=fig_dpi)

    ############### Making colorbar
    pc = plt.figure(figsize=(1,8))
    axc = pc.add_subplot(111)
    if   fig_unit in ['mm','mm/yr']: v_scale = 1000
    elif fig_unit in ['cm','cm/yr']: v_scale = 100
    elif fig_unit in ['m',  'm/yr']: v_scale = 1
    norm = mpl.colors.Normalize(vmin=Vmin*v_scale, vmax=Vmax*v_scale)
    clb  = mpl.colorbar.ColorbarBase(axc,cmap=ccmap,norm=norm, orientation='vertical')

    #clb.set_label(fig_unit)
    clb.set_label(cbar_label+' ['+fig_unit+']')
    clb.locator = mpl.ticker.MaxNLocator(nbins=cbar_bin_num)
    clb.update_ticks()

    pc.subplots_adjust(left=0.2,bottom=0.3,right=0.4,top=0.7)
    pc.patch.set_facecolor('white')
    pc.patch.set_alpha(0.7)
    pc.savefig('colorbar.png',bbox_inches='tight',facecolor=pc.get_facecolor(),dpi=300)

    ############## Generate KMZ file
    print 'generating kml file ...'
    try:     doc = KML.kml(KML.Folder(KML.name(atr['PROJECT_NAME'])))
    except:  doc = KML.kml(KML.Folder(KML.name('PySAR product')))
    slc = KML.GroundOverlay(KML.name(figName),KML.Icon(KML.href(figName)),\
                            KML.altitudeMode('clampToGround'),\
                            KML.LatLonBox(KML.north(str(North)),KML.south(str(South)),\
                                          KML.east( str(East)), KML.west( str(West))))
    doc.Folder.append(slc)

    #############################
    print 'adding colorscale ...'
    cb_rg = min(North-South, East-West)
    cb_N = (North+South)/2.0 + 0.5*0.5*cb_rg
    cb_W = East  + 0.1*cb_rg
    slc1 = KML.GroundOverlay(KML.name('colorbar'),KML.Icon(KML.href('colorbar.png')),\
                             KML.altitude('2000'),KML.altitudeMode('absolute'),\
                             KML.LatLonBox(KML.north(str(cb_N)),KML.south(str(cb_N-0.5*cb_rg)),\
                                           KML.west( str(cb_W)),KML.east( str(cb_W+0.14*cb_rg))))
    doc.Folder.append(slc1)

    #############################
    kmlstr = etree.tostring(doc, pretty_print=True) 
    kmlname = outName + '.kml'
    print 'writing '+kmlname
    kmlfile = open(kmlname,'w')
    kmlfile.write(kmlstr)
    kmlfile.close()

    kmzName = outName + '.kmz'
    print 'writing '+kmzName
    cmdKMZ = 'zip ' + kmzName +' '+ kmlname +' ' + figName + ' colorbar.png'
    os.system(cmdKMZ)

    cmdClean = 'rm '+kmlname;      print cmdClean;    os.system(cmdClean)
    cmdClean = 'rm '+figName;      print cmdClean;    os.system(cmdClean)
    cmdClean = 'rm colorbar.png';  print cmdClean;    os.system(cmdClean)


#######################################################
if __name__ == '__main__':
    main(sys.argv[1:])


