#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.2                            #
# Copyright(c) 2017, Zhang Yunjun                          #
# Author:  Yunjun Zhang, 2017 Mar 24                       #
# based on selectPairs.py written by Scott Baker at 2010   #
############################################################


import os
import sys
import argparse
import re

import matplotlib.pyplot as plt
import numpy as np

import pysar
import pysar._readfile as readfile
import pysar._datetime as ptime
import pysar._network as pnet


sar_sensor_list=['Ers','Env','Jers','Alos','Alos2','Tsx','Csk','Rsat','Rsat2','S1','Kmps5','G3']

#########################################################################
def project_name2sensor(projectName):
    if    re.search('Ers'    , projectName):  sensor = 'Ers'
    elif  re.search('Env'    , projectName):  sensor = 'Env'
    elif  re.search('Jers'   , projectName):  sensor = 'Jers'
    elif  re.search('Alos'   , projectName):  sensor = 'Alos'
    elif  re.search('Alos2'  , projectName):  sensor = 'Alos2' 
    elif  re.search('Tsx'    , projectName):  sensor = 'Tsx'
    elif  re.search('Tdm'    , projectName):  sensor = 'Tsx'
    elif  re.search('Csk'    , projectName):  sensor = 'Csk'
    elif  re.search('Rsat'   , projectName):  sensor = 'Rsat'
    elif  re.search('Rsat2'  , projectName):  sensor = 'Rsat2'
    elif  re.search('Sen'    , projectName):  sensor = 'S1'
    elif  re.search('Kmps5'  , projectName):  sensor = 'Kmps5'
    elif  re.search('Gaofen3', projectName):  sensor = 'G3'
    else: print 'satellite not found';  sensor = None
    return sensor


def read_template2inps(templateFile, inps=None):
    '''Read network options from template file into Namespace variable inps'''
    template_dict = readfile.read_template(templateFile)
    if not template_dict:
        print 'Empty template: '+templateFile
        return None
    keyList = template_dict.keys()

    if not inps:
        inps = cmdLineParse([''])

    # Read old template option in rsmas_insar
    if 'selectMethod'     in keyList:  inps.method          = template_dict['selectMethod']
    if 'perpBaseMax'      in keyList:  inps.perp_base_max   = template_dict['perpBaseMax']
    if 'lengthDayMax'     in keyList:  inps.temp_base_max   = template_dict['lengthDayMax']
    if 'lengthDayMin'     in keyList:  inps.temp_base_min   = template_dict['lengthDayMin']
    if 'seasonal'         in keyList:  inps.keep_seasonal   = template_dict['seasonal']
    if 'DopOverlapThresh' in keyList:  inps.dop_overlap_min = template_dict['DopOverlapThresh']

    # Read selectPairs.* options
    prefix='selectPairs.'
    if prefix+'selectMethod'     in keyList:  inps.method          = template_dict[prefix+'selectMethod']
    if prefix+'method'           in keyList:  inps.method          = template_dict[prefix+'method']
    if prefix+'perpBaseMax'      in keyList:  inps.perp_base_max   = template_dict[prefix+'perpBaseMax']
    if prefix+'lengthDayMax'     in keyList:  inps.temp_base_max   = template_dict[prefix+'lengthDayMax']
    if prefix+'lengthDayMin'     in keyList:  inps.temp_base_min   = template_dict[prefix+'lengthDayMin']
    if prefix+'seasonal'         in keyList:  inps.keep_seasonal   = template_dict[prefix+'seasonal']
    if prefix+'dopOverlapThresh' in keyList:  inps.dop_overlap_min = template_dict[prefix+'DopOverlapThresh']
    if prefix+'incrementNum'     in keyList:  inps.increment_num   = template_dict[prefix+'incrementNum']
    if prefix+'dayPerpList'      in keyList:  inps.temp_perp_list  = template_dict[prefix+'dayPerpList']
    if prefix+'referenceFile'    in keyList:  inps.reference_file  = template_dict[prefix+'referenceFile']    
    if prefix+'excludeDate' in keyList:
        ex_date_list = [i for i in template_dict[prefix+'excludeDate'].split(',')]
        inps.exclude_date = ptime.yymmdd(ex_date_list)

    # Read pysar.network.* options
    prefix='pysar.network.'
    if prefix+'method'        in keyList:  inps.method          = template_dict[prefix+'method']
    if prefix+'perpBaseMax'   in keyList:  inps.perp_base_max   = template_dict[prefix+'perpBaseMax']
    if prefix+'tempBaseMax'   in keyList:  inps.temp_base_max   = template_dict[prefix+'tempBaseMax']
    if prefix+'tempBaseMin'   in keyList:  inps.temp_base_min   = template_dict[prefix+'tempBaseMin']
    if prefix+'keepSeasonal'  in keyList:  inps.keep_seasonal   = template_dict[prefix+'keepSeasonal']
    if prefix+'dopOverlapMin' in keyList:  inps.dop_overlap_min = template_dict[prefix+'dopOverlapMin']
    if prefix+'referenceFile' in keyList:  inps.reference_file  = template_dict[prefix+'referenceFile']     
    if prefix+'incrementNum'  in keyList:  inps.increment_num   = template_dict[prefix+'incrementNum']
    if prefix+'tempPerpList'  in keyList:  inps.temp_perp_list  = template_dict[prefix+'tempPerpList']
    if prefix+'excludeDate' in keyList:
        ex_date_list = [i for i in template_dict[prefix+'excludeDate'].split(',')]
        inps.exclude_date = ptime.yymmdd(ex_date_list)

    return inps


#########################################################################
REFERENCE='''References:
  Berardino, P., G. Fornaro, R. Lanari, and E. Sansosti (2002), A new algorithm for surface deformation monitoring
    based on small baseline differential SAR interferograms, IEEE TGRS, 40(11), 2375-2383.
  Fattahi, H., and F. Amelung (2013), DEM Error Correction in InSAR Time Series, IEEE TGRS, 51(7), 4249-4259.
  Ferretti, A., C. Prati, and F. Rocca (2001), Permanent scatterers in SAR interferometry, IEEE TGRS, 39(1), 8-20.
  Pepe, A., and R. Lanari (2006), On the extension of the minimum cost flow algorithm for phase unwrapping
    of multitemporal differential SAR interferograms, IEEE TGRS, 44(9), 2374-2383.
  Perissin D., Wang T. (2012), Repeat-pass SAR interferometry with partially coherent targets. IEEE TGRS. 271-280
  Zebker, H. A., and J. Villasenor (1992), Decorrelation in interferometric radar echoes, IEEE TGRS, 30(5), 950-959.
  Zhao, W., (2015), Small deformation detected from InSAR time-series and their applications in geophysics, Doctoral
    dissertation, Univ. of Miami, Section 6.3.
'''

METHOD='''
all          - all possible pairs, pair number = N*(N-1)/2 where N is acquisition num 
               default, (Berardino et al., 2002, TGRS).
delaunay     - Delaunay triangulation (Fattahi and Amelung, 2013, TGRS).
               By default, temporal baseline is normalized using a maxPerpDiff/maxTempDiff
               ratio (Pepe and Lanari, 2006, TGRS);
               to disable this, use 'delaunay-noweight' option.
hierarchical - Select pairs in a hierarchical way using a list of temp/perp thresholds
               (Zhao, 2015, PhD Thesis)
               i.e. 16 days, 1600 m
                    32 days, 800  m
                    48 days, 600  m
                    64 days, 200  m
mst          - Minimum Spanning Tree (Perissin and Wang, 2012, TGRS).
               Find the MST based on the graph of temporal and perpendicular matrix.
sequential   - for each acquisition, select its incrementNum nearest neighbors in the past
               (Fattahi and Amelung, 2013, TGRS). 
               Pair number = N*m - m*(m-1)/2; it's N when m = 1.
               Designed for new satellites like Sentinel-1 and ALOS-2.
star         - Star-like/PS-like network/pairs, single common master interferogram
               (Ferretti et al., 2001, TGRS)
'''

EXAMPLE='''Examples:
  select_network.py KirishimaT246EnvD2.template
  select_network.py KirishimaT246EnvD2.template -b bl_list.txt
  select_network.py KirishimaT246EnvD2.template -b bl_list.txt -o $SC/KirishimaT246EnvD2/PROCESS/ifgram_list.txt
  select_network.py KirishimaT246EnvD2.template -r Pairs.list
  select_network.py KirishimaT246EnvD2.template -r Modified_LoadedData.h5
'''

TEMPLATE='''Template:
pysar.network.method        = all              # all,hierarchical,sequential,mst,delaunay,star
pysar.network.perpBaseMax   = 500              # max perpendicular baseline
pysar.network.tempBaseMax   = 365              # max      temporal baseline
pysar.network.tempBaseMin   = 0                # min       emporal baseline
pysar.network.keepSeasonal  = yes              # keep pairs with seasonal temporal baseline
pysar.network.dopOverlapMin = 15               # min dopploer overlap percentage

pysar.network.referenceFile = ifgram_list.txt  # reference pairs list file
pysar.network.referenceFile = unwrapIfgram.h5  # reference HDF5 file with pairs info
pysar.network.excludeDate   = 080520,100726    # exclude dates for pairs selection
pysar.network.incrementNum  = 2                # for sequential method, pairs num per new acquisition
pysar.network.tempPerpList  = 16,1600;32,800;48,600;64,200  # for hierarchical method, list of max temp/perp baseline
'''

def cmdLineParse():
    parser = argparse.ArgumentParser(description='Select Interferometric Network / Pairs.',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=REFERENCE+'\n'+TEMPLATE+'\n'+EXAMPLE)

    parser.add_argument('template_file', help='template file with options')
    parser.add_argument('-b', dest='baseline_file', \
                        help='baseline list file of all SLCs, e.g.'+pnet.BASELINE_LIST_FILE)
    parser.add_argument('-o','--outfile',\
                        help='Output list file for network, ifgram_list.txt by default.')
    parser.add_argument('--display', dest='disp_fig', action='store_true', help='display network ploting result')

    # Method
    method = parser.add_argument_group('Methods to generate the initial network')
    method.add_argument('--method', default='all',\
                        help='network type with info on temp/perp baseline and doppler centroid frequency.'+METHOD)
    method.add_argument('-r', dest='reference_file',\
                        help='Reference hdf5 / list file with network information. e.g.\n'+\
                             'unwrapIfgram.h5\n'+\
                             'ifgram_list.txt with content as below:'+pnet.IFGRAM_LIST_FILE+\
                             '\nIt could also be generated using plot_network.py --list option, e.g.\n'+\
                             'plot_network.py unwrapIfgram.h5 --list\n\n')
    method.add_argument('--ex','--exclude', dest='exclude_date', nargs='*', \
                        help='date(s) excluded for network selection, e.g. -ex 060713 070831')
    method.add_argument('--increment-num', dest='increment_num', type=int, default=3,\
                        help='number of new pairs for each new acquisition, for sequential method')
    method.add_argument('--temp-perp-list', dest='temp_perp_list', default='16,1600;32,800;48,600;64,200',\
                        help='list of max temp/perp baselines, for hierarchical method, e.g.\n'+\
                             '--temp-perp-list 16,1600;32,800;48,600;64,200')
    method.add_argument('--no-norm', dest='norm', action='store_false',\
                        help='do not normalize temp/perp baseline, for delaunay method')
    method.add_argument('--sensor', help='Name of sensor, choose from the list below:\n'+str(sar_sensor_list))

    # Thresholds
    threshold = parser.add_argument_group('Thresholds to filter the initial network')
    threshold.add_argument('--nothreshold', dest='threshold', action='store_false', \
                           help='do not remove pairs using min/max temp/perp baseline and dop\n'+\
                                'auto applied this option when --reference-file is specified.')
    threshold.add_argument('--dop-overlap-min', dest='dop_overlap_min', type=float, default=15.0,\
                           help='min doppler overlap percentage')
    threshold.add_argument('--bperp-max', dest='perp_base_max', type=float, default=500.0, \
                           help='max perpendicular spatial baseline in meters')
    threshold.add_argument('--btemp-min', dest='temp_base_min', type=float, default=0.0, \
                           help='min temporal baseline in days')
    threshold.add_argument('--btemp-max', dest='temp_base_max', type=float, default=365.0, \
                           help='max temporal baseline in days')
    threshold.add_argument('--noseasonal', dest='keep_seasonal', action='store_false',\
                           help='do not keep seasonal pairs, i.e. pairs in same/adjcent month within 3 years.')

    inps = parser.parse_args()
    try:    inps.reference_file = glob.glob(inps.reference_file)[0]
    except: inps.reference_file = None

    return inps


#########################################################################
def main(argv):
    
    # Read inputs
    inps = cmdLineParse()
    inps = read_template2inps(inps.template_file, inps)
    inps.temp_perp_list = [[float(j) for j in i.split(',')] for i in inps.temp_perp_list.split(';')]
    project_name = os.path.splitext(os.path.basename(inps.template_file))[0]
    print 'project name: '+project_name
    if not inps.sensor:
        inps.sensor = project_name2sensor(project_name)
 
    # Pair selection from reference
    if inps.reference_file:
        print 'Use pairs info from reference file: '+inps.reference_file
        date12_list = pnet.get_date12_list(inps.reference_file)

    # Pair selection from temp/perp/dop baseline info
    else:
        # Auto path setting for Miami user
        if not inps.baseline_file and pysar.miami_path and 'SCRATCHDIR' in os.environ:
            if pysar.miami_path and 'SCRATCHDIR' in os.environ:
                try:    inps.baseline_file = glob.glob(os.getenv('SCRATCHDIR')+'/'+projectName+'/SLC/bl_list.txt')
                except: inps.baseline_file = None
        if not inps.baseline_file:
            raise Exception('ERROR: No baseline file found!')
        
        # Read baseline list file: bl_list.txt
        if inps.exclude_date:
            print 'exclude dates: '+str(inps.exclude_date)
        date8_list, pbase_list, dop_list = pnet.read_baseline_file(inps.baseline_file, inps.exclude_date)[0:3]
        date6_list = ptime.yymmdd(date8_list)
        tbase_list = ptime.date_list2tbase(date8_list)[0]

        # Initial network using input methods
        inps.method = inps.method.lower().replace('-','_')
        print 'select method: '+inps.method

        if   inps.method == 'all':           date12_list = pnet.select_pairs_all(date6_list)
        elif inps.method == 'delaunay':      date12_list = pnet.select_pairs_delaunay(date6_list, pbase_list, inps.norm)
        elif inps.method in ['star','ps']:   date12_list = pnet.select_pairs_star(date6_list)
        elif inps.method.startswith('seq'):  date12_list = pnet.select_pairs_sequential(date6_list, inps.increment_num)
        elif inps.method.startswith('hierar'):
            date12_list = pnet.select_pairs_hierarchical(date6_list, pbase_list, inps.temp_perp_list)
        elif inps.method in ['mst', 'minimum_spanning_tree', 'min_spanning_tree']:
            date12_list = pnet.select_pairs_mst(date6_list, pbase_list)
        else:
            raise Exception('Unrecoganized select method: '+inps.method)
        print 'initial number of interferograms: '+str(len(date12_list))

        # Filter pairs (optional) using temp/perp/doppler baseline threshold
        if inps.threshold:
            # Temporal baseline
            date12_list = pnet.threshold_temporal_baseline(date12_list, inps.temp_base_max,\
                                                           inps.keep_seasonal, inps.temp_base_min)
            print 'number of interferograms after filtering of <%d, %d> days in temporal baseline: %d'\
                  % (inps.temp_base_min, inps.temp_base_max, len(date12_list))
            if inps.keep_seasonal:
                print '\tkeep seasonal pairs, i.e. pairs with temporal baseline == N*years +/- one month'

            # Perpendicular spatial baseline
            date12_list = pnet.threshold_perp_baseline(date12_list, date6_list, pbase_list, inps.perp_base_max)
            print 'number of interferograms after filtering of max %d meters in perpendicular baseline: %d'\
                  % (inps.perp_base_max, len(date12_list))
            
            # Doppler Overlap Percentage
            if inps.sensor:
                bandwidth_az = pnet.azimuth_bandwidth(inps.sensor)
                date12_list = pnet.threshold_doppler_overlap(date12_list, date6_list, dop_list,\
                                                             bandwidth_az, inps.dop_overlap_min/100.0)
                print 'number of interferograms after filtering of min '+str(inps.dop_overlap_min)+'%'+\
                      ' overlap in azimuth Doppler frequency: '+str(len(date12_list))

    # Write ifgram_list.txt
    if not date12_list:
        print 'WARNING: No interferogram selected!'
        return None
    
    # date12_list to date_list
    m_dates = [date12.split('-')[0] for date12 in date12_list]
    s_dates = [date12.split('-')[1] for date12 in date12_list]
    print 'number of acquisitions   selected: '+str(len(list(set(m_dates + s_dates))))
    print 'number of interferograms selected: '+str(len(date12_list))
    
    # Write txt file
    if not inps.outfile:
        try:    inps.work_dir = os.path.dirname(os.path.abspath(inps.reference_file))
        except: inps.work_dir = os.path.dirname(os.path.abspath(inps.baseline_file))
        inps.outfile = inps.work_dir+'/ifgram_list.txt'
    print 'writing >>> '+inps.outfile
    np.savetxt(inps.outfile, date12_list, fmt='%s')

    # Plot network info
    if not inps.disp_fig:
        plt.switch_backend('Agg')

    out_fig_name = 'BperpHist.pdf'
    print 'plotting baseline history in temp/perp baseline domain to file: '+out_fig_name
    fig2, ax2 = plt.subplots()
    ax2 = pnet.plot_perp_baseline_hist(ax2, date8_list, pbase_list)
    plt.savefig(os.path.dirname(inps.outfile)+'/'+out_fig_name, bbox_inches='tight')
        
    out_fig_name = 'Network.pdf'
    print 'plotting network in temp/perp baseline domain to file: '+out_fig_name
    fig1, ax1 = plt.subplots()
    ax1 = pnet.plot_network(ax1, date12_list, date8_list, pbase_list)
    plt.savefig(os.path.dirname(inps.outfile)+'/'+out_fig_name, bbox_inches='tight')
    
    if inps.disp_fig:
        plt.show()

    return inps.outfile

###########################################################################
if __name__ == '__main__':
    main(sys.argv[1:])

 
