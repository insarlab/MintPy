#!/usr/bin/env python3
############################################################
# Program is part of PySAR v2.0                            #
# Copyright(c) 2017, Zhang Yunjun                          #
# Author:  Zhang Yunjun, 2017 Mar 24                       #
############################################################
# Based on selectPairs.py written by Scott Baker at 2010
#


import os
import sys
import argparse
import re
import glob
import datetime
import inspect

import matplotlib.pyplot as plt
import numpy as np

import pysar
import pysar.utils.datetime as ptime
import pysar.utils.readfile as readfile
import pysar.utils.network as pnet
import pysar.utils.plot as pp


sar_sensor_list=['Ers','Env','Jers','Alos','Alos2','Tsx','Csk','Rsat','Rsat2','Sen','Kmps5','G3']

#########################################################################
def log(msg):
    '''Log function writen by Falk'''
    f = open('log','a')
    callingFunction = os.path.basename(inspect.stack()[1][1])
    dateStr = datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%dT%H:%M:%S')
    string = dateStr+" * "+msg
    print(string)
    f.write(string+"\n")
    f.close()

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
    elif  re.search('Sen'    , projectName):  sensor = 'Sen'
    elif  re.search('Kmps5'  , projectName):  sensor = 'Kmps5'
    elif  re.search('Gaofen3', projectName):  sensor = 'G3'
    else: print('satellite not found');  sensor = None
    return sensor


def read_template2inps(templateFile, inps=None):
    '''Read network options from template file into Namespace variable inps'''
    if not inps:
        inps = cmdLineParse()

    ##Read template file
    template = readfile.read_template(templateFile)
    if not template:
        print('Empty template: '+templateFile)
        return None
    prefix = 'select.network.'

    ##Extra keys
    #extra_key_list = ['masterDate','startDate','endDate']
    #for extra_key in extra_key_list:
    #    if extra_key in template.keys():
    #        template[prefix+extra_key] = template[extra_key]

    #Check option prefix
    for i in ['selectPairs.','selectNetwork.']:
        if any(i in key for key in template.keys()):
            print('\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
            print('WARNING: un-supported option prefix detected: {}'.format(i))
            print("         Use select.network. instead")
            print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')

    if all(prefix not in key for key in template.keys()):
        print('\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        print('ERROR: no valid input option deteced in template file!')
        print('Check the template below for supported options:')
        print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
        print(TEMPLATE)
        sys.exit(-1)


    ##Read template dict into inps namespace
    key = prefix+'method'
    if key in template.keys():
        value = template[key]
        if value == 'auto':
            inps.method = 'all'
        else:
            inps.method = value

    key = prefix+'referenceFile'
    if key in template.keys():
        value = template[key]
        if value in ['auto','no']:
            inps.reference_file = None
        else:
            inps.reference_file = value

    key = prefix+'perpBaseMax'
    if key in template.keys():
        value = template[key]
        if value == 'auto':
            inps.perp_base_max = 500.0
        elif value == 'no':
            inps.perp_base_max = 1e5
        else:
            inps.perp_base_max = float(value)

    key = prefix+'tempBaseMax'
    if key in template.keys():
        value = template[key]
        if value == 'auto':
            inps.temp_base_max = 1800.0
        elif value == 'no':
            inps.temp_base_max = 3.65e5
        else:
            inps.temp_base_max = float(value)

    key = prefix+'tempBaseMin'
    if key in template.keys():
        value = template[key]
        if value in ['auto','no']:
            inps.temp_base_min = 0.0
        else:
            inps.temp_base_min = float(value)

    key = prefix+'keepSeasonal'
    if key in template.keys():
        value = template[key]
        if value in ['auto','no']:
            inps.keep_seasonal = False
        else:
            inps.keep_seasonal = True

    key = prefix+'dopOverlapMin'
    if key in template.keys():
        value = template[key]
        if value == 'auto':
            inps.dop_overlap_min = 15.0
        elif value == 'no':
            inps.dop_overlap_min = 0.0
        else:
            inps.dop_overlap_min = float(value)

    key = 'PLATFORM'
    if key in template.keys() and not inps.sensor:
        inps.sensor = template[key]

    key = 'COH_COLOR_JUMP'
    if key in template.keys():
        inps.coh_thres = float(template[key])

    key = prefix+'masterDate'
    if key in template.keys():
        value = template[key]
        if value in ['auto','no']:
            inps.m_date = None
        else:
            inps.m_date = ptime.yymmdd(value)

    key = prefix+'startDate'
    if key in template.keys():
        value = template[key]
        if value in ['auto','no']:
            inps.start_date = None
        else:
            inps.start_date = ptime.yyyymmdd(value)

    key = prefix+'endDate'
    if key in template.keys():
        value = template[key]
        if value in ['auto','no']:
            inps.end_date = None
        else:
            inps.end_date = ptime.yyyymmdd(value)

    key = prefix+'excludeDate'
    if key in template.keys():
        value = template[key]
        if value in ['auto','no']:
            inps.exclude_date = []
        else:
            inps.exclude_date = ptime.yyyymmdd([i for i in value.split(',')])

    key = prefix+'incrementNum'
    if key in template.keys():
        value = template[key]
        if value in ['auto']:
            inps.increment_num = 3
        else:
            inps.increment_num = int(value)

    key = prefix+'tempPerpList'
    if key in template.keys():
        value = template[key]
        if value in ['auto']:
            inps.temp_perp_list = '16,1600;32,800;48,600;64,200'
        else:
            inps.temp_perp_list = value
    if isinstance(inps.temp_perp_list, str):
        inps.temp_perp_list = [[float(j) for j in i.split(',')] for i in inps.temp_perp_list.split(';')]

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

EXAMPLE='''Examples:
  select_network.py KirishimaT246EnvD2.template
  select_network.py KirishimaT246EnvD2.template -b bl_list.txt
  select_network.py KirishimaT246EnvD2.template -b bl_list.txt -o $SC/KirishimaT246EnvD2/PROCESS/ifgram_list.txt
  select_network.py KirishimaT246EnvD2.template -r Pairs.list
  select_network.py KirishimaT246EnvD2.template -r Modified_LoadedData.h5
'''

TEMPLATE='''Template:
## Select network (interferogram combination) in two steps
## 1) select initial network using method / reference File
## 2) filter network/pairs using temp/perp baseline, doppler overlap threshold, etc.
## selection method includes:
##     all          - all possible pairs, pair number = N*(N-1)/2 where N is acquisition num 
##                    default, (Berardino et al., 2002, TGRS).
##     delaunay     - Delaunay triangulation (Fattahi and Amelung, 2013, TGRS).
##                    By default, temporal baseline is normalized using a maxPerpDiff/maxTempDiff
##                    ratio (Pepe and Lanari, 2006, TGRS), use 'delaunay-noweight' to disable normalization. 
##     hierarchical - Select pairs in a hierarchical way using a list of temp/perp thresholds
##                    select.network.tempPerpList
##                    (Zhao, 2015, PhD Thesis)
##                    i.e. 16 days, 1600 m
##                         32 days, 800  m
##                         48 days, 600  m
##                         64 days, 200  m
##     mst          - Minimum Spanning Tree (Perissin and Wang, 2012, TGRS).
##                    Find the MST based on the graph of temporal and perpendicular matrix.
##     sequential   - for each acquisition, select its incrementNum nearest neighbors in the past
##                    (Fattahi and Amelung, 2013, TGRS). 
##                    Pair number = N*m - m*(m-1)/2; it's N when m = 1.
##                    Designed for new satellites like Sentinel-1 and ALOS-2.
##     star / ps    - Star-like/PS-like network/pairs, single common master interferogram
##                    (Ferretti et al., 2001, TGRS)
select.network.method        = auto  #[all / hierarchical / sequential / mst / delaunay / star], auto for all
select.network.referenceFile = auto  #[fname / no], auto for no, HDF5 or text file with pairs info

select.network.perpBaseMax   = auto  #[1-inf / no], auto for 500., max perpendicular spatial baseline
select.network.tempBaseMax   = auto  #[1-inf / no], auto for 1800., max temporal baseline
select.network.tempBaseMin   = auto  #[1-inf], auto for 0.,   min temporal baseline
select.network.keepSeasonal  = auto  #[yes / no], auto for no, keep pairs with seasonal temporal baseline
select.network.dopOverlapMin = auto  #[1-inf / no], auto for 15., min dopploer overlap percentage

select.network.masterDate    = auto  #[100102 / no], auto for no, master date for star/ps network and reference interferogram
select.network.startDate     = auto  #[070101 / no], auto for no, date in YYMMDD or YYYYMMDD format
select.network.endDate       = auto  #[110101 / no], auto for no
select.network.excludeDate   = auto  #[080520,100726 / no], auto for no, exclude dates for pairs selection

select.network.incrementNum  = auto  #[1-inf], auto for 3, for sequential method, pairs num per new acquisition
select.network.tempPerpList  = auto  #[btemp1,bperp1;...], auto for '16,1600;32,800;48,600;64,200'; max temp/perp baseline
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
    parser.add_argument('--show-fig', dest='disp_fig', action='store_true', help='display network ploting result')

    # Method
    method = parser.add_argument_group('Methods to generate the initial network')
    method.add_argument('--method', default='all',\
                        help='network type with info on temp/perp baseline and doppler centroid frequency.')
    method.add_argument('-r', dest='reference_file', default=None,\
                        help='Reference hdf5 / list file with network information. e.g.\n'+\
                             'unwrapIfgram.h5\n'+\
                             'ifgram_list.txt with content as below:'+pnet.IFGRAM_LIST_FILE+\
                             '\nIt could also be generated using plot_network.py --list option, e.g.\n'+\
                             'plot_network.py unwrapIfgram.h5 --list\n\n')
    method.add_argument('--exclude','--ex', dest='exclude_date', nargs='*', \
                        help='date(s) excluded for network selection, e.g. -ex 060713 070831')
    method.add_argument('--start-date', dest='start_date', type=str, help='start/min date of network')
    method.add_argument('--end-date', dest='end_date', type=str, help='end/max date of network')
    method.add_argument('--increment-num', dest='increment_num', type=int, default=3,\
                        help='number of new pairs for each new acquisition, for sequential method')
    method.add_argument('--temp-perp-list', dest='temp_perp_list', default='16,1600;32,800;48,600;64,200',\
                        help='list of max temp/perp baselines, for hierarchical method, e.g.\n'+\
                             '--temp-perp-list 16,1600;32,800;48,600;64,200')
    method.add_argument('--no-norm', dest='norm', action='store_false',\
                        help='do not normalize temp/perp baseline, for delaunay method')
    method.add_argument('--sensor', help='Name of sensor, choose from the list below:\n'+str(sar_sensor_list))
    method.add_argument('--master-date', dest='m_date', help='Master date in YYMMDD or YYYYMMDD format, for star/ps method')

    # Thresholds
    threshold = parser.add_argument_group('Thresholds to filter the initial network')
    threshold.add_argument('--nothreshold', dest='threshold', action='store_false', \
                           help='do not remove pairs using min/max temp/perp baseline and dop\n'+\
                                'auto applied this option when --reference-file is specified.')
    threshold.add_argument('--dop-overlap-min', dest='dop_overlap_min', type=float, default=15.0,\
                           help='min doppler overlap percentage')
    threshold.add_argument('--bperp-max', dest='perp_base_max', type=float, default=500.0,\
                           help='max perpendicular spatial baseline in meters')
    threshold.add_argument('--btemp-min', dest='temp_base_min', type=float, default=0.0, \
                           help='min temporal baseline in days')
    threshold.add_argument('--btemp-max', dest='temp_base_max', type=float, default=1800.0,\
                           help='max temporal baseline in days')
    threshold.add_argument('--keep-seasonal', dest='keep_seasonal', action='store_true',\
                           help='keep seasonal pairs, even they are out of temporal baseline limit\n'+\
                                'i.e. pairs in same/adjcent month within 3 years.')

    parser.add_argument('--inc-angle', dest='inc_angle', type=float, help='Center incidence angle in degrees.')
    inps = parser.parse_args()
    try:    inps.reference_file = glob.glob(inps.reference_file)[0]
    except: inps.reference_file = None
    if inps.temp_perp_list:
        inps.temp_perp_list = [[float(j) for j in i.split(',')] for i in inps.temp_perp_list.split(';')]

    return inps


#########################################################################
def main(argv):
    
    # Read inputs
    inps = cmdLineParse()
    inps = read_template2inps(inps.template_file, inps)
    log(os.path.basename(sys.argv[0])+' '+inps.template_file)

    project_name = os.path.splitext(os.path.basename(inps.template_file))[0]
    print('project name: '+project_name)
    if not inps.sensor:
        inps.sensor = project_name2sensor(project_name)
 
    # Auto path setting for Miami user
    if not inps.baseline_file and pysar.miami_path and 'SCRATCHDIR' in os.environ:
        if pysar.miami_path and 'SCRATCHDIR' in os.environ:
            try:    inps.baseline_file = glob.glob(os.getenv('SCRATCHDIR')+'/'+project_name+'/SLC/bl_list.txt')[0]
            except: inps.baseline_file = None

    # Pair selection from reference
    if inps.reference_file:
        print('Use pairs info from reference file: '+inps.reference_file)
        date12_list = pnet.get_date12_list(inps.reference_file)
        date12_list = [i.replace('_','-') for i in date12_list]

        if inps.baseline_file:
            date8_list, pbase_list, dop_list = pnet.read_baseline_file(inps.baseline_file)[0:3]
            date6_list = ptime.yymmdd(date8_list)
            tbase_list = ptime.date_list2tbase(date8_list)[0]

    # Pair selection from temp/perp/dop baseline info
    else:
        if not inps.baseline_file:
            raise Exception('ERROR: No baseline file found!')

        # Check start/end/exclude date
        date8_list = pnet.read_baseline_file(inps.baseline_file)[0]
        inps.exclude_date = ptime.yyyymmdd(inps.exclude_date)
        if not inps.exclude_date:
            inps.exclude_date = []
        else:
            print('input exclude dates: '+str(inps.exclude_date))
        if inps.start_date:
            print('input start date: '+inps.start_date)
            inps.exclude_date += [i for i in date8_list if float(i) < float(ptime.yyyymmdd(inps.start_date))]
            inps.exclude_date = sorted(inps.exclude_date)
        if inps.end_date:
            print('input end   date: '+inps.end_date)
            inps.exclude_date += [i for i in date8_list if float(i) > float(ptime.yyyymmdd(inps.end_date))]
            inps.exclude_date = sorted(inps.exclude_date)
        if inps.exclude_date:
            print('exclude    dates: ')
            print(inps.exclude_date)

        # Read baseline list file: bl_list.txt
        inps.exclude_date = ptime.yymmdd(inps.exclude_date)
        date8_list, pbase_list, dop_list = pnet.read_baseline_file(inps.baseline_file, inps.exclude_date)[0:3]
        date6_list = ptime.yymmdd(date8_list)
        tbase_list = ptime.date_list2tbase(date8_list)[0]

        # Initial network using input methods
        inps.method = inps.method.lower().replace('-','_')
        if inps.method in ['star','ps']:       inps.method = 'star'
        elif inps.method.startswith('seq'):    inps.method = 'sequential'
        elif inps.method.startswith('hierar'): inps.method = 'hierarchical'
        elif inps.method in ['mst','min_spanning_tree','minimum_spanning_tree']:  inps.method = 'mst'
        print('select method: '+inps.method)

        if   inps.method == 'all':
            date12_list = pnet.select_pairs_all(date6_list)
        elif inps.method == 'delaunay':
            date12_list = pnet.select_pairs_delaunay(date6_list, pbase_list, inps.norm)
        elif inps.method == 'star':
            date12_list = pnet.select_pairs_star(date6_list)
        elif inps.method == 'sequential':
            date12_list = pnet.select_pairs_sequential(date6_list, inps.increment_num)
        elif inps.method == 'hierarchical':
            date12_list = pnet.select_pairs_hierarchical(date6_list, pbase_list, inps.temp_perp_list)
        elif inps.method == 'mst':
            date12_list = pnet.select_pairs_mst(date6_list, pbase_list)
        else:
            raise Exception('Unrecoganized select method: '+inps.method)
        print('initial number of interferograms: '+str(len(date12_list)))

        # Filter pairs (optional) using temp/perp/doppler baseline threshold
        if inps.method in ['star','hierarchical','mst']:
            inps.threshold = False
        if inps.threshold:
            # Temporal baseline
            date12_list = pnet.threshold_temporal_baseline(date12_list, inps.temp_base_max,\
                                                           inps.keep_seasonal, inps.temp_base_min)
            print('number of interferograms after filtering of <%d, %d> days in temporal baseline: %d'\
                  % (inps.temp_base_min, inps.temp_base_max, len(date12_list)))
            if inps.keep_seasonal:
                print('\tkeep seasonal pairs, i.e. pairs with temporal baseline == N*years +/- one month')

            # Perpendicular spatial baseline
            date12_list = pnet.threshold_perp_baseline(date12_list, date6_list, pbase_list, inps.perp_base_max)
            print('number of interferograms after filtering of max %d meters in perpendicular baseline: %d'\
                  % (inps.perp_base_max, len(date12_list)))

            # Doppler Overlap Percentage
            if inps.sensor:
                bandwidth_az = pnet.azimuth_bandwidth(inps.sensor)
                date12_list = pnet.threshold_doppler_overlap(date12_list, date6_list, dop_list,\
                                                             bandwidth_az, inps.dop_overlap_min/100.0)
                print('number of interferograms after filtering of min '+str(inps.dop_overlap_min)+'%'+\
                      ' overlap in azimuth Doppler frequency: '+str(len(date12_list)))

    # Write ifgram_list.txt
    if not date12_list:
        print('WARNING: No interferogram selected!')
        return None

    # date12_list to date_list
    m_dates = [date12.replace('_','-').split('-')[0] for date12 in date12_list]
    s_dates = [date12.replace('_','-').split('-')[1] for date12 in date12_list]
    try: print('number of acquisitions   input   : '+str(len(date6_list)))
    except: pass
    print('number of acquisitions   selected: '+str(len(list(set(m_dates + s_dates)))))
    print('number of interferograms selected: '+str(len(date12_list)))

    # Output directory/filename
    if not inps.outfile:
        if pysar.miami_path and 'SCRATCHDIR' in os.environ:
            inps.out_dir = os.getenv('SCRATCHDIR')+'/'+project_name+'/PROCESS'
        else:
            try:    inps.out_dir = os.path.dirname(os.path.abspath(inps.reference_file))
            except: inps.out_dir = os.path.dirname(os.path.abspath(inps.baseline_file))
        inps.outfile = inps.out_dir+'/ifgram_list.txt'
    inps.outfile = os.path.abspath(inps.outfile)
    inps.out_dir = os.path.dirname(inps.outfile)
    if not os.path.isdir(inps.out_dir):
        os.makedirs(inps.out_dir)

    print('writing >>> '+inps.outfile)
    if not inps.baseline_file:
        np.savetxt(inps.outfile, date12_list, fmt='%s')
        return inps.outfile

    ## Calculate Bperp, Btemp and predicted coherence
    ifgram_num = len(date12_list)
    ifgram_pbase_list = []
    ifgram_tbase_list = []

    for i in range(ifgram_num):
        m_date, s_date = date12_list[i].split('-')
        m_idx = date6_list.index(m_date)
        s_idx = date6_list.index(s_date)
        pbase = pbase_list[s_idx] - pbase_list[m_idx]
        tbase = tbase_list[s_idx] - tbase_list[m_idx]
        ifgram_pbase_list.append(pbase)
        ifgram_tbase_list.append(tbase)

    try:
        inps.coherence_list = pnet.simulate_coherence(date12_list, inps.baseline_file, sensor=inps.sensor).flatten().tolist()
        inps.cbar_label = 'Simulated coherence'
    except:
        inps.coherence_list = None

    ##### Write txt file
    fl = open(inps.outfile, 'w')
    fl.write('#Interferograms configuration generated by select_network.py\n')
    fl.write('#   Date12      Btemp(days)    Bperp(m)    sim_coherence\n')
    for i in range(len(date12_list)):
        line = '%s   %6.0f         %6.1f' % (date12_list[i], ifgram_tbase_list[i], ifgram_pbase_list[i])
        if inps.coherence_list:
            line += '       %1.4f' % (inps.coherence_list[i])
        fl.write(line+'\n')
    fl.close()


    ##### Plot network info
    if not inps.disp_fig:
        plt.switch_backend('Agg')

    out_fig_name = 'BperpHistory.pdf'
    print('plotting baseline history in temp/perp baseline domain to file: '+out_fig_name)
    fig2, ax2 = plt.subplots()
    ax2 = pp.plot_perp_baseline_hist(ax2, date8_list, pbase_list)
    plt.savefig(inps.out_dir+'/'+out_fig_name, bbox_inches='tight')

    out_fig_name = 'Network.pdf'
    print('plotting network / pairs  in temp/perp baseline domain to file: '+out_fig_name)
    fig1, ax1 = plt.subplots()
    ax1 = pnet.plot_network(ax1, date12_list, date8_list, pbase_list, plot_dict=vars(inps), print_msg=False)
    plt.savefig(inps.out_dir+'/'+out_fig_name, bbox_inches='tight')

    out_fig_name = 'CoherenceMatrix.pdf'
    if inps.coherence_list:
        print('plotting predicted coherence matrix to file: '+out_fig_name)
        fig3, ax3 = plt.subplots()
        ax3 = pp.plot_coherence_matrix(ax3, date12_list, inps.coherence_list, plot_dict=vars(inps))
        plt.savefig(inps.out_dir+'/'+out_fig_name, bbox_inches='tight')

    if inps.disp_fig:
        plt.show()

    return inps.outfile

###########################################################################
if __name__ == '__main__':
    main(sys.argv[1:])

 
