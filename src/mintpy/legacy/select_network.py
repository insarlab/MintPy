#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2017                               #
############################################################


import argparse
import datetime
import glob
import inspect
import os
import sys

import matplotlib.pyplot as plt

import mintpy
from mintpy.defaults.template import get_template_content
from mintpy.objects import ifgramStack, sensor
from mintpy.utils import (
    network as pnet,
    plot as pp,
    ptime,
    readfile,
    utils as ut,
)

SENSOR_NAMES = [i.capitalize() for i in sensor.SENSOR_NAMES]


#########################################################################
REFERENCE = """references:
  Berardino, P., G. Fornaro, R. Lanari, and E. Sansosti (2002), A new algorithm for surface deformation monitoring
    based on small baseline differential SAR interferograms, IEEE TGRS, 40(11), 2375-2383.
  Fattahi, H., and F. Amelung (2013), DEM Error Correction in InSAR Time Series, IEEE TGRS, 51(7), 4249-4259.
  Ferretti, A., C. Prati, and F. Rocca (2001), Permanent scatterers in SAR interferometry, IEEE TGRS, 39(1), 8-20.
  Pepe, A., and R. Lanari (2006), On the extension of the minimum cost flow algorithm for phase unwrapping
    of multitemporal differential SAR interferograms, IEEE TGRS, 44(9), 2374-2383.
  Perissin D., Wang T. (2012), Repeat-pass SAR interferometry with partially coherent targets. IEEE TGRS. 271-280
  Yunjun, Z., H. Fattahi, and F. Amelung (2019), Small baseline InSAR time series analysis: Unwrapping error
    correction and noise reduction, Computers & Geosciences, 133, 104331, doi:10.1016/j.cageo.2019.104331.
  Zebker, H. A., and J. Villasenor (1992), Decorrelation in interferometric radar echoes, IEEE TGRS, 30(5), 950-959.
  Zhao, W., (2017), Small deformation detected from InSAR time-series and their applications in geophysics, Doctoral
    dissertation, Univ. of Miami, Section 6.3.
"""

EXAMPLE = """examples:
  select_network.py KirishimaT246EnvD2.template
  select_network.py KirishimaT246EnvD2.template -b bl_list.txt
  select_network.py KirishimaT246EnvD2.template -b bl_list.txt -o $SC/KirishimaT246EnvD2/PROCESS/ifgram_list.txt
  select_network.py KirishimaT246EnvD2.template -r Pairs.list
  select_network.py KirishimaT246EnvD2.template -r Modified_LoadedData.h5
"""

TEMPLATE = get_template_content('modify_network')


def create_parser():
    parser = argparse.ArgumentParser(description='Select Interferometric Network / Pairs.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=REFERENCE+'\n'+TEMPLATE+'\n'+EXAMPLE)

    parser.add_argument('template_file', help='template file with options')
    parser.add_argument('-b', dest='baseline_file', default='bl_list.txt',
                        help='baseline list file of all SLCs, e.g.'+pnet.BASELINE_LIST_FILE)
    parser.add_argument('-o', '--outfile',
                        help='Output list file for network, ifgram_list.txt by default.')

    # Figure
    fig = parser.add_argument_group('Figure settings')
    fig.add_argument('--figsize', dest='fig_size', type=float, nargs=2,
                     help='figure size in inches - width and length')
    fig.add_argument('--ms', '--markersize', dest='markersize',
                     type=int, default=16, help='marker size in points')
    fig.add_argument('--fs', '--fontsize', type=int,
                     default=12, help='font size in points')
    fig.add_argument('--show-fig', dest='disp_fig', action='store_true',
                     help='display network plotting result')
    fig.add_argument('--figext', dest='figext', default='.pdf',
                     help='file extension to be saved.')
    fig.add_argument('--dpi', dest='figdpi', type=int, default=150,
                     help='Figure dpi to be saved.')
    fig.add_argument('--coh-threshold', dest='coh_thres', type=float, default=0.4,
                     help='Coherence threshold for colormap jump')
    fig.add_argument('--notitle', dest='disp_title', action='store_false',
                     help='Do not display figure title.')
    fig.add_argument('--number', dest='number', type=str,
                     help='number mark to be plot at the corner of figure.')

    # Method
    method = parser.add_argument_group('Methods to generate the initial network')
    method.add_argument('--method', default='all',
                        help='network type with info on temp/perp baseline and doppler centroid frequency.')
    method.add_argument('-r', dest='referenceFile', default=None,
                        help='Reference hdf5 / list file with network information. e.g.\n' +
                             'ifgramStack.h5\n' +
                             'ifgram_list.txt with content as below:'+pnet.IFGRAM_LIST_FILE +
                             '\nIt could also be generated using plot_network.py --list option, e.g.\n' +
                             'info.py ifgramStack.h5 --date --show kept > date12_list.txt\n\n')
    method.add_argument('--exclude', '--ex', dest='excludeDate', nargs='*', default=[],
                        help='date(s) excluded for network selection, e.g. -ex 060713 070831')
    method.add_argument('--start-date', dest='startDate',
                        type=str, help='start/min date of network')
    method.add_argument('--end-date', dest='endDate',
                        type=str, help='end/max date of network')
    method.add_argument('--conn-num', dest='connNum', type=int, default=3,
                        help='number of new pairs for each new acquisition, for sequential method')
    method.add_argument('--temp-perp-list', dest='tempPerpList', default='16,1600;32,800;48,600;64,200',
                        help='list of max temp/perp baselines, for hierarchical method, e.g.\n' +
                             '--temp-perp-list 16,1600;32,800;48,600;64,200')
    method.add_argument('--no-norm', dest='norm', action='store_false',
                        help='do not normalize temp/perp baseline, for delaunay method')
    method.add_argument('--sensor', help='Name of sensor, choose from the list below:\n'+str(SENSOR_NAMES))
    method.add_argument('--reference-date', dest='referenceDate',
                        help='reference date in YYMMDD or YYYYMMDD format, for star/ps method')

    # Thresholds
    threshold = parser.add_argument_group('Thresholds to filter the initial network')
    threshold.add_argument('--nothreshold', dest='threshold', action='store_false',
                           help='do not remove pairs using min/max temp/perp baseline and dop\n' +
                                'auto applied this option when --reference-file is specified.')
    threshold.add_argument('--dop-overlap-min', dest='dopOverlapMin', type=float, default=0.0,
                           help='min doppler overlap percentage')
    threshold.add_argument('--bperp-max', dest='perpBaseMax', type=float, default=1e5,
                           help='max perpendicular spatial baseline in meters')
    threshold.add_argument('--btemp-min', dest='tempBaseMin', type=float, default=0.0,
                           help='min temporal baseline in days')
    threshold.add_argument('--btemp-max', dest='tempBaseMax', type=float, default=3.65e5,
                           help='max temporal baseline in days')
    threshold.add_argument('--keep-seasonal', dest='keepSeasonal', action='store_true',
                           help='keep seasonal pairs, even they are out of temporal baseline limit\n' +
                                'i.e. pairs in same/adjacent month within 3 years.')

    parser.add_argument('--inc-angle', dest='inc_angle',
                        type=float, help='Center incidence angle in degrees.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    if not os.path.isfile(inps.baseline_file):
        inps.baseline_file = None

    try:
        inps.referenceFile = glob.glob(inps.referenceFile)[0]
    except:
        inps.referenceFile = None

    return inps


#########################################################################
def log(msg):
    """Log function written by Falk"""
    f = open('log', 'a')
    callingFunction = os.path.basename(inspect.stack()[1][1])
    dateStr = datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%dT%H:%M:%S')
    string = dateStr+" * "+msg
    print(string)
    f.write(string+"\n")
    f.close()


def read_template2inps(templateFile, inps=None):
    """Read network options from template file into Namespace variable inps"""
    if not inps:
        inps = cmd_line_parse()
    inpsDict = vars(inps)

    # Read template file
    template = readfile.read_template(templateFile)
    auto_file = os.path.join(os.path.dirname(mintpy.__file__), 'defaults/selectNetwork.cfg')
    template = ut.check_template_auto_value(template, auto_file=auto_file)
    if not template:
        log('Empty template: '+templateFile)
        return None

    prefix = 'selectNetwork.'
    # Check obsolete option prefix
    for i in ['selectPairs.', 'select.network.']:
        if any(i in key for key in template.keys()):
            msg = f'obsolete option prefix detected: {i}\n'
            msg += f'Use {prefix} instead'
            raise Exception(msg)
    if all(prefix not in key for key in template.keys()):
        msg = 'no valid input option detected in template file!\n'
        msg += 'Check the template below for supported options:\n'
        msg += TEMPLATE
        raise Exception(msg)

    # convert template into inpsDict
    keyList = [i for i in list(inpsDict.keys()) if prefix+i in template.keys()]
    for key in keyList:
        value = template[prefix+key]
        # bool
        if key in ['keepSeasonal']:
            inpsDict[key] = value
        elif value:
            # str
            if key in ['method', 'referenceFile', 'tempPerpList']:
                inpsDict[key] = value
            # date in YYYYMMDD
            elif key in ['referenceDate', 'startDate', 'endDate']:
                inpsDict[key] = ptime.yyyymmdd(value)
            # list of dates in YYYYMMDD
            elif key in ['excludeDate']:
                inps.excludeDate = ptime.yyyymmdd([i.strip() for i in value.split(',')])
            # float
            elif key in ['perpBaseMax', 'tempBaseMax', 'tempBaseMin', 'dopOverlapMin']:
                inpsDict[key] = float(value)
            # int
            elif key in ['connNum']:
                inpsDict[key] = int(value)

    # read tempPerpList from str
    if isinstance(inps.tempPerpList, str):
        inps.tempPerpList = [[float(j) for j in i.split(',')]
                              for i in inps.tempPerpList.split(';')]

    # Initial network using input methods
    inps.method = inps.method.lower().replace('-', '_')
    if inps.method in ['star', 'ps']:
        inps.method = 'star'
    elif inps.method.startswith('seq'):
        inps.method = 'sequential'
    elif inps.method.startswith('hierar'):
        inps.method = 'hierarchical'
    elif inps.method in ['mst', 'min_spanning_tree', 'minimum_spanning_tree']:
        inps.method = 'mst'
    elif inps.method in ['all', 'sb']:
        inps.method = 'all'

    # for coherence prediction
    key = 'PLATFORM'
    if key in template.keys() and not inps.sensor:
        inps.sensor = template[key]

    key = 'COH_COLOR_JUMP'
    if key in template.keys():
        inps.coh_thres = float(template[key])

    # project name and sensor
    project_name = os.path.splitext(os.path.basename(inps.template_file))[0]
    log('project name: '+project_name)
    if not inps.sensor:
        inps.sensor = sensor.project_name2sensor_name(project_name)[0]

    # Output directory/filename
    if not inps.outfile:
        if 'SCRATCHDIR' in os.environ:
            inps.out_dir = os.getenv('SCRATCHDIR')+'/'+project_name+'/PROCESS'
        else:
            try:
                inps.out_dir = os.path.dirname(os.path.abspath(inps.referenceFile))
            except:
                inps.out_dir = os.path.dirname(os.path.abspath(inps.baseline_file))
        inps.outfile = inps.out_dir+'/ifgram_list.txt'

    # Auto path of bl_list.txt file (for Miami user)
    if not inps.baseline_file and 'SCRATCHDIR' in os.environ:
        bl_file = os.path.join(os.getenv('SCRATCHDIR'), f'{project_name}/SLC/bl_list.txt')
        if os.path.isfile(bl_file):
            inps.baseline_file = bl_file

    if not inps.referenceFile and not inps.baseline_file:
        raise Exception('No baseline file or reference file found! At least one is required.')

    return inps


def read_baseline_info(baseline_file, reference_file):
    """Read date, bperp and/or DOP info
    Parameters: baseline_file : str, path of bl_list.txt file
                reference_file : str, path of ifgramStack.h5 file
    Returns:    date_list : list of str in YYMMDD format
                tbase_list : list of int in days
                pbase_list : list of float in meter
                dop_list : None, list of 1D array in size of (3,)
    """
    dop_list = None
    if baseline_file:
        date_list, pbase_list, dop_list = pnet.read_baseline_file(baseline_file)[0:3]
        date_list = ptime.yymmdd(date_list)
        tbase_list = ptime.date_list2tbase(date_list)[0]

    elif reference_file:
        obj = ifgramStack(reference_file)
        date12_list_all = obj.get_date12_list(dropIfgram=False)
        date12_list_all = ptime.yymmdd_date12(date12_list_all)
        m_dates = [i.split('-')[0] for i in date12_list_all]
        s_dates = [i.split('-')[1] for i in date12_list_all]
        date_list = sorted(list(set(m_dates + s_dates)))
        tbase_list = ptime.date_list2tbase(date_list)[0]

        pbase_list = obj.get_perp_baseline_timeseries(dropIfgram=False).tolist()
    return date_list, tbase_list, pbase_list, dop_list


def read_exclude_date(inps, date_list_all):
    inps.excludeDate = ptime.yyyymmdd(inps.excludeDate)
    if not inps.excludeDate:
        inps.excludeDate = []
    else:
        log(f'input exclude dates: {inps.excludeDate}')

    if inps.startDate:
        log('input start date: '+inps.startDate)
        inps.excludeDate += [i for i in date_list_all
                             if float(i) < float(ptime.yyyymmdd(inps.startDate))]
        inps.excludeDate = sorted(inps.excludeDate)

    if inps.endDate:
        log('input end   date: '+inps.endDate)
        inps.excludeDate += [i for i in date_list_all
                             if float(i) > float(ptime.yyyymmdd(inps.endDate))]
        inps.excludeDate = sorted(inps.excludeDate)

    if inps.excludeDate:
        log(f'exclude    dates: ({len(inps.excludeDate)})\n{inps.excludeDate}')
    return inps.excludeDate


def select_network_candidate(inps):
    date_list, tbase_list, pbase_list, dop_list = read_baseline_info(baseline_file=inps.baseline_file,
                                                                     reference_file=inps.referenceFile)

    # Pair selection from reference
    if inps.referenceFile:
        log(f'select initial network from reference file: {inps.referenceFile}')
        stack_obj = ifgramStack(inps.referenceFile)
        date12_list = stack_obj.get_date12_list(dropIfgram=True)
        date12_list = ptime.yymmdd_date12(date12_list)

    # Pais selection from method
    elif inps.baseline_file:
        log(f'select initial network with method: {inps.method}')
        if inps.method == 'all':
            date12_list = pnet.select_pairs_all(date_list)
        elif inps.method == 'delaunay':
            date12_list = pnet.select_pairs_delaunay(date_list, pbase_list, inps.norm)
        elif inps.method == 'star':
            date12_list = pnet.select_pairs_star(date_list)
        elif inps.method == 'sequential':
            date12_list = pnet.select_pairs_sequential(date_list, inps.connNum)
        elif inps.method == 'hierarchical':
            date12_list = pnet.select_pairs_hierarchical(date_list, pbase_list, inps.tempPerpList)
        elif inps.method == 'mst':
            date12_list = pnet.select_pairs_mst(date_list, pbase_list)
        else:
            raise Exception('Unrecoganized select method: '+inps.method)

    log(f'initial number of interferograms: {len(date12_list)}')
    inps.date12_list = date12_list
    inps.date_list = date_list
    inps.tbase_list = tbase_list
    inps.pbase_list = pbase_list
    inps.dop_list = dop_list
    return inps


def prune_network(date12_list, inps):
    """Pruning network candidates based on temp/perp baseline and DOP overlap"""

    # Turn off pruning for some methods
    if inps.referenceFile or inps.method in ['star', 'hierarchical', 'mst']:
        inps.threshold = False

    if inps.threshold:
        # Temporal baseline
        date12_list = pnet.threshold_temporal_baseline(date12_list,
                                                       inps.tempBaseMax,
                                                       inps.keepSeasonal,
                                                       inps.tempBaseMin)
        log('number of interferograms with temporal baseline in <{}, {}> days: {}'.format(inps.tempBaseMin,
                                                                                            inps.tempBaseMax,
                                                                                            len(date12_list)))
        if inps.keepSeasonal:
            log('\tkeep seasonal pairs, i.e. pairs with temporal baseline == N*years +/- one month')

        # Perpendicular spatial baseline
        date12_list = pnet.threshold_perp_baseline(date12_list,
                                                   inps.date_list,
                                                   inps.pbase_list,
                                                   inps.perpBaseMax)
        log('number of interferograms with perp baseline <= {} meters: {}'.format(inps.perpBaseMax,
                                                                                    len(date12_list)))

        # Doppler Overlap Percentage
        if inps.sensor and inps.dop_list:
            bandwidth_az = sensor.SENSOR_DICT[inps.sensor.lower()]['doppler_bandwidth']
            date12_list = pnet.threshold_doppler_overlap(date12_list,
                                                         inps.date_list,
                                                         inps.dop_list,
                                                         bandwidth_az,
                                                         inps.dopOverlapMin/100.0)
            log('number of interferograms with azimuth Doppler frequency overlap >= {}%: {}'.format(inps.dopOverlapMin,
                                                                                                      len(date12_list)))

    if not date12_list:
        raise Exception('No interferogram left after pruning!')

    # print stat info
    m_dates = [date12.replace('_', '-').split('-')[0] for date12 in date12_list]
    s_dates = [date12.replace('_', '-').split('-')[1] for date12 in date12_list]
    log(f'number of acquisitions   input   : {len(inps.date_list)}')
    log(f'number of acquisitions   selected: {len(list(set(m_dates + s_dates)))}')
    log(f'number of interferograms selected: {len(date12_list)}')

    return date12_list


def write_ifgram_list(inps):
    # Output directory/filename
    inps.outfile = os.path.abspath(inps.outfile)
    inps.out_dir = os.path.dirname(inps.outfile)
    if not os.path.isdir(inps.out_dir):
        os.makedirs(inps.out_dir)

    # Calculate ifgram's temp/perp baseline
    ifgram_num = len(inps.date12_list)
    ifgram_pbase_list = []
    ifgram_tbase_list = []
    for i in range(ifgram_num):
        m_date, s_date = inps.date12_list[i].split('-')
        m_idx = inps.date_list.index(m_date)
        s_idx = inps.date_list.index(s_date)
        pbase = inps.pbase_list[s_idx] - inps.pbase_list[m_idx]
        tbase = inps.tbase_list[s_idx] - inps.tbase_list[m_idx]
        ifgram_pbase_list.append(pbase)
        ifgram_tbase_list.append(tbase)

    # calculate ifgram's predicted coherence
    try:
        inps.cohList = pnet.simulate_coherence(inps.date12_list,
                                                      inps.baseline_file,
                                                      sensor_name=inps.sensor).flatten().tolist()
    except:
        inps.cohList = None

    # Write txt file
    f = open(inps.outfile, 'w')
    f.write('#Interferograms configuration generated by select_network.py\n')
    f.write('#   Date12      Btemp(days)    Bperp(m)    sim_coherence\n')
    for i in range(ifgram_num):
        line = '{}   {:6.0f}         {:6.1f}'.format(inps.date12_list[i],
                                                     ifgram_tbase_list[i],
                                                     ifgram_pbase_list[i])
        if inps.cohList:
            line += f'       {inps.cohList[i]:1.4f}'
        f.write(line+'\n')
    f.close()
    log(f'write network/pairs info into file: {inps.outfile}')
    return inps.outfile


def plot_network_info(inps):
    if not inps.disp_fig:
        plt.switch_backend('Agg')

    out_fig_name = os.path.join(inps.out_dir, f'network{inps.figext}')
    log('plot network / pairs to file: '+os.path.basename(out_fig_name))
    fig1, ax1 = plt.subplots(figsize=inps.fig_size)
    ax1 = pp.plot_network(ax1,
                          inps.date12_list,
                          inps.date_list,
                          inps.pbase_list,
                          p_dict=vars(inps),
                          print_msg=False)
    plt.savefig(out_fig_name, bbox_inches='tight', dpi=inps.figdpi)

    out_fig_name = os.path.join(inps.out_dir, f'bperpHistory{inps.figext}')
    log('plot baseline history to file: '+os.path.basename(out_fig_name))
    fig2, ax2 = plt.subplots(figsize=inps.fig_size)
    ax2 = pp.plot_perp_baseline_hist(ax2,
                                     inps.date_list,
                                     inps.pbase_list)
    plt.savefig(out_fig_name, bbox_inches='tight', dpi=inps.figdpi)

    out_fig_name = os.path.join(inps.out_dir, f'coherenceMatrix{inps.figext}')
    if inps.cohList:
        log('plot predicted coherence matrix to file: '+os.path.basename(out_fig_name))
        fig3, ax3 = plt.subplots(figsize=inps.fig_size)
        ax3 = pp.plot_coherence_matrix(ax3,
                                       inps.date12_list,
                                       inps.cohList,
                                       p_dict=vars(inps))[0]
        plt.savefig(out_fig_name, bbox_inches='tight', dpi=inps.figdpi)

    if inps.disp_fig:
        plt.show()
    return


#########################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    log(os.path.basename(sys.argv[0])+' '+inps.template_file)

    inps = read_template2inps(inps.template_file, inps)

    inps = select_network_candidate(inps)

    inps.date12_list = prune_network(date12_list=inps.date12_list, inps=inps)

    inps.outfile = write_ifgram_list(inps)

    plot_network_info(inps)

    return inps.outfile


###########################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
