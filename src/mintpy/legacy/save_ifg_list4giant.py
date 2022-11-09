#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Jan 2018                           #
############################################################


import argparse
import os

from mintpy.objects import ifgramStack
from mintpy.utils import ptime, readfile, utils as ut

##################################################################################################
EXAMPLE = """example:
  save_ifg_list4giant.py  filt_*.unw
  save_ifg_list4giant.py  inputs/ifgramStack.h5  --sensor SEN
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Prepare ifg.list file for GIAnT.\n',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)
    parser.add_argument('file', nargs='+', help='Interferogram file(s)')
    parser.add_argument('--sensor', '--mission',
                        dest='sensor', help='Sensor name of SAR data')
    parser.add_argument('-o', '--output', dest='outfile',
                        default='ifg.list', help='Output list file')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    inps.file = ut.get_file_list(inps.file, abspath=True)
    inps.outfile = os.path.abspath(inps.outfile)
    return inps


def get_mission_name(meta_dict):
    """Get mission name in UNAVCO InSAR Archive format from attribute mission/PLATFORM
    Input:  meta_dict : dict, attributes
    Output: mission   : string, mission name in standard UNAVCO format.
    """
    mission = None

    if 'mission' in meta_dict.keys():
        value = meta_dict['mission'].lower()
    elif 'PLATFORM' in meta_dict.keys():
        value = meta_dict['PLATFORM'].lower()
    else:
        print('No PLATFORM nor mission attribute found, can not identify mission name.')
        print('return None')
        return mission

    # Convert to UNAVCO Mission name
    ## ERS, ENV, S1, RS1, RS2, CSK, TSX, JERS, ALOS, ALOS2
    if value.startswith('ers'):
        mission = 'ERS'
    elif value.startswith(('env', 'asar')):
        mission = 'ENV'
    elif value.startswith(('s1', 'sen')):
        mission = 'S1'
    elif value.startswith(('rs', 'rsat', 'radarsat')):
        mission = 'RS'
        if value.endswith('1'):
            mission += '1'
        else:
            mission += '2'
    elif value.startswith(('csk', 'cos')):
        mission = 'CSK'
    elif value.startswith(('tsx', 'tdx', 'terra', 'tandem')):
        mission = 'TSX'
    elif value.startswith('jers'):
        mission = 'JERS'
    elif value.startswith(('alos', 'palsar')):
        if value.endswith('2'):
            mission = 'ALOS2'
        else:
            mission = 'ALOS'
    else:
        print('Un-recognized PLATFORM attribute: '+value)
        print('return None')
    return mission


def get_giant_ifg_list(fnames):
    m_date_list = []
    s_date_list = []
    pbase_list = []

    ext = os.path.splitext(fnames[0])[1]
    if ext == '.h5':
        obj = ifgramStack(fnames[0])
        obj.open()
        m_date_list = obj.mDates[obj.dropIfgram].tolist()
        s_date_list = obj.sDates[obj.dropIfgram].tolist()
        pbase_list = obj.pbaseIfgram[obj.dropIfgram].tolist()

    else:
        ifgramNum = len(fnames)
        print('Number of interferograms: %d' % (ifgramNum))
        for fname in fnames:
            atr = readfile.read_attribute(fname)
            m_date, s_date = ptime.yymmdd(atr['DATE12'].split('-'))
            pbase = (float(atr['P_BASELINE_TOP_HDR']) +
                     float(atr['P_BASELINE_BOTTOM_HDR'])) / 2.
            m_date_list.append(m_date)
            s_date_list.append(s_date)
            pbase_list.append(pbase)
    return m_date_list, s_date_list, pbase_list


##################################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    atr = readfile.read_attribute(inps.file[0])

    if not inps.sensor:
        inps.sensor = get_mission_name(atr)
    print('Sensor name: %s' % (inps.sensor))

    m_date_list, s_date_list, pbase_list = get_giant_ifg_list(inps.file)

    # Write to text file
    with open(inps.outfile, 'w') as f:
        for i in range(len(pbase_list)):
            f.write('{} {}     {:<10.2f}   {}\n'.format(m_date_list[i],
                                                        s_date_list[i],
                                                        pbase_list[i],
                                                        inps.sensor))
    print('write to %s' % (inps.outfile))
    return


###################################################################################################
if __name__ == '__main__':
    main()
