#!/usr/bin/env python3
# Author: Zhang Yunjun, 2018-Jan-27

import os
import sys
import re
import argparse

import h5py

import pysar.utils.datetime as ptime
import pysar.utils.readfile as readfile
import pysar.utils.utils as ut


def get_mission_name(meta_dict):
    '''Get mission name in UNAVCO InSAR Archive format from attribute mission/PLATFORM
    Input:  meta_dict : dict, attributes
    Output: mission   : string, mission name in standard UNAVCO format.
    '''
    mission = None

    if 'mission' in meta_dict.keys():
        value = meta_dict['mission'].lower()
    elif 'PLATFORM' in meta_dict.keys():
        value = meta_dict['PLATFORM'].lower()
    else:
        print('No PLATFORM nor mission attribute found, can not identify mission name.')
        print('return None')
        return mission

    ## Convert to UNAVCO Mission name
    ## ERS, ENV, S1, RS1, RS2, CSK, TSX, JERS, ALOS, ALOS2
    if value.startswith('ers'):
        mission = 'ERS'
    elif value.startswith(('env','asar')):
        mission = 'ENV'
    elif value.startswith(('s1','sen')):
        mission = 'S1'
    elif value.startswith(('rs','rsat','radarsat')):
        mission = 'RS'
        if value.endswith('1'):
            mission += '1'
        else:
            mission += '2'
    elif value.startswith(('csk','cos')):
        mission = 'CSK'
    elif value.startswith(('tsx','tdx','terra','tandem')):
        mission = 'TSX'
    elif value.startswith('jers'):
        mission = 'JERS'
    elif value.startswith(('alos','palsar')):
        if value.endswith('2'):
            mission = 'ALOS2'
        else:
            mission = 'ALOS'
    else:
        print('Un-recognized PLATFORM attribute: '+value)
        print('return None')
    return mission


##################################################################################################
EXAMPLE='''example:
  prep_giant_ifg_list.py  filt_*.unw
  prep_giant_ifg_list.py  unwrapIfgram.h5  --sensor SEN
'''

def cmdLineParse():
    parser = argparse.ArgumentParser(description='Prepare ifg.list file for GIAnT.\n',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)
    parser.add_argument('file', nargs='+', help='Interferogram file(s)')
    parser.add_argument('--sensor','--mission', dest='sensor', help='Sensor name of SAR data')
    parser.add_argument('-o','--output',dest='outfile', default='ifg.list', help='Output list file')
    inps = parser.parse_args()
    return inps


##################################################################################################
def main(argv):
    inps = cmdLineParse()
    inps.outfile = os.path.abspath(inps.outfile)
    atr = readfile.read_attribute(inps.file[0])
    k = atr['FILE_TYPE']

    if not inps.sensor:
        inps.sensor = get_mission_name(atr)
    print('Sensor name: %s' % (inps.sensor))

    m_date_list = []
    s_date_list = []
    bperp_list = []

    inps.file = ut.get_file_list(inps.file, abspath=True)
    if os.path.splitext(inps.file[0])[1] not in ['.h5','.he5']:
        ifgramNum = len(inps.file)
        print('Number of interferograms: %d' % (ifgramNum))
        for fname in inps.file:
            try:    date12 = str(re.findall('\d{8}[-_]\d{8}', os.path.basename(fname))[0]).replace('_','-')
            except: date12 = str(re.findall('\d{6}[-_]\d{6}', os.path.basename(fname))[0]).replace('_','-')
            m_date, s_date = date12.split('-')
            bperp = readfile.read_attribute(fname)['P_BASELINE_TOP_HDR']
            m_date_list.append(m_date)
            s_date_list.append(s_date)
            bperp_list.append(bperp)

    else:
        h5 = h5py.File(inps.file[0],'r')
        ifgram_list = ut.check_drop_ifgram(h5)
        date12_list = ptime.list_ifgram2date12(ifgram_list)
        m_date_list = [date12.split('-')[0] for date12 in date12_list]
        s_date_list = [date12.split('-')[1] for date12 in date12_list]
        for ifgram in ifgram_list:
           bperp = h5[k][ifgram].attrs['P_BASELINE_TOP_HDR']
           bperp_list.append(bperp)
        ifgramNum = len(ifgram_list)

    fout = '{0} {1}     {2:<15}   {3}\n'
    fl = open(inps.outfile, 'w')
    for i in range(ifgramNum):
        fl.write(fout.format(m_date_list[i], s_date_list[i], bperp_list[i], inps.sensor))
    fl.close()
    print('write to %s' % (inps.outfile))
    return


###################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])


