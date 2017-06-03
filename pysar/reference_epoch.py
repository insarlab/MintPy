#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################

import sys
import os
import argparse

import h5py
import numpy as np
import matplotlib as mpl; mpl.use('Agg')
import matplotlib.pyplot as plt

import pysar._readfile as readfile
import pysar._datetime as ptime
import pysar._network as pnet
import pysar._pysar_utilities as ut


##################################################################
def ref_date_attribute(atr_in, ref_date, date_list):
    '''Update attribute dictionary for reference date'''
    ref_date = ptime.yyyymmdd(ref_date)
    date_list = ptime.yyyymmdd(date_list)
    ref_index = date_list.index(ref_date)

    atr = dict()
    for key, value in atr_in.iteritems():
        atr[key] = str(value)

    # Update ref_date
    atr['ref_date'] = ref_date
    print 'update ref_date'

    # Update Bperp time series
    try:
        pbase = np.array([float(i) for i in atr['P_BASELINE_TIMESERIES'].split()])
        pbase -= pbase[ref_index]
        atr['P_BASELINE_TIMESERIES'] = str(pbase.tolist()).translate(None,'[],')
        print 'update P_BASELINE_TIMESERIES'
    except:
        pass

    try:
        pbase_top    = np.array([float(i) for i in atr['P_BASELINE_TOP_TIMESERIES'].split()])
        pbase_bottom = np.array([float(i) for i in atr['P_BASELINE_BOTTOM_TIMESERIES'].split()])
        pbase_top    -= pbase_top[ref_index]
        pbase_bottom -= pbase_bottom[ref_index]
        atr['P_BASELINE_TOP_TIMESERIES']    = str(pbase_top.tolist()).translate(None,'[],')
        atr['P_BASELINE_BOTTOM_TIMESERIES'] = str(pbase_bottom.tolist()).translate(None,'[],')
        print 'update P_BASELINE_TOP/BOTTOM_TIMESERIES'
    except:
        pass

    return atr


def ref_date_file(inFile, ref_date, outFile=None):
    '''Change input file reference date to a different one.'''
    if not outFile:
        outFile = os.path.splitext(inFile)[0]+'_refDate.h5'

    # Input file type 
    atr = readfile.read_attribute(inFile)
    k = atr['FILE_TYPE']
    if not k in ['timeseries']:
        print 'Input file is '+k+', only timeseries is supported.'
        return None

    # Input reference date
    h5 = h5py.File(inFile, 'r')
    date_list = sorted(h5[k].keys())
    date_num = len(date_list)

    ref_date = ptime.yyyymmdd(ref_date)
    print 'input reference date: '+ref_date
    if not ref_date in date_list:
        print 'Input reference date was not found!\nAll dates available: '+str(date_list)
        return None

    # Referencing in time
    ref_data = h5[k].get(ref_date)[:]

    print 'writing >>> '+outFile
    h5out = h5py.File(outFile,'w')
    group = h5out.create_group(k)
    prog_bar = ptime.progress_bar(maxValue=date_num)
    for i in range(date_num):
        date = date_list[i]
        data = h5[k].get(date)[:]
        dset = group.create_dataset(date, data=data-ref_data, compression='gzip')
        prog_bar.update(i+1, suffix=date)
    prog_bar.close()
    h5.close()

    ## Update attributes
    atr = ref_date_attribute(atr, ref_date, date_list)
    for key,value in atr.iteritems():
        group.attrs[key] = value
    h5out.close()

    return outFile


def read_template2inps(templateFile, inps=None):
    '''Update inps with options from templateFile'''
    if not inps:
        inps = cmdLineParse()

    template = readfile.read_template(templateFile)
    key_list = template.keys()

    key = 'pysar.reference.date'
    if key in key_list:
        inps.ref_date = template[key]

    prefix = 'pysar.residualStd.'
    key = prefix+'maskFile'
    if key in key_list:
        value = template[key]
        if value == 'auto':
            inps.mask_file = 'maskTempCoh.h5'
        elif value == 'no':
            inps.mask_file = None
        else:
            inps.mask_file = value

    key = prefix+'ramp'
    if key in key_list:
        value = template[key]
        if value == 'auto':
            inps.ramp_type = 'quadratic'
        else:
            inps.ramp_type = 'no'

    return inps


##################################################################
TEMPLATE='''
## 8.1 Residual Standard Deviation (RSD)
## calculate the deramped standard deviation (STD) for each epoch of timeseries residual from DEM error inversion
## To get rid of long wavelength component in space, a ramp is removed for each epoch.
pysar.residualStd.maskFile        = auto  #[file name / no], auto for maskTempCoh.h5, mask for ramp estimation
pysar.residualStd.ramp            = auto  #[quadratic / plane / no], auto for quadratic

## 9. Reference in Time
## reference all timeseries to one date in time
## auto - choose date with minimum residual STD using value from step 8.1
## no   - do not change reference date, keep the defaut one (1st date usually) and skip this step
pysar.reference.date = auto   #[auto / txtFile / 20090214 / no]
'''

EXAMPLE='''example:
  reference_epoch.py timeseries_ECMWF_demErr.h5  --ref-date 20050107
  reference_epoch.py timeseries_ECMWF_demErr.h5  --ref-date auto
  reference_epoch.py timeseries_ECMWF_demErr.h5  --template KujuAlosAT422F650.template
'''

def cmdLineParse():
    parser = argparse.ArgumentParser(description='Change reference date of timeseries.',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)

    parser.add_argument('timeseries_file', help='Timeseries File')
    parser.add_argument('--ref-date', dest='ref_date', default='auto',\
                        help='reference date or method, default: auto. e.g.\n'+\
                             '20101120\n'+\
                             'reference_date.txt - text file with date in YYYYMMDD format in it\n'+\
                             'auto               - choose date with min residual standard deviation')
    parser.add_argument('--template', dest='template_file',\
                        help='template file with options below:\n'+TEMPLATE+'\n')

    auto = parser.add_argument_group('Auto referencing in time based on max phase residual coherence')
    auto.add_argument('--residual-file', dest='resid_file',\
                      help='timeseries of phase residual file from DEM error inversion.\n'+\
                           'Deramped Residual Standard Deviation')
    auto.add_argument('--deramp', dest='ramp_type', default='quadratic',\
                      help='ramp type to remove for each epoch from phase residual\n'+\
                           'default: quadratic\n'+\
                           'no - do not remove ramp')
    auto.add_argument('--mask', dest='mask_file', default='maskTempCoh.h5',\
                      help='mask file used for ramp estimation\n'+'default: maskTempCoh.h5')
    parser.add_argument('-o','--outfile', help='Output file name.')

    inps = parser.parse_args()
    return inps


##################################################################
def main(argv):
    inps = cmdLineParse()
    if inps.template_file:
        inps = read_template2inps(inps.template_file)

    if inps.ref_date == 'no':
        print 'No reference date input, skip this step.'
        return inps.timeseries_file

    elif inps.ref_date.lower() in ['auto']:
        print '------------------------------------------------------------'
        print 'auto choose reference date based on minimum Residual Standard Deviation (RSD)'
        if not inps.resid_file:
            inps.resid_file = os.path.splitext(inps.timeseries_file)[0]+'InvResid.h5'
        std_list, date_list = ut.get_residual_std(inps.resid_file, inps.mask_file, inps.ramp_type)
        ref_idx = np.argmin(std_list)
        inps.ref_date = date_list[ref_idx]
        print 'date with minimum residual std: %s - %.4f' % (inps.ref_date, std_list[ref_idx])

        #txtFile = 'reference_date.txt'
        #f = open(txtFile, 'w')
        #f.write(inps.ref_date+'\n')
        #f.close()
        #print 'save date to file: '+txtFile
        print '------------------------------------------------------------'

    elif os.path.isfile(inps.ref_date):
        print 'read reference date from file: '+inps.ref_date
        inps.ref_date = ptime.read_date_list(inps.ref_date)[0]

    # Referencing input file
    inps.outfile = ref_date_file(inps.timeseries_file, inps.ref_date, inps.outfile)
    return inps.outfile


##################################################################
if __name__ == '__main__':
    main(sys.argv[1:])


