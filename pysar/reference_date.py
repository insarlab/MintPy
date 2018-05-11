#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2013, Zhang Yunjun, Heresh Fattahi          #
# Author:  Zhang Yunjun, Heresh Fattahi                    #
############################################################

import os
import shutil
import argparse
import numpy as np
from pysar.utils import readfile, ptime, utils as ut
from pysar.objects import timeseries


##################################################################
TEMPLATE = """
## 4.1 Phase Residual Root Mean Square
## calculate the deramped Root Mean Square (RMS) for each epoch of timeseries residual from DEM error inversion
## To get rid of long wavelength component in space, a ramp is removed for each epoch.
## Recommendation: quadratic for whole image, plane for local/small area
pysar.residualRms.maskFile        = auto  #[filename / no], auto for maskTempCoh.h5, mask for ramp estimation
pysar.residualRms.ramp            = auto  #[quadratic / plane / no], auto for quadratic
pysar.residualRms.threshold       = auto  #[0.0-inf], auto for 0.02, minimum RMS in meter for exclude date(s)

## 4.2 Select Reference Date
## reference all timeseries to one date in time
## minRMS - choose date with minimum residual RMS using value from step 8.1
## no     - do not change the default reference date (1st date)
pysar.reference.date = auto   #[reference_date.txt / 20090214 / minRMS / no], auto for minRMS
"""

EXAMPLE = """example:
  reference_date.py timeseries_ECMWF_demErr.h5  --ref-date 20050107
  reference_date.py timeseries_ECMWF_demErr.h5  --ref-date minRMS
  reference_date.py timeseries_ECMWF_demErr.h5  --template KujuAlosAT422F650.template
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Change reference date of timeseries.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('timeseries_file', help='Timeseries File')
    parser.add_argument('-r', '--ref-date', dest='refDate', default='minRMS',
                        help='reference date or method, default: auto. e.g.\n' +
                             '20101120\n' +
                             'reference_date.txt - text file with date in YYYYMMDD format in it\n' +
                             'minRMS             - choose date with min residual standard deviation')
    parser.add_argument('-t', '--template', dest='template_file',
                        help='template file with options below:\n' + TEMPLATE + '\n')

    auto = parser.add_argument_group('Auto referencing in time based on max phase residual coherence')
    auto.add_argument('--residual-file', dest='resid_file', default='timeseriesResidual.h5',
                      help='timeseries of phase residual file from DEM error inversion.\n' +
                           'Deramped Residual RMS')
    auto.add_argument('--deramp', dest='ramp_type', default='quadratic',
                      help='ramp type to remove for each epoch from phase residual\n' +
                           'default: quadratic\n' +
                           'no - do not remove ramp')
    auto.add_argument('--mask', dest='maskFile', default='maskTempCoh.h5',
                      help='mask file used for ramp estimation\n' + 'default: maskTempCoh.h5')
    parser.add_argument('-o', '--outfile', help='Output file name.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


def read_template2inps(templateFile, inps=None):
    """Update inps with options from templateFile"""
    if not inps:
        inps = cmd_line_parse()
    template = readfile.read_template(templateFile)
    template = ut.check_template_auto_value(template)

    key = 'pysar.reference.date'
    if key in template.keys() and template[key]:
        inps.refDate = template[key]

    key = 'pysar.residualStd.maskFile'
    if key in template.keys() and template[key]:
        inps.maskFile = template[key]

    key = 'pysar.residualStd.ramp'
    if key in template.keys() and template[key]:
        inps.ramp_type = template[key]
    return inps


def read_ref_date(inps):
    if not inps.refDate:
        print('No reference date input, skip this step.')
        return inps.timeseries_file

    elif inps.refDate.lower() == 'minrms':
        print('-'*50)
        print('auto choose reference date based on minimum residual RMS')
        rms_list, date_list = ut.get_residual_rms(inps.resid_file,
                                                  inps.maskFile,
                                                  inps.ramp_type)[0:2]
        ref_idx = np.argmin(rms_list)
        inps.refDate = date_list[ref_idx]
        print('date with minimum residual RMS: %s - %.4f' % (inps.refDate, rms_list[ref_idx]))
        print('-'*50)

    elif os.path.isfile(inps.refDate):
        print('read reference date from file: ' + inps.refDate)
        inps.refDate = ptime.read_date_list(inps.refDate)[0]
    return inps.refDate


##################################################################
def ref_date_file(inFile, refDate, outFile=None):
    """Change input file reference date to a different one."""
    if not outFile:
        outFile = os.path.splitext(inFile)[0] + '_refDate.h5'
    refDate = ptime.yyyymmdd(refDate)
    print('input reference date: ' + refDate)

    # Input file info
    obj = timeseries(inFile)
    obj.open()
    atr = dict(obj.metadata)
    if atr['FILE_TYPE'] != 'timeseries':
        print('ERROR: input file is {}, only timeseries is supported.'.format(atr['FILE_TYPE']))
        return None
    if refDate not in obj.dateList:
        print('ERROR: Input reference date was not found!\nAll dates available: {}'.format(obj.dateList))
        return None
    if refDate == atr['REF_DATE']:
        print('Same reference date chosen as existing reference date.')
        print('Copy {} to {}'.format(inFile, outFile))
        shutil.copy2(inFile, outFile)
        return outFile

    # Referencing in time
    data = obj.read()
    data -= np.tile(data[obj.dateList.index(refDate), :, :].reshape(1, data.shape[1], data.shape[2]),
                    (data.shape[0], 1, 1))
    atr['REF_DATE'] = refDate

    if not outFile:
        outFile = '{}_refDate.h5'.format(os.path.splitext(inFile)[0])
    outObj = timeseries(outFile)
    outObj.write2hdf5(data, refFile=inFile, metadata=atr)
    return outFile


##################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    if inps.template_file:
        inps = read_template2inps(inps.template_file)

    inps.refDate = read_ref_date(inps)

    inps.outfile = ref_date_file(inps.timeseries_file,
                                 inps.refDate,
                                 inps.outfile)
    return inps.outfile


##################################################################
if __name__ == '__main__':
    main()
