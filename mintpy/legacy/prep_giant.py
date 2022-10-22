#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Jul 2018                           #
############################################################


import argparse
import os
import sys

from lxml import objectify

from mintpy.objects import sensor
from mintpy.utils import readfile, utils as ut

key_giant2mintpy = {
    'xmin':'SUBSET_XMIN', 'xmax':'SUBSET_XMAX',
    'ymin':'SUBSET_YMIN', 'ymax':'SUBSET_YMAX',
}


##################################################################################################
EXAMPLE = """example:
  prep_giant.py  LS-PARAMS.h5
  prep_giant.py  TS-PARAMS.h5
  prep_giant.py  NSBAS-PARAMS.h5
  prep_giant.py  RAW-STACK.h5
  prep_giant.py  PROC-STACK.h5
  prep_giant.py  LS-PARAMS.h5 -x ../data.xml ../sbas.xml ../mints.xml
  prep_giant.py  LS-PARAMS.h5 -x ../data.xml ../sbas.xml ../mints.xml ../filt_fine.unw.rsc
"""


def create_parser():
    parser = argparse.ArgumentParser(description='Prepare attributes for GIAnT timeseries file.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', help='GIAnT timeseries file')
    parser.add_argument('-x','--xml', nargs='+', dest='xml_file',
                        help='XML file with data setting info.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    if not inps.xml_file:
        inps.xml_file = auto_xml_file4giant(inps.file)
    if not inps.xml_file:
        parser.print_usage()
        raise SystemExit('ERROR: no xml file found.')

    return inps


def auto_xml_file4giant(fname):
    file_list = [os.path.join(os.path.dirname(fname), f'../{i}')
                 for i in ['data.xml',
                           'sbas.xml',
                           'mints.xml',
                           'filt_fine.unw.rsc']]
    file_list = [i for i in file_list if os.path.isfile(i)]
    return file_list


def read_giant_xml(fname):
    odict = {}
    root = objectify.parse(fname).getroot()

    if root.find('master') is not None:
        comp = root['master']
        for key in ['wavelength', 'incidence']:
            odict[key] = comp[key].value

    if root.find('subimage') is not None:
        comp = root['subimage']
        for key in ['width', 'length',
                    'xmin', 'xmax',
                    'ymin', 'ymax',
                    'rxmin', 'rxmax',
                    'rymin', 'rymax']:
            odict[key] = comp[key].value

        odict = readfile.standardize_metadata(odict, standardKeys=key_giant2mintpy)
        odict['REF_Y'] = int((int(odict['rymin']) +
                              int(odict['rymax'])) / 2. + 0.5)
        odict['REF_X'] = int((int(odict['rxmin']) +
                              int(odict['rxmax'])) / 2. + 0.5)

    if root.find('proc/masterdate') is not None:
        odict['REF_DATE'] = root['proc']['masterdate'].value
    return odict


def prepare_metadata4giant(fname, meta_files=None):
    """Extract metadata from xml files for GIAnT time-series file."""
    # check xml files
    if not meta_files:
        meta_files = auto_xml_file4giant(fname)
    if not meta_files:
        raise FileNotFoundError("no xml file found.")

    # extract metadata from xml files
    rsc_files = [i for i in meta_files if i.endswith('.rsc')]
    xml_files = [i for i in meta_files if i.endswith('.xml')]
    xml_dict = {}
    for rsc_file in rsc_files:
        print(f'reading {rsc_file}')
        rsc_dict = readfile.read_roipac_rsc(rsc_file)
        for key in ['length', 'LENGTH', 'FILE_LENGTH', 'width', 'WIDTH']:
            if key in rsc_dict.keys():
                rsc_dict.pop(key)
        xml_dict.update(rsc_dict)
    for xml_file in xml_files:
        print(f'reading {xml_file}')
        xml_dict.update(read_giant_xml(xml_file))

    if not xml_dict:
        raise ValueError('No metadata found in file: '+xml_file)

    # standardize metadata names
    xml_dict = readfile.standardize_metadata(xml_dict)

    # project name
    sensor_name, project_name = sensor.project_name2sensor_name(os.path.abspath(fname))
    if sensor_name:
        xml_dict['PLATFORM'] = sensor_name
    if project_name:
        xml_dict['PROJECT_NAME'] = project_name
        if sensor_name in project_name:
            tmp = project_name.split(sensor_name)[1][0]
            if tmp == 'A':
                xml_dict['ORBIT_DIRECTION'] = 'ASCENDING'
            else:
                xml_dict['ORBIT_DIRECTION'] = 'DESCENDING'

    # update GIAnT HDF5 file
    fname = ut.add_attribute(fname, xml_dict, print_msg=True)
    return fname


##################################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    prepare_metadata4giant(inps.file, inps.xml_file)
    return


###################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
