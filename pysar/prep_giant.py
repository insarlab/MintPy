#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2018, Zhang Yunjun                          #
# Author:  Zhang Yunjun, Jul 2018                          #
############################################################


import os
import argparse
from lxml import objectify
from pysar.utils import readfile, utils as ut
from pysar.objects import sensor


key_giant2pysar = {'xmin':'SUBSET_XMIN', 'xmax':'SUBSET_XMAX',
                   'ymin':'SUBSET_YMIN', 'ymax':'SUBSET_YMAX',
                  }


##################################################################################################
EXAMPLE = """example:
  prep_giant.py  LS-PARAMS.h5 -x ../data.xml
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
        print('error: no xml file found.')
        parser.print_usage()
        raise SystemExit()

    return inps


def auto_xml_file4giant(fname):
    file_list = [os.path.join(os.path.dirname(fname), '../{}'.format(i))
                 for i in ['data.xml',
                           'sbas.xml',
                           'mints.xml']]
    file_list = [i for i in file_list if os.path.isfile(i)]
    return file_list


def prepare_metadata4giant(fname, xml_files=None):
    """Extract metadata from xml files for GIAnT time-series file."""
    # check xml files
    if not xml_files:
        xml_files = auto_xml_file4giant(fname)
    if not xml_files:
        raise FileNotFoundError("no xml file found.")

    # extract metadata from xml files
    xml_dict = {}
    for xml_file in xml_files:
        print('reading {}'.format(xml_file))
        root = objectify.parse(xml_file).getroot()

        if root.find('master') is not None:
            comp = root['master']
            for key in ['wavelength', 'incidence']:
                xml_dict[key] = comp[key].value

        if root.find('subimage') is not None:
            comp = root['subimage']
            for key in ['width', 'length',
                        'xmin', 'xmax',
                        'ymin', 'ymax',
                        'rxmin', 'rxmax',
                        'rymin', 'rymax']:
                xml_dict[key] = comp[key].value

            xml_dict = readfile.standardize_metadata(xml_dict, standardKeys=key_giant2pysar)
            xml_dict['REF_Y'] = int((int(xml_dict['rymin']) +
                                     int(xml_dict['rymax'])) / 2. + 0.5)
            xml_dict['REF_X'] = int((int(xml_dict['rxmin']) +
                                     int(xml_dict['rxmax'])) / 2. + 0.5)

        if root.find('proc/masterdate') is not None:
            xml_dict['REF_DATE'] = root['proc']['masterdate'].value

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
    main()

