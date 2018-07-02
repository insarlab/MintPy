#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2018, Zhang Yunjun                          #
# Author:  Zhang Yunjun, Jul 2018                          #
############################################################


import os
import argparse
import h5py
from lxml import objectify
from pysar.utils import readfile


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
        file_dir = os.path.dirname(inps.file)
        file_list = ['data.xml', 'sbas.xml', 'mints.xml']
        file_list = [os.path.join(file_dir, '../{}'.format(i)) for i in file_list]
        inps.xml_file = [i for i in file_list if os.path.isfile(i)]
    if not inps.xml_file:
        print('error: no xml file found.')
        parser.print_usage()
        raise SystemExit()

    return inps

def prepare_metadata4giant(fname, xml_files):
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

    xml_dict = readfile.standardize_metadata(xml_dict)

    print('open {} with +r mode'.format(fname))
    f = h5py.File(fname, 'r+')
    for key, value in xml_dict.items():
        f.attrs[key] = str(value)
    f.close()
    print('finished updating {}'.format(fname))
    return fname


##################################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    prepare_metadata4giant(inps.file, inps.xml_file)
    return


###################################################################################################
if __name__ == '__main__':
    main()

