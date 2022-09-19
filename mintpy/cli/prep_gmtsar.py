############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Aug 2022                      #
############################################################


import os
import sys
import glob
from mintpy.utils.arg_utils import create_argument_parser


#########################################################################
EXAMPLE = """example:
  prep_gmtsar.py StHelensEnvDT156.txt
"""

def create_parser(subparsers=None):
    """Command line parser."""
    synopsis = 'Prepare GMTSAR metadata files.'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('template_file', type=str, help='MintPy template file for GMTSAR products.')
    parser.add_argument('--mintpy-dir', dest='mintpy_dir', default='./',
                        help='MintPy directory (default: %(default)s).')
    parser.add_argument('--force', dest='update_mode', action='store_false',
                        help='Force to overwrite all .rsc metadata files.')
    return parser


def cmd_line_parse(iargs = None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    inps.template_file = os.path.abspath(inps.template_file)
    inps.mintpy_dir = os.path.expanduser(inps.mintpy_dir)
    inps.mintpy_dir = os.path.abspath(inps.mintpy_dir)
    return inps


#########################################################################
def main(iargs=None):
    from mintpy.utils import readfile
    from mintpy.prep_gmtsar import extract_gmtsar_metadata, prepare_geometry, prepare_stack

    inps = cmd_line_parse(iargs)

    # read file path from template file
    template = readfile.read_template(inps.template_file)
    inps.unw_files = sorted(glob.glob(template['mintpy.load.unwFile']))
    inps.cor_files = sorted(glob.glob(template['mintpy.load.corFile']))
    inps.dem_file = glob.glob(template['mintpy.load.demFile'])[0]

    # extract common metadata
    rsc_file = os.path.join(inps.mintpy_dir, 'inputs/data.rsc')
    meta = extract_gmtsar_metadata(unw_file=inps.unw_files[0],
                                   template_file=inps.template_file,
                                   rsc_file=rsc_file,
                                   update_mode=inps.update_mode)

    # prepare metadata for geometry files
    prepare_geometry([inps.dem_file], meta=meta, update_mode=inps.update_mode)

    # prepare metadata for interferogram files
    prepare_stack(inps.unw_files, meta=meta, update_mode=inps.update_mode)

    print('Done.')


#########################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
