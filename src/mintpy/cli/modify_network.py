#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Zhang Yunjun, Aug 2022        #
############################################################


import os
import sys

from mintpy.defaults.template import get_template_content
from mintpy.utils.arg_utils import create_argument_parser

###############################  Usage  ################################
TEMPLATE = get_template_content('modify_network')

REFERENCE = """reference:
  Yunjun, Z., Fattahi, H. and Amelung, F. (2019), Small baseline InSAR time series analysis:
    Unwrapping error correction and noise reduction, Computers & Geosciences, 133, 104331,
    doi:10.1016/j.cageo.2019.104331.
  Chaussard, E., Bürgmann, R., Fattahi, H., Nadeau, R. M., Taira, T., Johnson, C. W. and Johanson, I.
    (2015), Potential for larger earthquakes in the East San Francisco Bay Area due to the direct
    connection between the Hayward and Calaveras Faults, Geophysical Research Letters, 42(8),
    2734-2741, doi:10.1002/2015GL063575.
  Kang, Y., Lu, Z., Zhao, C., Xu, Y., Kim, J. W., & Gallegos, A. J. (2021).InSAR monitoring
    of creeping landslides in mountainous regions: A case study in Eldorado National Forest,
    California. Remote Sensing of Environment, 258, 112400. doi:10.1016/j.rse.2021.112400
"""

EXAMPLE = """example:
  modify_network.py inputs/ifgramStack.h5 -t smallbaselineApp.cfg
  modify_network.py inputs/ifgramStack.h5 --reset
  modify_network.py inputs/ifgramStack.h5 --manual
"""


def create_parser(subparsers=None):
    synopsis = 'Modify the network of interferograms'
    epilog = REFERENCE + '\n' + TEMPLATE + '\n' + EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', help='Files to modify/drop network, e.g. inputs/ifgramStack.h5.')
    parser.add_argument('-t', '--template', dest='template_file',
                        help='Template file with input options')
    parser.add_argument('--reset', action='store_true',
                        help='restore all interferograms in the file, by marking all dropIfgram=True')
    parser.add_argument('--noaux', dest='update_aux', action='store_false',
                        help='Do not update auxiliary files, e.g.\n' +
                             'maskConnComp.h5 or avgSpatialCoh.h5 from ifgramStack.h5')

    # 1. temp/perp baseline, num of conn., dates, pair index, etc.
    parser.add_argument('--max-tbase', dest='tempBaseMax',
                        type=float, help='max temporal baseline in days')
    parser.add_argument('--max-pbase', dest='perpBaseMax',
                        type=float, help='max perpendicular baseline in meters')
    parser.add_argument('--max-conn-num', dest='connNumMax', type=int,
                        help='max number of connections/neighbors per acquisition')
    parser.add_argument('-r', '--reference', dest='referenceFile',
                        help='Reference hdf5 / list file with network information.\n'
                             'i.e. ifgramStack.h5, date12_list.txt')
    parser.add_argument('--exclude-ifg-index', dest='excludeIfgIndex', nargs='*',
                        help='index of interferograms to remove/drop.\n1 as the first')
    parser.add_argument('--exclude-date', dest='excludeDate', nargs='*',
                        help='date(s) to remove/drop, all interferograms included date(s) will be removed')
    parser.add_argument('--start-date', '--min-date', dest='startDate',
                        help='remove/drop interferograms with date earlier than start-date in YYMMDD or YYYYMMDD format')
    parser.add_argument('--end-date', '--max-date', dest='endDate',
                        help='remove/drop interferograms with date later than end-date in YYMMDD or YYYYMMDD format')

    # 2. coherence-based network
    cohBased = parser.add_argument_group('Data-driven network modification', 'Drop/modify network based on data')
    # 2.1 coherence-based
    cohBased.add_argument('--coherence-based', dest='coherenceBased', action='store_true',
                          help='Enable coherence-based network modification (default: %(default)s).')
    cohBased.add_argument('--min-coherence', dest='minCoherence', type=float, default=0.7,
                          help='Minimum coherence value (default: %(default)s).')
    # 2.2 area-ratio-based
    cohBased.add_argument('--area-ratio-based', dest='areaRatioBased', action='store_true',
                          help='Enable area ratio-based network modification (default: %(default)s).')
    cohBased.add_argument('--min-area-ratio', dest='minAreaRatio', type=float, default=0.75,
                          help='Minimum area ratio value (default: %(default)s).')
    # common parameters
    cohBased.add_argument('--no-mst', dest='keepMinSpanTree', action='store_false',
                          help='Do not keep interferograms in Min Span Tree network based on inversed mean coherene')
    cohBased.add_argument('--mask', dest='maskFile', default='waterMask.h5',
                          help='Mask file used to calculate the spatial coherence '
                               '(default: waterMask.h5 or None)')
    cohBased.add_argument('--aoi-yx', dest='aoiYX', type=int, nargs=4, metavar=('X0', 'Y0', 'X1', 'Y1'), default=None,
                          help='AOI in row/column range for coherence calculation (default: %(default)s).')
    cohBased.add_argument('--aoi-lalo', dest='aoiLALO', type=float, nargs=4, metavar=('W', 'S', 'E', 'N'), default=None,
                          help='AOI in lat/lon range for coherence calculation (default: %(default)s).')
    cohBased.add_argument('--lookup', dest='lookupFile',
                          help='Lookup table/mapping transformation file for geo/radar coordinate conversion.\n' +
                               'Needed for mask AOI in lalo')

    # 3. manual selection
    manual = parser.add_argument_group('Manual Network', 'Manually select/drop/modify network')
    manual.add_argument('--manual', action='store_true',
                        help='display network to manually choose line/interferogram to remove')
    return parser


def cmd_line_parse(iargs=None):
    # parse
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # import
    from mintpy.utils import readfile, utils as ut

    # check: --mask option
    if not os.path.isfile(inps.maskFile):
        inps.maskFile = None

    # check: --exclude-ifg-index option (convert input index to continuous index list)
    inps.excludeIfgIndex = read_input_index_list(inps.excludeIfgIndex, stackFile=inps.file)

    # check: -t / --template option
    if inps.template_file:
        inps = read_template2inps(inps.template_file, inps)

    # check: input arguments (required at least one)
    required_args = [
        inps.referenceFile, inps.tempBaseMax, inps.perpBaseMax, inps.connNumMax,
        inps.excludeIfgIndex, inps.excludeDate, inps.coherenceBased, inps.areaRatioBased,
        inps.startDate, inps.endDate, inps.reset, inps.manual,
    ]
    if all(not i for i in required_args + [inps.template_file]):
        msg = 'No input option found to remove interferogram, exit.\n'
        msg += 'To manually modify network, please use --manual option '
        raise Exception(msg)

    # default: --lookup option
    if not inps.lookupFile:
        inps.lookupFile = ut.get_lookup_file()

    # check: --aoi-lalo option (not for radar-coded products without lookup files)
    if inps.aoiLALO:
        atr = readfile.read_attribute(inps.file)
        if not inps.lookupFile and 'Y_FIRST' not in atr.keys():
            msg = 'WARNING: Can NOT use --aoi-lalo option for files in radar coordinates '
            msg += 'without lookup file. Ignore this option and continue.'
            print(msg)
            inps.aoiLALO = None

    # default: turn --reset ON if:
    # 1) no input options found to drop ifgram AND
    # 2) there is template input
    if inps.template_file and all(not i for i in required_args):
        print('No input option found to remove interferogram')
        print('Keep all interferograms by enable --reset option')
        inps.reset = True

    return inps


def read_template2inps(template_file, inps):
    """Read input template options into Namespace inps"""
    print('read options from template file: '+os.path.basename(template_file))

    from mintpy.utils import ptime, readfile, utils as ut

    iDict = vars(inps)
    template = readfile.read_template(inps.template_file, skip_chars=['[', ']'])
    template = ut.check_template_auto_value(template)

    # Update inps if key existed in template file
    prefix = 'mintpy.network.'
    key_list = [i for i in list(iDict.keys()) if prefix+i in template.keys()]
    for key in key_list:
        value = template[prefix+key]
        if key in ['coherenceBased', 'areaRatioBased', 'keepMinSpanTree']:
            iDict[key] = value

        elif value:
            if key in ['minCoherence', 'minAreaRatio', 'tempBaseMax', 'perpBaseMax']:
                iDict[key] = float(value)
            elif key in ['connNumMax']:
                iDict[key] = int(value)
            elif key in ['maskFile', 'referenceFile']:
                iDict[key] = value

            elif key == 'aoiYX':
                tmp = [i.strip() for i in value.split(',')]
                sub_y = sorted(int(i.strip()) for i in tmp[0].split(':'))
                sub_x = sorted(int(i.strip()) for i in tmp[1].split(':'))
                inps.aoiYX = (sub_x[0], sub_y[0], sub_x[1], sub_y[1])

            elif key == 'aoiLALO':
                tmp = [i.strip() for i in value.split(',')]
                sub_lat = sorted(float(i.strip()) for i in tmp[0].split(':'))
                sub_lon = sorted(float(i.strip()) for i in tmp[1].split(':'))
                inps.aoiLALO = (sub_lon[0], sub_lat[1], sub_lon[1], sub_lat[0])

            elif key in ['startDate', 'endDate']:
                iDict[key] = ptime.yyyymmdd(value)
            elif key == 'excludeDate':
                iDict[key] = ptime.yyyymmdd(value.split(','))
            elif key == 'excludeIfgIndex':
                iDict[key] += value.split(',')
                iDict[key] = read_input_index_list(iDict[key], stackFile=inps.file)

    return inps


def read_input_index_list(idxList, stackFile=None):
    """Read ['2','3:5','10'] into ['2','3','4','5','10']"""
    from mintpy.objects import ifgramStack

    if not idxList:
        return []

    idxListOut = []
    for idx in idxList:
        c = sorted(int(i) for i in idx.split(':'))
        if len(c) == 2:
            idxListOut += list(range(c[0], c[1]+1))
        elif len(c) == 1:
            idxListOut.append(c[0])
        else:
            print('Unrecoganized input: '+idx)
    idxListOut = sorted(set(idxListOut))

    # remove index not existing in the input ifgram stack file
    if stackFile:
        obj = ifgramStack(stackFile)
        obj.open(print_msg=False)
        idxListOut = [i for i in idxListOut if i < obj.numIfgram]

    return idxListOut


#########################  Main Function  ##############################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.modify_network import modify_network

    # run
    modify_network(inps)


########################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
