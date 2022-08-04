#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2015                               #
############################################################


import os
import sys
import numpy as np
from mintpy.utils import readfile, writefile
from mintpy.utils.arg_utils import create_argument_parser


#######################################################################################
EXAMPLE = """example:
  image_math.py  velocity.h5            '+'  0.5
  image_math.py  geo_080212_101120.cor  '-'  0.2
  image_math.py  timeseries.h5          '*'  1.5
  image_math.py  velocity.h5            '/'  2.0
  image_math.py  velocity.h5            '^'  2.0
"""


def create_parser(subparsers=None):
    synopsis = 'Basic Mathmatic Operation of file'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('file', help='input file')
    parser.add_argument('-o', '--output', dest='outfile',
                        help='output file name.')
    parser.add_argument('operator', choices=[
                        '+', '-', '*', '/', '^'], help='mathmatical operator')
    parser.add_argument('operand', metavar='VALUE', type=float,
                        help='value to be operated with input file')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


#######################################################################################
def data_operation(data, operator, operand):
    """Mathmatic operation of 2D matrix"""
    if operator == '+':
        data_out = data + operand
    elif operator == '-':
        data_out = data - operand
    elif operator == '*':
        data_out = data * operand
    elif operator == '/':
        data_out = data * (1.0/operand)
    elif operator == '^':
        data_out = data**operand
    data_out = np.array(data_out, dtype=np.float32)
    return data_out


def file_operation(fname, operator, operand, out_file=None):
    """Mathmathic operation of file"""

    # Basic Info
    atr = readfile.read_attribute(fname)
    k = atr['FILE_TYPE']
    print('input is '+k+' file: '+fname)
    print('operation: file %s %f' % (operator, operand))

    # default output filename
    if not out_file:
        if operator in ['+', 'plus',  'add',      'addition']:
            suffix = 'plus'
        elif operator in ['-', 'minus', 'substract', 'substraction']:
            suffix = 'minus'
        elif operator in ['*', 'times', 'multiply', 'multiplication']:
            suffix = 'multiply'
        elif operator in ['/', 'obelus', 'divide',   'division']:
            suffix = 'divide'
        elif operator in ['^', 'pow', 'power']:
            suffix = 'pow'
        out_file = '{}_{}{}{}'.format(os.path.splitext(fname)[0], suffix,
                                      str(operand), os.path.splitext(fname)[1])

    atr = readfile.read_attribute(fname)
    dsNames = readfile.get_dataset_list(fname)
    dsDict = {}
    for dsName in dsNames:
        data = readfile.read(fname, datasetName=dsName)[0]
        data = data_operation(data, operator, operand)
        dsDict[dsName] = data
    writefile.write(dsDict, out_file=out_file, metadata=atr, ref_file=fname)
    return out_file


#######################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    file_operation(inps.file, inps.operator, inps.operand, inps.outfile)

    print('Done.')
    return


#######################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
