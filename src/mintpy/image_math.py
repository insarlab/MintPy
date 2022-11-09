############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2015                               #
############################################################


import os

import numpy as np

from mintpy.utils import readfile, writefile


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
    print(f'operation: file {operator} {operand:f}')

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
