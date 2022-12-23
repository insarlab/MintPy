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
    """Mathmathic operation of file

    Parameters: fname     - str, path to the input file
                operator  - str,   math operator, e.g., + - * / etc.
                operand   - float, math operand
                out_file  - str, path to the output file
    Returns:    out_file  - str, path to the output file
    """

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
        fbase, fext = os.path.splitext(fname)
        out_file = f'{fbase}_{suffix}{operand}{fext}'

    # basic info
    atr = readfile.read_attribute(fname)
    print(f'input is {atr["PROCESSOR"]} {atr["FILE_TYPE"]} file: {fname}')
    print(f'operation: file {operator} {operand:f}')

    ds_dict = {}
    ds_names = readfile.get_dataset_list(fname)
    for ds_name in ds_names:
        # read
        data = readfile.read(fname, datasetName=ds_name)[0]

        # apply math operation
        data = data_operation(data, operator, operand)

        # save
        ds_dict[ds_name] = data

    # write
    writefile.write(ds_dict, out_file=out_file, metadata=atr, ref_file=fname)

    return out_file
