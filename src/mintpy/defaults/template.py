"""Utilities to grab template content for processing steps."""
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Mar 2020                           #
############################################################
# Recommend usage:
#   from mintpy.defaults.template import STEP_LIST, get_template_content


import os
import re

STEP_LIST4OFFSET = [
    'load_data',
    'modify_network',
    'invert_network',
    'deramp',
    'velocity',
    'geocode',
    'google_earth',
]

STEP_LIST = [
    'load_data',
    'modify_network',
    'reference_point',
    'quick_overview',
    'correct_unwrap_error',
    'invert_network',
    'correct_LOD',
    'correct_SET',
    'correct_troposphere',
    'deramp',
    'correct_topography',
    'residual_RMS',
    'reference_date',
    'velocity',
    'geocode',
    'google_earth',
    'hdfeos5',
]


def get_template_content(step_name, template_file=None, indentation=2, header_footer=True):
    """Grab the related template content for each step
    To avoid duplication in each utility script.

    Parameters: step_name     - str, step name
                template_file - str, path of the template file
    Returns:    step_content  - str, comments and options of the step
    """
    import mintpy

    # check
    if step_name not in STEP_LIST:
        raise ValueError(f'input step name "{step_name}" not found! STEP_LIST={STEP_LIST}')

    # read template file into a list of strings
    if template_file is None:
        template_file = os.path.join(os.path.dirname(mintpy.__file__), 'defaults/smallbaselineApp.cfg')
    lines = open(template_file).readlines()
    lines = [line.strip(' ') for line in lines]

    # get starting line index
    pattern = r"^##########[ ]?\d{0,2}\.*\d{0,2} " + step_name
    inds = [lines.index(line) for line in lines if re.match(pattern, line)]
    if len(inds) > 0:
        ind0 = inds[0] + 1
    else:
        raise ValueError(f'pattern "{pattern}" is not found in file: {template_file}!')

    # get ending line index: the next line with ten of #
    ind1 = -1
    for i in range(ind0+1, len(lines)):
        if lines[i].startswith('########## '):
            ind1 = i
            break

    # merge the related list of strings into one string
    step_content = ''.join(' '*indentation + line for line in lines[ind0:ind1])
    step_content = step_content.rstrip()

    if header_footer:
        step_content = 'template options:\n' + step_content + '\n'
    return step_content
