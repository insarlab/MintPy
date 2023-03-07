#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Zhang Yunjun, Aug 2022        #
############################################################


import datetime
import os
import sys

import mintpy
from mintpy.defaults.template import STEP_LIST
from mintpy.utils.arg_utils import create_argument_parser

##########################################################################
STEP_HELP = """Command line options for steps processing with names are chosen from the following list:

{}
{}
{}

In order to use either --start or --dostep, it is necessary that a
previous run was done using one of the steps options to process at least
through the step immediately preceding the starting step of the current run.
""".format(STEP_LIST[0:5], STEP_LIST[5:11], STEP_LIST[11:])

REFERENCE = """reference:
  Yunjun, Z., H. Fattahi, and F. Amelung (2019), Small baseline InSAR time series analysis:
    Unwrapping error correction and noise reduction, Computers & Geosciences, 133, 104331,
    doi:10.1016/j.cageo.2019.104331.
"""

EXAMPLE = """example:
  smallbaselineApp.py                         # run with default template 'smallbaselineApp.cfg'
  smallbaselineApp.py <custom_template>       # run with default and custom templates
  smallbaselineApp.py -h / --help             # help
  smallbaselineApp.py -H                      # print    default template options
  smallbaselineApp.py -g                      # generate default template if it does not exist
  smallbaselineApp.py -g <custom_template>    # generate/update default template based on custom template
  smallbaselineApp.py --plot                  # plot results w/o run [to populate the 'pic' folder after failed runs]

  # step processing with --start/stop/dostep options
  smallbaselineApp.py GalapagosSenDT128.txt --dostep velocity  # run step 'velocity' only
  smallbaselineApp.py GalapagosSenDT128.txt --end load_data    # end run after step 'load_data'
"""


def create_parser(subparsers=None):
    synopsis = 'Routine Time Series Analysis for Small Baseline InSAR Stack'
    epilog = REFERENCE + '\n' + EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('customTemplateFile', nargs='?',
                        help='custom template with option settings.\n' +
                             "ignored if the default smallbaselineApp.cfg is input.")
    parser.add_argument('--dir', '--work-dir', dest='workDir', default='./',
                        help='work directory, (default: %(default)s).')

    parser.add_argument('-g', dest='generate_template', action='store_true',
                        help='generate default template (if it does not exist) and exit.')
    parser.add_argument('-H', dest='print_template', action='store_true',
                        help='print the default template file and exit.')
    parser.add_argument('-v','--version', action='store_true', help='print software version and exit')

    parser.add_argument('--plot', dest='plot', action='store_true',
                        help='plot results [only] without running smallbaselineApp.')

    step = parser.add_argument_group('steps processing (start/end/dostep)', STEP_HELP)
    step.add_argument('--start', dest='startStep', metavar='STEP', default=STEP_LIST[0],
                      help='start processing at the named step (default: %(default)s).')
    step.add_argument('--end','--stop', dest='endStep', metavar='STEP',  default=STEP_LIST[-1],
                      help='end processing at the named step (default: %(default)s)')
    step.add_argument('--dostep', dest='doStep', metavar='STEP',
                      help='run processing at the named step only')

    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    # parse
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    # save argv (to check the manually specified arguments)
    # use iargs        for python call
    # use sys.argv[1:] for command line call
    inps.argv = iargs if iargs else sys.argv[1:]

    # check
    template_file = os.path.join(os.path.dirname(mintpy.__file__), 'defaults/smallbaselineApp.cfg')

    # check: -H option (print default template)
    if inps.print_template:
        with open(template_file) as f:
            lines = f.read()
        try:
            # syntax highlight via rich
            from rich.console import Console
            from rich.syntax import Syntax
            console = Console()
            console.print(Syntax(lines, "cfg", background_color='default'))
        except ImportError:
            print(lines)
        sys.exit(0)

    # check: -v option (print software version)
    if inps.version:
        print(mintpy.version.version_description)
        sys.exit(0)

    # check: existence of input template files
    if (not inps.customTemplateFile
            and not os.path.isfile(os.path.basename(template_file))
            and not inps.generate_template):
        parser.print_usage()
        print(EXAMPLE)

        msg = "no template file found! It requires:"
        msg += "\n  1) input a custom template file, OR"
        msg += "\n  2) there is a default template 'smallbaselineApp.cfg' in current directory."
        raise SystemExit(f'ERROR: {msg}')

    # check: custom input template file
    if inps.customTemplateFile:
        # check the existence
        if not os.path.isfile(inps.customTemplateFile):
            raise FileNotFoundError(inps.customTemplateFile)

        inps.customTemplateFile = os.path.abspath(inps.customTemplateFile)

        # ignore if it is the default template file (smallbaselineApp.cfg)
        if os.path.basename(inps.customTemplateFile) == os.path.basename(template_file):
            inps.customTemplateFile = None

    # check: --plot option
    if inps.argv == ['--plot']:
        plot_only = True
        print('plot smallbaselineApp results without run.')
    else:
        plot_only = False

    # check: --start/end/dostep options
    inps.runSteps = read_inps2run_steps(inps, step_list=STEP_LIST, plot_only=plot_only)

    return inps


def read_inps2run_steps(inps, step_list, plot_only=False):
    """read/get run_steps from input arguments."""
    # check: if start/end/do step input is valid
    for key in ['startStep', 'endStep', 'doStep']:
        value = vars(inps)[key]
        if value and value not in step_list:
            msg = f'Input step not found: {value}'
            msg += f'\nAvailable steps: {step_list}'
            raise ValueError(msg)

    # check: ignore --start/end input if --dostep is specified
    if inps.doStep:
        inps.startStep = inps.doStep
        inps.endStep = inps.doStep

    # get list of steps to run
    idx0 = step_list.index(inps.startStep)
    idx1 = step_list.index(inps.endStep)
    if idx0 > idx1:
        msg = f'start step "{inps.startStep}" is CAN NOT be after the end step "{inps.endStep}"'
        raise ValueError(msg)
    run_steps = step_list[idx0:idx1+1]

    # empty the step list if:
    # a) -g OR
    # b) iargs == ['--plot']
    if inps.generate_template or plot_only:
        run_steps = []

    # print mssage - processing steps
    if len(run_steps) > 0:
        # for single step - compact version info
        if len(run_steps) == 1:
            print(mintpy.version.version_description)
        else:
            print(mintpy.version.logo)

        print(f'--RUN-at-{datetime.datetime.now()}--')
        print(f'Current directory: {os.getcwd()}')
        print(f'Run routine processing with {os.path.basename(__file__)} on steps: {run_steps}')
        print(f'Remaining steps: {step_list[idx0+1:]}')
    print('-'*50)

    return run_steps


##########################################################################
def main(iargs=None):
    # parse
    inps = cmd_line_parse(iargs)

    # import
    from mintpy.smallbaselineApp import run_smallbaselineApp

    # run
    run_smallbaselineApp(inps)


###########################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
