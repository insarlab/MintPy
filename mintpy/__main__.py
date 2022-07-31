"""Command line interface for MintPy.

The Miami INsar Time-series software in PYthon (MintPy as /mint pai/)
is an open-source package for Interferometric Synthetic Aperture Radar
(InSAR) time series analysis.

It reads the stack of interferograms (coregistered and unwrapped) in
ISCE, ARIA, FRInGE, HyP3, GMTSAR, SNAP, GAMMA or ROI_PAC format, and
produces three dimensional (2D in space and 1D in time) ground surface
displacement in line-of-sight direction.

It includes a routine time series analysis (`smallbaselineApp.py`) and
some independent toolbox.

This is research code provided to you "as is" with NO WARRANTIES OF
CORRECTNESS. Use at your own risk.
"""

# PYTHON_ARGCOMPLETE_OK

import sys
import logging
import argparse

try:
    from os import EX_OK
except ImportError:
    EX_OK = 0
EX_FAILURE = 1
EX_INTERRUPT = 130


from . import __version__

PROG = __package__
LOGFMT = '%(asctime)s %(levelname)-8s -- %(message)s'


def _autocomplete(parser):
    try:
        import argcomplete
    except ImportError:
        pass
    else:
        argcomplete.autocomplete(parser)


# processing
def get_smallbaseline_parser(subparsers=None):
    from . import smallbaselineApp
    parser = smallbaselineApp.create_parser(subparsers)
    parser.set_defaults(func=smallbaselineApp.main)
    return parser


def get_geocode_parser(subparsers=None):
    from . import geocode
    parser = geocode.create_parser(subparsers)
    parser.set_defaults(func=geocode.main)
    return parser


def get_multilook_parser(subparsers=None):
    from . import multilook
    parser = multilook.create_parser(subparsers)
    parser.set_defaults(func=multilook.main)
    return parser


def get_spatial_average_parser(subparsers=None):
    from . import spatial_average
    parser = spatial_average.create_parser(subparsers)
    parser.set_defaults(func=spatial_average.main)
    return parser


def get_spatial_filter_parser(subparsers=None):
    from . import spatial_filter
    parser = spatial_filter.create_parser(subparsers)
    parser.set_defaults(func=spatial_filter.main)
    return parser


def get_temporal_average_parser(subparsers=None):
    from . import temporal_average
    parser = temporal_average.create_parser(subparsers)
    parser.set_defaults(func=temporal_average.main)
    return parser


def get_temporal_derivative_parser(subparsers=None):
    from . import temporal_derivative
    parser = temporal_derivative.create_parser(subparsers)
    parser.set_defaults(func=temporal_derivative.main)
    return parser


def get_temporal_filter_parser(subparsers=None):
    from . import temporal_filter
    parser = temporal_filter.create_parser(subparsers)
    parser.set_defaults(func=temporal_filter.main)
    return parser


# pre-processing
def get_prep_aria_parser(subparsers=None):
    from . import prep_aria
    parser = prep_aria.create_parser(subparsers)
    parser.set_defaults(func=prep_aria.main)
    return parser


def get_prep_cosicorr_parser(subparsers=None):
    from . import prep_cosicorr
    parser = prep_cosicorr.create_parser(subparsers)
    parser.set_defaults(func=prep_cosicorr.main)
    return parser


def get_prep_fringe_parser(subparsers=None):
    from . import prep_fringe
    parser = prep_fringe.create_parser(subparsers)
    parser.set_defaults(func=prep_fringe.main)
    return parser


def get_parser():
    """Instantiate the command line argument parser."""
    parser = argparse.ArgumentParser(prog=PROG, description=__doc__)
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}" 
    )

    # Sub-command management
    sp = parser.add_subparsers(title="sub-commands", dest='func')

    # processing
    get_smallbaseline_parser(sp)
    get_geocode_parser(sp)
    get_multilook_parser(sp)
    get_spatial_average_parser(sp)
    get_spatial_filter_parser(sp)
    get_temporal_average_parser(sp)
    get_temporal_derivative_parser(sp)
    get_temporal_filter_parser(sp)

    # pre-processing
    get_prep_aria_parser(sp)
    get_prep_cosicorr_parser(sp)
    get_prep_fringe_parser(sp)

    _autocomplete(parser)

    return parser


def main(*argv):
    """Main CLI interface."""
    # setup logging
    logging.basicConfig(format=LOGFMT)
    logging.captureWarnings(True)
    log = logging.getLogger(PROG)

    # parse cmd line arguments
    parser = get_parser()
    args = parser.parse_args(argv if argv else None)

    # execute main tasks
    exit_code = EX_OK
    try:
        return args.func(sys.argv[2:])
    except Exception as exc:
        log.critical(
            "unexpected exception caught: {!r} {}".format(
                type(exc).__name__, exc)
        )
        log.debug("stacktrace:", exc_info=True)
        exit_code = EX_FAILURE
    except KeyboardInterrupt:
        log.warning("Keyboard interrupt received: exit the program")
        exit_code = EX_INTERRUPT

    return exit_code


if __name__ == "__main__":
    sys.exit(main())
