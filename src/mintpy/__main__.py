############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Aug 2022                      #
############################################################


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

import argparse
import logging
import sys

try:
    from os import EX_OK
except ImportError:
    EX_OK = 0
EX_FAILURE = 1
EX_INTERRUPT = 130


from mintpy import __version__

PROG = __package__
LOGFMT = '%(asctime)s %(levelname)-8s -- %(message)s'
DEFAULT_LOGLEVEL = "WARNING"

EXAMPLE = """example:
  mintpy --version
  mintpy smalbaselineApp FernandinaSenDT128.txt
  mintpy tsview timeseries_ERA5_ramp_demErr.h5
  mintpy view velocity.h5

  # debug mode with complete trace back info
  mintpy --debug view velocity.h5
"""


################################################################################################
def _autocomplete(parser):
    try:
        import argcomplete
    except ImportError:
        pass
    else:
        argcomplete.autocomplete(parser)


#######################################  Sub-Parsers  ##########################################
# in alphabetical order
# A-G
def get_add_parser(subparsers=None):
    from mintpy.cli import add
    parser = add.create_parser(subparsers)
    parser.set_defaults(func=add.main)
    return parser


def get_asc_desc2horz_vert_parser(subparsers=None):
    from mintpy.cli import asc_desc2horz_vert
    parser = asc_desc2horz_vert.create_parser(subparsers)
    parser.set_defaults(func=asc_desc2horz_vert.main)
    return parser


def get_closure_phase_bias_parser(subparsers=None):
    from mintpy.cli import closure_phase_bias
    parser = closure_phase_bias.create_parser(subparsers)
    parser.set_defaults(func=closure_phase_bias.main)
    return parser


def get_dem_error_parser(subparsers=None):
    from mintpy.cli import dem_error
    parser = dem_error.create_parser(subparsers)
    parser.set_defaults(func=dem_error.main)
    return parser


def get_dem_gsi_parser(subparsers=None):
    from mintpy.cli import dem_gsi
    parser = dem_gsi.create_parser(subparsers)
    parser.set_defaults(func=dem_gsi.main)
    return parser


def get_diff_parser(subparsers=None):
    from mintpy.cli import diff
    parser = diff.create_parser(subparsers)
    parser.set_defaults(func=diff.main)
    return parser


def get_generate_mask_parser(subparsers=None):
    from mintpy.cli import generate_mask
    parser = generate_mask.create_parser(subparsers)
    parser.set_defaults(func=generate_mask.main)
    return parser


def get_geocode_parser(subparsers=None):
    from mintpy.cli import geocode
    parser = geocode.create_parser(subparsers)
    parser.set_defaults(func=geocode.main)
    return parser


# H-N
def get_ifgram_inversion_parser(subparsers=None):
    from mintpy.cli import ifgram_inversion
    parser = ifgram_inversion.create_parser(subparsers)
    parser.set_defaults(func=ifgram_inversion.main)
    return parser


def get_image_math_parser(subparsers=None):
    from mintpy.cli import image_math
    parser = image_math.create_parser(subparsers)
    parser.set_defaults(func=image_math.main)
    return parser


def get_image_stitch_parser(subparsers=None):
    from mintpy.cli import image_stitch
    parser = image_stitch.create_parser(subparsers)
    parser.set_defaults(func=image_stitch.main)
    return parser


def get_info_parser(subparsers=None):
    from mintpy.cli import info
    parser = info.create_parser(subparsers)
    parser.set_defaults(func=info.main)
    return parser


def get_iono_tec_parser(subparsers=None):
    from mintpy.cli import iono_tec
    parser = iono_tec.create_parser(subparsers)
    parser.set_defaults(func=iono_tec.main)
    return parser


def get_mask_parser(subparsers=None):
    from mintpy.cli import mask
    parser = mask.create_parser(subparsers)
    parser.set_defaults(func=mask.main)
    return parser


def get_multilook_parser(subparsers=None):
    from mintpy.cli import multilook
    parser = multilook.create_parser(subparsers)
    parser.set_defaults(func=multilook.main)
    return parser


def get_load_data_parser(subparsers=None):
    from mintpy.cli import load_data
    parser = load_data.create_parser(subparsers)
    parser.set_defaults(func=load_data.main)
    return parser


def get_load_gbis_parser(subparsers=None):
    from mintpy.cli import load_gbis
    parser = load_gbis.create_parser(subparsers)
    parser.set_defaults(func=load_gbis.main)
    return parser


def get_local_oscilator_drift_parser(subparsers=None):
    from mintpy.cli import local_oscilator_drift
    parser = local_oscilator_drift.create_parser(subparsers)
    parser.set_defaults(func=local_oscilator_drift.main)
    return parser


def get_lookup_geo2radar_parser(subparsers=None):
    from mintpy.cli import lookup_geo2radar
    parser = lookup_geo2radar.create_parser(subparsers)
    parser.set_defaults(func=lookup_geo2radar.main)
    return parser


def get_modify_network_parser(subparsers=None):
    from mintpy.cli import modify_network
    parser = modify_network.create_parser(subparsers)
    parser.set_defaults(func=modify_network.main)
    return parser


# O-Q
def get_plate_motion_parser(subparsers=None):
    from mintpy.cli import plate_motion
    parser = plate_motion.create_parser(subparsers)
    parser.set_defaults(func=plate_motion.main)
    return parser


def get_plot_coherence_matrix_parser(subparsers=None):
    from mintpy.cli import plot_coherence_matrix
    parser = plot_coherence_matrix.create_parser(subparsers)
    parser.set_defaults(func=plot_coherence_matrix.main)
    return parser


def get_plot_network_parser(subparsers=None):
    from mintpy.cli import plot_network
    parser = plot_network.create_parser(subparsers)
    parser.set_defaults(func=plot_network.main)
    return parser


def get_plot_transection_parser(subparsers=None):
    from mintpy.cli import plot_transection
    parser = plot_transection.create_parser(subparsers)
    parser.set_defaults(func=plot_transection.main)
    return parser


def get_prep_aria_parser(subparsers=None):
    from mintpy.cli import prep_aria
    parser = prep_aria.create_parser(subparsers)
    parser.set_defaults(func=prep_aria.main)
    return parser


def get_prep_cosicorr_parser(subparsers=None):
    from mintpy.cli import prep_cosicorr
    parser = prep_cosicorr.create_parser(subparsers)
    parser.set_defaults(func=prep_cosicorr.main)
    return parser


def get_prep_fringe_parser(subparsers=None):
    from mintpy.cli import prep_fringe
    parser = prep_fringe.create_parser(subparsers)
    parser.set_defaults(func=prep_fringe.main)
    return parser


def get_prep_gamma_parser(subparsers=None):
    from mintpy.cli import prep_gamma
    parser = prep_gamma.create_parser(subparsers)
    parser.set_defaults(func=prep_gamma.main)
    return parser


def get_prep_gmtsar_parser(subparsers=None):
    from mintpy.cli import prep_gmtsar
    parser = prep_gmtsar.create_parser(subparsers)
    parser.set_defaults(func=prep_gmtsar.main)
    return parser


def get_prep_hyp3_parser(subparsers=None):
    from mintpy.cli import prep_hyp3
    parser = prep_hyp3.create_parser(subparsers)
    parser.set_defaults(func=prep_hyp3.main)
    return parser


def get_prep_isce_parser(subparsers=None):
    from mintpy.cli import prep_isce
    parser = prep_isce.create_parser(subparsers)
    parser.set_defaults(func=prep_isce.main)
    return parser


def get_prep_nisar_parser(subparsers=None):
    from mintpy.cli import prep_nisar
    parser = prep_nisar.create_parser(subparsers)
    parser.set_defaults(func=prep_nisar.main)
    return parser


def get_prep_roipac_parser(subparsers=None):
    from mintpy.cli import prep_roipac
    parser = prep_roipac.create_parser(subparsers)
    parser.set_defaults(func=prep_roipac.main)
    return parser


def get_prep_snap_parser(subparsers=None):
    from mintpy.cli import prep_snap
    parser = prep_snap.create_parser(subparsers)
    parser.set_defaults(func=prep_snap.main)
    return parser


# R-T
def get_reference_date_parser(subparsers=None):
    from mintpy.cli import reference_date
    parser = reference_date.create_parser(subparsers)
    parser.set_defaults(func=reference_date.main)
    return parser


def get_reference_point_parser(subparsers=None):
    from mintpy.cli import reference_point
    parser = reference_point.create_parser(subparsers)
    parser.set_defaults(func=reference_point.main)
    return parser


def get_remove_hdf5_dset(subparsers=None):
    from mintpy.cli import remove_hdf5_dset
    parser = remove_hdf5_dset.create_parser(subparsers)
    parser.set_defaults(func=remove_hdf5_dset.main)
    return parser


def get_remove_ramp_parser(subparsers=None):
    from mintpy.cli import remove_ramp
    parser = remove_ramp.create_parser(subparsers)
    parser.set_defaults(func=remove_ramp.main)
    return parser


def get_s1ab_range_bias_parser(subparsers=None):
    from mintpy.cli import s1ab_range_bias
    parser = s1ab_range_bias.create_parser(subparsers)
    parser.set_defaults(func=s1ab_range_bias.main)
    return parser


def get_save_gbis_parser(subparsers=None):
    from mintpy.cli import save_gbis
    parser = save_gbis.create_parser(subparsers)
    parser.set_defaults(func=save_gbis.main)
    return parser


def get_save_gdal_parser(subparsers=None):
    from mintpy.cli import save_gdal
    parser = save_gdal.create_parser(subparsers)
    parser.set_defaults(func=save_gdal.main)
    return parser


def get_save_gmt_parser(subparsers=None):
    from mintpy.cli import save_gmt
    parser = save_gmt.create_parser(subparsers)
    parser.set_defaults(func=save_gmt.main)
    return parser


def get_save_hdfeos5_parser(subparsers=None):
    from mintpy.cli import save_hdfeos5
    parser = save_hdfeos5.create_parser(subparsers)
    parser.set_defaults(func=save_hdfeos5.main)
    return parser


def get_save_kite_parser(subparsers=None):
    from mintpy.cli import save_kite
    parser = save_kite.create_parser(subparsers)
    parser.set_defaults(func=save_kite.main)
    return parser


def get_save_kmz_timeseries_parser(subparsers=None):
    from mintpy.cli import save_kmz_timeseries
    parser = save_kmz_timeseries.create_parser(subparsers)
    parser.set_defaults(func=save_kmz_timeseries.main)
    return parser


def get_save_kmz_parser(subparsers=None):
    from mintpy.cli import save_kmz
    parser = save_kmz.create_parser(subparsers)
    parser.set_defaults(func=save_kmz.main)
    return parser


def get_save_qgis_parser(subparsers=None):
    from mintpy.cli import save_qgis
    parser = save_qgis.create_parser(subparsers)
    parser.set_defaults(func=save_qgis.main)
    return parser


def get_save_roipac_parser(subparsers=None):
    from mintpy.cli import save_roipac
    parser = save_roipac.create_parser(subparsers)
    parser.set_defaults(func=save_roipac.main)
    return parser


def get_smallbaselineApp_parser(subparsers=None):
    from mintpy.cli import smallbaselineApp
    parser = smallbaselineApp.create_parser(subparsers)
    parser.set_defaults(func=smallbaselineApp.main)
    return parser


def get_solid_earth_tides_parser(subparsers=None):
    from mintpy.cli import solid_earth_tides
    parser = solid_earth_tides.create_parser(subparsers)
    parser.set_defaults(func=solid_earth_tides.main)
    return parser


def get_spatial_average_parser(subparsers=None):
    from mintpy.cli import spatial_average
    parser = spatial_average.create_parser(subparsers)
    parser.set_defaults(func=spatial_average.main)
    return parser


def get_spatial_filter_parser(subparsers=None):
    from mintpy.cli import spatial_filter
    parser = spatial_filter.create_parser(subparsers)
    parser.set_defaults(func=spatial_filter.main)
    return parser


def get_subset_parser(subparsers=None):
    from mintpy.cli import subset
    parser = subset.create_parser(subparsers)
    parser.set_defaults(func=subset.main)
    return parser


def get_temporal_average_parser(subparsers=None):
    from mintpy.cli import temporal_average
    parser = temporal_average.create_parser(subparsers)
    parser.set_defaults(func=temporal_average.main)
    return parser


def get_temporal_derivative_parser(subparsers=None):
    from mintpy.cli import temporal_derivative
    parser = temporal_derivative.create_parser(subparsers)
    parser.set_defaults(func=temporal_derivative.main)
    return parser


def get_temporal_filter_parser(subparsers=None):
    from mintpy.cli import temporal_filter
    parser = temporal_filter.create_parser(subparsers)
    parser.set_defaults(func=temporal_filter.main)
    return parser


def get_timeseries_rms_parser(subparsers=None):
    from mintpy.cli import timeseries_rms
    parser = timeseries_rms.create_parser(subparsers)
    parser.set_defaults(func=timeseries_rms.main)
    return parser


def get_timeseries2velocity_parser(subparsers=None):
    from mintpy.cli import timeseries2velocity
    parser = timeseries2velocity.create_parser(subparsers)
    parser.set_defaults(func=timeseries2velocity.main)
    return parser


def get_tropo_gacos_parser(subparsers=None):
    from mintpy.cli import tropo_gacos
    parser = tropo_gacos.create_parser(subparsers)
    parser.set_defaults(func=tropo_gacos.main)
    return parser


def get_tropo_phase_elevation_parser(subparsers=None):
    from mintpy.cli import tropo_phase_elevation
    parser = tropo_phase_elevation.create_parser(subparsers)
    parser.set_defaults(func=tropo_phase_elevation.main)
    return parser


def get_tropo_pyaps3_parser(subparsers=None):
    from mintpy.cli import tropo_pyaps3
    parser = tropo_pyaps3.create_parser(subparsers)
    parser.set_defaults(func=tropo_pyaps3.main)
    return parser


def get_tsview_parser(subparsers=None):
    from mintpy.cli import tsview
    parser = tsview.create_parser(subparsers)
    parser.set_defaults(func=tsview.main)
    return parser


# U-Z
def get_unwrap_error_bridging_parser(subparsers=None):
    from mintpy.cli import unwrap_error_bridging
    parser = unwrap_error_bridging.create_parser(subparsers)
    parser.set_defaults(func=unwrap_error_bridging.main)
    return parser


def get_unwrap_error_phase_closure_parser(subparsers=None):
    from mintpy.cli import unwrap_error_phase_closure
    parser = unwrap_error_phase_closure.create_parser(subparsers)
    parser.set_defaults(func=unwrap_error_phase_closure.main)
    return parser


def get_view_parser(subparsers=None):
    from mintpy.cli import view
    parser = view.create_parser(subparsers)
    parser.set_defaults(func=view.main)
    return parser


#######################################  Main Parser  ##########################################
def _add_logging_control_args(parser, default_loglevel=DEFAULT_LOGLEVEL):
    """Add command line options for logging control."""
    loglevels = [logging.getLevelName(level) for level in range(10, 60, 10)]

    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--loglevel",
        default=default_loglevel,
        choices=loglevels,
        help="logging level (default: %(default)s)",
    )
    group.add_argument(
        "--quiet",
        dest="loglevel",
        action="store_const",
        const="ERROR",
        help="suppress standard output messages, "
        "only errors are printed to screen",
    )
    group.add_argument(
        "--verbose",
        dest="loglevel",
        action="store_const",
        const="INFO",
        help="print verbose output messages",
    )
    group.add_argument(
        "--debug",
        dest="loglevel",
        action="store_const",
        const="DEBUG",
        help="print debug messages",
    )


def get_parser():
    """Instantiate the command line argument parser."""
    parser = argparse.ArgumentParser(
        prog=PROG,
        description=__doc__,
        epilog=EXAMPLE,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("-v","--version", action="version", version=f"{__version__}")

    # Sub-command management
    sp = parser.add_subparsers(title="sub-commands", dest='func', required=True, metavar='')

    _add_logging_control_args(parser)

    # workflow
    get_smallbaselineApp_parser(sp)

    # standard processing
    get_asc_desc2horz_vert_parser(sp)
    get_geocode_parser(sp)
    get_ifgram_inversion_parser(sp)
    get_mask_parser(sp)
    get_modify_network_parser(sp)
    get_multilook_parser(sp)
    get_reference_date_parser(sp)
    get_reference_point_parser(sp)
    get_spatial_average_parser(sp)
    get_spatial_filter_parser(sp)
    get_temporal_average_parser(sp)
    get_temporal_derivative_parser(sp)
    get_temporal_filter_parser(sp)
    get_timeseries_rms_parser(sp)
    get_timeseries2velocity_parser(sp)

    # image operations
    get_add_parser(sp)
    get_diff_parser(sp)
    get_image_math_parser(sp)
    get_image_stitch_parser(sp)
    get_subset_parser(sp)

    # noise reduction / error correction
    get_closure_phase_bias_parser(sp)
    get_dem_error_parser(sp)
    get_iono_tec_parser(sp)
    get_local_oscilator_drift_parser(sp)
    get_plate_motion_parser(sp)
    get_remove_ramp_parser(sp)
    get_s1ab_range_bias_parser(sp)
    get_solid_earth_tides_parser(sp)
    get_tropo_gacos_parser(sp)
    get_tropo_phase_elevation_parser(sp)
    get_tropo_pyaps3_parser(sp)
    get_unwrap_error_bridging_parser(sp)
    get_unwrap_error_phase_closure_parser(sp)

    # misc
    get_dem_gsi_parser(sp)
    get_generate_mask_parser(sp)
    try:
        get_lookup_geo2radar_parser(sp)
    except ImportError:
        pass

    # pre-processing
    get_prep_aria_parser(sp)
    get_prep_cosicorr_parser(sp)
    get_prep_fringe_parser(sp)
    get_prep_gamma_parser(sp)
    get_prep_gmtsar_parser(sp)
    get_prep_hyp3_parser(sp)
    get_prep_isce_parser(sp)
    get_prep_nisar_parser(sp)
    get_prep_roipac_parser(sp)
    get_prep_snap_parser(sp)

    # I/O
    get_load_data_parser(sp)
    get_load_gbis_parser(sp)
    get_remove_hdf5_dset(sp)
    get_save_gbis_parser(sp)
    get_save_gdal_parser(sp)
    get_save_gmt_parser(sp)
    get_save_hdfeos5_parser(sp)
    get_save_kite_parser(sp)
    get_save_kmz_timeseries_parser(sp)
    get_save_kmz_parser(sp)
    get_save_qgis_parser(sp)
    get_save_roipac_parser(sp)

    # visualization
    get_info_parser(sp)
    # get_multi_transect_parser(sp)
    get_plot_coherence_matrix_parser(sp)
    get_plot_network_parser(sp)
    get_plot_transection_parser(sp)
    get_tsview_parser(sp)
    get_view_parser(sp)

    _autocomplete(parser)

    return parser


################################################################################################
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
        # NOTE: use the root logger to set the logging level
        logging.getLogger().setLevel(args.loglevel)
        log.debug("args: %s", args)

        exit_code = args.func(sys.argv[2:])
    except Exception as exc:
        log.critical(
            "unexpected exception caught: {!r} {}".format(
                type(exc).__name__, exc)
        )
        log.exception("stacktrace:")  #, exc_info=True)
        exit_code = EX_FAILURE
    except KeyboardInterrupt:
        log.warning("Keyboard interrupt received: exit the program")
        exit_code = EX_INTERRUPT

    return exit_code


################################################################################################
if __name__ == "__main__":
    sys.exit(main())
