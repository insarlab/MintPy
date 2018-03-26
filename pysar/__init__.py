#! /usr/bin/env python2

########################################################################
miami_path = True    # Package-wide variable, Auto setting for University of Miami
                     # change it to False if you are not using the file structure of University of Miami
parallel_num = 8     # max core number used in parallel processing
figsize_single_min = 6.0       # default min size in inch, for single plot
figsize_single_max = 12.0      # default min size in inch, for single plot
figsize_multi = [15.0, 8.0]    # default size in inch, for multiple subplots


###################### Do not change below this line ###################
from . import _datetime
from . import _gmt
from . import _readfile
from . import _writefile

from . import _network
from . import _remove_surface
from . import _pysar_utilities

from . import subset
from . import mask
from . import multilook

from . import view
#import tsviewer

from . import add
from . import asc_desc
from . import baseline_error
from . import baseline_trop
#import correlation_with_dem
from . import dem_error
from . import diff
#import generate_mask
#import geocode
#import ifgram_closure
#import ifgram_inversion
#import ifgram_reconstruction
#import ifgram_simulation
from . import image_math
from . import incidence_angle
from . import info
from . import insar_vs_gps
#import l1
from . import load_data
#import load_dem
from . import lod
#import look_angle
from . import match
#import temporal_average
#import spatial_average
from . import modify_network
#import multi_transect
#import plot_network
#import pysarApp
#import quality_map
#import reference_epoch
#import remove_plane
#import save_gmt
#import save_kml
#import save_unw
#import save_mat
#import save_hdfeos5
from . import seed_data
#import spatial_filter
from . import sum_epochs
#import temporal_filter
#import temporal_coherence
#import temporal_derivative
from . import timeseries2velocity
from . import transect
from . import tropcor_phase_elevation
#import tropcor_pyaps
from . import unwrap_error
from . import view_gui
from . import ts_gui
