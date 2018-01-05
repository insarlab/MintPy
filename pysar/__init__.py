#! /usr/bin/env python2

########################################################################
miami_path = True    # Package-wide variable, Auto setting for University of Miami
                     # change it to False if you are not using the file structure of University of Miami
parallel_num = 8     # max core number used in parallel processing
figsize_single_min = 6.0       # default min size in inch, for single plot
figsize_single_max = 12.0      # default min size in inch, for single plot
figsize_multi = [15.0, 8.0]    # default size in inch, for multiple subplots


###################### Do not change below this line ###################
import _datetime
import _gmt
import _readfile
import _writefile

import _network
import _remove_surface
import _pysar_utilities

import subset
import mask
import multilook

import view
#import tsviewer

import add
import asc_desc
import baseline_error
import baseline_trop
#import correlation_with_dem
import dem_error
import diff
#import generate_mask
#import geocode
#import ifgram_closure
#import ifgram_inversion
#import ifgram_reconstruction
#import ifgram_simulation
import image_math
import incidence_angle
import info
import insar_vs_gps
#import l1
import load_data
#import load_dem
import lod
#import look_angle
import match
#import temporal_average
#import spatial_average
import modify_network
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
import seed_data
#import spatial_filter
import sum_epochs
#import temporal_filter
#import temporal_coherence
#import temporal_derivative
import timeseries2velocity
import transect
import tropcor_phase_elevation
#import tropcor_pyaps
import unwrap_error
