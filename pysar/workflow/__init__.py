## Dynamic import for modules used in routine workflows, i.e. pysarApp
## Recommended usage:
##     import pysar
##     import pysar.workflow


from pathlib import Path
import importlib


# expose the following modules
__all__ = [
    'dem_error',
    'diff',
    'generate_mask',
    'geocode',
    'ifgram_inversion',
    'load_data',
    'local_oscilator_drift',
    'modify_network',
    'plot_network',
    'reference_date',
    'reference_point',
    'remove_ramp',
    'save_hdfeos5',
    'save_kmz',
    'temporal_average',
    'timeseries2velocity',
    'timeseries_rms',
    'tropo_phase_elevation',
    'version',
]

root_module = Path(__file__).parent.parent.name   #pysar
for module in __all__:
    importlib.import_module(root_module + '.' + module)
