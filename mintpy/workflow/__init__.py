## Dynamic import for modules used in routine workflows, i.e. smallbaselineApp
## Recommended usage:
##     import mintpy
##     import mintpy.workflow


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
    'tropo_gacos',
    'unwrap_error_bridging',
    'unwrap_error_phase_closure',
    'version',
    'view',
]

root_module = Path(__file__).parent.parent.name   #mintpy
for module in __all__:
    importlib.import_module(root_module + '.' + module)
