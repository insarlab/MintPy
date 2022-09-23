## Dynamic import for modules used in routine workflows, i.e. smallbaselineApp
## Recommended usage:
##     import mintpy
##     import mintpy.workflow


import importlib
from pathlib import Path

# this fixes the UnboundLocalError: local variable 'mintpy' referenced before assignment
import mintpy

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
    'solid_earth_tides',
    'temporal_average',
    'timeseries2velocity',
    'timeseries_rms',
    'tropo_gacos',
    'tropo_phase_elevation',
    'tropo_pyaps3',
    'unwrap_error_bridging',
    'unwrap_error_phase_closure',
    'view',
]

root_module = Path(__file__).parent.parent.name   #mintpy
for module in __all__:
    importlib.import_module(f'{root_module}.cli.{module}')
