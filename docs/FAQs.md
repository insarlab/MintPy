## Frequently Asked Questions

#### 1. What's the sign convention of the line-of-sight data?

For line-of-sight (LOS) phase in the unit of radians, i.e. 'unwrapPhase' dataset in `ifgramStack.h5` file, positive value represents motion away from the satellite. We assume the "date1_date2" format for the interferogram with "date1" being the earlier acquisition.

For LOS displacement in the unit of meters, i.e. 'timeseries' dataset in `timeseries.h5` file positive value represents motion toward the satellite (uplift for pure vertical motion).

#### 2. How to prepare the input for MintPy if I am using InSAR software rather than ISCE stack processors and ARIA-tools?




