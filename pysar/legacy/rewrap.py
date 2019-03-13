#!/usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2013-2018, Heresh Fattahi                   #
# Author:  Heresh Fattahi                                  #
############################################################


import os
import sys
import h5py
import numpy as np
from pysar.utils import readfile, ptime


############################################################
USAGE = """
usage: rewrap.py  ifgram_file   [cycle_unit]

Re-wrap unwraped interferograms to wrapped interferograms.

example:
  rewrap.py  interferograms.h5
  rewrap.py  timeseries_ECMWF_demErr.h5
"""

def usage():
    print(USAGE)
    return


def rewrap(unw, cycle=2*np.pi):
    rewrapped = unw - np.round(unw/cycle) * cycle
    return rewrapped


def main(argv):
    try:
        file = argv[0]
    except:
        usage()
        sys.exit(1)

    outfile = os.path.splitext(file)[0]+'_wrap'+os.path.splitext(file)[1]
    one_cycle = 2*np.pi
    one_cycle = 0.05

    atr = readfile.read_attribute(file)
    k = atr['FILE_TYPE']

    if k in ['interferograms', 'coherence', 'wrapped', 'timeseries']:
        h5 = h5py.File(file, 'r')
        epochList = sorted(h5[k].keys())
        epoch_num = len(epochList)
        prog_bar = ptime.progressBar(maxValue=epoch_num)

        print('writing >>> '+outfile)
        h5out = h5py.File(outfile, 'w')
        group = h5out.create_group(k)

        if k in ['interferograms', 'coherence', 'wrapped']:
            date12_list = ptime.list_ifgram2date12(epochList)
            print('number of interferograms: '+str(len(epochList)))
            for i in range(epoch_num):
                epoch = epochList[i]
                data = h5[k][epoch].get(epoch)[:]

                data_wrap = rewrap(data)

                gg = group.create_group(epoch)
                dset = gg.create_dataset(epoch, data=data_wrap)
                for key, value in h5[k][epoch].attrs.items():
                    gg.attrs[key] = value
                prog_bar.update(i+1, suffix=date12_list[i])

        elif k == 'timeseries':
            print('number of acquisitions: '+str(len(epochList)))
            for i in range(epoch_num):
                epoch = epochList[i]
                data = h5[k].get(epoch)[:]

                data_wrap = rewrap(data, one_cycle)

                dset = group.create_dataset(epoch, data=data_wrap)
                prog_bar.update(i+1, suffix=epoch)
            for key, value in h5[k].attrs.items():
                group.attrs[key] = value

        h5.close()
        h5out.close()
        prog_bar.close()

    print('Done.')
    return outfile


if __name__ == '__main__':
    main(sys.argv[1:])
