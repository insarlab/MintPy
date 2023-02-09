############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2017                               #
############################################################


import os
import time

import numpy as np

from mintpy.objects.resample import resample
from mintpy.utils import attribute as attr, readfile, utils as ut, writefile


############################################################################################
def auto_output_filename(in_file, inps):
    if len(inps.file) == 1 and inps.outfile:
        return inps.outfile

    fbase, fext = os.path.splitext(os.path.basename(in_file))
    prefix = 'geo_' if inps.radar2geo else 'rdr_'
    suffix = inps.dset if inps.dset else fbase
    out_file = f'{prefix}{suffix}{fext}'

    if inps.out_dir:
        if not os.path.isdir(inps.out_dir):
            os.makedirs(inps.out_dir)
            print(f'create directory: {inps.out_dir}')
        out_file = os.path.join(inps.out_dir, out_file)

    return out_file


def run_geocode(inps):
    """geocode all input files"""
    start_time = time.time()

    # feed the largest file for resample object initiation
    ind_max = np.argmax([os.path.getsize(i) for i in inps.file])

    # prepare geometry for geocoding
    kwargs = dict(interp_method=inps.interpMethod,
                  fill_value=inps.fillValue,
                  nprocs=inps.nprocs,
                  max_memory=inps.maxMemory,
                  software=inps.software,
                  print_msg=True)
    if inps.latFile and inps.lonFile:
        kwargs['lat_file'] = inps.latFile
        kwargs['lon_file'] = inps.lonFile
    res_obj = resample(lut_file=inps.lookupFile,
                       src_file=inps.file[ind_max],
                       SNWE=inps.SNWE,
                       lalo_step=inps.laloStep,
                       **kwargs)
    res_obj.open()
    res_obj.prepare()

    # resample input files one by one
    for infile in inps.file:
        print('-' * 50+f'\nresampling file: {infile}')
        atr = readfile.read_attribute(infile, datasetName=inps.dset)
        outfile = auto_output_filename(infile, inps)

        # update_mode
        if inps.updateMode:
            print('update mode: ON')
            if ut.run_or_skip(outfile, in_file=[infile, inps.lookupFile]) == 'skip':
                continue

        ## prepare output
        # update metadata
        if inps.radar2geo:
            atr = attr.update_attribute4radar2geo(atr, res_obj=res_obj)
        else:
            atr = attr.update_attribute4geo2radar(atr, res_obj=res_obj)

        # instantiate output file
        hdf5_file = os.path.splitext(outfile)[1] in ['.h5', '.he5']
        if hdf5_file:
            # grab metadata from input file: dataset compression and UNIT
            compression = readfile.get_hdf5_compression(infile)
            ds_unit_dict = readfile.get_hdf5_dataset_attrs(infile, key='UNIT')
            # initiate output file
            writefile.layout_hdf5(
                outfile,
                metadata=atr,
                ds_unit_dict=ds_unit_dict,
                ref_file=infile,
                compression=compression,
            )
        else:
            dsDict = dict()

        ## run
        dsNames = readfile.get_dataset_list(infile, datasetName=inps.dset)
        maxDigit = max(len(i) for i in dsNames)
        for dsName in dsNames:

            if not hdf5_file:
                dsDict[dsName] = np.zeros((res_obj.length, res_obj.width), dtype=atr['DATA_TYPE'])

            # loop for block-by-block IO
            for i in range(res_obj.num_box):
                src_box = res_obj.src_box_list[i]
                dest_box = res_obj.dest_box_list[i]

                # read
                print('-'*50 + f'{i+1}/{res_obj.num_box}')
                print('reading {d:<{w}} in block {b} from {f} ...'.format(
                    d=dsName, w=maxDigit, b=src_box, f=os.path.basename(infile)))

                data = readfile.read(infile,
                                     datasetName=dsName,
                                     box=src_box,
                                     print_msg=False)[0]

                # resample
                data = res_obj.run_resample(src_data=data, box_ind=i)

                # write / save block data
                if data.ndim == 3:
                    block = [0, data.shape[0],
                             dest_box[1], dest_box[3],
                             dest_box[0], dest_box[2]]
                else:
                    block = [dest_box[1], dest_box[3],
                             dest_box[0], dest_box[2]]

                if hdf5_file:
                    print(f'write data in block {block} to file: {outfile}')
                    writefile.write_hdf5_block(outfile,
                                               data=data,
                                               datasetName=dsName,
                                               block=block,
                                               print_msg=False)
                else:
                    dsDict[dsName][block[0]:block[1],
                                   block[2]:block[3]] = data

            # for binary file: ensure same data type
            if not hdf5_file:
                dsDict[dsName] = np.array(dsDict[dsName], dtype=data.dtype)

        # write binary file
        if not hdf5_file:
            atr['BANDS'] = len(dsDict.keys())
            writefile.write(dsDict, out_file=outfile, metadata=atr, ref_file=infile)

            # create ISCE XML and GDAL VRT file if using ISCE lookup table file
            if inps.latFile and inps.lonFile:
                writefile.write_isce_xml(atr, fname=outfile)

    # used time
    m, s = divmod(time.time()-start_time, 60)
    print(f'time used: {m:02.0f} mins {s:02.1f} secs.\n')

    return outfile
