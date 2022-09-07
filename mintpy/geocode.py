############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2017                               #
############################################################


import os
import sys
import time
import numpy as np

from mintpy.objects.resample import resample
from mintpy.utils import (
    readfile,
    writefile,
    utils as ut,
    attribute as attr,
)


######################################################################################
def read_template2inps(template_file, inps):
    """Read input template options into Namespace inps"""
    print('read input option from template file: ' + template_file)

    inps_dict = vars(inps)
    template = readfile.read_template(template_file, skip_chars=['[', ']'])
    template = ut.check_template_auto_value(template)

    prefix = 'mintpy.geocode.'
    key_list = [i for i in list(inps_dict.keys()) if prefix + i in template.keys()]
    for key in key_list:
        value = template[prefix + key]
        if value:
            if key in ['SNWE', 'laloStep']:
                inps_dict[key] = [float(i) for i in value.split(',')]
            elif key in ['interpMethod']:
                inps_dict[key] = value
            elif key == 'fillValue':
                if 'nan' in value.lower():
                    inps_dict[key] = np.nan
                else:
                    inps_dict[key] = float(value)

    # computing configurations
    key = 'mintpy.compute.maxMemory'
    if key in template.keys() and template[key]:
        inps.maxMemory = float(template[key])

    return inps


############################################################################################
def check_num_processor(nprocs):
    """Check number of processors
    Note by Yunjun, 2019-05-02:
    1. conda install pyresample will install pykdtree and openmp, but it seems not working:
        geocode.py is getting slower with more processors
            Test on a TS HDF5 file in size of (241, 2267, 2390)
            Memory: up to 10GB
            Run time: 2.5 mins for nproc=1, 3 mins for nproc=4
    2. macports seems to have minor speedup when more processors
    Thus, default number of processors is set to 1; although the capability of using multiple
    processors is written here.
    """
    if not nprocs:
        #OMP_NUM_THREADS is defined in environment variable for OpenMP
        if 'OMP_NUM_THREADS' in os.environ:
            nprocs = int(os.getenv('OMP_NUM_THREADS'))
        else:
            nprocs = int(os.cpu_count() / 2)
    nprocs = min(os.cpu_count(), nprocs)
    print('number of processor to be used: {}'.format(nprocs))
    return nprocs


def auto_output_filename(infile, inps):
    if len(inps.file) == 1 and inps.outfile:
        return inps.outfile

    if inps.radar2geo:
        prefix = 'geo_'
    else:
        prefix = 'rdr_'

    if inps.dset:
        outfile = '{}{}.h5'.format(prefix, inps.dset)
    else:
        outfile = '{}{}'.format(prefix, os.path.basename(infile))

    if inps.out_dir:
        if not os.path.isdir(inps.out_dir):
            os.makedirs(inps.out_dir)
            print('create directory: {}'.format(inps.out_dir))
        outfile = os.path.join(inps.out_dir, outfile)
    return outfile


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
        print('-' * 50+'\nresampling file: {}'.format(infile))
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
        file_is_hdf5 = os.path.splitext(outfile)[1] in ['.h5', '.he5']
        if file_is_hdf5:
            compression = readfile.get_hdf5_compression(infile)
            writefile.layout_hdf5(outfile, metadata=atr, ref_file=infile, compression=compression)
        else:
            dsDict = dict()

        ## run
        dsNames = readfile.get_dataset_list(infile, datasetName=inps.dset)
        maxDigit = max([len(i) for i in dsNames])
        for dsName in dsNames:

            if not file_is_hdf5:
                dsDict[dsName] = np.zeros((res_obj.length, res_obj.width))

            # loop for block-by-block IO
            for i in range(res_obj.num_box):
                src_box = res_obj.src_box_list[i]
                dest_box = res_obj.dest_box_list[i]

                # read
                print('-'*50+'\nreading {d:<{w}} in block {b} from {f} ...'.format(
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

                if file_is_hdf5:
                    print('write data in block {} to file: {}'.format(block, outfile))
                    writefile.write_hdf5_block(outfile,
                                               data=data,
                                               datasetName=dsName,
                                               block=block,
                                               print_msg=False)
                else:
                    dsDict[dsName][block[0]:block[1],
                                   block[2]:block[3]] = data

            # for binary file: ensure same data type
            if not file_is_hdf5:
                dsDict[dsName] = np.array(dsDict[dsName], dtype=data.dtype)

        # write binary file
        if not file_is_hdf5:
            atr['BANDS'] = len(dsDict.keys())
            writefile.write(dsDict, out_file=outfile, metadata=atr, ref_file=infile)

            # create ISCE XML and GDAL VRT file if using ISCE lookup table file
            if inps.latFile and inps.lonFile:
                writefile.write_isce_xml(atr, fname=outfile)

    m, s = divmod(time.time()-start_time, 60)
    print('time used: {:02.0f} mins {:02.1f} secs.\n'.format(m, s))
    return outfile
