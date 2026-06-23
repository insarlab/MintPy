############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, 2017                               #
############################################################


# os 用来处理路径、创建输出目录、判断文件大小等。
import os
# time 用来统计地理编码处理耗时。
import time

# numpy 用来处理数组，例如创建输出数组、判断最大文件索引。
import numpy as np

# resample 是 MintPy 的重采样对象，负责根据查找表把数据从一个坐标系重采样到另一个坐标系。
from mintpy.objects.resample import resample
# attribute(attr) 负责更新文件元数据；readfile/writefile 负责读写数据；ut 是通用工具函数。
from mintpy.utils import attribute as attr, readfile, utils as ut, writefile


############################################################################################
def auto_output_filename(in_file, inps):
    # 自动生成输出文件名。
    # 如果只处理一个输入文件，并且用户明确指定了输出文件名，就直接使用用户指定的名字。
    if len(inps.file) == 1 and inps.outfile:
        return inps.outfile

    # os.path.basename() 取文件名；os.path.splitext() 把文件名拆成“主名”和“扩展名”。
    fbase, fext = os.path.splitext(os.path.basename(in_file))
    # radar2geo=True 表示雷达坐标转地理坐标，输出名前缀用 geo_；否则用 rdr_。
    prefix = 'geo_' if inps.radar2geo else 'rdr_'
    # 如果用户只处理某个数据集 dset，就用数据集名做输出名后缀；否则用原文件主名。
    suffix = inps.dset if inps.dset else fbase
    out_file = f'{prefix}{suffix}{fext}'

    if inps.out_dir:
        # 如果指定输出目录但目录不存在，就先创建目录。
        if not os.path.isdir(inps.out_dir):
            os.makedirs(inps.out_dir)
            print(f'create directory: {inps.out_dir}')
        # os.path.join() 把目录和文件名拼成完整路径。
        out_file = os.path.join(inps.out_dir, out_file)

    return out_file


def run_geocode(inps):
    """geocode all input files"""
    # geocode 的核心任务：把输入文件中的二维/三维栅格数据重采样到另一个坐标系。
    # 常见方向是 radar2geo：雷达坐标 -> 经纬度地理坐标。
    start_time = time.time()

    # feed the largest file for resample object initiation
    # 选择最大的输入文件初始化 resample 对象，通常可以代表最完整的数据尺寸和结构。
    # os.path.getsize(i) 返回文件大小；np.argmax() 返回最大值的索引。
    ind_max = np.argmax([os.path.getsize(i) for i in inps.file])

    # prepare geometry for geocoding
    # kwargs 是传给 resample 对象的配置字典，包括插值方法、填充值、并行数和内存限制。
    kwargs = dict(interp_method=inps.interpMethod,
                  fill_value=inps.fillValue,
                  nprocs=inps.nprocs,
                  max_memory=inps.maxMemory,
                  software=inps.software,
                  print_msg=True)
    if inps.latFile and inps.lonFile:
        # 有些数据用经纬度文件而不是 lookupFile 表示坐标转换关系。
        kwargs['lat_file'] = inps.latFile
        kwargs['lon_file'] = inps.lonFile
    # 创建 resample 对象。lut_file 是 lookup table，描述源坐标和目标坐标的对应关系。
    res_obj = resample(lut_file=inps.lookupFile,
                       src_file=inps.file[ind_max],
                       SNWE=inps.SNWE,
                       lalo_step=inps.laloStep,
                       **kwargs)
    # open() 读取查找表和源文件的基本信息；prepare() 预计算重采样需要的块和索引。
    res_obj.open()
    res_obj.prepare()

    # resample input files one by one
    for infile in inps.file:
        print('-' * 50+f'\nresampling file: {infile}')
        # 读取输入文件属性；datasetName=inps.dset 表示可只读取某个指定数据集的属性。
        atr = readfile.read_attribute(infile, datasetName=inps.dset)
        outfile = auto_output_filename(infile, inps)

        # update_mode
        if inps.updateMode:
            print('update mode: ON')
            # 如果输出文件已存在且比输入文件/查找表新，就跳过，避免重复计算。
            if ut.run_or_skip(outfile, in_file=[infile, inps.lookupFile]) == 'skip':
                continue

        ## prepare output
        # update metadata
        if inps.radar2geo:
            # 雷达坐标转地理坐标时，需要把元数据中的坐标信息更新成经纬度网格。
            atr = attr.update_attribute4radar2geo(atr, res_obj=res_obj)
        else:
            # 地理坐标转雷达坐标时，需要把元数据更新成雷达坐标尺寸和坐标定义。
            atr = attr.update_attribute4geo2radar(atr, res_obj=res_obj)

        # instantiate output file
        # 判断输出是否是 HDF5 文件；.h5 和 .he5 需要先创建文件结构，再分块写数据。
        hdf5_file = os.path.splitext(outfile)[1] in ['.h5', '.he5']
        if hdf5_file:
            # grab metadata from input file: dataset compression and UNIT
            # 保留原 HDF5 文件的数据压缩方式和每个数据集的 UNIT 属性。
            compression = readfile.get_hdf5_compression(infile)
            ds_unit_dict = readfile.get_hdf5_dataset_attrs(infile, key='UNIT')
            # initiate output file
            # layout_hdf5() 先创建空的 HDF5 文件和数据集布局，后面再按块写入数据。
            writefile.layout_hdf5(
                outfile,
                metadata=atr,
                ds_unit_dict=ds_unit_dict,
                ref_file=infile,
                compression=compression,
            )
        else:
            # 非 HDF5 输出先用字典在内存里收集所有数据集，最后统一写出。
            dsDict = dict()

        ## run
        # get_dataset_list() 获取要处理的数据集列表；如果指定 inps.dset，就只返回相关数据集。
        dsNames = readfile.get_dataset_list(infile, datasetName=inps.dset)
        maxDigit = max(len(i) for i in dsNames)
        for dsName in dsNames:

            if not hdf5_file:
                # 非 HDF5 文件需要先创建目标大小的空数组。
                dsDict[dsName] = np.zeros((res_obj.length, res_obj.width), dtype=atr['DATA_TYPE'])

            # loop for block-by-block IO
            for i in range(res_obj.num_box):
                # 为了节省内存，大文件不会一次性读入，而是按 res_obj 预先划分的小块逐块处理。
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
                # run_resample() 对当前数据块执行坐标转换和插值。
                data = res_obj.run_resample(src_data=data, box_ind=i)

                # write / save block data
                if data.ndim == 3:
                    # 三维数据通常是 time/date x row x col，因此 block 要包含第 0 维范围。
                    block = [0, data.shape[0],
                             dest_box[1], dest_box[3],
                             dest_box[0], dest_box[2]]
                else:
                    # 二维数据只需要行列范围。
                    block = [dest_box[1], dest_box[3],
                             dest_box[0], dest_box[2]]

                if hdf5_file:
                    print(f'write data in block {block} to file: {outfile}')
                    # HDF5 文件直接把当前块写入目标文件对应位置。
                    writefile.write_hdf5_block(outfile,
                                               data=data,
                                               datasetName=dsName,
                                               block=block,
                                               print_msg=False)
                else:
                    # 非 HDF5 文件先把当前块放进内存数组对应位置。
                    dsDict[dsName][block[0]:block[1],
                                   block[2]:block[3]] = data

            # for binary file: ensure same data type
            if not hdf5_file:
                # 保持输出数组数据类型和实际重采样结果一致。
                dsDict[dsName] = np.array(dsDict[dsName], dtype=data.dtype)

        # write binary file
        if not hdf5_file:
            # 非 HDF5 格式最后一次性写出，并附带更新后的元数据。
            atr['BANDS'] = len(dsDict.keys())
            writefile.write(dsDict, out_file=outfile, metadata=atr, ref_file=infile)

            # create ISCE XML and GDAL VRT file if using ISCE lookup table file
            if inps.latFile and inps.lonFile:
                # ISCE 生态常用 XML/VRT 描述二进制文件，这里一并生成。
                writefile.write_isce_xml(atr, fname=outfile)

    # used time
    # 统计所有输入文件地理编码总耗时。
    m, s = divmod(time.time()-start_time, 60)
    print(f'time used: {m:02.0f} mins {s:02.1f} secs.\n')

    return outfile
