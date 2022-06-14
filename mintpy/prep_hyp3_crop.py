#!/usr/bin/env python3
############################################################
# Copyright (c) 2021,         VEG                          #
# Author: Xie Lei, Wang Jiageng,  Jun 2022                 #
# A special thank to Benjy Marks                           #
############################################################
import argparse
import os
import shutil
import glob
from zipfile import ZipFile
import numpy as np
import pathlib
from osgeo import gdal
from osgeo import osr
from hyp3lib.cutGeotiffs import cutFiles
import mintpy.prep_hyp3


def gdal_clip(source_path, target_extent, output_path=None, target_format='GTiff', flag_prompt=False):
    """
    Generate and run gdal command for clip img
    :param source_path: source path need to be clipped
    :param target_extent: e_min north_min e_max north_max
    :param output_path: output path, can be None
    :param target_format: Default GTiff
    :param flag_prompt: Default False
    :return:
    """
    # gdalwarp cmd
    # e.g., gdalwarp -te 78.25 30.22 78.67 30.64 -of GTiff input output

    if output_path is None:
        p = pathlib.Path(os.path.abspath(source_path))
        output_path = os.path.join(p.parent, '%s_clip%s' % (p.stem, p.suffix))

    cmd = 'gdalwarp -te %s %s %s %s -of %s %s %s' % (
        target_extent[0], target_extent[1], target_extent[2], target_extent[3],
        target_format, source_path, output_path)
    if flag_prompt:
        print(':: GDAL cmd: %s' % cmd)
        print('-' * 50)

    os.system(cmd)
    return True


def read_img(filename):
    dataset = gdal.Open(filename)

    im_width = dataset.RasterXSize
    im_height = dataset.RasterYSize

    im_geotrans = dataset.GetGeoTransform()
    im_proj = dataset.GetProjection()
    im_data = dataset.ReadAsArray(0, 0, im_width, im_height)

    minx = im_geotrans[0]
    maxy = im_geotrans[3]
    maxx = minx + im_geotrans[1] * dataset.RasterXSize
    miny = maxy + im_geotrans[5] * dataset.RasterYSize

    im_extent = [minx, maxx, miny, maxy]
    del dataset
    return im_proj, im_geotrans, im_data, im_extent


def hyp3_reprojection(UTM_path, latlon_path):
    gdal.AllRegister()
    src_data = gdal.Open(UTM_path)
    srcSRS_wkt = src_data.GetProjection()
    srcSRS = osr.SpatialReference()
    srcSRS.ImportFromWkt(srcSRS_wkt)
    src_width = src_data.RasterXSize
    src_height = src_data.RasterYSize
    src_count = src_data.RasterCount
    src_trans = src_data.GetGeoTransform()
    OriginLX_src = src_trans[0]
    OriginTY_src = src_trans[3]
    pixl_w_src = src_trans[1]
    pixl_h_src = src_trans[5]
    OriginRX_src = OriginLX_src + pixl_w_src * src_width
    OriginBY_src = OriginTY_src + pixl_h_src * src_height
    driver = gdal.GetDriverByName("GTiff")
    driver.Register()
    dst_data = driver.Create(latlon_path, src_width, src_height, 1, gdal.GDT_Float32)
    dstSRS = osr.SpatialReference()
    dstSRS.ImportFromEPSG(4326)
    dstSRS.ImportFromProj4("proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
    ct = osr.CoordinateTransformation(srcSRS, dstSRS)
    OriginLX_dst, OriginTY_dst, temp = ct.TransformPoint(OriginLX_src, OriginTY_src)
    OriginRX_dst, OriginBY_dst, temp1 = ct.TransformPoint(OriginRX_src, OriginBY_src)
    pixl_w_dst = (OriginRX_dst - OriginLX_dst) / src_width
    pixl_h_dst = (OriginBY_dst - OriginTY_dst) / src_height
    dst_trans = [OriginLX_dst, pixl_w_dst, 0, OriginTY_dst, 0, pixl_h_dst]

    dstSRS_wkt = dstSRS.ExportToWkt()
    dst_data.SetGeoTransform(dst_trans)
    dst_data.SetProjection(dstSRS_wkt)
    data = src_data.GetRasterBand(1).ReadAsArray()
    dst_data.WriteArray(data)
    dst_data.FlushCache()


def unzip_crop(wk_dir: str, geo_box):
    #
    # lower left
    lonlat_ll = [geo_box[2], geo_box[1]]
    # up right
    lonlat_ur = [geo_box[3], geo_box[0]]

    os.chdir(wk_dir)
    file_list = glob.glob('*_*')
    file_list.sort()
    target_dir = os.path.join(pathlib.Path(os.path.abspath(
        wk_dir)).parent, pathlib.Path(os.path.abspath(wk_dir)).stem + '_clip')

    if os.path.isdir(target_dir):
        shutil.rmtree(target_dir)
    os.mkdir(target_dir)
    target_extent = (lonlat_ll[0], lonlat_ll[1], lonlat_ur[0], lonlat_ur[1])
    # Format of target_extent S,W,N,E
    # Preview of the target area
    # preview_img_path = glob.glob('%s/*.tif' % file_list[0])[0]
    # _, _, im_data, im_extent = read_img(preview_img_path)

    print(" Wait")
    # plot but blocked. can be gedited
    # used by Jupyter-Notebook to show the cut region
    """
    fig, ax = plt.subplots(figsize=[7, 7])
    ax.tick_params(direction='in', labelsize=10, length=7)
    ax.ticklabel_format(style='plain')
    ax.add_patch(Rectangle((target_extent[0:2]),
                           np.abs(target_extent[2] - target_extent[0]),
                           np.abs(target_extent[3] - target_extent[1]),
                           fill=False,
                           edgecolor='indianred',
                           lw=3))

    ax.imshow(im_data, cmap='gist_earth', extent=im_extent)
    """
    # conduct the clip operation
    for i in range(len(file_list)):
        tmp_tif_file_list = glob.glob('%s/*.tif' % file_list[i])
        tmp_target_dir = os.path.join(target_dir, file_list[i])
        if os.path.isdir(tmp_target_dir):
            shutil.rmtree(tmp_target_dir)
        os.mkdir(tmp_target_dir)
        for tmp_file in tmp_tif_file_list:
            p = pathlib.Path(os.path.abspath(tmp_file))
            tmp_output_path = os.path.join(
                tmp_target_dir, '%s_clip%s' % (p.stem, p.suffix))
            # Transfer the reference from WGS84 UTM into WGS84
            p_out = os.path.abspath(tmp_output_path)[:-9] + '_latlon.tif'
            hyp3_reprojection(os.path.abspath(p), p_out)
            _, _, im_data, im_extent = read_img(p_out)
            im_extent = crop_region_banpick(target_extent, im_extent)
            if target_extent == (-181, -91, 181, 91):
                p_rename = os.path.abspath(p_out)[:-4] + '.tif'
                os.rename(p_out, p_rename)
                # leave the generate of *_clip.tif below
            else:
                # transfer *.clip.tif  into *_clip.tif
                # to meet the requirement of prep_hyp3
                p_rename = os.path.abspath(p_out)[:-4] + '_clip.tif'
                gdal_clip(p_out, im_extent, output_path=p_rename)
                os.remove(p_out)
        # copy the txt file into
        txt_name = file_list[i] + '/' + file_list[i] + '.txt'
        txt_name1 = file_list[i] + '.txt'
        txt_path0 = os.path.join(wk_dir, txt_name)
        txt_path1 = os.path.join(tmp_target_dir, txt_name1)
        shutil.copy(txt_path0, txt_path1)
        print(':: [%s] Transferd and Clipped files in: %s' % (i + 1, os.path.abspath(file_list[i])))

    # do the clip for transformed but uncliped via hyp3lib.cutFiles
    files = []
    wk1_dir = wk_dir + '_clip'

    print("current path:  " + wk_dir)
    if os.path.exists(wk1_dir):
        filetypes = ['unw_phase', 'corr', 'dem', 'lv_theta', 'water_mask']
        for root, dir, file_list in os.walk(wk1_dir):
            for file in file_list:
                for f in filetypes:
                    if f + '_latlon.tif' in file and f + '_latlon.tif.xml' not in file:
                        full_path = os.path.join(root, file)  # Get the full path to the file
                        files.append(full_path)
        print("There are " + str(len(files)) + " satisfied files")
        # Generate
        cutFiles(files)
    # generate _clip.tif.rsc
    for i in files:
        if target_extent == (-181, -91, 181, 91):
            mintpy.prep_hyp3.main(files)
        else:
            clip_files = ['']
            clip_files.append(i[:-4] + '_clip.tif')
            mintpy.prep_hyp3.main(clip_files)


def unzip_files(work_dir, target_postfix_list):
    # return the path that contains unziped files

    os.chdir(work_dir)
    file_list = glob.glob('*.zip')
    file_list.sort()

    primary_date_list = [file_list[i][5:13] for i in range(len(file_list))]
    secondary_date_list = [file_list[i][21:29] for i in range(len(file_list))]
    extract_name_list = ['%s_%s' % (primary_date_list[i], secondary_date_list[i])
                         for i in range(len(file_list))]

    target_dir = os.path.join(os.path.dirname(pathlib.Path(work_dir)), 'unzip')
    if os.path.isdir(target_dir):
        shutil.rmtree(target_dir)
    os.mkdir(target_dir)
    for i in range(len(file_list)):
        print(':: [%s] Unpack file: %s' % (i + 1, file_list[i]))
        # make a new directory for the target unzip files
        current_target_path = os.path.join(target_dir, extract_name_list[i])
        if os.path.isdir(current_target_path):
            shutil.rmtree(current_target_path)
        os.mkdir(current_target_path)
        # unzip target files
        with ZipFile(file_list[i], 'r') as zipObject:
            post_idx = 0  # index for postfix
            zip_file_list = zipObject.namelist()
            for tmp_file in zip_file_list:
                if tmp_file.endswith(tuple(target_postfix_list)):
                    # get tmp_post name
                    # get index -> get name
                    tmp_post = np.array([tmp_file.find(tmp)
                                         for tmp in target_postfix_list])
                    tmp_post = tmp_post[tmp_post != -1][0]
                    tmp_post = tmp_file[tmp_post:]
                    # unzip step
                    zipObject.extract(tmp_file, target_dir)
                    # copy and rename target file
                    copy_path = current_target_path
                    copy_path = os.path.join(copy_path, '%s_%s' % (
                        extract_name_list[i], tmp_post))
                    source_path = glob.glob('%s/S1*/*%s' %
                                            (target_dir, tmp_post))[0]
                    shutil.copy(source_path, copy_path)

            # delete tmp file
            tmp_unzip_path = glob.glob('%s/S1*' % target_dir)[0]
            shutil.rmtree(tmp_unzip_path)
    return target_dir


def crop_region_banpick(geobox, imgbox):
    region = [0.0, 0.0, 0.0, 0.0]
    if geobox[0] == -181:
        region[0] = imgbox[2]
    if geobox[1] == -91:
        region[1] = imgbox[0]
    if geobox[2] == 181:
        region[2] = imgbox[3]
    if geobox[3] == 91:
        region[3] = imgbox[1]
    return region


def ASF_crop(work_dir, zip_state, geobox, crop_state):
    # directory for zip file
    # example: wkdir = '/media/lei/disk1/hyp3'
    # example:wkdir = '/media/lei/disk1/unzip'
    # target directory
    # example: target_dir = '../unzip'

    # target file to be extracted, another options,
    # geobox represent the n,s,w,e
    # e.g., wrapped_phase.tif, inc_map.tif, water_mask.tif,
    if geobox[:] == [91, -91, -181, 181] and crop_state == 1:
        print("the crop state code is crop but the region has not been input")
        exit()
    elif crop_state == 0:
        geobox = [91, -91, -181, 181]
    # adjusted geobox before, the bounding should no need to be adjusted after.
    target_postfix_list = ['corr.tif', 'unw_phase.tif', 'dem.tif', 'lv_theta.tif', 'water_mask.tif']

    if zip_state == 1:
        work_dir = unzip_files(work_dir, target_postfix_list)
        unzip_crop(wk_dir=work_dir, geo_box=geobox)
    elif zip_state == 0:
        unzip_crop(wk_dir=work_dir, geo_box=geobox)
    else:
        exit()

    return 0


def main():
    str_h = """Example usage: 
    ASF_crop.py --i /Examples/asf_zip/ --o /Examples/asf_unzip/ --z 0 
    ASF_crop.py --i /Examples/asf_zip/ --o /Examples/asf_unzip/ --c 1 --e -91.5  --n 1.5 
    ASF_crop.py --i /Examples/asf_zip/ --o /Examples/asf_unzip/ --c 1 --e -91.5 --w -90 --n 1.5 --s -0.75  
    
    # rely on asf-hyp3-lib
    # add a filter for --sd start-date
    # add a filter for --ed end-date
    # add a filter for area co-coverage, although this may be covered by frame and path.
    # add a if-else function for invaild area recognise.
    # to avoid some edge set by 0.0, we use impossible value for N,S,W,E. This may be improved soon.
    
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description=str_h)
    parser.add_argument('--i', type=str, help="The input folder path that contains all the zipped or unzipped datas")
    parser.add_argument('--o', type=str, default='', help="The output folder path that contains all the croped files")
    parser.add_argument('--z', type=int, default=1, help="0 = unzipped, 1 = zipped")
    parser.add_argument('--c', type=int, default=0, help="0 = do not crop, 1 = do the crop")
    parser.add_argument('--n', type=float, default=91, help="north latitude")
    parser.add_argument('--s', type=float, default=-91, help="south latitude")
    parser.add_argument('--w', type=float, default=-181, help="west latitude")
    parser.add_argument('--e', type=float, default=181, help="east latitude")
    args = parser.parse_args()

    # edge = n,s,w,e

    try:
        ASF_crop(work_dir=args.i, zip_state=args.z, geobox=[args.n, args.s, args.w, args.e], crop_state=args.c)
    except:
        print("Something wrong")
    return 0


if __name__ == '__main__':
    main()
