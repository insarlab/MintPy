#! /usr/bin/env python2
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Bhuvan Varugu, Zhang Yunjun, 2018                #
############################################################


import os
import sys
import argparse
import h5py
import numpy as np
from scipy.interpolate import griddata, RegularGridInterpolator as RGI
from mintpy.utils import ptime, readfile, utils as ut, network as pnet


##########################################################
def read_params(filename):
    fil=open(filename+'.rsc','r')
    line= fil.readline()
    rscdict = {}
    while line:
        llist = line.split()
        if len(llist)>0:
            rscdict[llist[0]] = llist[1]
        line = fil.readline()
    fil.close()
    nx=np.int(rscdict['WIDTH']);ny=np.int(rscdict['FILE_LENGTH'])
    lat=np.zeros((4,1));lon=np.zeros((4,1));
    lat[0]=np.float(rscdict['Y_FIRST']);lat[1]=np.float(rscdict['Y_FIRST'])
    lat[2]=np.float(rscdict['Y_FIRST'])+(ny-1)*np.float(rscdict['Y_STEP']);
    lat[3]=np.float(rscdict['Y_FIRST'])+(ny-1)*np.float(rscdict['Y_STEP'])
    lon[0]=np.float(rscdict['X_FIRST']);lon[1]=np.float(rscdict['X_FIRST'])+(nx-1)*np.float(rscdict['X_STEP']);
    lon[2]=np.float(rscdict['X_FIRST']);lon[3]=np.float(rscdict['X_FIRST'])+(nx-1)*np.float(rscdict['X_STEP']);
    return lat,lon,nx,ny


def get_delay(delay_file,atr,lookup_file,cinc):
    [lat,lon,nx,ny]= read_params(delay_file)
    date1=np.fromfile(delay_file, dtype=np.float32, sep=(""))
    data=date1.reshape(ny,nx)
    lats,step=np.linspace(lat[2],lat[0],num=ny,endpoint=True,retstep=True);lats=np.asarray(lats)
    lons,step=np.linspace(lon[0],lon[1],num=nx,endpoint=True,retstep=True);lons=np.asarray(lons)
    pts_old=tuple((lats,lons));del step
    RGI_func=RGI(pts_old,data,method='linear',bounds_error=False)
    latd=np.zeros((4,1));lond=np.zeros((4,1));
    rx=np.int(atr['ref_x']);ry=np.int(atr['ref_y'])
    if 'X_FIRST' in list(atr.keys()):
        nxd=np.int(atr['WIDTH']);nyd=np.int(atr['FILE_LENGTH'])
        latd[0]=np.float(atr['Y_FIRST']);latd[2]=np.float(atr['Y_FIRST'])+(nyd-1)*np.float(atr['Y_STEP']);
        lond[0]=np.float(atr['X_FIRST']);lond[1]=np.float(atr['X_FIRST'])+(nxd-1)*np.float(atr['X_STEP']);
        xarr,step = np.linspace(lond[0],lond[1],num=nxd,endpoint=True,retstep=True);
        yarr,step = np.linspace(latd[2],latd[0],num=nyd,endpoint=True,retstep=True);
        lons3=np.tile(xarr,nyd);lats3=np.repeat(yarr,nxd);
        pts_new = np.hstack((lats3.reshape(-1,1), lons3.reshape(-1,1)));
        delay_geo=RGI_func(pts_new)
        delay_geo=delay_geo.reshape(1858,1820);
        rlat=np.float(atr['ref_lat']);rlon=np.float(atr['ref_lon'])
        delay_geo=delay_geo-delay_geo[ry,rx]
        delay_geo=delay_geo/cinc;
        delay=delay_geo; del delay_geo, pts_new
        
    else:
        if atr['processor']=='isce':
            print('interferogram processor is isce')
            len_rdr = int(atr['FILE_LENGTH']);wid_rdr = int(atr['WIDTH']);
            lons=np.tile(lons,ny);lats=np.repeat(lats,nx);
            pts_old=np.hstack((lats.reshape(-1,1), lons.reshape(-1,1)));
            atr_lut = readfile.read_attribute(lookup_file)
            lon_lut = readfile.read(lookup_file, datasetName='longitude')[0]
            lat_lut = readfile.read(lookup_file, datasetName='latitude')[0]
            pts_new=np.hstack((lat_lut.reshape(-1,1), lon_lut.reshape(-1,1)));
            delay_rdr=griddata(pts_old, data.flatten(), pts_new, method='linear')
            delay_rdr=delay_rdr.reshape(len_rdr,wid_rdr)
            delay_rdr=delay_rdr-delay_rdr[ry,rx];
            delay_rdr=delay_rdr/cinc;
            delay=delay_rdr; del delay_rdr, pts_new
            
        else:
            print('interferogram processor is roi_pac')
            [latd,lond,nxd,nyd] = read_params(lookup_file); 
            xarr,step = np.linspace(lond[0],lond[1],num=nxd,endpoint=True,retstep=True);
            yarr,step = np.linspace(latd[2],latd[0],num=nyd,endpoint=True,retstep=True);
            lons3=np.tile(xarr,nyd);lats3=np.repeat(yarr,nxd);
            pts_new = np.hstack((lats3.reshape(-1,1), lons3.reshape(-1,1)));
            delay_geo=RGI_func(pts_new)
            delay_geo=delay_geo.reshape(nyd,nxd)
            len_rdr = int(atr['FILE_LENGTH']);wid_rdr = int(atr['WIDTH']);
            yy, xx = np.mgrid[0:len_rdr:len_rdr*1j, 0:wid_rdr:wid_rdr*1j]
            yx_rdr = np.hstack((yy.reshape(-1,1), xx.reshape(-1,1)))
            atr_lut = readfile.read_attribute(lookup_file)
            rg = readfile.read(lookup_file, datasetName='range')[0]
            az = readfile.read(lookup_file, datasetName='azimuth')[0]
            idx = (az>0.0)*(az<=len_rdr)*(rg>0.0)*(rg<=wid_rdr)
            pts_new2 = np.hstack((az[idx].reshape(-1,1), rg[idx].reshape(-1,1)))
            delay_geo=delay_geo[idx]
            delay_rdr=griddata(pts_new2, delay_geo, yx_rdr, method='linear')
            delay_rdr=delay_rdr.reshape(len_rdr,wid_rdr)
            delay_rdr=delay_rdr-delay_rdr[ry,rx];
            delay_rdr=delay_rdr/cinc;
            delay=delay_rdr; del delay_rdr, pts_new, pts_new2
        
    return delay


###############################################################
EXAMPLE='''example:
  tropo_gacos.py timeseries.h5 -l geomap_*rlks.trans -i incidenceAngle.h5
'''
TEMPLATE='''
mintpy.troposphericDelay.method        = GACOS   #[pyaps, height_correlation,GACOS] 
'''

def cmdLineParse():
    parser = argparse.ArgumentParser(description='Tropospheric correction using GACOS delays\n',\
                                     formatter_class=argparse.RawTextHelpFormatter,\
                                     epilog=EXAMPLE)

    parser.add_argument(dest='timeseries_file', nargs='?', help='timeseries HDF5 file, i.e. timeseries.h5')
    parser.add_argument('-i', dest='inc_angle',\
                        help='a file containing all incidence angles, or a number representing for the whole image.')
    parser.add_argument('-l', dest='lookup_file',\
                        help='a file containing all information to tranfer from radar to geo coordinates.')
    parser.add_argument('--GACOS-dir', dest='GACOS_dir', \
                        help='directory to downloaded GACOS delays data, i.e. ./../WEATHER/GACOS\n'+\
                             'use directory of input timeseries_file if not specified.')
    parser.add_argument('--date-list', dest='date_list_file',\
                        help='Read the first column of text file as list of date to download data\n'+\
                             'in YYYYMMDD or YYMMDD format')
    parser.add_argument('--ref-yx', dest='ref_yx', type=int, nargs=2, help='reference pixel in y/x')
    parser.add_argument('--template', dest='template_file',\
                        help='template file with input options below:\n'+TEMPLATE)
    parser.add_argument('-o', dest='out_file', help='Output file name for trospheric corrected timeseries.')

    inps = parser.parse_args()

    # Correcting TIMESERIES or DOWNLOAD DATA ONLY, required one of them
    if not inps.timeseries_file and not inps.download:
        parser.print_help()
        sys.exit(1)
    return inps

def main(argv):
    inps = cmdLineParse()
    
    if inps.timeseries_file:
        inps.timeseries_file=ut.get_file_list([inps.timeseries_file])[0]
        atr=readfile.read_attribute(inps.timeseries_file)
        k = atr['FILE_TYPE']
        if 'ref_y' not in list(atr.keys()) and inps.ref_yx:
            print('No reference info found in input file, use input ref_yx: '+str(inps.ref_yx))
            atr['ref_y'] = inps.ref_yx[0]
            atr['ref_x'] = inps.ref_yx[1]

    #****reading incidence angle file***/
    if os.path.isfile(inps.inc_angle):
        inps.inc_angle=readfile.read(inps.inc_angle, datasetName='incidenceAngle')[0]
        inps.inc_angle=np.nan_to_num(inps.inc_angle)
    else:
        inps.inps.inc_angle = float(inps.inc_angle)
        print('incidence angle: '+str(inps.inc_angle))
    cinc=np.cos(inps.inc_angle*np.pi/180.0);

    #****look up file****/
    if inps.lookup_file:
        inps.lookup_file = ut.get_file_list([inps.lookup_file])[0] #'geomap_32rlks_tight.trans'

    #****GACOS****/
    delay_source = 'GACOS'
    # Get weather directory
    if not inps.GACOS_dir:
        if inps.timeseries_file:
            inps.GACOS_dir = os.path.dirname(os.path.abspath(inps.timeseries_file))+'/../WEATHER/GACOS'
        elif inps.lookup_file:
            inps.GACOS_dir = os.path.dirname(os.path.abspath(inps.lookup_file))+'/../WEATHER/GACOS'
        else:
            inps.GACOS_dir = os.path.abspath(os.getcwd())
    
    print('Store weather data into directory: '+inps.GACOS_dir)
        
    #source_dir=os.path.dirname(os.path.abspath('timeseries_file'))+'/Agung/GACOS/data';print source_dir
    #os.makedirs(GACOS_dir)  -----------------------------------------------add part to copy/download weather data------#
    #----get date list-----#
    if not inps.date_list_file:
        print('read date list info from: '+inps.timeseries_file)
        h5=h5py.File(inps.timeseries_file,'r')
        if 'timeseries' in list(h5.keys()):
            date_list=sorted(h5[k].keys())
        elif k in ['interferograms','coherence','wrapped']:
            ifgram_list = sorted(h5[k].keys())
            date12_list = pnet.get_date12_list(inps.timeseries_file)
            m_dates = [i.split('-')[0] for i in date12_list]
            s_dates = [i.split('-')[1] for i in date12_list]
            date_list = ptime.yyyymmdd(sorted(list(set(m_dates + s_dates))))
        else:
            raise ValueError('Un-support input file type:'+k)
        h5.close()
    else:
        date_list = ptime.yyyymmdd(np.loadtxt(inps.date_list_file, dtype=str, usecols=(0,)).tolist())
        print('read date list info from: '+inps.date_list_file)

    #****cheacking availability of delays****/
    print('checking availability of delays')
    delay_file_list=[]
    for d in date_list:
        if   delay_source == 'GACOS':  delay_file = inps.GACOS_dir+'/'+d+'.ztd';
        delay_file_list.append(delay_file)
    delay_file_existed = ut.get_file_list(delay_file_list)

    if len(delay_file_existed)==len(date_list):
        print('no missing files')
    else:
        print('no. of date files found:', len(delay_file_existed));
        print('no. of dates:', len(date_list))

    #*****Calculating delays***/
    print('calculating delays')

    length=int(atr['FILE_LENGTH'])
    width=int(atr['WIDTH'])
    #initialise delay files
    date_num=len(date_list)
    trop_ts=np.zeros((date_num, length, width), np.float32)

    #reading wrf files for each epoch and getting delay
    for i in range(date_num):
        delay_file=delay_file_existed[i]
        date=date_list[i]
        print('calculating delay for date',date)
        trop_ts[i] =get_delay(delay_file,atr,inps.lookup_file,cinc)  


    print('Delays Calculated')
    # Convert relative phase delay on reference date
    try:    ref_date = atr['ref_date']
    except: ref_date = date_list[0]
    print('convert to relative phase delay with reference date: '+ref_date)
    ref_idx = date_list.index(ref_date)
    trop_ts -= np.tile(trop_ts[ref_idx,:,:], (date_num, 1, 1))

    ## Write tropospheric delay to HDF5
    tropFile = 'GACOSdelays'+'.h5'
    print('writing >>> %s' % (tropFile))
    h5trop = h5py.File(tropFile, 'w')
    group_trop = h5trop.create_group('timeseries')
    print('number of acquisitions: '+str(date_num))
    prog_bar = ptime.progress_bar(maxValue=date_num)
    for i in range(date_num):
        date = date_list[i]
        group_trop.create_dataset(date, data=trop_ts[i], compression='gzip')
        prog_bar.update(i+1, suffix=date)
    prog_bar.close()

    # Write Attributes
    for key,value in atr.items():
        group_trop.attrs[key] = value
    h5trop.close()

    ## Write corrected Time series to HDF5
    if k == 'timeseries':
        if not inps.out_file:
            inps.out_file = os.path.splitext(inps.timeseries_file)[0]+'_'+'GACOS'+'.h5'
        print('writing trop corrected timeseries file %s' % (inps.out_file))
        h5ts = h5py.File(inps.timeseries_file, 'r')
        h5tsCor = h5py.File(inps.out_file, 'w')
        group_tsCor = h5tsCor.create_group('timeseries')
        print('number of acquisitions: '+str(date_num))
        prog_bar = ptime.progress_bar(maxValue=date_num)
        for i in range(date_num):
            date = date_list[i];print(date)
            ts = h5ts['timeseries'].get(date)[:]
            group_tsCor.create_dataset(date, data=ts-trop_ts[i], compression='gzip')
            prog_bar.update(i+1, suffix=date)
        prog_bar.close()
        h5ts.close()
        # Write Attributes
        for key,value in atr.items():
            group_tsCor.attrs[key] = value
        h5tsCor.close()
        print('delays written to %s' % (inps.out_file))

    print('finished')
    return inps.out_file
    
###############################################################
if __name__ == '__main__':
    main(sys.argv[1:])



