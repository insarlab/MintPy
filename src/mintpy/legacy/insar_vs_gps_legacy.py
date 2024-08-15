#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Heresh Fattahi, 2013                             #
############################################################


import getopt
import os
import sys

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator


def readGPSfile(gpsFile, gps_source):
    if gps_source in ['cmm4', 'CMM4']:

        gpsData = np.loadtxt(gpsFile, usecols=(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
        Stations = np.loadtxt(gpsFile, dtype=str, usecols=(0, 1))[:, 0]

        St = []
        Lon = []
        Lat = []
        Ve = []
        Vn = []
        Se = []
        Sn = []
        for i in range(gpsData.shape[0]):
            if 'GPS' in Stations[i]:
                Lat.append(gpsData[i, 0])
                Lon.append(gpsData[i, 1])
                Ve.append(gpsData[i, 2])
                Se.append(gpsData[i, 3])
                Vn.append(gpsData[i, 4])
                Sn.append(gpsData[i, 5])
                St.append(Stations[i])
        Up = np.zeros(St.shape)
        Sup = np.zeros(St.shape)
    elif gps_source == 'mintpy':

        #gpsData = np.loadtxt(gpsFile,usecols = (1,2,3,4,5,6,7,8,9))
        gpsData = np.loadtxt(gpsFile, usecols=(1, 2, 3, 4, 5, 6))
        Stations = np.loadtxt(gpsFile, dtype=str, usecols=(0, 1))[:, 0]

        St = []
        Lon = []
        Lat = []
        Ve = []
        Vn = []
        Se = []
        Sn = []

        for i in range(gpsData.shape[0]):
            Lon.append(gpsData[i, 0]-360)
            Lat.append(gpsData[i, 1])
            Ve.append(gpsData[i, 2])
            Vn.append(gpsData[i, 3])
            Se.append(gpsData[i, 4])
            Sn.append(gpsData[i, 5])
            St.append(Stations[i])
        Up = np.zeros(St.shape)
        Sup = np.zeros(St.shape)
    elif gps_source in ['usgs', 'USGS']:
        gpsData_Hz = np.loadtxt(gpsFile, usecols=(0, 1, 2, 3, 4, 5, 6))
        gpsData_up = np.loadtxt(gpsFile, usecols=(8, 9))
        gpsData = np.hstack((gpsData_Hz, gpsData_up))
        Stations = np.loadtxt(gpsFile, dtype=str, usecols=(7, 8))[:, 0]

        St = []
        Lon = []
        Lat = []
        Ve = []
        Vn = []
        Se = []
        Sn = []
        Up = []
        Sup = []

        for i in range(gpsData.shape[0]):
            Lat.append(gpsData[i, 0])
            Lon.append(gpsData[i, 1])
            Vn.append(gpsData[i, 2])
            Ve.append(gpsData[i, 3])
            Sn.append(gpsData[i, 4])
            Se.append(gpsData[i, 5])
            Up.append(gpsData[i, 6])
            Sup.append(gpsData[i, 7])
            St.append(Stations[i])

    return list(St), Lat, Lon, Ve, Se, Vn, Sn, Up, Sup


def nearest(x, tbase, xstep):
    # """ find nearest neighbour """
    dist = np.sqrt((tbase - x)**2)
    if min(dist) <= np.abs(xstep):
        indx = dist == min(dist)
    else:
        indx = []
    return indx


def find_row_column(Lon, Lat, lon, lat, lon_step, lat_step):
    ################################################
    # finding row and column numbers of the GPS point

    idx = nearest(Lon, lon, lon_step)
    idy = nearest(Lat, lat, lat_step)
    if idx != [] and idy != []:
        IDX = np.where(idx == True)[0][0]
        IDY = np.where(idy == True)[0][0]
    else:
        IDX = np.nan
        IDY = np.nan
    return IDY, IDX


################################################
def usage():
    print("""
   ************************************************************************************************

   Compares InSAR and GPS velocities. Option G can be used to specify the mode of Comparison.
   Out put is InSARvsGPS.png

   Usage:

   insar_vs_gps.py  -v InSARvelocity.h5 -g GPS velocity file -r Name of the reference GPS station -l ststion list to be compared

   -l :  if station list is not specified, all stations in GPS velocity file are compared with InSAR

   -m min value of the x and y axis of the plot
   -M max value of the x and y axis of the plot
   -r reference GPS station
   -s second velocity map
   -S source of the GPS data: (usgs,cmm4,mintpy)
      see documentation for more information
   -I incidence angle (if not given average look angle is used instead)
   -H Heading angle (if not given then the program reads it from the attributes of the velocity file)

   -G GPS components to be used to compare with InSAR: los_3D , los_Hz , los_Up , gps_Up. [default is los_3D]
       los_3D: to project three gps components to LOS
       los_Hz: to project horizontal gps components to LOS
       los_Up: to project only vertical gps component to LOS
       gps_Up: uses vertical GPS to compare with InSAR (InSAR LOS will be projected to Up)

   -A annotate the GPS station name on the plot [yes] or no
   -C annotation font color [default is green]
   -x annotation offset from the point in x direction [default=0]
   -y annotation offset from the point in y direction [default=0]
   -u to plot 1, 2 or 3 sigma uncertainty [default=1]
   -B marker size [default = 15]
   Example:

    insar_vs_gps.py  -v geo_InSAR_velocity.h5  -g gpsVelocity.txt -S usgs  -r BEMT -l 'HNPS,DHLG,SLMS,USGC'
    insar_vs_gps.py  -v geo_InSAR_velocity.h5  -g gpsVelocity.txt -S cmm4 -r BEMT -l 'HNPS,DHLG,SLMS,USGC,ROCH,MONP,SIO3,IVCO,TMAP,BMHL,BILL,OGHS'
    insar_vs_gps.py  -v geo_InSAR_velocity.h5  -g gpsVelocity.txt -S usgs -r BEMT
    insar_vs_gps.py  -v geo_InSAR_velocity.h5  -g gpsVelocity.txt -S usgs -r BEMT -c geo_temporal_coherence.h5 -t 0.95
    insar_vs_gps.py  -v geo_InSAR_velocity.h5  -g gpsVelocity.txt -S usgs -r BEMT -c geo_temporal_coherence.h5 -t 0.95  -l 'HNPS,DHLG,SLMS,USGC'

    insar_vs_gps.py  -v geo_velocity_New_masked.h5 -g usgs_velocities_NAfixed.txt -r BEMT -S usgs -A yes -C green -x 0 -y 0.5 -H 193.0 -I 23.0

   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   When user wants to compare two InSAR velocities at GPS stations:

   Example: (to comaper geo_InSAR_velocity.h5 with simulated_velocity.h5 in GPS stations)

    insar_vs_gps.py  -v geo_InSAR_velocity.h5  -g gpsVelocity.txt  -r BEMT -s simulated_velocity.h5

   ************************************************************************************************
    """)


def main(argv):

    annotation = 'yes'
    ann_x = 0
    ann_y = 0
    annotation_Color = 'green'
    disp_velocity = 'yes'
    GPS_InSAR_dif_thr = 1
    gps_comp = 'los_3D'
    uncertainty_fac = 1.0
    MarkerSize = 5
    try:
        opts, args = getopt.getopt(
            argv, "v:r:g:G:l:c:t:m:M:s:S:A:B:C:x:y:I:H:u:")

    except getopt.GetoptError:
        usage()
        sys.exit(1)

    for opt, arg in opts:

        if opt == '-v':
            velocityFile = arg
        elif opt == '-s':
            velocityFile2 = arg
        elif opt == '-g':
            gpsFile = arg
        elif opt == '-r':
            refStation = arg
        elif opt == '-l':
            stationsList = arg.split(',')
        elif opt == '-c':
            coherenceFile = arg
        elif opt == '-t':
            thr = float(arg)
        elif opt == '-m':
            minV = float(arg)
        elif opt == '-M':
            maxV = float(arg)
        elif opt == '-S':
            gps_source = arg
        elif opt == '-A':
            annotation = arg
        elif opt == '-C':
            annotation_Color = arg
        elif opt == '-x':
            ann_x = float(arg)
        elif opt == '-y':
            ann_y = float(arg)
        elif opt == '-I':
            theta = float(arg)
        elif opt == '-H':
            heading = float(arg)
        elif opt == '-G':
            gps_comp = arg
        elif opt == '-u':
            uncertainty_fac = float(arg)
        elif opt == '-B':
            MarkerSize = float(arg)

    try:
        velocityFile
        gpsFile
        refStation
    except:
        usage()
        sys.exit(1)

    try:
        thr
    except:
        thr = 0.9

    h5file = h5py.File(velocityFile, 'r')
    dset = h5file['velocity'].get('velocity')
    insarData = dset[0:dset.shape[0], 0:dset.shape[1]]
    k = list(h5file.keys())

    try:
        h5file2 = h5py.File(velocityFile2, 'r')
        dset2 = h5file2['velocity'].get('velocity')
        insarData2 = dset2[0:dset2.shape[0], 0:dset2.shape[1]]
    except:
        print('')

    ullon = float(h5file[k[0]].attrs['X_FIRST'])
    ullat = float(h5file[k[0]].attrs['Y_FIRST'])
    lon_step = float(h5file[k[0]].attrs['X_STEP'])
    lat_step = float(h5file[k[0]].attrs['Y_STEP'])
    #lon_unit = h5file[k[0]].attrs['Y_UNIT']
    #lat_unit = h5file[k[0]].attrs['X_UNIT']

    Length, Width = np.shape(insarData)

    lllat = ullat+Length*lat_step
    urlon = ullon+Width*lon_step
    lat = np.arange(ullat, lllat, lat_step)
    lon = np.arange(ullon, urlon, lon_step)
    #################################################################################################
    # finding the raw an column of the reference gps station and referencing insar data to this pixel
    Stations, Lat, Lon, Ve, Se, Vn, Sn, Vu, Su = readGPSfile(
        gpsFile, gps_source)
    idxRef = Stations.index(refStation)
    IDYref, IDXref = find_row_column(
        Lon[idxRef], Lat[idxRef], lon, lat, lon_step, lat_step)

    #############################################
    #  Stations, gpsData = redGPSfile(gpsFile)
    #  idxRef=Stations.index(refStation)
    #  Lat,Lon,Vn,Ve,Sn,Se,Corr,Vu,Su = gpsData[idxRef,:]
    #  IDYref,IDXref=find_row_column(Lon,Lat,lon,lat,lon_step,lat_step)
    ###################################################

    if (not np.isnan(IDYref)) and (not np.isnan(IDXref)):
        print('')
        print('-----------------------------------------------------------------------')
        print('referencing InSAR data to the GPS station at : ' +
              str(IDYref) + ' , ' + str(IDXref))
        if not np.isnan(insarData[IDYref][IDXref]):
            insarData = insarData - insarData[IDYref][IDXref]
        else:
            print("""
%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       WARNING: nan value for InSAR data at the reference pixel!
                reference station should be a pixel with valid value in InSAR data.

                please select another GPS station as the reference station.
%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  """)
            sys.exit(1)
    else:
        print('WARNING:')
        print('Reference GPS station is out of the area covered by InSAR data')
        print('please select another GPS station as the reference station.')
        sys.exit(1)

    #######################################################################################

    try:
        stationsList
    except:
        stationsList = Stations
        # stationsList.remove(refStation)

    try:
        print('incidence angle = ' + str(theta))
    except:
        print('using look angle from the velocity file. For more precise results input the incidence angle using option -I.')
        look_n = float(h5file['velocity'].attrs['LOOK_REF1'])
        look_f = float(h5file['velocity'].attrs['LOOK_REF2'])
        theta = (look_n+look_f)/2.
        print('incidence angle = ' + str(theta))

    try:
        print('Heading angle = '+str(heading))
    except:
        heading = float(h5file['velocity'].attrs['HEADING'])
        if heading < 0:
            heading = heading+360

    theta = theta*np.pi/180.0
    heading = heading*np.pi/180.0

    if gps_comp in ['los_3D', 'LOS_3D', 'los_3d']:
        unitVec = [np.cos(heading)*np.sin(theta), -np.sin(theta)
                   * np.sin(heading), -np.cos(theta)]
        gps_comp_txt = ' projecting three gps components to LOS'
    elif gps_comp in ['los_Hz', 'LOS_HZ', 'los_hz', 'los_HZ', 'LOS_hz']:
        unitVec = [np.cos(heading)*np.sin(theta), -
                   np.sin(theta)*np.sin(heading), 0]
        gps_comp_txt = ' projecting horizontal gps components to LOS'
    elif gps_comp in ['LOS_UP', 'los_Up', 'los_up', 'LOS_up']:
        unitVec = [0, 0, -np.cos(theta)]
        gps_comp_txt = ' projecting vertical gps components to LOS'
    elif gps_comp in ['gps_up', 'GPS_UP', 'GPS_Up', 'gps_Up']:
        unitVec = [0, 0, 1]
        gps_comp_txt = ' comparing veryical gps with InSAR'
        print('-------------------------')
        print('Projecting InSAR to vertical')
        insarData = -insarData/np.cos(theta)
    print('-------------------------')
    print('unit vector for :' + gps_comp_txt)
    print(unitVec)
    print('-------------------------')
    gpsLOS_ref = (unitVec[0] * Ve[idxRef]
                  + unitVec[1] * Vn[idxRef]
                  + unitVec[2] * Vu[idxRef])
    Sr = ((unitVec[0]**2) * Se[idxRef]**2
          + (unitVec[1]**2) * Sn[idxRef]**2
          + (unitVec[2]**2) * Su[idxRef]**2)**0.5

    print('######################################################################')
    try:
        h5coh = h5py.File(coherenceFile)
        kh5coh = list(h5coh.keys())
        dset = h5coh[kh5coh[0]].get(kh5coh[0])
        Coh = dset[0:dset.shape[0], 0:dset.shape[1]]
    except:
        print('No information about the coherence of the points')

    InSAR = []
    GPS = []
    InSAR1 = []
    GPS1 = []
    InSAR2 = []
    GPS2 = []
    coherence = []
    GPSx = []
    GPSy = []
    GPSx1 = []
    GPSy1 = []
    GPSx2 = []
    GPSy2 = []
    GPS_station = []
    GPS_std = []
    for st in stationsList:

        try:
            idx = Stations.index(st)
            #Lat,Lon,Vn,Ve,Sn,Se,Corr,Vu,Su = gpsData[idx,:]
            gpsLOS = unitVec[0]*Ve[idx]+unitVec[1]*Vn[idx]+unitVec[2]*Vu[idx]
            Sg = ((unitVec[0]**2)*Se[idx]**2+(unitVec[1]**2)
                  * Sn[idx]**2+(unitVec[2]**2)*Su[idx]**2)**0.5
            S = (Sg**2+Sr**2)**0.5
            gpsLOS = gpsLOS-gpsLOS_ref
            IDY, IDX = find_row_column(
                Lon[idx], Lat[idx], lon, lat, lon_step, lat_step)
            insar_velocity = -insarData[IDY][IDX]

            try:
                gpsLOS = insarData2[IDY][IDX]-insarData2[IDYref][IDXref]
                gpsLOS = -1000.0*gpsLOS
            except:
                InSAR_GPS_Copmarison = 'yes'

            if not np.isnan(insarData[IDY][IDX]):
                print('%%%%%%%%%%%%%%%%%%%%')
                print(st)
                print('GPS: ' + str(gpsLOS) + '  +/- '+str(S))
                print('INSAR: ' + str(-insarData[IDY][IDX]*1000.0))
                try:
                    print('Coherence: ' + str(Coh[IDY][IDX]))
                    coherence.append(Coh[IDY][IDX])
                    if Coh[IDY][IDX] > thr:
                        InSAR1.append(-insarData[IDY][IDX]*1000.0)
                        GPS1.append(gpsLOS)
                    else:
                        InSAR2.append(-insarData[IDY][IDX]*1000.0)
                        GPS2.append(gpsLOS)

                except:
                    print('No information about the coherence is available!')

                InSAR.append(-insarData[IDY][IDX]*1000.0)
                GPS.append(gpsLOS)
                GPS_station.append(st)
                GPSx.append(IDX)
                GPSy.append(IDY)
                GPS_std.append(S)
                if np.abs(gpsLOS+insarData[IDY][IDX]*1000.0) < GPS_InSAR_dif_thr:
                    GPSx1.append(IDX)
                    GPSy1.append(IDY)
                else:
                    GPSx2.append(IDX)
                    GPSy2.append(IDY)
        except:
            NoInSAR = 'yes'

    InSAR = np.array(InSAR)
    GPS = np.array(GPS)
    GPS_std = np.array(GPS_std)
    lt = len(InSAR)
    SAD = np.sum(np.abs(InSAR-GPS), 0)/lt
    C1 = np.zeros([2, len(InSAR)])
    C1[0][:] = InSAR
    C1[1][:] = GPS
    Cor = np.corrcoef(C1)[0][1]
    print('++++++++++++++++++++++++++++++++++++++++++++++')
    print('Comparison summary:')
    print('')
    print('AAD (average absolute difference)= '+str(SAD) + ' [mm/yr]')
    print('Correlation = '+str(Cor))
    print('')
    print('++++++++++++++++++++++++++++++++++++++++++++++')
    ###############################################################
    try:
        minV
        maxV
    except:
        minV = np.min([InSAR, GPS])
        maxV = np.max([InSAR, GPS])

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.errorbar(GPS, InSAR, yerr=0.0, xerr=uncertainty_fac *
                GPS_std, fmt='ko', ms=MarkerSize)
    ax.plot([minV-3, maxV+3], [minV-3, maxV+3], 'k--')
    ax.set_ylabel('InSAR [mm/yr]', fontsize=26)
    ax.set_xlabel('GPS LOS [mm/yr]', fontsize=26)
    ax.set_ylim(minV-3, maxV+3)
    ax.set_xlim(minV-3, maxV+3)

    if annotation in ['yes', 'y', 'Y', 'Yes', 'YES']:
        for i in range(len(GPS)):
            ax.annotate(GPS_station[i], xy=(GPS[i], InSAR[i]), xytext=(
                GPS[i]+ann_x, InSAR[i]+ann_y), color=annotation_Color)

    majorLocator = MultipleLocator(5)
    ax.yaxis.set_major_locator(majorLocator)
    minorLocator = MultipleLocator(1)
    ax.yaxis.set_minor_locator(minorLocator)
    ax.xaxis.set_minor_locator(minorLocator)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(26)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(26)
    plt.tick_params(which='major', length=15, width=2)
    plt.tick_params(which='minor', length=6, width=2)

    figName = 'InSARvsGPS_errorbar.png'
    plt.savefig(figName)
    ###############################################################

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(GPS, InSAR, 'ko', ms=MarkerSize)
    ax.plot([minV-3, maxV+3], [minV-3, maxV+3], 'k--')
    # ax.plot([-10,20],[-10,20],'k--')
    ax.set_ylabel('InSAR [mm/yr]', fontsize=26)
    ax.set_xlabel('GPS LOS [mm/yr]', fontsize=26)
    ax.set_ylim(minV-3, maxV+3)
    ax.set_xlim(minV-3, maxV+3)
    # ax.set_ylim(-10,15)
    # ax.set_xlim(-10,15)
    if annotation in ['yes', 'y', 'Y', 'Yes', 'YES']:
        for i in range(len(GPS)):
            ax.annotate(GPS_station[i], xy=(GPS[i], InSAR[i]), xytext=(
                GPS[i]+ann_x, InSAR[i]+ann_y), color=annotation_Color)

    majorLocator = MultipleLocator(5)
    ax.yaxis.set_major_locator(majorLocator)
    minorLocator = MultipleLocator(1)
    ax.yaxis.set_minor_locator(minorLocator)
    ax.xaxis.set_minor_locator(minorLocator)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(26)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(26)
    plt.tick_params(which='major', length=15, width=2)
    plt.tick_params(which='minor', length=6, width=2)

    figName = 'InSARvsGPS.png'
    # plt.savefig(figName,pad_inches=0.0)
    plt.savefig(figName)

    ######################################################

    try:
        Coh
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.errorbar(GPS1, InSAR1, yerr=1.0, xerr=1.0, fmt='o')
        ax.errorbar(GPS2, InSAR2, yerr=1.0, xerr=1.0, fmt='^')
        ax.plot([minV-3, maxV+3], [minV-3, maxV+3], '--')
        ax.set_ylabel('InSAR [mm/yr]', fontsize=26)
        ax.set_xlabel('GPS LOS [mm/yr]', fontsize=26)
        ax.set_ylim(minV-3, maxV+3)
        ax.set_xlim(minV-3, maxV+3)

    except:
        print('')
    plt.show()


if __name__ == '__main__':
    main(sys.argv[1:])
