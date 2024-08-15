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


def usage():
    print("""
*****************************************************************************************
   Generating multiple profiles(each profile includes seeveral transects [specified by -n])
   perpendicular to a Fault . Fault is a path specified by lat and lon coordinates.

   Usage:
       -n number of transects used to generate one profile
       -d distance [in pixel] between individual transects to generate one profile
       -F a txt file including the fault coordinates (first column lon , second column: lat)
       -p flip profile left-right (yes or no) [default: no]
       -u flip up - down [default: no]
       -g gps_file (if exists)
       -S source of GPS velocities (usgs,cmm4,mintpy)
       -G gps stations to compare with InSAR  (all,insar,profile)
          "all": all gps stations is projected to the profile
          "insar": same as all but limited to the area covered by insar
          "profile": only those gps stations which are in the profile area]

       -x lower bound to display in x direction
       -X higher bound to display in x direction
       -l lower bound to display in y direction
       -h higher bound to display in y direction
       -I display InSAR velocity [on] or off
       -A display Average InSAR velocity [on] or off
       -U display Standard deviation of the InSAR velocity [on] or off
       -E Export the generated transect to a matlab file [off] or on
       -W Length of a profile
       -D Distance between two consequent average profile

   Example:
       multi_transect.py -f geo_velocity_masked.h5 -n 50 -d 1 -W 10 -D 2 -F Chaman_fault.txt

********************************************************************************************
    """)
    return


def dms2d(Coord):
    d, m, s = Coord.split(' ')
    d = float(d)
    m = float(m)
    s = float(s)
    d = d+m/60.+s/3600.
    return d


def gps_to_LOS(Ve, Vn, theta, heading):

    unitVec = [ np.cos(heading) * np.sin(theta),
               -np.sin(theta) * np.sin(heading),
                np.cos(theta)]

    gpsLOS = unitVec[0]*Ve+unitVec[1]*Vn

    return gpsLOS


def check_st_in_box(x, y, x0, y0, x1, y1, X0, Y0, X1, Y1):

    m1 = float(y1-y0)/float(x1-x0)
    c1 = float(y0-m1*x0)

    m2 = float(Y1-Y0)/float(X1-X0)
    c2 = float(Y0-m2*X0)

    m3 = float(y0-Y0)/float(x0-X0)
    c3 = float(Y0-m3*X0)

    m4 = float(y1-Y1)/float(x1-X1)
    c4 = float(Y1-m4*X1)

    yy1 = m1*x+c1
    yy2 = m2*x+c2

    yy = [yy1, yy2]
    yy.sort()

    xx3 = (y-c3)/m3
    xx4 = (y-c4)/m4
    xx = [xx3, xx4]
    xx.sort()
    if y >= yy[0] and y <= yy[1] and x >= xx[0] and x <= xx[1]:
        Check_result = 'True'
    else:
        Check_result = 'False'

    return Check_result


def check_st_in_box2(x, y, x0, y0, x1, y1, X0, Y0, X1, Y1):

    m1, c1 = line(x0, y0, x1, y1)
    m2, c2 = line(X0, Y0, X1, Y1)
    m3, c3 = line(x0, y0, X0, Y0)
    m4, c4 = line(x1, y1, X1, Y1)

    D1 = np.sqrt((x1-x0)**2+(y1-y0)**2)
    D2 = np.sqrt((X0-x0)**2+(Y0-y0)**2)

    d1 = dist_point_from_line(m1, c1, x, y, 1, 1)
    d2 = dist_point_from_line(m2, c2, x, y, 1, 1)

    d3 = dist_point_from_line(m3, c3, x, y, 1, 1)
    d4 = dist_point_from_line(m4, c4, x, y, 1, 1)

    if np.round(d1+d2) == np.round(D2) and np.round(d3+d4) == np.round(D1):
        Check_result = 'True'
    else:
        Check_result = 'False'

    return Check_result


def line(x0, y0, x1, y1):
    m = float(y1-y0)/float(x1-x0)
    c = float(y0-m*x0)
    return m, c


def dist_point_from_line(m, c, x, y, dx, dy):
    # finds the distance of a point at x ,y xoordinate
    # from a line with Y =  mX +c

    d = np.sqrt((((x+m*y-m*c)/(m**2+1)-x)*dx)**2 +
                ((m*(x+m*y-m*c)/(m**2+1)+c-y)*dy)**2)
    #a=m;b=-1;
    #d=np.abs(a*x+b*y+c)/np.sqrt(a**2+b**2)
    return d


def get_intersect(m, c, x, y):

    xp = (x+m*y-m*c)/(m**2+1)
    yp = m*(x+m*y-m*c)/(m**2+1)+c
    return xp, yp


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

    elif gps_source == 'mintpy':

        gpsData = np.loadtxt(gpsFile, usecols=(1, 2, 3, 4, 5, 6))  # ,7,8,9))
        Stations = np.loadtxt(gpsFile, dtype=str, usecols=(0, 1))[:, 0]

        St = []
        Lon = []
        Lat = []
        Ve = []
        Vn = []
        Se = []
        Sn = []

        for i in range(gpsData.shape[0]):
            Lon.append(gpsData[i, 0])
            Lat.append(gpsData[i, 1])
            Ve.append(gpsData[i, 2])
            Vn.append(gpsData[i, 3])
            Se.append(gpsData[i, 4])
            Sn.append(gpsData[i, 5])
            St.append(Stations[i])

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

        for i in range(gpsData.shape[0]):
            Lat.append(gpsData[i, 0])
            Lon.append(gpsData[i, 1])
            Vn.append(gpsData[i, 2])
            Ve.append(gpsData[i, 3])
            Sn.append(gpsData[i, 4])
            Se.append(gpsData[i, 5])
            St.append(Stations[i])

    return list(St), Lat, Lon, Ve, Se, Vn, Sn


def redGPSfile(gpsFile):
    gpsData_Hz = np.loadtxt(gpsFile, usecols=(0, 1, 2, 3, 4, 5, 6))
    gpsData_up = np.loadtxt(gpsFile, usecols=(8, 9))
    gpsData = np.hstack((gpsData_Hz, gpsData_up))
    Stations = np.loadtxt(gpsFile, dtype=str, usecols=(7, 8))[:, 0]
    return list(Stations), gpsData


def redGPSfile_cmm4(gpsFile):
    #gpsData = np.loadtxt(gpsFile,usecols = (1,2,3,4,5,6,7,8,9,10))
    gpsData = np.loadtxt(gpsFile, usecols=(1, 2, 3, 4, 5, 6, 7, 8, 9))
    #gpsData_up = np.loadtxt(gpsFile,usecols = (8,9))
    #gpsData=np.hstack((gpsData_Hz,gpsData_up))
    Stations = np.loadtxt(gpsFile, dtype=str, usecols=(0, 1))[:, 0]
    return list(Stations), gpsData


def nearest(x, tbase, xstep):
    """ find nearest neighbour """
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


def get_lat_lon(h5file):
    k = list(h5file.keys())
    Length = float(h5file[k[0]].attrs['LENGTH'])
    Width = float(h5file[k[0]].attrs['WIDTH'])
    ullon = float(h5file[k[0]].attrs['X_FIRST'])
    ullat = float(h5file[k[0]].attrs['Y_FIRST'])
    lon_step = float(h5file[k[0]].attrs['X_STEP'])
    lat_step = float(h5file[k[0]].attrs['Y_STEP'])

    #Length,Width = np.shape(insarData)

    lllat = ullat+(Length-1)*lat_step
    urlon = ullon+(Width-1)*lon_step
    lat = np.linspace(ullat, lllat, Length)
    lon = np.linspace(ullon, urlon, Width)

    lon_all = np.tile(lon, (Length, 1))
    lat_all = np.tile(lat, (Width, 1)).T

    #lat=np.arange(ullat,lllat,lat_step)
    #lon=np.arange(ullon,urlon,lon_step)
    return lat, lon, lat_step, lon_step, lat_all, lon_all


def nanmean(data, **args):
    return np.ma.filled(np.ma.masked_array(data, np.isnan(data)).mean(**args), fill_value=np.nan)


def nanstd(data, **args):
    return np.ma.filled(np.ma.masked_array(data, np.isnan(data)).std(**args), fill_value=np.nan)


def get_transect(z, x0, y0, x1, y1):

    length = int(np.hypot(x1-x0, y1-y0))
    x, y = np.linspace(x0, x1, length), np.linspace(y0, y1, length)
    zi = z[y.astype(int), x.astype(int)]
    return zi


def get_start_end_point(Xf0, Yf0, Xf1, Yf1, L, dx, dy):

    L = L*1000/dx
    m, c = line(Xf0, Yf0, Xf1, Yf1)
    mp = -1./m

    a = np.arctan(1./mp)
    xs = L*np.sin(a)+Xf0
    ys = L*np.cos(a)+Yf0

    xe = -L*np.sin(a)+Xf0
    ye = -L*np.cos(a)+Yf0

    return int(ys), int(xs), int(ye), int(xe)


def point_with_distance_from_line(Xf0, Yf0, Xf1, Yf1, L):
    m, c = line(Xf0, Yf0, Xf1, Yf1)
    mp = -1./m

    a = np.arctan(1./mp)
    x = L*np.sin(a)+Xf0
    y = L*np.cos(a)+Yf0
    return x, y


def point_on_line_with_distance_from_beginning(Xf0, Yf0, Xf1, Yf1, L):
    m, c = line(Xf0, Yf0, Xf1, Yf1)

    a = np.arctan(1./m)
    x = L*np.sin(a)+Xf0
    y = L*np.cos(a)+Yf0
    return x, y


def read_fault_coords(Fault_coord_file, Dp):

    Coords = np.loadtxt(Fault_coord_file)
    Dp = Dp*180.0/np.pi/6375.0
    Np = Coords.shape[0]

    Fault_lon = []
    Fault_lat = []
    Fault_lon.append(Coords[0, 0])
    Fault_lat.append(Coords[0, 1])

    for i in range(Np-1):
        # Fault_lon.append(Coords[i,0])
        # Fault_lat.append(Coords[i,1])
        Xf0 = Coords[i, 0]
        Yf0 = Coords[i, 1]
        Xf1 = Coords[i+1, 0]
        Yf1 = Coords[i+1, 1]
        DY = Yf1-Yf0
        DX = Xf1-Xf0
        Fault_Segment = np.hypot(DX, DY)

        if Fault_Segment <= Dp:

            Fault_lon.append(Coords[i+1, 0])
            Fault_lat.append(Coords[i+1, 1])

        else:
            L = Dp
            ii = 0
            while Fault_Segment > L:
                ii = ii+1
                x, y = point_on_line_with_distance_from_beginning(
                    Xf0, Yf0, Xf1, Yf1, Dp*ii)
                Dx = y-Yf0
                Dy = x-Xf0
                L = np.hypot(Dx, Dy)
                if L <= Fault_Segment:
                    Fault_lon.append(x)
                    Fault_lat.append(y)
                print(ii)
    return Fault_lon, Fault_lat


#####################################################################
def main(argv=None):
    ntrans = 1
    save_to_mat = 'off'
    flip_profile = 'no'
    which_gps = 'all'
    flip_updown = 'yes'
    incidence_file = 'incidence_file'
    display_InSAR = 'on'
    display_Average = 'on'
    display_Standard_deviation = 'on'

    try:
        opts, args = getopt.getopt(
            argv, "f:s:e:n:d:g:l:h:r:L:F:p:u:G:S:i:I:A:U:E:D:W:x:X:")

    except getopt.GetoptError:
        usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-f':
            velocityFile = arg
        elif opt == '-s':
            pnt1 = arg.split(',')
            y0 = int(pnt1[0])
            x0 = int(pnt1[1])
        elif opt == '-e':
            pnt2 = arg.split(',')
            y1 = int(pnt2[0])
            x1 = int(pnt2[1])
        elif opt == '-n':
            ntrans = int(arg)
        elif opt == '-d':
            dp = float(arg)
        elif opt == '-g':
            gpsFile = arg
        elif opt == '-r':
            refStation = arg
        elif opt == '-i':
            incidence_file = arg
        elif opt == '-L':
            stationsList = arg.split(',')
        elif opt == '-F':
            Fault_coord_file = arg
        elif opt == '-p':
            flip_profile = arg
        elif opt == '-u':
            flip_updown = arg
            print(flip_updown)
        elif opt == '-G':
            which_gps = arg
        elif opt == '-S':
            gps_source = arg
        elif opt == '-l':
            lbound = float(arg)
        elif opt == '-I':
            display_InSAR = arg
        elif opt == '-A':
            display_Average = arg
        elif opt == '-U':
            display_Standard_deviation = arg
        elif opt == '-E':
            save_to_mat = arg
        elif opt == '-h':
            hbound = float(arg)
        elif opt == '-D':
            Dp = float(arg)
        elif opt == '-W':
            profile_Length = float(arg)
        elif opt == '-x':
            x_lbound = float(arg)
        elif opt == '-X':
            x_hbound = float(arg)

    try:
        h5file = h5py.File(velocityFile, 'r')
    except:
        usage()
        sys.exit(1)

    k = list(h5file.keys())
    dset = h5file[k[0]].get(k[0])
    z = dset[0:dset.shape[0], 0:dset.shape[1]]
    dx = float(h5file[k[0]].attrs['X_STEP'])*6375000.0*np.pi/180.0
    dy = float(h5file[k[0]].attrs['Y_STEP'])*6375000.0*np.pi/180.0

    #############################################################################

    try:
        lat, lon, lat_step, lon_step, lat_all, lon_all = get_lat_lon(h5file)
    except:
        print('radar coordinate')

    Fault_lon, Fault_lat = read_fault_coords(Fault_coord_file, Dp)

    # Fault_lon=[66.40968453947265,66.36000186563085,66.31103920134248]
    # Fault_lat=[30.59405079532564,30.51565960186412,30.43928430936202]

    Num_profiles = len(Fault_lon)-1
    print('*********************************************')
    print('*********************************************')
    print('Number of profiles to be generated: '+str(Num_profiles))
    print('*********************************************')
    print('*********************************************')

    for Np in range(Num_profiles):
        FaultCoords = [Fault_lat[Np], Fault_lon[Np],
                       Fault_lat[Np+1], Fault_lon[Np+1]]
        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        print('')
        print('Profile '+str(Np) + ' [of total '+str(Num_profiles)+']')
        print('')

        try:
            #  Lat0 = dms2d(FaultCoords[0]); Lon0 = dms2d(FaultCoords[1])
            #  Lat1 = dms2d(FaultCoords[2]); Lon1 = dms2d(FaultCoords[3])

            Lat0 = FaultCoords[0]
            Lon0 = FaultCoords[1]
            Lat1 = FaultCoords[2]
            Lon1 = FaultCoords[3]
            Length, Width = np.shape(z)
            Yf0, Xf0 = find_row_column(
                Lon0, Lat0, lon, lat, lon_step, lat_step)
            Yf1, Xf1 = find_row_column(
                Lon1, Lat1, lon, lat, lon_step, lat_step)

            print('*********************************************')
            print(' Fault Coordinates:')
            print('   --------------------------  ')
            print('    Lat          Lon')
            print(str(Lat0) + ' , ' + str(Lon0))
            print(str(Lat1) + ' , ' + str(Lon1))
            print('   --------------------------  ')
            print('    row          column')
            print(str(Yf0) + ' , ' + str(Xf0))
            print(str(Yf1) + ' , ' + str(Xf1))
            print('*********************************************')
            #mf=float(Yf1-Yf0)/float((Xf1-Xf0))  # slope of the fault line
            #cf=float(Yf0-mf*Xf0)   # intercept of the fault line
            #df0=dist_point_from_line(mf,cf,x0,y0,1,1) #distance of the profile start point from the Fault line
            #df1=dist_point_from_line(mf,cf,x1,y1,1,1) #distance of the profile end point from the Fault line

            #mp=-1./mf  # slope of profile which is perpendicualr to the fault line
            # correcting the end point of the profile to be on a line perpendicular to the Fault
            #x1=int((df0+df1)/np.sqrt(1+mp**2)+x0)
            #y1=int(mp*(x1-x0)+y0)

        except:
            print('*********************************************')
            print('No information about the Fault coordinates!')
            print('*********************************************')

        #############################################################################
        y0, x0, y1, x1 = get_start_end_point(
            Xf0, Yf0, Xf1, Yf1, profile_Length, dx, dy)

        try:
            x0
            y0
            x1
            y1
        except:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.imshow(z)
            try:
                ax.plot([Xf0, Xf1], [Yf0, Yf1], 'k-')
            except:
                print('Fault line is not specified')

            xc = []
            yc = []
            print('please click on start and end point of the desired profile')

            def onclick(event):
                if event.button == 1:
                    print('click')
                    xc.append(int(event.xdata))
                    yc.append(int(event.ydata))
            cid = fig.canvas.mpl_connect('button_press_event', onclick)
            plt.show()
            x0 = xc[0]
            x1 = xc[1]
            y0 = yc[0]
            y1 = yc[1]

        ##############################################################################
        print('******************************************************')
        print('First profile coordinates:')
        print('Start point:  y = '+str(y0) + ',x = ' + str(x0))
        print('End point:   y = ' + str(y1) + '  , x = '+str(x1))
        print('')
        print(str(y0) + ',' + str(x0))
        print(str(y1) + ',' + str(x1))
        print('******************************************************')
        length = int(np.hypot(x1-x0, y1-y0))
        x, y = np.linspace(x0, x1, length), np.linspace(y0, y1, length)
        zi = z[y.astype(int), x.astype(int)]
        try:
            lat_transect = lat_all[y.astype(int), x.astype(int)]
            lon_transect = lon_all[y.astype(int), x.astype(int)]
        except:
            lat_transect = 'Nan'
            lon_transect = 'Nan'

        # zi=get_transect(z,x0,y0,x1,y1)

        try:
            dx = float(h5file[k[0]].attrs['X_STEP'])*6375000.0*np.pi/180.0
            dy = float(h5file[k[0]].attrs['Y_STEP'])*6375000.0*np.pi/180.0
            DX = (x-x0)*dx
            DY = (y-y0)*dy
            D = np.hypot(DX, DY)
            print('geo coordinate:')
            print('profile length = ' + str(D[-1]/1000.0) + ' km')
            #df0_km=dist_point_from_line(mf,cf,x0,y0,dx,dy)
        except:
            dx = float(h5file[k[0]].attrs['RANGE_PIXEL_SIZE'])
            dy = float(h5file[k[0]].attrs['AZIMUTH_PIXEL_SIZE'])
            DX = (x-x0)*dx
            DY = (y-y0)*dy
            D = np.hypot(DX, DY)
            print('radar coordinate:')
            print('profile length = ' + str(D[-1]/1000.0) + ' km')
            #df0_km=dist_point_from_line(mf,cf,x0,y0,dx,dy)

        try:
            mf, cf = line(Xf0, Yf0, Xf1, Yf1)
            df0_km = dist_point_from_line(mf, cf, x0, y0, dx, dy)
        except:
            print('Fault line is not specified')

        transect = np.zeros([len(D), ntrans])
        transect[:, 0] = zi
        XX0 = []
        XX1 = []
        YY0 = []
        YY1 = []
        XX0.append(x0)
        XX1.append(x1)
        YY0.append(y0)
        YY1.append(y1)

        if ntrans > 1:

            m = float(y1-y0)/float(x1-x0)
            c = float(y0-m*x0)
            m1 = -1.0/m
            try:
                dp
            except:
                dp = 1.0
            if lat_transect == 'Nan':
                for i in range(1, ntrans):

                    X0 = i*dp/np.sqrt(1+m1**2)+x0
                    Y0 = m1*(X0-x0)+y0
                    X1 = i*dp/np.sqrt(1+m1**2)+x1
                    Y1 = m1*(X1-x1)+y1
                    zi = get_transect(z, X0, Y0, X1, Y1)
                    transect[:, i] = zi
                    XX0.append(X0)
                    XX1.append(X1)
                    YY0.append(Y0)
                    YY1.append(Y1)
            else:
                transect_lat = np.zeros([len(D), ntrans])
                transect_lat[:, 0] = lat_transect
                transect_lon = np.zeros([len(D), ntrans])
                transect_lon[:, 0] = lon_transect

                for i in range(1, ntrans):

                    X0 = i*dp/np.sqrt(1+m1**2)+x0
                    Y0 = m1*(X0-x0)+y0
                    X1 = i*dp/np.sqrt(1+m1**2)+x1
                    Y1 = m1*(X1-x1)+y1
                    zi = get_transect(z, X0, Y0, X1, Y1)
                    lat_transect = get_transect(lat_all, X0, Y0, X1, Y1)
                    lon_transect = get_transect(lon_all, X0, Y0, X1, Y1)
                    transect[:, i] = zi
                    transect_lat[:, i] = lat_transect
                    transect_lon[:, i] = lon_transect
                    XX0.append(X0)
                    XX1.append(X1)
                    YY0.append(Y0)
                    YY1.append(Y1)

        #############################################
        try:
            m_prof_edge, c_prof_edge = line(XX0[0], YY0[0], XX0[-1], YY0[-1])
        except:
            print('Plotting one profile')

        ###############################################################################
        if flip_profile == 'yes':
            transect = np.flipud(transect)
            try:
                df0_km = np.max(D)-df0_km
            except:
                print('')

        print('******************************************************')
        try:
            gpsFile
        except:
            gpsFile = 'Nogps'
        print('GPS velocity file:')
        print(gpsFile)
        print('*******************************************************')
        if os.path.isfile(gpsFile):
            insarData = z
            del z
            Stations, Lat, Lon, Ve, Se, Vn, Sn = readGPSfile(gpsFile, gps_source)
            idxRef = Stations.index(refStation)
            Length, Width = np.shape(insarData)
            lat, lon, lat_step, lon_step = get_lat_lon(h5file)
            IDYref, IDXref = find_row_column(
                Lon[idxRef], Lat[idxRef], lon, lat, lon_step, lat_step)
            if (not np.isnan(IDYref)) and (not np.isnan(IDXref)):
                print('referencing InSAR data to the GPS station at : ' +
                      str(IDYref) + ' , ' + str(IDXref))
                if not np.isnan(insarData[IDYref][IDXref]):
                    transect = transect - insarData[IDYref][IDXref]
                    insarData = insarData - insarData[IDYref][IDXref]

                else:
                    print("""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      WARNING: nan value for InSAR data at the reference pixel!
               reference station should be a pixel with valid value in InSAR data.

               please select another GPS station as the reference station.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   """)
                    sys.exit(1)
            else:
                print('WARNING:')
                print('Reference GPS station is out of the area covered by InSAR data')
                print('please select another GPS station as the reference station.')
                sys.exit(1)

            try:
                stationsList
            except:
                stationsList = Stations

            # theta=23.0*np.pi/180.0
            if os.path.isfile(incidence_file):
                print('Using exact look angle for each pixel')
                h5file_theta = h5py.File(incidence_file, 'r')
                dset = h5file_theta['mask'].get('mask')
                theta = dset[0:dset.shape[0], 0:dset.shape[1]]
                theta = theta*np.pi/180.0
            else:
                print('Using average look angle')
                theta = np.ones(np.shape(insarData))*23.0*np.pi/180.0

            heading = 193.0*np.pi/180.0

            #unitVec=[-np.sin(theta)*np.sin(heading),-np.cos(heading)*np.sin(theta),-np.cos(theta)]
            # -np.cos(theta)]
            unitVec = [np.cos(heading)*np.sin(theta), -
                       np.sin(theta)*np.sin(heading), 0]

            #  [0.0806152480932643, 0.34918300221540616, -0.93358042649720174]
            # print unitVec
            # unitVec=[0.3,-0.09,0.9]
            # unitVec=[-0.3,0.09,-0.9]
            # unitVec=[-0.3,0.09,0]

            # print '*******************************************'
            # print 'unit vector to project GPS to InSAR LOS:'
            # print unitVec
            # print '*******************************************'
            # gpsLOS_ref=unitVec[0]*Ve[idxRef]+unitVec[1]*Vn[idxRef]#+unitVec[2]*Vu[idxRef]

            gpsLOS_ref = gps_to_LOS(
                Ve[idxRef], Vn[idxRef], theta[IDYref, IDXref], heading)
            print('%%%%%%^^^^^^^%%%%%%%%')
            print(gpsLOS_ref/1000.0)
            # insarData=insarData -gpsLOS_ref/1000.0
            # transect = transect -gpsLOS_ref/1000.0

            GPS = []
            GPS_station = []
            GPSx = []
            GPSy = []
            GPS_lat = []
            GPS_lon = []
            for st in stationsList:
                try:
                    idx = Stations.index(st)

                    #gpsLOS = unitVec[0]*Ve[idx]+unitVec[1]*Vn[idx]#+unitVec[2]*Vu[idx]

                    #gpsLOS = gps_to_LOS(Ve[idx],Vn[idx],theta[idx],heading)
                    #gpsLOS=gpsLOS-gpsLOS_ref

                    IDY, IDX = find_row_column(
                        Lon[idx], Lat[idx], lon, lat, lon_step, lat_step)
                    print(theta[IDY, IDX])
                    gpsLOS = gps_to_LOS(
                        Ve[idx], Vn[idx], theta[IDY, IDX], heading)
                    #gpsLOS = gpsLOS-gpsLOS_ref

                    if which_gps == 'all':
                        if theta[IDY, IDX] != 0.0:
                            GPS.append(gpsLOS-gpsLOS_ref)
                            GPS_station.append(st)
                            GPSx.append(IDX)
                            GPSy.append(IDY)
                            GPS_lat.append(Lat[idx])
                            GPS_lon.append(Lon[idx])
                    elif not np.isnan(insarData[IDY][IDX]):
                        if theta[IDY, IDX] != 0.0:
                            GPS.append(gpsLOS-gpsLOS_ref)
                            GPS_station.append(st)
                            GPSx.append(IDX)
                            GPSy.append(IDY)
                            GPS_lat.append(Lat[idx])
                            GPS_lon.append(Lon[idx])
                except:
                    NoInSAR = 'yes'

            DistGPS = []
            GPS_in_bound = []
            GPS_in_bound_st = []
            GPSxx = []
            GPSyy = []
            for i in range(len(GPS_station)):
                gx = GPSx[i]
                gy = GPSy[i]

                if which_gps in ['all', 'insar']:
                    check_result = 'True'
                else:
                    check_result = check_st_in_box(
                        gx, gy, x0, y0, x1, y1, X0, Y0, X1, Y1)

                if check_result == 'True':
                    check_st_in_box2(gx, gy, x0, y0, x1, y1, X0, Y0, X1, Y1)
                    GPS_in_bound_st.append(GPS_station[i])
                    GPS_in_bound.append(GPS[i])
                    GPSxx.append(GPSx[i])
                    GPSyy.append(GPSy[i])
                    # gy=y0+1
                    # gx=x0+1
                    # gxp,gyp=get_intersect(m,c,gx,gy)
                    # Dx=dx*(gx-gxp);Dy=dy*(gy-gyp)
                    # print gxp
                    # print gyp
                    # distance of GPS station from the first profile line
                    dg = dist_point_from_line(m, c, gx, gy, 1, 1)
                    # DistGPS.append(np.hypot(Dx,Dy))
                    # X0=dg/np.sqrt(1+m1**2)+x0
                    # Y0=m1*(X0-x0)+y0
                    # DistGPS.append(np.hypot(dx*(gx-X0), dy*(gy-Y0)))

                    DistGPS.append(dist_point_from_line(
                        m_prof_edge, c_prof_edge, GPSx[i], GPSy[i], dx, dy))

            print('****************************************************')
            print('GPS stations in the profile area:')
            print(GPS_in_bound_st)
            print('****************************************************')
            GPS_in_bound = np.array(GPS_in_bound)
            DistGPS = np.array(DistGPS)
            #axes[1].plot(DistGPS/1000.0, -1*GPS_in_bound/1000, 'bo')

        if gpsFile == 'Nogps':

            insarData = z
            GPSxx = []
            GPSyy = []
            GPSx = []
            GPSy = []
            GPS = []
            XX0[0] = x0
            XX1[0] = x1
            YY0[0] = y0
            YY1[0] = y1

        print('****************')
        print('flip up-down')
        print(flip_updown)

        if flip_updown == 'yes' and gpsFile != 'Nogps':
            print('Flipping up-down')
            transect = -1*transect
            GPS_in_bound = -1*GPS_in_bound
        elif flip_updown == 'yes':
            print('Flipping up-down')
            transect = -1*transect

        if flip_profile == 'yes' and gpsFile != 'Nogps':

            GPS = np.flipud(GPS)
            GPS_in_bound = np.flipud(GPS_in_bound)
            DistGPS = np.flipud(max(D)-DistGPS)

        fig, axes = plt.subplots(nrows=2)
        axes[0].imshow(insarData)
        for i in range(ntrans):
            axes[0].plot([XX0[i], XX1[i]], [YY0[i], YY1[i]], 'r-')

        axes[0].plot(GPSx, GPSy, 'b^')
        axes[0].plot(GPSxx, GPSyy, 'k^')
        if gpsFile != 'Nogps':
            axes[0].plot(IDXref, IDYref, 'r^')
        axes[0].axis('image')
        axes[1].plot(D/1000.0, transect, 'ko', ms=1)

        avgInSAR = np.array(nanmean(transect, axis=1))
        stdInSAR = np.array(nanstd(transect, axis=1))

        # std=np.std(transect,1)
        # axes[1].plot(D/1000.0, avgInSAR, 'r-')
        try:
            axes[1].plot(DistGPS/1000.0, -1*GPS_in_bound/1000, 'b^', ms=10)
        except:
            print('')
        # pl.fill_between(x, y-error, y+error,alpha=0.6, facecolor='0.20')
        # print transect
        #############################################################################

        fig2, axes2 = plt.subplots(nrows=1)
        axes2.imshow(insarData)
        # for i in range(ntrans):
        axes2.plot([XX0[0], XX1[0]], [YY0[0], YY1[0]], 'k-')
        axes2.plot([XX0[-1], XX1[-1]], [YY0[-1], YY1[-1]], 'k-')
        axes2.plot([XX0[0], XX0[-1]], [YY0[0], YY0[-1]], 'k-')
        axes2.plot([XX1[0], XX1[-1]], [YY1[0], YY1[-1]], 'k-')

        try:
            axes2.plot([Xf0, Xf1], [Yf0, Yf1], 'k-')
        except:
            FaultLine = 'None'

        axes2.plot(GPSx, GPSy, 'b^')
        axes2.plot(GPSxx, GPSyy, 'k^')
        if gpsFile != 'Nogps':
            axes2.plot(IDXref, IDYref, 'r^')
        axes2.axis('image')

        figName = 'transect_area_'+str(Np)+'.png'
        print('writing '+figName)
        plt.savefig(figName)

        #############################################################################
        fig = plt.figure()
        fig.set_size_inches(10, 4)
        ax = plt.Axes(fig, [0., 0., 1., 1.], )
        ax = fig.add_subplot(111)
        if display_InSAR in ['on', 'On', 'ON']:
            ax.plot(D/1000.0, transect*1000, 'o',
                    ms=1, mfc='Black', linewidth='0')


        ############################################################################
        # save the profile data:
        if save_to_mat in ['ON', 'on', 'On']:
            import scipy.io as sio
            matFile = 'transect'+str(Np)+'.mat'
            dataset = {}
            dataset['datavec'] = transect
            try:
                dataset['lat'] = transect_lat
                dataset['lon'] = transect_lon
            except:
                dataset['lat'] = 'Nan'
                dataset['lon'] = 'Nan'
            dataset['Unit'] = 'm'
            dataset['Distance_along_profile'] = D
            print('*****************************************')
            print('')
            print('writing transect to >>> '+matFile)
            sio.savemat(matFile, {'dataset': dataset})
            print('')
            print('*****************************************')

        #############################################################################
        if display_Standard_deviation in ['on', 'On', 'ON']:

            for i in np.arange(0.0, 1.01, 0.01):
                # ,color='#DCDCDC')#'LightGrey')
                ax.plot(D/1000.0, (avgInSAR-i*stdInSAR) *
                        1000, '-', color='#DCDCDC', alpha=0.5)
            for i in np.arange(0.0, 1.01, 0.01):
                ax.plot(D/1000.0, (avgInSAR+i*stdInSAR)*1000, '-',
                        color='#DCDCDC', alpha=0.5)  # 'LightGrey')
        #############################################################################
        if display_Average in ['on', 'On', 'ON']:
            ax.plot(D/1000.0, avgInSAR*1000, 'r-')
        ###########
        try:
            ax.plot(DistGPS/1000.0, -1*GPS_in_bound, '^', ms=10, mfc='Cyan')
        except:
            print('')
        ax.set_ylabel('LOS velocity [mm/yr]', fontsize=26)
        ax.set_xlabel('Distance along profile [km]', fontsize=26)

        ###################################################################
        # lower and higher bounds for displaying the profile

        try:
            lbound
            hbound
        except:
            lbound = np.nanmin(transect)*1000
            hbound = np.nanmax(transect)*1000

        ###################################################################
        # To plot the Fault location on the profile
        ax.plot([df0_km/1000.0, df0_km/1000.0], [lbound, hbound],
                '--', color='black', linewidth='2')

        ###################################################################

        try:
            ax.set_ylim(lbound, hbound)
        except:
            ylim = 'no'

        try:
            ax.set_xlim(x_lbound, x_hbound)
        except:
            xlim = 'no'


        ##########
        # Temporary To plot DEM
        # majorLocator = MultipleLocator(5)
        # ax.yaxis.set_major_locator(majorLocator)
        # minorLocator   = MultipleLocator(1)
        # ax.yaxis.set_minor_locator(minorLocator)

        # plt.tick_params(which='major', length=15,width=2)
        # plt.tick_params(which='minor', length=6,width=2)

        # try:
        #    for tick in ax.xaxis.get_major_ticks():
        #             tick.label.set_fontsize(26)
        #    for tick in ax.yaxis.get_major_ticks():
        #             tick.label.set_fontsize(26)
        #
        #    plt.tick_params(which='major', length=15,width=2)
        #    plt.tick_params(which='minor', length=6,width=2)
        # except:
        #    print 'couldn not fix the ticks! '

        figName = 'transect_'+str(Np)+'.png'
        print('writing '+figName)
        plt.savefig(figName)
        print('')
        print('________________________________')
        # plt.show()


#############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
