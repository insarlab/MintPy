############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Antonio Valentino, Aug 2022                      #
############################################################


import sys
from mintpy.utils.arg_utils import create_argument_parser


################################################################################
EXAMPLE = '''examples:
    lookup_geo2radar.py geometryGeo.h5 
    lookup_geo2radar.py geometryGeo.h5 -w geometryRadar.h5 
    lookup_geo2radar.py geometryGeo.h5 -w geometryRadar.h5 --parallel 4
'''

def create_parser(subparsers=None):
    synopsis = 'Convert lookup table from geo-coord (GAMMA, ROI_PAC) into radar-coord (ISCE)'
    epilog = EXAMPLE
    name = __name__.split('.')[-1]
    parser = create_argument_parser(
        name, synopsis=synopsis, description=synopsis, epilog=epilog, subparsers=subparsers)

    parser.add_argument('geometryGeo',help='geometryGeo file which includes geo-coordinates based lookup-table')
    parser.add_argument('-w','--write', dest='write', metavar='FILE', default = 'geometryRadar.h5',
                      help='update geometryRadar.h5 file by adding the radar-coordinates based lookup-table.')
    parser.add_argument('--parallel', dest='parallelNumb', type=int, metavar='NUM',default = 1,
                      help='Enable parallel processing and specify the the used processor number.[default: 1]')

    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    return inps


################################################################################        
def main(iargs=None):
    import numpy as np
    from ..utils import readfile
    from ..lookup_geo2radar import split_range, split_box, parallel_process, get_dataNames, write_h5

    inps = cmd_line_parse(iargs) 
    geom = inps.geometryGeo
    rangeCoord = readfile.read(geom,datasetName = 'rangeCoord')[0]
    azimuthCoord = readfile.read(geom,datasetName = 'azimuthCoord')[0]
    rangeCoord = rangeCoord.astype(np.float64)
    azimuthCoord = azimuthCoord.astype(np.float64)
    #CPX_lt =complex(rangeCoord + '+' + azimuthCoord+'j')
    #CPX_lt = rangeCoord  + 1j *azimuthCoord

    meta_geo = readfile.read_attribute(geom)
    post_Lat = meta_geo['Y_STEP']
    post_Lon = meta_geo['X_STEP']
    Corner_LAT = meta_geo['Y_FIRST']
    Corner_LON = meta_geo['X_FIRST']

    if inps.write:
        meta = readfile.read_attribute(inps.write)
    elif inps.reference:
        meta = readfile.read_attribute(inps.reference)
    else:
        print('write_file or the reference_file should be provided at least one.')
        sys.exit(1)

    WIDTH_geo  = int(meta_geo['WIDTH'])
    LENGTH_geo  = int(meta_geo['LENGTH'])

    x = np.arange(0,WIDTH_geo)
    y = np.arange(0,LENGTH_geo)
    xv, yv = np.meshgrid(x, y)

    LAT = float(Corner_LAT) + yv*float(post_Lat)
    LON = float(Corner_LON) + xv*float(post_Lon)
    LAT = LAT.flatten()
    LON = LON.flatten() 

    WIDTH  = int(meta['WIDTH'])
    LENGTH  = int(meta['LENGTH'])

    xx0 = rangeCoord.flatten()
    yy0 = azimuthCoord.flatten()

    zz01 = LAT.flatten()
    zz02 = LON.flatten()

    xx = xx0[xx0!=0]
    yy = yy0[xx0!=0]
    zz1 = zz01[xx0!=0] #lat 
    zz2 = zz02[xx0!=0] # lon

    #points = (xx,yy)
    #points = np.zeros((len(xx),2))
    #points[:,0] = xx
    #points[:,1] = yy

    x = np.arange(0,WIDTH)
    y = np.arange(0,LENGTH)
    grid_x, grid_y = np.meshgrid(x, y)

    row_sample = 10
    col_sample = 10

    list_row = split_range(LENGTH, row_sample)
    list_col = split_range(WIDTH, col_sample)

    split_grid_y = split_box(grid_y,row_sample,col_sample)
    split_grid_x = split_box(grid_x,row_sample,col_sample)

    data_parallel = []
    for i, (ay, ax) in enumerate(zip(split_grid_y, split_grid_x)):
        # extend the search area by 5 pixels
        max_ax = max(ax.flatten()) + 5 
        min_ax = min(ax.flatten()) - 5
        max_ay = max(ay.flatten()) + 5
        min_ay = min(ay.flatten()) - 5

        f0 = np.where((min_ax < xx) & (xx < max_ax) & (min_ay < yy) & (yy < max_ay))
        xx0 = xx[f0]
        yy0 = yy[f0]
        zz10 = zz1[f0]
        zz20 = zz2[f0]

        points0 = np.zeros((len(xx0),2))
        points0[:,0] = xx0
        points0[:,1] = yy0
        #print(split_grid_x[i].shape)

        data0 = (points0, zz10, zz20, ax, ay)
        data_parallel.append(data0)

    #grid_lat_all = []
    #grid_lon_all = []

    grid_lat = np.zeros((LENGTH,WIDTH), dtype=np.float32)
    grid_lon = np.zeros((LENGTH,WIDTH), dtype=np.float32)

    proNumb = inps.parallelNumb
    future = np.zeros((len(data_parallel),))
    future = list(future)
    future = parallel_process(data_parallel, function, n_jobs= proNumb, use_kwargs=False, front_num=1)

    for i in range(row_sample):
        for j in range(col_sample):
            k0 = i*col_sample + j
            kk = future[k0]
            y0 = min(list_row[i])
            y1 = max(list_row[i])
            x0 = min(list_col[j])
            x1 = max(list_col[j])
            #print(kk)
            try:
                lat0 = kk[0]
                lon0 = kk[1]
                grid_lat[y0:y1+1,x0:x1+1] = lat0
                grid_lon[y0:y1+1,x0:x1+1] = lon0
            except Exception as e:
                del e

    #grid_lat = griddata(points, zz1, (grid_x, grid_y), method='nearest')
    #grid_lon = griddata(points, zz2, (grid_x, grid_y), method='nearest')

    dataNames = get_dataNames(inps.write)
    datasetDict = dict()
    meta = readfile.read_attribute(inps.write)
    for k0 in dataNames:
        datasetDict[k0] = readfile.read(inps.write,datasetName = k0)[0]

    DEM = readfile.read(inps.write,datasetName = 'height')[0]
    grid_lat[DEM==0] = 0
    grid_lon[DEM==0] = 0
    grid_lat[grid_lat==0] = 'nan'
    grid_lon[grid_lon==0] = 'nan'
    datasetDict['latitude'] = grid_lat
    datasetDict['longitude'] = grid_lon
    write_h5(datasetDict, inps.write, metadata=meta, ref_file=None, compression=None)
    print('done.')


##############################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
