#!/usr/bin/env python3
############################################################
# Program is used for extract accumulated displacement of period
# Author: Lv Xiaoran                                       #
# Created: August 2019                                     #
############################################################

import os
import argparse

import mintpy
import mintpy.workflow  #dynamic import for modules used by pysarApp workflow
from mintpy.objects import sensor
from mintpy.utils import readfile, writefile
from mintpy.objects import timeseries

######################################################################################
EXAMPLE = """example:
  save_geodmod.py timeseries_ECMWF_demErr.h5 -b 34.2 35.2 45.0 46.3 -y 0.001 -x 0.001 -startDate 20171117 -endDate 20180603 -outdir $MODELDIR 
  save_geodmod.py $SCRATCHDIR/DarbandikhanSenAT72/PYSAR/INPUTS/ifgramStack.h5 20180603_20180615 -b 34.2 35.2 45.0 46.3 -y 0.001 -x 0.001 -t cor -o geo_20180603_20180615.cor -outdir $SCRATCHDIR/preGeodmod/
  save_geodmod.py $SCRATCHDIR/DarbandikhanSenAT72/PYSAR/temporalCoherence.h5 20171117_20180603 -b 34.2 35.2 45.0 46.3 -y 0.001 -x 0.001 -t cor -o geo_20171117_20180603.cor -outdir $SCRATCHDIR/preGeodmod/
  save_geodmod.py $SCRATCHDIR/DarbandikhanSenAT72/PYSAR/INPUTS/geometryRadar.h5 20171117_20180603 -b 34.2 35.2 45.0 46.3 -y 0.001 -x 0.001 -t dem -o srtm.dem -outdir $SCRATCHDIR/preGeodmod/

"""

def create_parser():
    parser = argparse.ArgumentParser(description='Prepare data for Geodmod software',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('file', nargs=1, help='ascending and descending timeseries files\n')
    #parser.add_argument('dset', nargs=1, help='date12 of timeseries to be converted')
    parser.add_argument('-b', '--bbox', dest='SNWE', type=float, nargs=4, metavar=('S', 'N', 'W', 'E'),
                        help='Bounding box of area to be geocoded.\n' +
                        'Include the uppler left corner of the first pixel' +
                        '    and the lower right corner of the last pixel')
    parser.add_argument('-y', '--latstep', dest='latStep', type=float,
                        help='output pixel size in degree in latitude.')
    parser.add_argument('-x', '--lonstep', dest='lonStep', type=float,
                        help='output pixel size in degree in longitude.')
    #parser.add_argument('-o', '--output', dest='outfile', nargs=1,
    #                    help='output file name')
    parser.add_argument('-startDate','--startDate',dest='StartDate',nargs=?,help='date1 of timeseires to be converted')
    parser.add_argument('-endDate','--endDate',dest='EndDate',nargs=?,help='date2 of timeseries to be converted')
    parser.add_argument('-outdir','--outdir',dest='outdir',nargs=1,
                        help='output directory')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)
    
    # default startDate and endDate
    atr = readfile.read_attribute(inps.file)
    if not inps.StartDate:
        inps.StartDate=atr['START_DATE']
    if not inps.EndDate:
        inps.EndDate=atr['END_DATE']
    #dateList = None;
    #if k in ['timeseries']:
    #    dateList = timeseries(inps.file).get_date_list()
    #if not inps.StartDate:
    #    if k in ['timeseries']:
    #       inps.StartDate=dateList[0]
    #    else:
    #        inps.StartDate=atr['START_DATE']
    #   #datevector = ptime.date_list2vector(dateList)[1]
    #if not inps.EndDate:
    #    if k in ['timeseries']:
    #        inps.EndDate=dateList[-1]
    #    else:
    #        inps.EndDate=atr['END_DATE']

    return inps
    
def format_args(arr, dst=""):
    """ parse list array(list item or values) to string """
    for k in arr:
        if isinstance(k,list):
        #if type(k) == list:
            dst = dst + " " + format_args(k)
        else:
            dst = dst + " " + str(k)
    return dst.strip()
   
def get_path_com(path):
    """ return directory(filepath), filename(shotname), extension """
    (filepath, tempfilename) = os.path.split(os.path.abspath(path))
    (filename, extension) = os.path.splitext(tempfilename)
    return filepath, filename, extension

def write_rsc_file(inps,in_file,out_file):
    """ write rsc file for Geodmod just estract several properities from rsc file"""
    # read file
    meta = readfile.read_roipac_rsc(in_file)
    # initiate dict
    rsc = dict()
    rsc['FILE_DIR'] = "".join(inps.outdir)
    rsc['FILE_LENGTH'] = meta["FILE_LENGTH"]
    rsc['WIDTH'] = meta["WIDTH"]
    rsc['XMIN'] = 0
    rsc['XMAX'] = int(meta["WIDTH"]) - 1
    rsc['YMIN'] = 0
    rsc['YMAX'] = int(meta["FILE_LENGTH"]) - 1
    rsc['X_FIRST'] = float(meta["X_FIRST"])
    rsc['Y_FIRST'] = float(meta["Y_FIRST"])
    rsc['X_STEP'] = float(meta["X_STEP"])
    rsc['Y_STEP'] = float(meta["Y_STEP"])
    rsc['X_UNIT'] = 'degrees'
    rsc['Y_UNIT'] = 'degrees'
    rsc['RLOOKS'] = meta["RLOOKS"]
    rsc['ALOOKS'] = meta["ALOOKS"]
    rsc['Z_OFFSET'] = 0
    rsc['Z_SCALE'] = 1
    rsc['PROJECTION'] = 'LATLON'
    rsc['DATE12'] = '111111-222222'
    # write rsc file
    writefile.write_roipac_rsc(rsc, out_file, print_msg=True)
    return out_file 

def processdata(inps):
    #use geocode.py and save_roipac.py to process data"
    atr_asc = inps.file[0]
    # find data path
    #os.chdir(os.path.dirname(atr_asc))
    
    if os.path.exists("".join(inps.outdir))=='False':
        os.mkdir("".join(inps.outdir))
   
    # process cor and dem dataset
    corname='temporalCoherence'
    cor_args = [corname, '-b',inps.SNWE, '-y',inps.latStep, '-x',inps.lonStep, '--outdir',"".join(inps.outdir)]
    print("geocode.py", cor_args)
    args_str = format_args(cor_args)
    mintpy.geocode.main(args_str.split())
    demname='geometryRadar'
    dem_args = [demname, '-b',inps.SNWE, '-y',inps.latStep, '-x',inps.lonStep, '--outdir',"".join(inps.outdir)]
    args_str = format_args(dem_args)
    mintpy.geocode.main(args_str.split())
   
    #unw file
    asc_args = [atr_asc, '-b',inps.SNWE, '-y',inps.latStep, '-x',inps.lonStep, '--outdir',"".join(inps.outdir)]
    print("geocode.py", asc_args)
    args_str = format_args(asc_args)
    mintpy.geocode.main(args_str.split())
        
    #save dataset of unw cor and dem
    os.chdir("".join(inps.outdir))
    filename, extension = get_path_com(atr_asc)[1:3]
    asct_args = ['geo_'+filename+extension, inps.StartDate,'_',inps.EndDate]
    print("save_roipac.py", asct_args)
    asct_str = format_args(asct_args)
    mintpy.save_roipac.main(asct_str.split())     

    asct_args = ['geo_temporalCoherence.h5', '-o', inps.StartDate,'_',inps.EndDate,'.cor']    
    print("save_roipac.py", asct_args)
    asct_str = format_args(asct_args)
    mintpy.save_roipac.main(asct_str.split())
    
    asct_args = ['geo_geometryRadar.h5', 'height', '-o', 'srtm.dem']
    print("save_roipac.py", asct_args)
    asct_str = format_args(asct_args)
    mintpy.save_roipac.main(asct_str.split())

######################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)
    
    #  geocode timeseries**.h5 file and get the deformation field of two time periods
    #  and geocode ifgramStack.h5 file and get the coherence of two time periods
    #  and geocode geometryRadar.h5 file and get the dem    
    processdata(inps)
    
    # rename *.rsc1 to *.rsc
    outfile = format_args(['srtm.dem' + '.rsc'])
    write_rsc_file(inps,outfile,format_args(['srtm.dem' +'.rsc1']))
    os.remove(outfile)
    print('rename *.rsc1 to *.rsc')
    os.rename(format_args(['srtm.dem' +'.rsc1']),outfile)

######################################################################################
if __name__ == '__main__':
    main()
