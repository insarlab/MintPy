#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 23:00:53 2021

@author: Marin Govorcin
"""

from mintpy.objects import ifgramStack
from mintpy.utils import arg_group, readfile, writefile, utils as ut
import numpy as np
import argparse
import sys

EXAMPLE = """example:
  # Get the coherence change for the event [volcano eruption] on 20170904 in the FernandinaSenDT128 test dataset
  # Coherence Change is calculated as:
       ccd = mean(coherenceStack_beforeEvent) - mean(coherenceStack_coveringEvent)

  coherenceCD.py inputs/ifgramStack.h5 -e 20170904 

  # Use option '--start' to define the first date of coherence stack before the event, if not used start date is the first date of timeseries
  # You might want to use this option if you want to calculate avg. pre-event stack coherence for defined "short" period (let say 1-2 month) before 
  # the event to avoid seasonal change in coherence

  coherenceCD.py inputs/ifgramStack.h5 -e 20170904 --start 20160505

  # use option --thresh to remove coherence datasets with average spatial coherence for the whole scene lower than defined threshold 
  # this option can be played with subset option

  coherenceCD.py inputs/ifgramStack.h5 -e 20170904 --start 20160505 -t 0.71
  ## 
 
"""

def create_parser():
    parser = argparse.ArgumentParser(description=f'Generate Coherence Change Detection map from Mintpy ifgram.h5 .',
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('file', type=str, help='stack of coherences, e.g. ifgramStack.h5')
    parser.add_argument('-e', '--event', dest='event', type=str, required=True,
                        help='date of event to calculate pre-event and co-event average coherence.\n'+
                             'e.g.:   date12 in YYYYMMDD')
    parser.add_argument('-s', '--start', dest='sdate', type=str,
                        help='start date to calculate pre-event average coherence.\n'+
                             'e.g.:   date12 in YYYYMMDD')
    parser.add_argument('-t', '--thresh', dest='min_coh', type=str,
                        help='use threshold to remove data based on avg. spatial coherence for the whole scene')
    parser.add_argument('-g', '--geom', dest='geom', type=str,
                        help='geometryRadar.h5, use with subset lat,lon option')
    parser.add_argument('-o', '--output', dest='outfile', type=str,
                        help='output filename')
    parser = arg_group.add_subset_argument(parser)
    return parser

def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    return inps

#########################################################################################################

def getCoherence(ifgram, event, sdate=None, ibox=None):
    # read the stack and get the date12 list
    obj = ifgramStack(ifgram)
    obj.open()

    date12List = obj.date12List
    
    #if the start date for pre-event coherence is not defined use the first acquistion date
    if sdate is None:
        sdate = obj.date12List[0]

    pre_eventIdx = []
    co_eventIdx = []
    
    for i in range(len(date12List)):
        if date12List[i].split('_')[0] < event and date12List[i].split('_')[0] > sdate and date12List[i].split('_')[1] < event:
            pre_eventIdx.append(i)
        elif date12List[i].split('_')[0] < event and date12List[i].split('_')[1] > event:
            co_eventIdx.append(i)
  
    #read coherence from the stack
    coherence = obj.read(datasetName='coherence',box=ibox)
    pre_coherence = coherence[pre_eventIdx,:,:]
    co_coherence = coherence[co_eventIdx,:,:]
    
    print('List of pre-event coherence data:\n{}'.format(np.array(date12List)[pre_eventIdx]))
    print('List of co-event coherence data:\n{}'.format(np.array(date12List)[co_eventIdx]))
    obj.close()
    
    return pre_coherence, co_coherence

def rm_coh(coh,thresh, pmsg = None):
    
    idx = []
    nidx =[]
    for i in range(len(coh)):
        avg_coh = np.nanmean(coh[i,:,:])
        if pmsg:
            print('{}: mean coh. {:.2f}'.format(i+1,avg_coh))
        if avg_coh > float(thresh):
            idx.append(i)
        else:
            nidx.append(i+1)
    print('Removed coherence data: {}'.format(nidx))
    print('#####################################\n')
    coh = coh[idx]
    
    return coh

def coherenceCD(coh1,coh2):
    # get average coherence of the coherence stack
    if len(coh1) > 1:
        mcoh1 = np.mean(np.array(coh1),axis=0)
    else:
        mcoh1 = coh1
        
    if len(coh2) > 1:
        mcoh2 = np.mean(np.array(coh2),axis=0)
    else:
        mcoh2 = coh2


    ccd = mcoh1 - mcoh2
    
    return mcoh1, mcoh2, ccd

def main(iargs=None):
    inps = cmd_line_parse(iargs)
    
    print('Read metadata from file: {}'.format(inps.file))
    attr = readfile.read_attribute(inps.file)

    #Extract subset if defined
    #need to open ifgramstack
    obj = ifgramStack(inps.file)
    obj.open()

    if inps.subset_lon and inps.subset_lat:
        lalo = (inps.subset_lon[0], inps.subset_lat[1],inps.subset_lon[1],inps.subset_lat[0])
        coord=ut.coordinate(obj.metadata, lookup_file=inps.geom)
        inps.pix_box = coord.bbox_geo2radar(lalo)
        print(inps.pix_box)
    else:
        inps.pix_box = None

    if inps.subset_x and inps.subset_y:
        pix = (inps.subset_x[0],inps.subset_y[1],inps.subset_x[1],inps.subset_y[0])
        coord=ut.coordinate(obj.metadata, lookup_file=inps.geom)
        inps.pix_box =coord.check_box_within_data_coverage(pix)
        print(inps.pix_box)
    else:
        inps.pix_box = None
 
    #Read data from
    print('Read Coherence from file: {}'.format(inps.file))
    print('Event date {}'.format(inps.event))
    pc, cc = getCoherence(inps.file,inps.event, inps.sdate,ibox=inps.pix_box)
    
    #Remove the data below the min average threshold
    if inps.min_coh:
        print('########## Pre-Coherence Stack ########\n')
        pc = rm_coh(pc, inps.min_coh, pmsg=True)
        print('########## Co-Coherence Stack ########\n')
        cc = rm_coh(cc, inps.min_coh, pmsg=True)
    
    # Coherence Change Detection
    mpc, mcc, ccd = coherenceCD(pc,cc)
    
    dsDict={}
    dsDict['pre_event-avg_coherence'] = mpc
    dsDict['co_event-avg_coherencee'] = mcc
    dsDict['CCD'] = ccd
    
    size =ccd.shape

    attr['FILE_TYPE'] = 'coherence'
    attr['PROCESOR'] = 'ISCE'
    attr['LENGTH'] = size[0]
    attr['WIDTH'] = size[1]

    if inps.outfile is None:
        inps.outfile = 'ccd.h5'

     
    print('Write CCD in the file {}'.format(inps.outfile))
    writefile.write(dsDict,out_file=inps.outfile,metadata=attr)

#########################################################################################################
if __name__ == "__main__":
    main(sys.argv[1:])
