#! /usr/bin/env python
############################################################
# Program is part of PySAR v1.0                            #
# Copyright(c) 2013, Heresh Fattahi                        #
# Author:  Heresh Fattahi                                  #
############################################################
# Yunjun, Jan 2016: add bonding points correction
# Yunjun, Jul 2016: add ramp removal step


import sys
import os
import getopt

import h5py
import numpy as np
from scipy.linalg import pinv

import pysar._pysar_utilities as ut
import pysar._readfile as readfile
import pysar._writefile as writefile
import pysar._remove_surface as rm


##########################################################################################
def phase_bonding(data,mask,x,y):
    ## Phase Jump Correction, using phase continuity on bridge/bonding points in each pair of patches.
    ## data : phase matrix need to be corrected
    ## mask : mask file marks different patches with different positive integers
    ## x/y  : array of bridge points, lied as: x_ref, x, x_ref, x
  
    ## loop based on number of bridges
    n_bridge = len(x)/2
    for i in range(1,n_bridge+1):
        p_ref = data[y[2*i-2],x[2*i-2]]
        p     = data[y[2*i-1],x[2*i-1]]
        n_jump = (abs(p-p_ref)+np.pi)//(2*np.pi)
        if not n_jump == 0:
            if p-p_ref >=0:  n_jump *= -1
            id = np.where(mask == mask[y[2*i-1],x[2*i-1]])
            data[id] = data[id] + n_jump*2*np.pi;
  
    return data

####################################################################################################
def usage():
    print '''
************************************************************************************
  Unwrapping Error Correction based on:
      1. Triangular consistency of interferograms (phase closure), or
      2. Phase continuity of close bonding points (spatial continuity)
  
  -------------------------------------------------------------------
  1. Correct unwrapping errors based on triangular consistency
      Based on phase closure of pairs circle (ab + bc + ca == 0), this method assumes
      a. abundance of network: for interferogram with unwrapping error, there is
         at least of one triangular connection to form a closed circle; with more
         closed circles comes better constrain.
      b. majority rightness: most of interferograms have to be right (no unwrapping
         error) to correct the wrong minority. And if most of interferograms have 
         unwrapping errors, then the minor right interferograms will turn into wrong.

  Usage:
      unwrap_error.py interferograms_file    [ mask_file ]
      unwrap_error.py -f interferograms_file [ -m mask_file -o output_file]

      -f : unwrapped interferograms, i.e. Seeded_LoadedData.h5
      -m : mask file to specify those pixels which user wants to correct for unwrapping errors.
      -o : output file name [default is interferogram_file_unwCor.h5]

  Examples:
      unwrap_error.py Seeded_LoadedData.h5 mask.h5
      unwrap_error.py -f Seeded_LoadedData.h5 -m mask.h5
      unwrap_error.py Seeded_LoadedData.h5


  -------------------------------------------------------------------
  2. Correct unwrapping errors based on bonding points
      This method assumes:
      a. no phase unwrapping error within each patch marked by mask file.
      b. the absolute phase difference of bonding points (usually close in space) is 
         smaller than one pi. Considering prevalent ramps in InSAR data might break
         this assumptioin for bonding points that are not very close, across a bay
         for example, we first estimate and remove a linear phase ramp, then applied
         phase continuity constrain, and add the removed ramp back at the end.
 
      Phase unwrapping error is corrected epoch by epoch, following the steps below:
      a. estimate and remove a linear phase ramp from unwrapped phase;
      b. following the pair order of bonding points, correct patch by patch marked
         by point's coordinate and mask file:
         1) use 1st point as reference, calculate integer N, add N*2pi to 2nd point's
            phase to make sure their absolute phase difference is smaller than pi.
         2) add N*2pi to all pixels in 2nd point's patch.
      c. add linear phase ramp estimated in step a back to the corrected phase in step b.

  Usage:
      unwrap_error.py -f interferograms_file -m mask_file -t template_file      [ -o output_file ]
      unwrap_error.py -f interferograms_file -m mask_file -x x_ref,x -y y_ref,y [ -o output_file ]

      -f : unwrapped interferogram(s), i.e. LoadedData.h5, .unw file
      -m : mask file to mark different patches that want to be corrected.
           Masked out area is marked with 0, patches/area needed to be corrected marked with 
           positive integers, i.e. 1, 2, 3, ...
      -o : output file name [default is interferogram_file_unwCor.h5/.unw]

      -x : reference and to-be-corrected patches' bridge points coordinate in x direction
      -y : reference and to-be-corrected patches' bridge points coordinate in y direction
           Example: x_ref1, x1, x_ref2, x2, ...
      -t : template file with unwrapError.bonding_point option given value as
               y_ref1,x_ref1,y1,x1,y_ref2,x_ref2,y2,x2 ...
           Example: pysar.unwrapError.yx = 283,1177,305,1247 

           Note: choose x/y_ref point in the patch that also have seed point, for consistency
                 in multiple images.
      --ramp         : ramp type, i.e. plane, quadratic
      --no-ramp-save : save corrected data with the ramp removed.

  Examples:
      unwrap_error.py -f Seeded_LoadedData.h5     -m Mask.h5 -t ShikokuT417F650_690AlosA.template
      unwrap_error.py -f Seeded_LoadedData.h5     -m Mask.h5 -x 283,305 -y 1177,1247
      unwrap_error.py -f Seeded_081018_090118.unw -m Mask_all.h5 -x 283,305 -y 1177,1247 --ramp quadratic

************************************************************************************
    '''


####################################################################################################
def main(argv):

    method    = 'triangular_consistency'    ## or 'bonding_point'
    ramp_type = 'plane'
    save_rampCor = 'yes'
    plot_bonding_points = 'yes'
  
    ##### Check Inputs
    if len(sys.argv)>2:
        try: opts, args = getopt.getopt(argv,'h:f:m:x:y:o:t:',['ramp=','no-ramp-save'])
        except getopt.GetoptError:  print 'Error while getting args';  usage(); sys.exit(1)
  
        for opt,arg in opts:
            if   opt in ['-h','--help']:    usage(); sys.exit()
            elif opt in '-f':    File     = arg
            elif opt in '-m':    maskFile = arg
            elif opt in '-o':    outName  = arg
            elif opt in '-x':    x = [int(i) for i in arg.split(',')];    method = 'bonding_point'
            elif opt in '-y':    y = [int(i) for i in arg.split(',')];    method = 'bonding_point'
            elif opt in '-t':    templateFile = arg
            elif opt in '--ramp'         :  ramp_type    = arg.lower()
            elif opt in '--no-ramp-save' :  save_rampCor = 'no'
  
    elif len(sys.argv)==2:
        if argv[0] in ['-h','--help']:    usage();  sys.exit()
        elif os.path.isfile(argv[0]):     File = argv[0];  maskFile = argv[1]
        else:    print 'Input file does not existed: '+argv[0];  sys.exit(1)
  
    else:  usage(); sys.exit(1)
  
    ##### Check template file
    try:
        templateFile
        templateContents = readfile.read_template(templateFile)
    except: pass
  
    try:
        yx = [int(i) for i in templateContents['pysar.unwrapError.yx'].split(',')]
        x = yx[1::2]
        y = yx[0::2]
        method = 'bonding_point'
    except: pass

    ##### Read Mask File 
    ## Priority:
    ## Input mask file > pysar.mask.file > existed Modified_Mask.h5 > existed Mask.h5
    try:       maskFile
    except:
        try:    maskFile = templateContents['pysar.mask.file']
        except:
            if   os.path.isfile('Modified_Mask.h5'):  maskFile = 'Modified_Mask.h5'
            elif os.path.isfile('Mask.h5'):           maskFile = 'Mask.h5'
            else: print 'No mask found!'; sys.exit(1)
    try:    Mask,Matr = readfile.read(maskFile);   print 'mask: '+maskFile
    except: print 'Can not open mask file: '+maskFile; sys.exit(1)
  
    ##### Output file name
    ext = os.path.splitext(File)[1]
    try:    outName
    except: outName = File.split('.')[0]+'_unwCor'+ext
  
    print '\n**************** Unwrapping Error Correction ******************'

    ####################  Triangular Consistency (Phase Closure)  ####################
    if method == 'triangular_consistency':
        print 'Phase unwrapping error correction using Triangular Consistency / Phase Closure'
  
        h5file=h5py.File(File)
        ifgramList = h5file['interferograms'].keys()
        sx = int(h5file['interferograms'][ifgramList[0]].attrs['WIDTH'])
        sy = int(h5file['interferograms'][ifgramList[0]].attrs['FILE_LENGTH'])
        curls,Triangles,C=ut.get_triangles(h5file)
        A,B = ut.design_matrix(h5file)   
        ligram,lv=np.shape(B)
        lcurls=np.shape(curls)[0]
        print 'Number of all triangles: '+  str(lcurls)
        print 'Number of interferograms: '+ str(ligram)
        #print curls
  
        curlfile='curls.h5'
        if not os.path.isfile(curlfile):
            ut.generate_curls(curlfile,h5file,Triangles,curls)
         
        thr=0.50
        curls=np.array(curls);   n1=curls[:,0];   n2=curls[:,1];   n3=curls[:,2]
  
        numPixels=sy*sx
        print 'reading interferograms...'   
        data = np.zeros((ligram,numPixels),np.float32)
        for ni in range(ligram):
            dset=h5file['interferograms'][ifgramList[ni]].get(ifgramList[ni])
            d = dset[0:dset.shape[0],0:dset.shape[1]]
            data[ni] = d.flatten(1)   
  
        print np.shape(data)
        print 'reading curls ...' 
        h5curl=h5py.File(curlfile)
        curlList=h5curl['interferograms'].keys()
        curlData = np.zeros((lcurls,numPixels),np.float32)
        for ni in range(lcurls):
            dset=h5curl['interferograms'][curlList[ni]].get(curlList[ni])
            d = dset[0:dset.shape[0],0:dset.shape[1]]
            curlData[ni] = d.flatten(1)
        pi=np.pi
        EstUnwrap=np.zeros((ligram,numPixels),np.float32)
  
        #try:
        #    maskFile=argv[1]
        #    h5Mask=h5py.File(maskFile)
        #    dset = h5Mask['mask'].get('mask')
        #    Mask=dset[0:dset.shape[0],0:dset.shape[1]]
        #except:
        #    dset = h5file['mask'].get('mask')
        #    Mask=dset[0:dset.shape[0],0:dset.shape[1]]
        
        Mask=Mask.flatten(1)

        for ni in range(numPixels):
            #dU = np.zeros([ligram,1])
            #print np.shape(dU)
            #print np.shape(data[:,ni])
  
            if Mask[ni]==1:
                dU = data[:,ni]
                #nan_ndx = dataPoint == 0.
                unwCurl = np.array(curlData[:,ni])
                #print unwCurl
  
                ind  = np.abs(unwCurl)>=thr;      N1 =n1[ind];      N2 =n2[ind];      N3 =n3[ind]
                indC = np.abs(unwCurl)< thr;      Nc1=n1[indC];     Nc2=n2[indC];     Nc3=n3[indC]
  
                N =np.hstack([N1, N2, N3]);       UniN =np.unique(N)
                Nc=np.hstack([Nc1,Nc2,Nc3]);      UniNc=np.unique(Nc)
  
                inter=list(set(UniNc) & set(UniN)) # intersetion
                UniNc= list(UniNc)
                for x in inter:
                    UniNc.remove(x)
  
                D=np.zeros([len(UniNc),ligram])
                for i in range(len(UniNc)):
                    D[i,UniNc[i]]=1
  
                AAA=np.vstack([-2*pi*C,D])
                #AAA1=np.hstack([AAA,np.zeros([AAA.shape[0],lv])])
                #AAA2=np.hstack([-2*pi*np.eye(ligram),B]) 
                #AAAA=np.vstack([AAA1,AAA2])
                AAAA=np.vstack([AAA,0.25*np.eye(ligram)])
  
                #print '************************'
                #print np.linalg.matrix_rank(C)
                #print np.linalg.matrix_rank(AAA) 
                #print np.linalg.matrix_rank(AAAA)
                #print '************************'
  
                #LLL=list(np.dot(C,dU)) + list(np.zeros(np.shape(UniNc)[0]))# + list(dU)
                #ind=np.isnan(AAA)
                #M1=pinv(AAA)      
                #M=np.dot(M1,LLL)
                #EstUnwrap[:,ni]=np.round(M[0:ligram])*2.0*np.pi
  
                ##########
                # with Tikhonov regularization:
                AAAA=np.vstack([AAA,0.25*np.eye(ligram)])
                LLL=list(np.dot(C,dU)) + list(np.zeros(np.shape(UniNc)[0])) + list(np.zeros(ligram))
                ind=np.isnan(AAAA)
                M1=pinv(AAAA)
                M=np.dot(M1,LLL)
                EstUnwrap[:,ni]=np.round(M[0:ligram])*2.0*np.pi
                #print M[0:ligram]
                #print np.round(M[0:ligram])
  
            else:
                EstUnwrap[:,ni]=np.zeros([ligram])
                if not np.remainder(ni,10000): print 'Processing point: %7d of %7d ' % (ni,numPixels)

        ##### Output
        dataCor = data+EstUnwrap
        unwCorFile=File.replace('.h5','')+'_unwCor.h5';  print 'writing >>> '+unwCorFile
        h5unwCor=h5py.File(unwCorFile,'w') 
        gg = h5unwCor.create_group('interferograms') 
        for i in range(ligram):
            group = gg.create_group(ifgramList[i])
            dset = group.create_dataset(ifgramList[i], data=np.reshape(dataCor[i,:],[sx,sy]).T, compression='gzip')
            for key, value in h5file['interferograms'][ifgramList[i]].attrs.iteritems():
                group.attrs[key] = value
  
        try:
            MASK=h5file['mask'].get('mask')
            gm = h5unwCor.create_group('mask')
            dset = gm.create_dataset('mask', data=MASK, compression='gzip')
        except: pass
  
        h5unwCor.close()
        h5file.close()
        h5curl.close() 


    ####################  Bonding Points (Spatial Continuity)  ####################
    elif method == 'bonding_point':
        print 'Phase unwrapping error correction using Bonding Points / Spatial Continuity'
  
        ##### Read Bridge Points Info
        try:
            x
            y
            if len(x) != len(y) or np.mod(len(x),2) != 0:
                print 'Wrong number of bridge points input: '+str(len(x))+' for x, '+str(len(y))+' for y'
                usage();  sys.exit(1)
        except: print 'Error in reading bridge points info!';  usage();  sys.exit(1)
        for i in range(0,len(x)):
            if Mask[y[i],x[i]] == 0:
                print '\nERROR: Connecting point ('+str(y[i])+','+str(x[i])+') is out of masked area! Select them again!\n'
                sys.exit(1)
  
        print 'Number of bonding point pairs: '+str(len(x)/2)
        print 'Bonding points coordinates:\nx: '+str(x)+'\ny: '+str(y)
  
        ## Plot Connecting Pair of Points
        if plot_bonding_points == 'yes':
            point_yx = ''
            line_yx  = ''
            n_bridge = len(x)/2
            for i in range(n_bridge):
                pair_yx = str(y[2*i])+','+str(x[2*i])+','+str(y[2*i+1])+','+str(x[2*i+1])
                if not i == n_bridge-1:
                    point_yx += pair_yx+','
                    line_yx  += pair_yx+';'
                else:
                    point_yx += pair_yx
                    line_yx  += pair_yx

            try:
                plot_cmd = 'view.py --point-yx="'+point_yx+'" --line-yx="'+line_yx+\
                           '" --nodisplay -o bonding_points.png -f '+maskFile
                print plot_cmd
                os.system(plot_cmd)
            except: pass


        ##### Ramp Info
        ramp_mask = Mask==1
        print 'estimate phase ramp during the correction'
        print 'ramp type: '+ramp_type
        if save_rampCor == 'yes':
            outName_ramp = os.path.basename(outName).split(ext)[0]+'_'+ramp_type+ext
  
        ########## PySAR ##########
        if ext == '.h5':
            ##### Read
            try:     h5file=h5py.File(File,'r')
            except:  print 'ERROR: Cannot open input file: '+File; sys.exit(1)
            k=h5file.keys()
            if 'interferograms' in k: k[0] = 'interferograms';  print 'Input file is '+k[0]
            else: print 'Input file - '+File+' - is not interferograms.';  usage();  sys.exit(1)
            igramList = h5file[k[0]].keys()
            igramList = sorted(igramList)
  
            #### Write
            h5out = h5py.File(outName,'w')
            gg = h5out.create_group(k[0])
            print 'writing >>> '+outName
  
            if save_rampCor == 'yes':
                h5out_ramp = h5py.File(outName_ramp,'w')
                gg_ramp = h5out_ramp.create_group(k[0])
                print 'writing >>> '+outName_ramp
  
            ##### Loop
            print 'Number of interferograms: '+str(len(igramList))
            for igram in igramList:
                print igram
                data = h5file[k[0]][igram].get(igram)[:]
  
                data_ramp,ramp = rm.remove_data_surface(data,ramp_mask,ramp_type)
                #ramp = data_ramp - data
                data_rampCor = phase_bonding(data_ramp,Mask,x,y)
                dataCor = data_rampCor - ramp
  
                group = gg.create_group(igram)
                dset = group.create_dataset(igram, data=dataCor, compression='gzip')
                for key, value in h5file[k[0]][igram].attrs.iteritems():
                    group.attrs[key]=value
  
                if save_rampCor == 'yes':
                    group_ramp = gg_ramp.create_group(igram)
                    dset = group_ramp.create_dataset(igram, data=data_rampCor, compression='gzip')
                    for key, value in h5file[k[0]][igram].attrs.iteritems():
                        group_ramp.attrs[key]=value
  
            try:
                mask = h5file['mask'].get('mask');
                gm = h5out.create_group('mask')
                dset = gm.create_dataset('mask', data=mask[0:mask.shape[0],0:mask.shape[1]], compression='gzip')
            except: print 'no mask group found.'
  
            h5file.close()
            h5out.close()
            if save_rampCor == 'yes':
                h5out_ramp.close()

        ########## ROI_PAC ##########
        elif ext == '.unw':
            print 'Input file is '+ext
            a,data,atr = readfile.read_float32(File);
  
            data_ramp,ramp = rm.remove_data_surface(data,ramp_mask,ramp_type)
            #ramp = data_ramp - data
            data_rampCor = phase_bonding(data_ramp,Mask,x,y)
            dataCor = data_rampCor - ramp
  
            writefile.write(dataCor, atr, outName)
            if save_rampCor == 'yes':
                writefile.write(data_rampCor,atr,outName_ramp)
  
        else: print 'Un-supported file type: '+ext;  usage();  sys.exit(1)



####################################################################################################
if __name__ == '__main__':
    main(sys.argv[1:])

