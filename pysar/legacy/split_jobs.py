#! /usr/bin/env python2
################################################################
# Author: Yunjun Zhang, 2016-05-23
#   Modified from split_igramJobs_done.py written by Heresh.
################################################################

import sys
import os
import getopt
import glob
import time
import datetime

from messageRsmas import Message as msg
from messageRsmas import log


#####################################################################
def Usage():
  print '''
  ****************************************************************

  Split batch command line into job files and submit them.

  Usage:
        split_jobs.py run_file
        split_jobs.py -f run_file [ -n core_number -l command_number -w walltime
                                    -r memory -q queue -e email -p projectID
                                    --parallel ]

        -e/--email    : email address to send notification when job finished
        -f/--file     : list file of command, i.e. run_raw2slc
        -h/--help     : help, usage

        -n/--number   : number of processors per job          [1          by default]
        -l/--line     : number of command    per job          [same as -n by default]
        -p/--project  : project name for the job              [insarlab   by default]
        -q/--queue    : job queue, i.e. general, parallel     [general    by default]
        -w/--walltime : maximum time         per job          [5:00       by default]
        -r/--memory   : maximum RAM memory in Mb per command  [3.7 GB     by default]
                        if N lines of command in one job, still put memory requirement
                        for one command line, script will calculate and submit N*memory
                        to CCS.

        --parallel    : enable parallel computing using subprocess.Popen

  Pegasus Job Queues on CCS:
         name       processors        memory         walltime
                      cores                         default/max
        general        15-          25 GB max      1 day   / 7 days
        parallel       16+          25 GB max      1 day   / 7 days
        bigmem         64 max      250 GB max      4 hours / 3 days
        interactive    15-         250 GB max      6 hours / 1 day

  Example:
        split_jobs.py run_raw2slc
        split_jobs.py -f run_raw2slc -n 1  -w 3:00 -r 3700 -e $NOTIFICATIONEMAIL
        split_jobs.py -f run_raw2slc -n 8  -w 5:00 -r 8000
        split_jobs.py -f run_raw2slc -n 8  -w 3:00 -r 8000 --parallel
        split_jobs.py -f run_raw2slc -n 16 -w 3:00 -r 8000 --parallel -q parallel
        split_jobs.py --file=run_raw2slc --number=1 --walltime=5:00

  ****************************************************************
  '''


#####################################################################
def main(argv):

  ##### Inputs
  coreNum    = 1
  memory     = 3700
  projectID  = 'insarlab'
  parallel   = 'no'
  queue      = 'general'
  walltime   = '5:00'
  waitJob    = True

  if len(sys.argv)>2:
      try:   opts, args = getopt.getopt(argv,'h:f:n:w:r:q:p:e:l:',\
                                            ['help','file=','number=','walltime=',\
                                             'memory=','queue=','project=','email=',\
                                             'parallel','line=','nowait'])
      except getopt.GetoptError:   Usage(); sys.exit(1)
      if opts==[]: Usage(); sys.exit(1)

      for opt,arg in opts:
          if   opt in ['-h','--help'    ]: Usage(); sys.exit()
          elif opt in ['-f','--file'    ]: run_file = arg
          elif opt in ['-n','--number'  ]: coreNum  = int(arg)
          elif opt in ['-w','--walltime']: walltime = arg
          elif opt in ['-r','--memory'  ]: memory   = int(arg)
          elif opt in ['-q','--queue'   ]: queue    = arg
          elif opt in ['-p','--project' ]: projectID= arg
          elif opt in ['-e','--email'   ]: email    = arg
          elif opt in ['-l','--line'    ]: lineNum  = int(arg)
          elif opt in ['--parallel'     ]: parallel = 'yes'
          elif opt in ['--nowait'       ]: waitJob  = False

  elif len(sys.argv)==2 and os.path.isfile(argv[0]):   run_file = argv[0]
  else:  Usage(); sys.exit(1)

  run_file = os.path.realpath(run_file)
  run_name = os.path.basename(run_file)
  workDir  = os.path.dirname(run_file)
  os.chdir(workDir)
  log(os.path.basename(sys.argv[0])+' '+run_name)

  #if int(coreNum) >= 16 and queue == 'general':
  #    msg('More than 16 processors is setted, change job queue to parallel.')
  #    queue = 'parallel'

  try:    lineNum
  except: lineNum = coreNum

  memory *= lineNum

  #if int(memory) >= 25000:
  #    msg('Setted memory exceed the maximum, change it to the limit [25 GB].')
  #    memory = '25000'

  ##### Clean
  prefix_in  = 'z_input_' +run_name+'.'
  prefix_out = 'z_output_'+run_name+'.'
  cleanCmd = 'rm '+prefix_in+'* '+prefix_out+'*'
  msg(cleanCmd)
  os.system(cleanCmd)

  ##### Input file info
  msg('Input file: '+run_file)
  with open(run_file) as f_run:
      cmdNum = sum(1 for _ in f_run)
  f_run.close()
  msg('Number of lines: '+str(cmdNum))

  ##### Split run_file according to the job/processor number
  jobNum = int((cmdNum-0.5)/lineNum) + 1
  digitNum = len(str(jobNum))
  splitCmd = 'split -a '+str(digitNum)+' -l '+str(lineNum)+' -d '+run_file+' '+prefix_in
  msg(splitCmd)
  os.system(splitCmd)
  z_inList = glob.glob(prefix_in+'*')

  ##### Job Info
  z_inNum = len(z_inList)
  msg('Number of jobs to submit: '+str(z_inNum))
  msg('Queue: '+queue)
  msg('Walltime: '+walltime)
  msg('Memory:'+str(memory))
  msg('Processors:'+str(coreNum))
  try:
      email
      msg('Email: '+email)
  except: pass

  ##### Write and Submit job file
  for z_in in z_inList:
      ##### Write job setting
      count = z_in.split(prefix_in)[1]
      job = z_in+'.sh'
      f = open(job,'w')
      f.write('#! /bin/bash')
      f.write('\n#BSUB -J '+run_name+'.'+count)
      f.write('\n#BSUB -P '+projectID)
      f.write('\n#BSUB -o '+prefix_out+count+'.%J.o')
      f.write('\n#BSUB -e '+prefix_out+count+'.%J.e')
      f.write('\n#BSUB -W '+walltime)
      f.write('\n#BSUB -q '+queue)
      f.write('\n#BSUB -n '+str(coreNum))

      try:
          f.write('\n#BSUB -R "rusage[mem='+str(memory)+']"')
      except: pass

      if queue == 'parallel':
          f.write('\n#BSUB -R "span[ptile=16]"')

      try:
          f.write('\n#BSUB -u '+email)
          f.write('\n#BSUB -N')
      except: pass

      ##### Write job excutable commands
      f.write('\ncd '+workDir+'\n')
      fz_in = open(z_in,'r')

      ## parallel computing
      if parallel == 'yes' and lineNum>1:
          parallelFile = z_in+'.popen.py'
          f_parallel = open(parallelFile,'w')
          f_parallel.write('#! /usr/bin/env python\n')
          f_parallel.write('\nfrom subprocess import Popen\n')
          for line in fz_in:
              f_parallel.write('\nPopen('+str(line.split())+')')
          f_parallel.close()
          os.system('chmod 755 '+parallelFile)
          f.write('./'+parallelFile)

      ## sequential computing
      else:
          for line in fz_in:
              f.write(line)

      fz_in.close()
      f.close()
 
      submitCmd = 'bsub < '+job;   msg('\n'+submitCmd);   os.system(submitCmd)
      cleanCmd  = 'rm '+z_in;                             os.system(cleanCmd)

  ##### Loop waiting until all jobs are done to exit script
  if not waitJob:
      msg('Job submission completed, exit without waiting.')
      return

  msg('\nSleeping until '+str(z_inNum)+' jobs are done for '+run_name)
  im = 0
  z_outNum = len(glob.glob(prefix_out+'*.o'))
  while z_outNum < z_inNum:
      time.sleep(1)
      if im%60 == 0:       #every minute
          msg('Current # of '+prefix_out+'*.o files in '+os.getcwd()+\
              ': <'+str(z_outNum)+'> out of <'+str(z_inNum)+'> after <'+str(im/60)+'> minutes')
      z_outNum = len(glob.glob(prefix_out+'*.o'))
      im = im + 1
  msg('All '+str(z_inNum)+' jobs are done for '+run_name)
  m,s = divmod(im,60); h,m = divmod(m,60)
  msg('Total used time: %02d hours %02d mins %02d secs'%(h,m,s))

  ##### move z_* files into dedicated directory
  jobDir = run_file+'_jobs'
  msg('Moving z_* files into '+jobDir)
  rmCmd = 'rm -rf '+jobDir; os.system(rmCmd)
  mkCmd = 'mkdir ' +jobDir; os.system(mkCmd)
  mvCmd = 'mv '+prefix_in+'* '+prefix_out+'* '+jobDir; os.system(mvCmd)
  msg('Finished at '+str(datetime.datetime.now()))


#####################################################################

if __name__ == '__main__':
  main(sys.argv[1:])

 
