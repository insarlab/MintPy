#!/usr/bin/env python3


import os
import time
import glob
import shutil
import argparse
import subprocess
from pysar.utils import readfile


#####################################################################################
cDict = {}
config = {};  config['walltime'] = '2:00';  config['memory'] = '2000';  cDict['run_1_master']               = config
config = {};  config['walltime'] = '1:00';  config['memory'] = '2000';  cDict['run_2_focus_split']          = config
config = {};  config['walltime'] = '2:00';  config['memory'] = '7000';  cDict['run_3_geo2rdr_coarseResamp'] = config
config = {};  config['walltime'] = '2:00';  config['memory'] = '2000';  cDict['run_4_refineSlaveTiming']    = config
config = {};  config['walltime'] = '1:00';  config['memory'] = '2000';  cDict['run_5_invertMisreg']         = config
config = {};  config['walltime'] = '1:00';  config['memory'] = '2000';  cDict['run_6_fineResamp']           = config
config = {};  config['walltime'] = '1:00';  config['memory'] = '1000';  cDict['run_7_grid_baseline']        = config
config = {};  config['walltime'] = '2:00';  config['memory'] = '7000';  cDict['run_8_igram']                = config
config = {};  config['walltime'] = '0:30';  config['memory'] = '2000';  cDict['run_9_maskIgram']            = config
config = {};  config['walltime'] = '1:00';  config['memory'] = '5000';  cDict['run_10_unwrap']              = config
config = {};  config['walltime'] = '0:30';  config['memory'] = '2000';  cDict['run_11_maskUnwrap']          = config
num_run_file = len(cDict.keys())

#####################################################################################
EXAMPLE = """example:
  cd $SCRATCHDIR/KirishimaAlosAT424F620_630
  submit_run_files4stripmap_stack.py
  submit_run_files4stripmap_stack.py --start 2 --end 7
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Submit the run files to jobs for ISCE/stripmapStack',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

    parser.add_argument('--start', dest='startNum', type=int, default=1,
                        help='Start submitting at named run number. Default: 1')
    parser.add_argument('--end', dest='endNum', type=int, default=num_run_file,
                        help='End submitting at named run number. Default: {}'.format(num_run_file))

    parser.add_argument('--bsub', action='store_true', help='submit this script as a job to generic queue.')
    parser.add_argument('-r','--memory', type=int, default=2000, help='memory for bsub. Default: 2000')
    parser.add_argument('-w','--walltime', type=str, default='5:00', help='wall time for bsub. Default: 5:00')
    parser.add_argument('-e','--email', dest='email', help='email address to send notification when all jobs finished.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    if any(not os.path.isdir(i) for i in ['configs', 'run_files']):
        msg = 'ERROR: NO configs or run_files folder found in the current directory!'
        raise SystemExit(msg)
    return inps


#####################################################################################
def run_job_submission4run_files(start_num=1, end_num=num_run_file, run_file_dir='./run_files'):

    c_dict = dict()
    for key, value in cDict.items():
        num = int(key.split('_')[1])
        if start_num <= num <= end_num:
            c_dict[key] = value

    cwd = os.getcwd()
    os.chdir(run_file_dir)
    for run_file in c_dict.keys():
        config = c_dict[run_file]
        cmd = 'split_jobs.py {f} -r {r} -w {w}'.format(f=run_file,
                                                       r=config['memory'],
                                                       w=config['walltime'])
        print(cmd)
        status = subprocess.Popen(cmd, shell=True).wait()
        if status is not 0:
            raise RuntimeError("Error in step {}".format(run_file))
    os.chdir(cwd)
    return status


#####################################################################################
def main(iargs=None):
    inps = cmd_line_parse()
    start_time = time.time()

    if not inps.bsub:
        run_job_submission4run_files(inps.startNum, inps.endNum)
    else:
        print('bsub option is ON')
        # write run_stripmap_stack
        run_file = 'run_stripmap_stack'
        with open(run_file, 'w') as f:
            f.write('submit_stripmap_stack.py\n')
        cmd = 'chmod +x {}'.format(run_file)
        print(cmd)
        os.system(cmd)

        # submit job
        cmd = 'split_jobs.py {} -r {} -w {}'.format(run_file, inps.memory, inps.walltime)
        if inps.email:
            cmd += ' -e {}'.format(inps.email)
        print(cmd)
        os.system(cmd)

    # Timing
    m, s = divmod(time.time()-start_time, 60)
    print('\nTotal time: {:02.0f} mins {:02.1f} secs'.format(m, s))
    return 

#####################################################################################
if __name__ == '__main__':
    main()
