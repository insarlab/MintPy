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
config = {};  config['walltime'] = '1:00';  config['memory'] = '2000';  cDict['run_1_master']               = config
config = {};  config['walltime'] = '0:30';  config['memory'] = '2000';  cDict['run_2_focus_split']          = config
config = {};  config['walltime'] = '1:00';  config['memory'] = '5500';  cDict['run_3_geo2rdr_coarseResamp'] = config
config = {};  config['walltime'] = '1:00';  config['memory'] = '2000';  cDict['run_4_refineSlaveTiming']    = config
config = {};  config['walltime'] = '0:30';  config['memory'] = '2000';  cDict['run_5_invertMisreg']         = config
config = {};  config['walltime'] = '0:30';  config['memory'] = '2000';  cDict['run_6_fineResamp']           = config
config = {};  config['walltime'] = '0:30';  config['memory'] = '2000';  cDict['run_7_grid_baseline']        = config
config = {};  config['walltime'] = '1:00';  config['memory'] = '6000';  cDict['run_8_igram']                = config


#####################################################################################
EXAMPLE = """example:
  submit_run_files4stripmap_stack.py ./run_files -e $NOTIFICATION
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Submit the run files to jobs for ISCE/stripmapStack',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=EXAMPLE)

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
def run_job_submission4run_files():
    for run_file in sorted(cDict.keys()):
        config = cDict[run_file]
        cmd = 'split_jobs.py {f} -r {r} -w {w}'.format(f=run_file,
                                                       r=config['memory'],
                                                       w=config['walltime'])
        print(cmd)
        status = subprocess.Popen(cmd, shell=True).wait()
        if status is not 0:
            raise RuntimeError("Error in step {}".format(run_file))
    return status


def get_multilook_number():
    config_igram_file = glob.glob(os.path.join(os.getcwd(), 'configs/config_igram_*'))[0]
    config = readfile.read_template(config_igram_file, delimiter=':')
    return config['alks'], config['rlks']


def multilook_geometry():
    print('copy run_multilook_geometry file to the current directory')
    run_file = os.path.expandvars('${PYSAR_HOME}/sh/run_multilook_geometry')
    shutil.copy2(run_file, os.getcwd())

    alks, rlks = get_multilook_number()

    cmd = './{} {} {}'.format(run_file, alks, rlks)
    print(cmd)
    os.system(cmd)
    return


#####################################################################################
def main(iargs=None):
    inps = cmd_line_parse()
    start_time = time.time()

    os.chdir(inps.run_file_dir)
    print('Go to directory', inps.run_file_dir)

    if not inps.bsub:
        run_job_submission4run_files()
    else:
        print('bsub option is ON')
        # write run_stripmap_stack
        run_file = 'run_stripmap_stack'
        with open(run_file, 'w') as f:
            f.write('submit_stripmap_stack.py {}\n'.format(inps.run_file_dir))
        cmd = 'chmod +x {}'.format(run_file)
        print(cmd)
        os.system(cmd)

        # submit job
        cmd = 'split_jobs.py {} -r {} -w {}'.format(run_file, inps.memory, inps.walltime)
        if inps.email:
            cmd += ' -e {}'.format(inps.email)
        print(cmd)
        os.system(cmd)

    # prepare multilooked geometry
    print('-'*50)
    multilook_geometry()

    # Timing
    m, s = divmod(time.time()-start_time, 60)
    print('\nTotal time: {:02.0f} mins {:02.1f} secs'.format(m, s))
    return 

#####################################################################################
if __name__ == '__main__':
    main()
