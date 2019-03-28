#!/usr/bin/env python3
# Author: Zhang Yunjun, 2019-03-19


import os
import argparse
import subprocess
import configparser
from pysar.objects import sensor
from pysar.utils import readfile, utils as ut


#####################################################################################
EXAMPLE = """example:
  cd $SCRATCHDIR/KirishimaAlosAT424F620_630
  process_isce_stack.py -t KirishimaAlosAT424F620_630.txt --bsub -e $NOTIFICATIONEMAIL
  process_isce_stack.py -t KirishimaAlosAT424F620_630.txt --start 2 --end 7
  process_isce_stack.py --reset
"""

TEMPLATE = """template:
##------------------------------- ISCE/stripmapStack OPTIONS ------------------##
isce.processor          = stripmapStack  #[stripmapStack, topsStack]
isce.ALOS.fbd2fbs       = yes
isce.demSNWE            = 31.1, 32.8, 130.1, 131.9  #[S, N, W, E] in degree
isce.demFile            = ${SCRATCHDIR}/KirishimaAlosAT424F620_630/DEM/gsi10m.dem
isce.azimuthLooks       = 20
isce.rangeLooks         = 8
isce.maxTempBaseline    = 1800
isce.maxPerpBaseline    = 1800
isce.masterDate         = 20080212
isce.unwrapMethod       = snaphu
isce.filtStrength       = 0.5
isce.applyWaterMask     = yes
"""

def create_parser():
    parser = argparse.ArgumentParser(description='HPC Wrapper for InSAR stack processor.',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=TEMPLATE+'\n'+EXAMPLE)

    parser.add_argument('-t','--template', dest='templateFile', type=str, help='template file')
    parser.add_argument('--start', dest='startNum', type=int, help='Start submitting at named run number.')
    parser.add_argument('--end', dest='endNum', type=int, help='End submitting at named run number.')
    parser.add_argument('--reset', action='store_true', help='clean the directory before re-run.')

    parser.add_argument('--bsub', action='store_true', help='submit this script as a job to generic queue.')
    parser.add_argument('-r','--memory', type=int, default=2000, help='memory for bsub. Default: 2000')
    parser.add_argument('-w','--walltime', type=str, default='8:00', help='wall time for bsub. Default: 8:00')
    parser.add_argument('-e','--email', dest='email', help='email address to send notification when all jobs finished.')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    if not inps.reset and not inps.templateFile:
        raise ValueError('ERROR: the following arguments are required: -t/--template')

    if inps.templateFile:
        inps.templateFile = os.path.abspath(inps.templateFile)
    #if any(not os.path.isdir(i) for i in ['configs', 'run_files']):
    #    msg = 'ERROR: NO configs or run_files folder found in the current directory!'
    #    raise SystemExit(msg)
    return inps


def read_inps2dict(inps):
    print('read options from template file: '+os.path.basename(inps.templateFile))
    template = readfile.read_template(inps.templateFile)
    template = ut.check_template_auto_value(template)

    iDict = vars(inps)
    key_prefix = 'isce.'
    key_list = [i.split(key_prefix)[1] for i in template.keys() if i.startswith(key_prefix)]
    for key in key_list:
        iDict[key] = template[key_prefix+key]
    return iDict


#####################################################################################
def reset_process_directory():
    cmd_str="""------ Copy and paste the following the command to reset the process direction ----
rm -r baselines/ configs/ coregSLC/ geom_master/ Igrams/ merged/ offsets/ refineSlaveTiming/ run_* SLC/
cd download
rm -rf 20* ALP*
mv ARCHIVED_FILES/* .
cd ..
    """
    print(cmd_str)
    return


def prepare_ALOS(iDict):
    # uncompress tar/zip files
    cmd = 'prepRawALOS.py -i ./download -o ./SLC -t "" '
    if iDict['ALOS.fbd2fbs']:
        cmd += ' --dual2single '
    print(cmd)
    status = subprocess.Popen(cmd, shell=True).wait()
    if status is not 0:
        raise RuntimeError('Error while runing prepRawALOS.py.')

    # submit job for run_unPack*
    run_file = 'run_unPackALOS'
    cmd = 'split_jobs.py run_unPackALOS -r 2000 -w 0:30'
    print(cmd)
    status = subprocess.Popen(cmd, shell=True).wait()
    if status is not 0:
        raise RuntimeError('Error while runing run_unPackALOS.')
    return status


def prepare_stack(iDict):
    cmd = ('stackStripMap.py -W interferogram -s ./SLC -d {d} -u {u} -f {f} '
           ' -t {t} -b {b} -a {a} -r {r} -m {m}').format(d=iDict['demFile'],
                                                         u=iDict['unwrapMethod'],
                                                         f=iDict['filtStrength'],
                                                         t=iDict['maxTempBaseline'],
                                                         b=iDict['maxPerpBaseline'],
                                                         a=iDict['azimuthLooks'],
                                                         r=iDict['rangeLooks'],
                                                         m=iDict['masterDate'])
    if iDict['applyWaterMask']:
        cmd += ' --applyWaterMask'
    print(cmd)
    status = subprocess.Popen(cmd, shell=True).wait()
    if status is not 0:
        raise RuntimeError('Error while runing stackStripMap.py.')
    return status


def run_stack(iDict, run_file_dir='./run_files'):
    # read job config file
    config_dir = os.path.expandvars('${PYSAR_HOME}/pysar/defaults')
    config_file = os.path.join(config_dir, 'job4{}.cfg'.format(iDict['processor']))
    if not os.path.isfile(config_file):
        raise ValueError('job config file NOT found, it should be: {}'.format(config_file))

    config = configparser.ConfigParser(delimiters='=')
    config.read(config_file)

    num_step = len(config.sections())
    print('number of steps: {}'.format(num_step))
    if not iDict['startNum']:
        iDict['startNum'] = 1
    if not iDict['endNum']:
        iDict['endNum'] = num_step

    # check ./run_files directory
    if any(not os.path.isdir(i) for i in ['configs', 'run_files']):
        msg = 'ERROR: NO configs or run_files folder found in the current directory!'
        raise SystemExit(msg)

    run_file_dir = os.path.join(os.getcwd(), 'run_files')
    cwd = os.getcwd()
    os.chdir(run_file_dir)

    for step_num in range(iDict['startNum'], iDict['endNum']+1):
        step_prefix = 'run_{:d}_'.format(step_num)
        step_name = [i for i in config.sections() if i.startswith(step_prefix)][0]

        cmd = 'split_jobs.py {f} -r {r} -w {w}'.format(f=step_name,
                                                       r=config[step_name]['memory'],
                                                       w=config[step_name]['walltime'])
        print(cmd)
        status = subprocess.Popen(cmd, shell=True).wait()
        if status is not 0:
            raise RuntimeError("Error in step {}".format(step_name))
    os.chdir(cwd)
    return status


def write_job_file(iDict):
    # command line
    cmd = 'process_isce_stack.py -t {}'.format(iDict['templateFile'])
    if iDict['startNum']:
        cmd += ' --start {} '.format(iDict['startNum'])
    if iDict['endNum']:
        cmd += ' --end {} '.format(iDict['endNum'])
    print('run the following command in bsub mode')
    print(cmd)

    # write job file
    job_dir = os.getcwd()
    job_file = os.path.join(job_dir, 'z_input_process_isce_stack.job')
    job_name = sensor.project_name2sensor_name(iDict['templateFile'])[1]
    with open(job_file, 'w') as f:
        f.write('#! /bin/tcsh\n')
        f.write('#BSUB -J {}\n'.format(job_name))
        f.write('#BSUB -P insarlab\n')
        f.write('#BSUB -o z_output_{}.%J.o\n'.format(job_name))
        f.write('#BSUB -e z_output_{}.%J.e\n'.format(job_name))
        f.write('#BSUB -W {}\n'.format(iDict['walltime']))
        f.write('#BSUB -q general\n')
        f.write('#BSUB -n 1\n')
        f.write('#BSUB -R "rusage[mem={}]"\n'.format(iDict['memory']))
        if iDict['email']:
            f.write('#BSUB -u {}\n'.format(iDict['email']))
            f.write('#BSUB -N\n')

        # write cd work directory
        f.write('\n')
        f.write('cd {}\n'.format(job_dir))
        f.write('{}\n'.format(cmd))
    print('finished writing job file: {}'.format(job_file))
    return job_file


#####################################################################################
def main(iargs=None):
    inps = cmd_line_parse()

    if inps.reset:
        reset_process_directory()
        return

    iDict = read_inps2dict(inps)

    job_file = write_job_file(iDict)
    if iDict['bsub']:
        cmd = 'bsub < {}'.format(job_file)
        print(cmd)
        os.system(cmd)
        return

    if not inps.startNum:
        prepare_ALOS(iDict)
        prepare_stack(iDict)

    run_stack(iDict)

    return

#####################################################################################
if __name__ == '__main__':
    main()
