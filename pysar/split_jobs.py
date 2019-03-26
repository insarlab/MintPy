#! /usr/bin/env python3
############################################################
# Program is part of PySAR                                 #
# Copyright(c) 2016-2019, Zhang Yunjun, Heresh Fattahi     #
# Author:  Zhang Yunjun, Heresh Fattahi                    #
############################################################


import os
import glob
import time
import argparse
import datetime
from pysar.objects import sensor


#####################################################################################
EXAMPLE = """example:
  split_jobs.py run_1_master -w 1:00 -r 2000 -e $NOTIFICATIONEMAIL
  split_jobs.py run_stripmap_stack -w 8:00 -r 200 -l 10 -e $NOTIFICATIONEMAIL
"""

CCS = """Pegasus Job Queues on UM/CCS
   name       processors        memory         walltime             description
                (cores)                      default / max
  general        15-          25 GB max      1 day   / 7 days    -R "span[ptime=16]"
  parallel       16+          25 GB max      1 day   / 7 days
  bigmem         64 max      250 GB max      4 hours / 5 days
  interactive    15-         250 GB max      6 hours / 1 day     max 1 job per user
"""

def create_parser():
    parser = argparse.ArgumentParser(description='Split and submit batch command lines into job files',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog=CCS+'\n'+EXAMPLE)

    parser.add_argument('run_file', help='run file with >1 command line')
    parser.add_argument('-w','--walltime', type=str, default='3:00',
                        help='maximum time per job. Default: 3:00')
    parser.add_argument('-r','--memory', type=int, default='3700',
                        help='maximum RAM memoery in MB per command line. Default: 3700\n'
                             'if N lines of commands in one job, still put memory for one command line,\n'
                             'this script will calculate and request N*memory to CCS.')
    parser.add_argument('-e','--email', dest='email', help='email address to send notification when job finished.')

    parser.add_argument('-q','--queue', dest='queue_name', default='general', choices={'general','parallel','bigmem'},
                        help='job queue to submit. Default: general')
    parser.add_argument('-p','--project', dest='project_name', default='insarlab',
                        help='project name for the job. Default: insarlab')

    parser.add_argument('-n','--number', dest='num_processor', type=int, default=1,
                        help='number of processors per job. Default: 1')
    parser.add_argument('-l','--line', dest='num_cmd', type=int,
                        help='number of command per job. Default: same as -n')
    return parser


def cmd_line_parse(iargs=None):
    parser = create_parser()
    inps = parser.parse_args(args=iargs)

    inps.run_file = os.path.realpath(inps.run_file)

    if not inps.num_cmd:
        inps.num_cmd = inps.num_processor

    inps.memory *= inps.num_cmd

    return inps


#####################################################################################
def split_run_file(run_file, num_line_per_job):
    # get number of lines in the input file
    with open(run_file, 'r') as f:
        for i, l in enumerate(f):
            pass
        num_line = i + 1
    print('number of lines: {}'.format(num_line))
    if num_line == 0:
        raise ValueError('input run file is empty.')

    # split run_file according to the num_processor
    num_job = int((num_line - 0.5) / num_line_per_job) + 1
    num_digit = len(str(num_job))  #length of suffixes

    # Usage: split [OPTION]... [INPUT [PREFIX]]
    prefix = 'z_input_{}.'.format(os.path.basename(run_file))
    cmd = 'split -a {a} -l {l} -d {i} {p}'.format(a=num_digit, l=num_line_per_job, i=run_file, p=prefix)
    print(cmd)
    os.system(cmd)

    # search all split files
    zfiles = glob.glob(prefix+'*')
    return zfiles


def wait4jobs2finish(run_file, num_job):
    job_name = os.path.basename(run_file)
    file_pattern = 'z_output_{}.*.o'.format(job_name)

    proj_name = sensor.project_name2sensor_name(os.getcwd())[1]

    print('-'*50)
    print('sleeping until {} jobs are done for {}'.format(num_job, job_name))
    t_sec = 0
    num_file = len(glob.glob(file_pattern))
    while num_file < num_job:
        time.sleep(1)           #wait one second

        msg = '# of '
        if proj_name:
            msg += '{}/'.format(proj_name)
        msg += '{} files: {} / {} after {} mins'.format(file_pattern, num_file, num_job, int(t_sec/60))

        if t_sec <= 600:
            if t_sec % 60 == 0:
                print(msg)
        else:
            if t_sec % 300 == 0:
                print(msg)

        num_file = len(glob.glob(file_pattern))
        t_sec += 1

    print('-'*50)
    print('ALL {} jobs are done for {}'.format(num_job, job_name))

    m, s = divmod(t_sec, 60); h, m = divmod(m, 60)
    print('Total used time: {:02d} hours {:02d} mins {:02d} secs'.format(h,m,s))
    return


def move_zfiles2job_folder(run_file):
    """move z_* files into dedicated directory"""
    # prepare job folder
    job_folder = run_file+'_jobs'
    print('moving z_* files into '+job_folder)

    cmd = 'rm -rf {}'.format(job_folder)
    os.system(cmd)

    cmd = 'mkdir {}'.format(job_folder)
    os.system(cmd)

    # move z_input/output_* files
    job_name = os.path.basename(run_file)
    cmd = 'mv z_input_{j}.* z_output_{j}.* {d}'.format(j=job_name, d=job_folder)
    os.system(cmd)
    return


def remove_zfiles(run_file):
    job_name = os.path.basename(run_file)
    file_pattern1 = 'z_input_{}.*'.format(job_name)
    file_pattern2 = 'z_output_{}.*'.format(job_name)
    for file_pattern in [file_pattern1, file_pattern2]:
        zfiles = glob.glob(file_pattern)
        if len(zfiles) > 0:
            cmd = 'rm {} {}'.format(file_pattern1, file_pattern2)
            print(cmd)
            os.system(cmd)
    return


#####################################################################################
def main(iargs=None):
    inps = cmd_line_parse(iargs)

    job_dir = os.path.dirname(inps.run_file)
    os.chdir(job_dir)
    print('Go to job directory', job_dir)

    job_name = os.path.basename(inps.run_file)
    print('input run file: {}'.format(inps.run_file))

    # Remove existing z_input/output_* files
    remove_zfiles(inps.run_file)

    # split run_file according to the num_processor
    zfiles = split_run_file(inps.run_file, num_line_per_job=inps.num_cmd)
    num_file = len(zfiles)

    # job setting
    print('number of jobs to submit: {}'.format(num_file))
    print('job queue: {}'.format(inps.queue_name))
    print('walltime: {}'.format(inps.walltime))
    print('memory: {}'.format(inps.memory))
    print('number of processors: {}'.format(inps.num_processor))
    print('email: {}'.format(inps.email))

    # write and submit job files
    for zfile in zfiles:
        # write job setting
        job_num = zfile.split('.')[-1]
        job_file = zfile+'.sh'
        f = open(job_file,'w')
        f.write('#! /bin/tcsh\n')
        f.write('#BSUB -J {}.{}\n'.format(job_name, job_num))
        f.write('#BSUB -P {}\n'.format(inps.project_name))
        f.write('#BSUB -o z_output_{}.{}.%J.o\n'.format(job_name, job_num))
        f.write('#BSUB -e z_output_{}.{}.%J.e\n'.format(job_name, job_num))
        f.write('#BSUB -W {}\n'.format(inps.walltime))
        f.write('#BSUB -q {}\n'.format(inps.queue_name))
        f.write('#BSUB -n {}\n'.format(inps.num_processor))
        f.write('#BSUB -R "rusage[mem={}]"\n'.format(inps.memory))

        if inps.queue_name == 'parallel':
            f.write('#BSUB -R "span[ptile=16]"\n')

        if inps.email:
            f.write('#BSUB -u {}\n'.format(inps.email))
            f.write('#BSUB -N\n')

        # write cd work directory
        f.write('\n')
        f.write('cd {}\n'.format(job_dir))

        # write job excutable commands
        with open(zfile, 'r') as fz:
            for line in fz:
                f.write(line)
        f.close()

        # submit job file
        cmd = 'bsub < '+job_file
        print(cmd)
        os.system(cmd)

        cmd  = 'rm {}'.format(zfile)
        os.system(cmd)

    wait4jobs2finish(inps.run_file, num_job=num_file)

    move_zfiles2job_folder(inps.run_file)

    print('finished at {}'.format(datetime.datetime.now()))
    return


#####################################################################################
if __name__ == '__main__':
    main()
