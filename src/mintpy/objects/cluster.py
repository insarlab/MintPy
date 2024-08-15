"""Class wrapped around Dask for parallel computing."""
#############################################################
# Program is part of MintPy                                 #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi          #
# Author: Joshua Zahner, David Grossman, Zhang Yunjun, 2020 #
#############################################################
# Recommend import:
#   from mintpy.objects import cluster


import glob
import os
import shutil
import time

# supported / tested clusters
CLUSTER_LIST = ['lsf', 'pbs', 'slurm', 'local']
NUM_THREADS_ENV_LIST = [
    'OMP_NUM_THREADS',         # openmp
    'OPENBLAS_NUM_THREADS',    # openblas
    'MKL_NUM_THREADS',         # mkl
    'NUMEXPR_NUM_THREADS',     # numexpr
    'VECLIB_MAXIMUM_THREADS',  # accelerate
]


############################## Utilities functions #########################################

def split_box2sub_boxes(box, num_split, dimension='x', print_msg=False):
    """Divide the input box into `num_split` different sub_boxes.

    :param box: [x0, y0, x1, y1]: list[int] of size 4
    :param num_split: int, the initial number of sub_boxes to split a box into
    :param dimension: str = 'y' or 'x', the dimension along which to split the boxes
    :return: sub_boxes: list(list(4 int)), the splited sub boxes
    :return: num_split: int, the final number of split sub_boxes
    """
    import numpy as np

    dimension = dimension.lower()
    if num_split <= 1:
        return [box], num_split

    # basic info
    x0, y0, x1, y1 = box
    length, width = y1 - y0, x1 - x0

    # calc step
    if dimension == 'y':
        dim_size = length
    else:
        dim_size = width
    step = int(np.ceil(dim_size / num_split))
    # condition: step >= 10
    step = max(step, 10)
    # update num_split based on the final step size
    num_split = int(np.ceil(dim_size / step))
    # if the last step is too small, merge it into the 2nd last one
    last_step = dim_size - step * (num_split - 1)
    if last_step < step * 0.05 or last_step < 5:
        num_split -= 1

    # get list of boxes
    sub_boxes = []
    for i in range(num_split):
        if dimension == 'y':
            r0 = y0 + step * i
            r1 = y0 + step * (i + 1)
            if i == num_split - 1:
                r1 = y1
            sub_boxes.append([x0, r0, x1, r1])

        else:
            c0 = x0 + step * i
            c1 = x0 + step * (i + 1)
            if i == num_split - 1:
                c1 = x1
            sub_boxes.append([c0, y0, c1, y1])

    if print_msg:
        print(f'split along {dimension} dimension ({dim_size:d}) into {num_split:d} boxes')
        print(f'    with each box up to {step:d} in {dimension} dimension')

    return sub_boxes, num_split


def set_num_threads(num_threads=None, print_msg=True):
    """limit/set the number of threads for all environmental variables to the given value
    and save/return the original value for backup purpose.
    Link: https://stackoverflow.com/questions/30791550

    Parameters: num_threads      - str, number of threads
                                   Set to None to return without changing env variables
    Returns:    num_threads_dict - dict, dictionary of the original number of threads
    """

    # grab the original number of threads
    if print_msg:
        print(f'save the original settings of {NUM_THREADS_ENV_LIST}')
    num_threads_dict = {}
    for key in NUM_THREADS_ENV_LIST:
        num_threads_dict[key] = os.environ.get(key, None)

    # change the env variables
    if num_threads:
        num_threads = str(num_threads)
        for key in NUM_THREADS_ENV_LIST:
            os.environ[key] = num_threads
            if print_msg:
                print(f'set {key} = {num_threads}')

    return num_threads_dict


def roll_back_num_threads(num_threads_dict, print_msg=True):
    """Set back the number of threads for all environmental variables."""
    if print_msg:
        print(f'roll back to the original settings of {NUM_THREADS_ENV_LIST}')
    for key, value in num_threads_dict.items():
        if key in os.environ.keys():
            if value is None:
                os.environ.pop(key)
                if print_msg:
                    print(f'remove env variable {key}')
            else:
                os.environ[key] = value
                if print_msg:
                    print(f'set {key} = {value}')
    return



############################## Beginning of DaskCluster class ##############################

class DaskCluster:
    """
    Generic dask cluster wrapper for parallel processing in blocks.

    This object takes in a computing function for one block in space.
    For the computing function:
        1. the output is always several matrices (or None) and one box.
        2. the number of matrices may vary for different applications/functions.
        3. all matrices will be in 2D in size of (len, wid) or 3D in size of (n, len, wid),
           thus, the last two dimension (in space) will be the same.
    This charateristics allows the automatic result collection without prior knowledge
        of the computing function, thus being a generic wrapper.

    Check ifgram_inversion.py as an example.

    """

    def __init__(self, cluster_type, num_worker, config_name=None, **kwargs):
        """Initiate object
        :param cluster_type: str, cluster to use (local, slurm, lsf, pbs)
        :param num_worker: str, number of workers to use
        :param config_name: str, the name of configuratino section
        :other param **kwargs: dask configuration parameters
                 e.g. config_name: str, the user specified config name to use
        """

        self.cluster_type = cluster_type.lower()
        self.num_worker = num_worker
        self.config_name = config_name
        self.cluster_kwargs = kwargs

        ## format input arguments
        # num_worker
        self.num_worker = self.format_num_worker(self.cluster_type, self.num_worker)

        # config_name
        self.format_config_name()
        self.cluster_kwargs['config_name'] = self.config_name

        ## printout message
        print(f"input Dask cluster type: {self.cluster_type}")
        if self.config_name is not None:
            print(f"input Dask config name: {self.config_name}")

        ## initial value
        self.cluster = None
        self.client = None


    def open(self):
        """Initiate the cluster"""

        # initiate the cluster object
        # Look at the ~/.config/dask/mintpy.yaml file for changing the Dask configuration defaults
        print('initiate Dask cluster')
        if self.cluster_type == 'local':
            from dask.distributed import LocalCluster

            # initiate cluster object
            self.cluster = LocalCluster()

        else:
            # for non-local cluster, import related dask module only when it's needed
            # because job_queue is not available on macports, which make sense
            import dask_jobqueue

            # initiate cluster object
            if self.cluster_type == 'lsf':
                self.cluster = dask_jobqueue.LSFCluster(**self.cluster_kwargs)

            elif self.cluster_type == 'pbs':
                self.cluster = dask_jobqueue.PBSCluster(**self.cluster_kwargs)

            elif self.cluster_type == 'slurm':
                self.cluster = dask_jobqueue.SLURMCluster(**self.cluster_kwargs)

            else:
                msg = f'un-recognized input cluster: {self.cluster_type}'
                msg += f'\nsupported clusters: {CLUSTER_LIST}'
                raise ValueError(msg)

            # show dask cluster job script for reference
            print("\n", self.cluster.job_script())
            # for debug
            debug_mode = False
            if debug_mode:
                with open('dask_command_run_from_python.txt', 'w') as f:
                    f.write(self.cluster.job_script() + '\n')


    def run(self, func, func_data, results):
        """Wrapper function encapsulating submit_workers and compile_workers.

        For a generic result collection without prior knowledge of the computing function,
        we assume that the output of "func" is: several 2D or 3D matrices + a box

        :param func: function, a python function to run in parallel
        :param func_data: dict, a dictionary of the argument to pass to the function
        :param results: list[numpy.ndarray], arrays of the appropriate structure representing
               the final output of processed box (need to be in the same order as the function passed in
               submit_workers returns in)
        :return: results: tuple(numpy.ndarray), the processed results of the box
        """
        from dask.distributed import Client

        # split the primary box into sub boxes for workers AND
        # update the number of workers based on split result
        box = func_data["box"]
        sub_boxes, self.num_worker = split_box2sub_boxes(
            box,
            num_split=self.num_worker,
            dimension='x',
            print_msg=False,
        )
        print(f'split patch into {self.num_worker} sub boxes in x direction for workers to process')

        # start a bunch of workers from the cluster
        print(f'scale Dask cluster to {self.num_worker} workers')
        self.cluster.scale(self.num_worker)

        print('initiate Dask client')
        self.client = Client(self.cluster)
        self.client.get_versions(check=True)

        # submit job for each worker
        futures, submission_time = self.submit_job(func, func_data, sub_boxes)

        # assemble results from all workers
        results = self.collect_result(futures, results, box, submission_time)

        return results


    def submit_job(self, func, func_data, sub_boxes):
        """Submit dask workers to the networking client that run the specified function (func)
        on the specified data (func_data). Each dask worker is in charge of a small subbox of the main box.

        :param func: function, a python function to run in parallel
        :param func_data: dict, a dictionary of the argument to pass to the function
        :param sub_boxes: list(np.ndarray), list of boxes to be computed in parallel

        :return futures: list(dask.Future), list of futures representing future dask worker calculations
        :return submission_time: time, the time of submission of the dask workers (used to determine worker
                runtimes as a performance diagnostic)
        """

        submission_time = time.time()
        futures = []
        for i, sub_box in enumerate(sub_boxes):
            print(f'submit a job to the worker for sub box {i}: {sub_box}')
            func_data['box'] = sub_box

            # David: I haven't played with fussing with `retries`, however sometimes a future fails
            # on a worker for an unknown reason. retrying will save the whole process from failing.
            # TODO:  I don't know what to do if a future fails > 3 times. I don't think an error is
            # thrown in that case, therefore I don't know how to recognize when this happens.
            future = self.client.submit(func, **func_data, retries=3)
            futures.append(future)

        return futures, submission_time


    def collect_result(self, futures, results, box, submission_time):
        """Compile results from completed workers and recompiles their sub outputs into the output
        for the complete box being worked on.
        :param futures: list(dask.Future), list of futures representing future dask worker calculations
        :param results: list[numpy.ndarray], arrays of the appropriate structure representing
               the final output of processed box (need to be in the same order as the function passed in
               submit_workers returns in)
        :param box: numpy.ndarray, the initial complete box being processed
        :param submission_time: time, the time of submission of the dask workers (used to determine worker
               runtimes as a performance diagnostic)
        :return: results: tuple(numpy.ndarray), the processed results of the box
        """
        from dask.distributed import as_completed

        num_future = 0
        for future, sub_results in as_completed(futures, with_results=True):

            # message
            num_future += 1
            sub_t = time.time() - submission_time
            print(f"\nFUTURE #{num_future} complete. Time used: {sub_t:.0f} seconds")

            # catch result - sub_box
            # and convert the absolute sub_box into local col/row start/end relative to the primary box
            # to assemble the result from each worker
            sub_box = sub_results[-1]
            x0, y0, x1, y1 = sub_box
            x0 -= box[0]
            x1 -= box[0]
            y0 -= box[1]
            y1 -= box[1]

            # catch result - matrices
            # and loop across all of the returned data to rebuild complete box
            for i, sub_result in enumerate(sub_results[:-1]):
                if sub_result is not None:
                    num_dim = sub_result.ndim
                    if num_dim == 4:
                        results[i][:, :,
                                   y0:y1, x0:x1] = sub_result
                    elif num_dim == 3:
                        results[i][:,
                                   y0:y1, x0:x1] = sub_result
                    elif num_dim == 2:
                        results[i][y0:y1, x0:x1] = sub_result
                    else:
                        msg = f"worker result has unexpected dimension: {num_dim}"
                        msg += '\nit should be either 2 or 3 or 4!'
                        raise Exception(msg)

        return results


    def close(self):
        """Close connections to dask client and cluster and moves dask output/error files. """

        # close client before cluster -> less likely to have the CancelledError
        # https://github.com/dask/distributed/issues/2273
        self.client.close()
        print('close dask client')

        self.cluster.close()
        print('close dask cluster')

        # move *.o/.e files produced by dask in stdout/stderr
        self.move_dask_stdout_stderr_files()


    ##### Utilities functions

    @staticmethod
    def format_num_worker(cluster_type, num_worker):
        """Format dask num_worker.
        :param cluster_type: str
        :param num_worker: str, number of workers to use
        :return: num_worker: int, number of workers to use
        """

        if cluster_type == 'local':
            num_core = os.cpu_count()

            # all / percentage --> num_core
            msg = f'numWorker = {num_worker}'
            if num_worker == 'all':
                ## divide by the number of threads per core [for Linux only]
                #import subprocess
                #from mintpy.utils import utils0 as ut0
                #if ut0.which('lscpu') is not None:
                #    # get the number of threads per core
                #    # link: https://stackoverflow.com/questions/62652951
                #    ps = subprocess.run(['lscpu'], capture_output=True, text=True).stdout.split('\n')
                #    ns = [p.split(':')[1].strip() for p in ps if p.startswith('Thread(s) per core:')]
                #    if len(ns) > 0:
                #        num_thread = int(ns[0])
                #        num_core = int(num_core / num_thread)

                # set num_worker to the number of cores
                num_worker = str(num_core)
                print(f'translate {msg} to {num_worker}')

            elif num_worker.endswith('%'):
                num_worker = int(num_core * float(num_worker[:-1]) / 100)
                print(f'translate {msg} to {num_worker}')
                if num_worker < 1 or num_worker >= num_core:
                    raise ValueError('Invalid numWorker percentage!')

            # str --> int
            num_worker = int(num_worker)

            # if num_worker > num_core,
            # then we assume that the user is not aware of the available resources
            # and use max(num_core/2, 1) instead to be conservative.
            if num_worker > num_core:
                print(f'\nWARNING: input number of worker: {num_worker} > available cores: {num_core}')
                num_worker = max(int(num_core / 2), 1)
                print(f'change number of worker to {num_worker} and continue\n')

        else:
            if num_worker == 'all':
                msg = f'numWorker = all is NOT supported for cluster type: {cluster_type}'
                raise ValueError(msg)
            num_worker = int(num_worker)

        return num_worker


    def format_config_name(self):
        """Format dask config_name property based on presence or absence of user specified config name.

        :return: config_name: str, the config_name formatted as follows:
                 - the user specified config name if its exists in $DASK_CONFIG/dask.yaml
                 - the default cluster_type config in $DASK_CONFIG/dask.yaml
        """
        import dask

        # config_name is not needed for local cluster
        if self.cluster_type == 'local':
            self.config_name = None
            return self.config_name

        # translate config_name = None to config_name = cluster_type
        if self.config_name is None:
            print('input config name is None, thus use the default (same as cluster type)')
            self.config_name = self.cluster_type

        # if config_name not found, use cluster_type as defined in minpy.dask
        config_names = list(dask.config.get('jobqueue').keys())
        if self.config_name not in config_names:
            config_location = dask.config.get('config')
            msg = f'Dask configuration "{self.config_name}" was not found in {config_location}'
            msg += f'\nFalling back to default config name: "{self.cluster_type}"'
            print(msg)
            self.config_name = self.cluster_type

        return self.config_name


    def move_dask_stdout_stderr_files(self):
        """Move *o and *e files produced by dask into stdout and sderr directory"""

        stdout_files = glob.glob('*.o')
        stderr_files = glob.glob('*.e')
        job_files = glob.glob('dask_command_run_from_python.txt*')

        if len(stdout_files + stderr_files + job_files) == 0:
            return

        stdout_folder = 'stdout_dask'
        stderr_folder = 'stderr_dask'
        for std_dir in [stdout_folder, stderr_folder]:
            if os.path.isdir(std_dir):
                shutil.rmtree(std_dir)
            os.mkdir(std_dir)

        for item in stdout_files + job_files:
            shutil.move(item, stdout_folder)

        for item in stderr_files:
            shutil.move(item, stderr_folder)

############################## End of DaskCluster class ####################################
