#!/usr/bin/env python3
#############################################################
# Program is part of MintPy                                 #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi          #
# Author: Joshua Zahner, David Grossman, Zhang Yunjun, 2020 #
#############################################################
# Recommend import:
#     from mintpy.objects.cluster import DaskCluster


import os
import glob
import shutil
import time
try:
    import dask
    from dask.distributed import LocalCluster, Client, as_completed
except ImportError:
    raise ImportError('Cannot import dask!')


# supported / tested clusters
CLUSTER_LIST = ['lsf', 'pbs', 'slurm', 'local']


class DaskCluster:
    """
    Generic dask cluster wrapper for parallel processing blocks in space or in time

    Check ifgram_inversion.py as an example.
    """

    def __init__(self, cluster_type, num_worker, config_name=None, **kwargs):
        """Initiate object
        :param cluster_type: str, cluster to use (local, slurm, lsf, pbs)
        :param num_worker: int, number of workers to use
        :other param **kwargs: dask configuration parameters
                 e.g. config_name: str, the user specified config name to use
        """

        self.cluster_type = cluster_type.lower()
        self.num_worker = num_worker
        self.config_name = config_name
        self.kwargs = kwargs

        # format input config name
        self.format_config_name()

        # printout message
        print("input Dask cluster type: {}".format(self.cluster_type))
        if self.config_name is not None:
            print("input Dask config name: {}".format(self.config_name))

        # intitial value
        self.cluster = None
        self.client = None

        return


    def open(self):
        """Initiate and scale the cluster"""

        # initiate the cluster object
        # Look at the ~/.config/dask/mintpy.yaml file for changing the Dask configuration defaults
        print('initiate Dask cluster')
        if self.cluster_type == 'local':
            self.cluster = LocalCluster()

        else:
            # for non-local cluster, import related dask module only when it's needed
            # because job_queue is not available on macports, which make sense
            try:
                import dask_jobqueue as jobqueue
            except ImportError:
                raise ImportError('Cannot import dask_jobqueue!')

            self.kwargs['config_name'] = self.config_name

            # initiate cluster object
            if self.cluster_type == 'lsf':
                self.cluster = jobqueue.LSFCluster(**self.kwargs)

            elif self.cluster_type == 'pbs':
                self.cluster = jobqueue.PBSCluster(**self.kwargs)

            elif self.cluster_type == 'slurm':
                self.cluster = jobqueue.SLURMCluster(**self.kwargs)

            else:
                msg = 'un-recognized input cluster: {}'.format(self.cluster_type)
                msg += '\nsupported clusters: {}'.format(CLUSTER_LIST)
                raise ValueError(msg)

            # show dask cluster job script for reference
            print("\n", self.cluster.job_script())
            # for debug
            debug_mode = False
            if debug_mode:
                with open('dask_command_run_from_python.txt', 'w') as f:
                    f.write(self.cluster.job_script() + '\n')

        # This line submits num_worker jobs to the cluster to start a bunch of workers
        # In tests on Pegasus `general` queue in Jan 2019, no more than 40 workers could RUN
        # at once (other user's jobs gained higher priority in the general at that point)
        print('scale the cluster to {} workers'.format(self.num_worker))
        self.cluster.scale(self.num_worker)

        return


    def run(self, func, func_data, results):
        """Wrapper function encapsulating submit_workers and compile_workers.

        For a generic result collection without prior knowledge of the computing function,
        we assume that the output of "func" is: several 2D or 3D matrices + a box

        :param func: function, a python function to run in parallel
        :param func_data: dict, a dictionary of the argument to pass to the function
        :param results: list[numpy.nd.array], arrays of the appropriate structure representing
               the final output of processed box (need to be in the same order as the function passed in
               submit_workers returns in)
        :return: results: tuple(numpy.nd.arrays), the processed results of the box
        """

        # This line needs to be in a function or in a `if __name__ == "__main__":` block. If it is in no function
        # or "main" block, each worker will try to create its own client (which is bad) when loading the module
        print('initiate Dask client')
        self.client = Client(self.cluster)

        # split the primary box into sub boxes for each worker
        self.box = func_data["box"]
        self.sub_boxes = self.split_box2sub_boxes(self.box, num_split=self.num_worker, dimension='x')
        print('split patch into {} sub boxes in x direction for workers to process'.format(len(self.sub_boxes)))

        # submit job for each worker
        self.submit_job(func, func_data)

        # assemble results from all workers
        return self.collect_result(results)


    def submit_job(self, func, func_data):
        """Submit dask workers to the networking client that run the specified function (func)
        on the specified data (func_data). Each dask worker is in charge of a small subbox of the main box.

        :param func: function, a python function to run in parallel
        :param func_data: dict, a dictionary of the argument to pass to the function
        """

        self.start_time_sub = time.time()
        self.futures = []
        for i, sub_box in enumerate(self.sub_boxes):
            print('submit a job to the worker for sub box {}: {}'.format(i, sub_box))
            func_data['box'] = sub_box

            # David: I haven't played with fussing with `retries`, however sometimes a future fails
            # on a worker for an unknown reason. retrying will save the whole process from failing.
            # TODO:  I don't know what to do if a future fails > 3 times. I don't think an error is
            # thrown in that case, therefore I don't know how to recognize when this happens.
            future = self.client.submit(func, **func_data, retries=3)
            self.futures.append(future)

        return


    def collect_result(self, results):
        """Compile results from completed workers and recompiles their sub outputs into the output
        for the complete box being worked on.

        :param results: list[numpy.nd.array], arrays of the appropriate structure representing
               the final output of processed box (need to be in the same order as the function passed in
               submit_workers returns in)
        :return: results: tuple(numpy.nd.arrays), the processed results of the box
        """

        num_future = 0
        for future, sub_results in as_completed(self.futures, with_results=True):

            # message
            num_future += 1
            sub_t = time.time() - self.start_time_sub
            print("FUTURE #{} complete. Time used: {:.0f} seconds".format(num_future, sub_t))

            # catch result - sub_box
            # and convert the abosulte sub_box into local col/row start/end relative to the primary box
            # to assemble the result from each worker
            sub_box = sub_results[-1]
            x0, y0, x1, y1 = sub_box
            x0 -= self.box[0]
            x1 -= self.box[0]
            y0 -= self.box[1]
            y1 -= self.box[1]

            # catch result - matrices
            # and loop across all of the returned data to rebuild complete box
            for i, sub_result in enumerate(sub_results[:-1]):
                num_dim = sub_result.ndim

                if num_dim == 3:
                    results[i][:, y0:y1, x0:x1] = sub_result

                elif num_dim == 2:
                    results[i][y0:y1, x0:x1] = sub_result

                else:
                    msg = "worker result has unexpected dimension: {}".format(num_dim)
                    msg += '\nit should be either 2 or 3!'
                    raise Exception(msg)

        return results


    def close(self):
        """Close connections to dask client and cluster and moves dask output/error files. """

        self.cluster.close()
        print('close dask cluster')

        self.client.close()
        print('close dask client')

        # move *.o/.e files produced by dask in stdout/stderr
        self.move_dask_stdout_stderr_files()
        return


    ##### Utilities functions

    @staticmethod
    def split_box2sub_boxes(box, num_split, dimension='x'):
        """Further divides the box size into `num_split` different sub_boxes.
        Note that this is different from `split2boxes()`, whic splits based on chunk_size (memory-based).

        :param box: [x0, y0, x1, y1]: list[int] of size 4
        :param num_split: int, the number of sub_boxes to split a box into
        :param dimension: str = 'y' or 'x', the dimension along which to split the boxes
        """

        x0, y0, x1, y1 = box
        length, width = y1 - y0, x1 - x0

        sub_boxes = []
        if dimension == 'y':
            for i in range(num_split):
                start = (i * length) // num_split + y0
                end = ((i + 1) * length) // num_split + y0
                sub_boxes.append([x0, start, x1, end])

        else:
            for i in range(num_split):
                start = (i * width) // num_split + x0
                end = ((i + 1) * width) // num_split + x0
                sub_boxes.append([start, y0, end, y1])

        return sub_boxes


    def format_config_name(self):
        """Formats dask config_name property based on presence or absence of user specified config name.

        :param config_name: str, the user specified config name to use
        :param cluster_type: str, the type of HPC cluster being used (slurm, lsf, pbs)
        :return: config_name: str, the config_name formatted as follows:
                 - the user specified config name if its exists in $DASK_CONFIG/dask.yaml
                 - the default cluster_type config in $DASK_CONFIG/dask.yaml
        """
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
            msg = 'Dask configuration "{}" was not found in {}'.format(self.config_name, config_location)
            msg += '\nFalling back to default config name: "{}"'.format(self.cluster_type)
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

        return

