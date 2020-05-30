import glob
import os
import shutil
import time
import dask
from dask.distributed import LocalCluster, Client, as_completed


class DaskCluster:

    def __init__(self, cluster_type, num_workers, **kwargs):
        """Generic dask cluster wrapper"""

        self.cluster = None
        self.client = None
        self.num_workers = num_workers

        # Properly format cluster type for consistency
        cluster_type = cluster_type.lower()
        cluster_list = ['lsf', 'pbs', 'slurm', 'local']
        if cluster_type not in cluster_list:
            msg = "Cluster type '{}' not supported".format(cluster_type)
            msg += '\nSupported cluster types: {}'.format(cluster_list)
            raise ValueError(msg)
        print("Dask cluster type: {}".format(cluster_type))


        # for local cluster, NO need to do the extra configuration
        if cluster_type == 'local':
            self.cluster = LocalCluster()
        else:

            # import related dask module only
            # because job_queue is not available on macports, which make sense
            try:
                import dask_jobqueue as jobqueue
            except ImportError:
                raise ImportError('Cannot import dask_jobqueue!')

            # check input config name
            if 'config_name' in kwargs.keys():
                kwargs['config_name'] = self.format_config_name(kwargs['config_name'], cluster_type)
            print("Dask config name: {}".format(kwargs['config_name']))

            # initiate cluster object
            if cluster_type == 'lsf':
                self.cluster = jobqueue.LSFCluster(**kwargs)
            elif cluster_type == 'pbs':
                self.cluster = jobqueue.PBSCluster(**kwargs)
            elif cluster_type == 'slurm':
                self.cluster = jobqueue.SLURMCluster(**kwargs)

            print("\n", self.cluster.job_script())
            debug_mode=False
            if debug_mode:
                self.write_job_script()

        self.cluster.scale(num_workers)

    def write_job_script(self):
        """ Writes the dask cluster job script to a file for reference. """
        # Print and write job command file for HPC cluster types
        with open('dask_command_run_from_python.txt', 'w') as f:
            f.write(self.cluster.job_script() + '\n')

    @staticmethod
    def format_config_name(config_name, cluster_type):
        """ Formats dask config_name property based on presence or absence of user specified config name.
        :param config_name: str, the user specified config name to use
        :param cluster_type: str, the type of HPC cluster being used (slurm, lsf, pbs)
        :return: config_name: str, the config_name formatted as follows:
                 - the user specified config name if its exists in $DASK_CONFIG/dask.yaml
                 - the default cluster_type config in $DASK_CONFIG/dask.yaml
        """

        # check if config_name exists in ~/.config/dask/*.yaml
        config_names = list(dask.config.get('jobqueue').keys())
        if config_name not in config_names:
            config_location = dask.config.get('config')
            msg = 'Dask configuration "{}" was not found in {}'.format(config_name, config_location)
            msg += '\nFalling back to default config name: "{}"'.format(cluster_type)
            print(msg)

            # due to the pre-set in dask.yaml, default config_name is the same as cluster_type
            config_name = cluster_type

        return config_name

    def split_box2sub_boxes(self, box, num_split, dimension='x'):
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

    def submit_workers(self, func, func_data):
        """ Submits dask workers to the networking client that run the specified function (func)
        on the specified data (func_data). Each dask worker is in charge of a small subbox of the main box.

        :param func: function, a python function to run in parallel
        :param func_data: dict, a dictionary of the argument to pass to the function
        :return: futures: [Future], a list of dask futures containing the workers and their results
                 start_time_sub: time, the initial starting time for the futures (used for runtime profiling)
                 box: [x0, y0, x1, y1]: list[int], the complete box being worked on
        """

        print('Initiating dask client')
        self.client = Client(self.cluster)

        # split the primary box into sub boxes for each worker
        box = func_data["box"]
        sub_boxes = self.split_box2sub_boxes(box, num_split=self.num_workers, dimension='x')
        print('Split patch into {} sub boxes in x direction for workers to process'.format(len(sub_boxes)))

        # submit jobs for each worker
        start_time_sub = time.time()
        futures = []
        for i, sub_box in enumerate(sub_boxes):
            print('Submit job to workers for sub box {}: {}'.format(i, sub_box))
            func_data['box'] = sub_box

            # David: I haven't played with fussing with `retries`, however sometimes a future fails
            # on a worker for an unknown reason. retrying will save the whole process from failing.
            # TODO:  I don't know what to do if a future fails > 3 times. I don't think an error is
            # thrown in that case, therefore I don't know how to recognize when this happens.
            future = self.client.submit(func, **func_data, retries=3)
            futures.append(future)

        return futures, start_time_sub, box

    def compile_workers(self, futures, start_time_sub, box, master_result_boxes):
        """ Compiles results from completed workers and recompiles their sub outputs into the output
        for the complete box being worked on.

        :param futures: [Future], the to-be-completed dask futures
        :param start_time_sub: time, the starting time of the submitted futures
        :param box: [x0, y0, x1, y1]: list[int], the dimensions of the complete box
        :param master_result_boxes: list[numpy.nd.array], arrays of the appropriate structure representing
               the final output of processed box (need to be in the same order as the function passed in
               submit_workers returns in)
        :return: master_result_boxes: tuple(numpy.nd.arrays), the processed results of the box
        """

        i_future = 0
        for future, result in as_completed(futures, with_results=True):

            # message
            i_future += 1
            sub_t = time.time() - start_time_sub
            print("FUTURE #{} complete. Time used: {:.0f} seconds".format(i_future, sub_t))

            # catch result
            result_list = list(result)

            # the box always needs to be the final return item
            sub_box = result_list.pop()

            # Loop across all of the returned data in order to rebuild complete box
            for i, subresult in enumerate(result_list):

                # convert the abosulte sub_box into local col/row start/end relative to the primary box
                # to assemble the result from each worker
                x0, y0, x1, y1 = sub_box
                x0 -= box[0]
                x1 -= box[0]
                y0 -= box[1]
                y1 -= box[1]

                master_result_box = master_result_boxes[i]
                dim = subresult.ndim  # number of dimensions of dataset (should be 2 or 3)
                if dim == 3:
                    master_result_box[:, y0:y1, x0:x1] = subresult
                elif dim == 2:
                    master_result_box[y0:y1, x0:x1] = subresult
                else:
                    raise Exception("subresult has unexpected dimension {}".format(subresult.ndim))

        return tuple(master_result_boxes)

    def run(self, func, func_data, master_result_boxes):
        """ Wrapper function encapsulating submit_workers and compile_workers.

        :param func: function, a python function to run in parallel
        :param func_data: dict, a dictionary of the argument to pass to the function
        :param master_result_boxes: list[numpy.nd.array], arrays of the appropriate structure representing
               the final output of processed box (need to be in the same order as the function passed in
               submit_workers returns in)
        :return: master_result_boxes: tuple(numpy.nd.arrays), the processed results of the box
        """
        futures, start_time_sub, box = self.submit_workers(func, func_data)
        return self.compile_workers(futures, start_time_sub, box, master_result_boxes)

    def shutdown(self):
        """ Closes connnection to cluster and client objects. """
        self.cluster.close()
        self.client.close()

    def move_dask_stdout_stderr_files(self):
        """ move  *o and *e files produced by dask into stdout and sderr directory """

        stdout_files = glob.glob('*.o')
        stderr_files = glob.glob('*.e')
        job_files = glob.glob('dask_command_run_from_python.txt*')

        if len(stdout_files + stderr_files + job_files) == 0:
            return

        stdout_folder = 'stdout_ifgram_inversion_dask'
        stderr_folder = 'stderr_ifgram_inversion_dask'
        for std_dir in [stdout_folder, stderr_folder]:
            if os.path.isdir(std_dir):
                shutil.rmtree(std_dir)
            os.mkdir(std_dir)

        for item in stdout_files + job_files:
            shutil.move(item, stdout_folder)

        for item in stderr_files:
            shutil.move(item, stderr_folder)

        return

    def cleanup(self):
        """ Closes connections to dask client and cluster and moves dask output/error files. """
        self.shutdown()
        self.move_dask_stdout_stderr_files()
