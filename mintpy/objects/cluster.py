import time
import dask
from dask.distributed import LocalCluster, Client, as_completed


class DaskCluster:

    def __init__(self, cluster_type, num_workers, box, write_job_script=True, **kwargs):
        """Generic dask cluster wrapper"""

        self.cluster = None
        self.client = None
        self.box = box
        self.num_workers = num_workers

        # Properly format cluster type for consistency
        cluster_type = cluster_type.lower()
        cluster_list = ['lsf', 'pbs', 'slurm', 'local']
        if cluster_type not in cluster_list:
            msg = "Cluster type '{}' not supported".format(cluster_type)
            msg += '\nSupported cluster types: {}'.format(cluster_list)
            raise ValueError(msg)
        print("Dask cluster type: {}".format(cluster_type))

        # import related dask module only
        # because job_queue is not available on macports, which make sense
        if cluster_type is not 'local':
            try:
                import dask_jobqueue as jobqueue
            except ImportError:
                raise ImportError('Cannot import dask_jobqueue!')

        # for local cluster, NO need to do the extra configuration
        if cluster_type == 'local':
            self.cluster = LocalCluster()
            return

        # check input config name
        if 'config_name' in kwargs.keys():
            kwargs['config_name'] = self.format_config_name(kwargs['config_name'], cluster_type)
        print("Dask config name: {}".format(kwargs['config_name']))

        # check walltime format for each cluster type
        if 'walltime' in kwargs.keys():
            kwargs['walltime'] = self.format_walltime(kwargs["walltime"], cluster_type)
        print('Dask worker walltime: {}'.format(kwargs['walltime']))

        # initiate cluster object
        if cluster_type == 'lsf':
            self.cluster = jobqueue.LSFCluster(**kwargs)
        elif cluster_type == 'pbs':
            self.cluster = jobqueue.PBSCluster(**kwargs)
        elif cluster_type == 'slurm':
            self.cluster = jobqueue.SLURMCluster(**kwargs)

        self.cluster.scale(num_workers)

        if write_job_script:
            self.write_job_script()

    def write_job_script(self):
        # Print and write job command file for HPC cluster types
        print("JOB COMMAND CALLED FROM PYTHON:\n\n", self.cluster.job_script())
        with open('dask_command_run_from_python.txt', 'w') as f:
            f.write(self.cluster.job_script() + '\n')

    @staticmethod
    def format_config_name(config_name, cluster_type):

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

    @staticmethod
    def format_walltime(walltime, cluster_type):
        """format the walltime str for different clusters
        HH:MM:SS - pbs / slurm
        HH:MM    - lsf
        """

        num_digit = len(walltime)
        if cluster_type in ['pbs', 'slurm'] and num_digit != 8:
            if num_digit == 5:
                walltime += ":00"
            else:
                raise ValueError("input walltime ({}) not in HH:MM:SS or HH:MM format.".format(walltime))

        elif cluster_type in ['lsf'] and num_digit != 5:
            if num_digit == 8:
                walltime = walltime[:5]
            else:
                raise ValueError("input walltime ({}) not in HH:MM:SS or HH:MM format.".format(walltime))

        return walltime

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

        # This line needs to be in a function or in a `if __name__ == "__main__":` block. If it is in no function
        # or "main" block, each worker will try to create its own client (which is bad) when loading the module
        print('Initiating dask client')
        self.client = Client(self.cluster)

        # split the primary box into sub boxes for each worker
        sub_boxes = self.split_box2sub_boxes(self.box, num_split=self.num_workers, dimension='x')
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

        return futures, start_time_sub