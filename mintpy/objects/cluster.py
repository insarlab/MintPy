import time
import dask
from dask.distributed import LocalCluster, Client


class DaskCluster:

    def __init__(self, cluster_type, write_job_script=True, **kwargs):
        """Generic dask cluster wrapper"""

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
