from dask_jobqueue import LSFCluster, PBSCluster, SLURMCluster


class Cluster:

    def __init__(self, type, **kwargs):

        self.cluster = None

        print("Using cluster type: {}".format(type))

        if type == 'LSF':
            self.cluster = LSFCluster(**kwargs)
        elif type == 'PBS':
            self.cluster = PBSCluster(**kwargs)
        elif type == 'SLURM':
            self.cluster = SLURMCluster(**kwargs)

