from dask_jobqueue import LSFCluster, PBSCluster, SLURMCluster


def get_cluster(type, **kwargs):
    print("Using cluster type: {}".format(type))

    if type == 'LSF':
        cluster = LSFCluster(**kwargs)
    elif type == 'PBS':
        cluster = PBSCluster(**kwargs)
    elif type == 'SLURM':
        cluster = SLURMCluster(**kwargs)

    return cluster

