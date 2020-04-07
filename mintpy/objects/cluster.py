from dask_jobqueue import LSFCluster, PBSCluster, SLURMCluster


def get_cluster(cluster_type, **kwargs):
    print("Using cluster type: {}".format(type))
    print("Using config name: {}".format(kwargs['config_name']))

    cluster_type = cluster_type.lower()

    if cluster_type == 'lsf':
        cluster = LSFCluster(**kwargs)
    elif cluster_type == 'pbs':
        cluster = PBSCluster(**kwargs)
    elif cluster_type == 'slurm':
        cluster = SLURMCluster(**kwargs)

    return cluster

