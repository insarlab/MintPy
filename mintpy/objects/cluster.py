try:
    from dask_jobqueue import LSFCluster, PBSCluster, SLURMCluster
except ImportError:
    raise ImportError('Cannot import dask_jobqueue!')


def get_cluster(cluster_type, **kwargs):

    # check input cluster type
    cluster_type = cluster_type.lower()
    cluster_list = ['lsf','pbs','slurm']
    if cluster_type not in cluster_list:
        msg = "Cluster type '{}' not supported".format(cluster_type)
        msg += '\nsupported cluster types: {}'.format(cluster_list)
        raise ValueError(msg)

    print("Using cluster type: {}".format(cluster_type))
    print("Using config name: {}".format(kwargs['config_name']))

    # check walltime format for each cluster type
    kwargs['walltime'] = format_walltime_for_cluster(cluster_type, kwargs["walltime"])

    # initiate cluster object
    if cluster_type == 'lsf':
        cluster = LSFCluster(**kwargs)
    elif cluster_type == 'pbs':
        cluster = PBSCluster(**kwargs)
    elif cluster_type == 'slurm':
        cluster = SLURMCluster(**kwargs)

    return cluster


def format_walltime_for_cluster(cluster, wall):
    wall_len = len(wall)
    formatted_wall = wall
    if cluster in ['pbs', 'slurm'] and wall_len != 8:
        if wall_len == 5:
            formatted_wall = wall + ":00"
        else:
            raise ValueError("Walltime not in HH:MM:SS or HH:MM format.")
    elif cluster in ['lsf'] and wall_len != 5:
        if wall_len == 8:
            formatted_wall = wall[:5]
        else:
            raise ValueError("Walltime not in HH:MM:SS or HH:MM format.")

    return formatted_wall
