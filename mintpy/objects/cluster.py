from dask_jobqueue import LSFCluster, PBSCluster, SLURMCluster


def get_cluster(cluster_type, **kwargs):
    print("Using cluster type: {}".format(cluster_type))
    print("Using config name: {}".format(kwargs['config_name']))

    cluster_type = cluster_type.lower()
    kwargs['walltime'] = format_walltime_for_cluster(cluster_type, kwargs["walltime"])

    if cluster_type == 'lsf':
        cluster = LSFCluster(**kwargs)
    elif cluster_type == 'pbs':
        cluster = PBSCluster(**kwargs)
    elif cluster_type == 'slurm':
        cluster = SLURMCluster(**kwargs)
    else:
        raise ValueError("Cluster type '{}' not supported".format(cluster_type))

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
