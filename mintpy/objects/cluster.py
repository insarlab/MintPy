#!/usr/bin/env python3
############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Joshua Zahner, March 2020                        #
############################################################


try:
    import dask
    from dask_jobqueue import LSFCluster, PBSCluster, SLURMCluster
    from dask.distributed import LocalCluster
except ImportError:
    raise ImportError('Cannot import dask or dask_jobqueue!')


def get_cluster(cluster_type, **kwargs):
    """Generic dask cluster wrapper"""

    # check input cluster type
    cluster_type = cluster_type.lower()
    cluster_list = ['lsf','pbs','slurm','local']
    if cluster_type not in cluster_list:
        msg = "Cluster type '{}' not supported".format(cluster_type)
        msg += '\nsupported cluster types: {}'.format(cluster_list)
        raise ValueError(msg)
    print("Dask cluster type: {}".format(cluster_type))

    # check input config name
    if 'config_name' in kwargs.keys():
        kwargs['config_name'] = check_config_name(kwargs['config_name'], cluster_type)
    print("Dask config name: {}".format(kwargs['config_name']))

    # check walltime format for each cluster type
    if 'walltime' in kwargs.keys():
        kwargs['walltime'] = check_walltime_format(kwargs["walltime"], cluster_type)
    print('Dask worker walltime: {}'.format(kwargs['walltime']))

    # initiate cluster object
    if cluster_type == 'lsf':
        cluster = LSFCluster(**kwargs)
    elif cluster_type == 'pbs':
        cluster = PBSCluster(**kwargs)
    elif cluster_type == 'slurm':
        cluster = SLURMCluster(**kwargs)
    else:
        cluster = LocalCluster()

    return cluster


def check_config_name(config_name, cluster_type):
    # due to the pre-set in mintpy.yaml, default config_name is the same as cluster_type
    if not config_name:
        config_name = cluster_type

    # check if config_name exists in ~/.config/dask/*.yaml
    config_names = list(dask.config.get('jobqueue').keys())
    if config_name not in config_names:
        msg = 'Dask configuration "{}" was not found in ~/.config/dask/*.yaml'.format(config_name)
        msg += '\nFall back to default config name: "{}"'.format(cluster_type)
        print(msg)
        config_name = cluster_type

    return config_name


def check_walltime_format(walltime, cluster_type):
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
