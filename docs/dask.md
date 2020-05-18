# Run and configure `dask` and `dask_jobqueue` #

`dask` and `dask_jobqueue`, maintained by the [dask organization](https://dask.org/), provide the ability to run python code in parallel across multiple machines or across multiple computing cores on a single machine. The base `dask` library provides wrappers around two common python libraries (`numpy` and `pandas`) to perform common computations, such as data transformations and matrix operations, in parallel. The `dask_jobqueue` library adds the ability to submit data and python commands to an High Performance Cluster (HPC) job scheduler, such as `LSF` and `SLURM`. Internally, dask handles waiting for each task to complete (since jobs could finish in a different order they were submitted) and properly merging the outputs back together. 

`MintPy` uses dask for parallel processing (currently implemeted only in `ifgram_inversion.py`) in two ways: 

+ locally on a single machine with multiple CPU cores via `Dask.distributed.LocalCluster`. 
+ on a distributed HPC via job schedulera defined in `dask_jobqueue`. Various configuration parameters is specified in a dask configuration file, i.e. `~/.config/dask/dask.yaml`, to determine the amount of computing resources that allocated to each dask “worker” by the job scheduler.

The required options and recommended best practices for each cluster/scheduler differ slightly and are covered in the following section.

## via `LocalCluster` ##

Recommended if you are running `MintPy` on a local machine with multiple available cores, or if you have access to an HPC cluster but wish to allocate only a single node's worth of resources towards MintPy but still take advantage of multiple cores. For that, simply set the following options in `smallbaselineApp.cfg`:

```cfg
mintpy.networkInversion.parallel   = yes
mintpy.networkInversion.cluster    = local
mintpy.networkInversion.numWorkers = auto   #auto for 4 (local) or 40 (non-local), set to "all" to use all available cores.
```

"mintpy.networkInversion.numWorkers = all" will allocate `multiprocessing.cpu_count()` number of workers to the dask computation.

### Test on Stampede2 ###

To show the run time improvement, we test three datasets (Galapagos, Fernandina, and Kilauea) with different number of cores on a compute node in the [Stampede2 cluster's skx-normal queue](https://portal.tacc.utexas.edu/user-guides/stampede2#overview-skxcomputenodes). Results are as below:
 
| Property              | Fernandina Dataset | Galapagos Dataset | Kilauea Dataset            |
|-----------------------|--------------------|-------------------|----------------------------|
| Input File Size       | 0.623 GB           | 0.233 GB          | 15.0 GB                    |
| Dask Configuration    | None               | None              | None                       |
| Cluster Memory Used   | 16 GB              | 16 GB             | 16 GB (64 GB for 1-2 cores)|
| Chunk Size Used       | 100e6 (1 GB)       | 100e6 (1 GB)      | 100e6 (1 GB)               | 

<p align="center">
  <img src="https://github.com/insarlab/MintPy-tutorial/blob/master/docs/dask_local_cluster_performance.png">
</p>

## via `dask_jobqueue` on HPC ##

The dask_jobqueue cluster objects (`LSFCluster`, `PBSCluster`, `SLURMCluster`, etc.) accepts configuration parameters in two ways: (a) passed directly to the object within the code, or (b) specified in a YAML file and be sourced by dask at runtime. MintPy uses option (b). Any `*.yaml` file in the `~/.config/dask/` directory will be sourced and used by dask regardless of its name. `dask.yaml`,  `distributed.yaml` and `jobqueue.yaml` file will be created in this directory by default during dask installation. 

We provide a MintPy-specific YAML file `$MINTPY_HOME/mintpy/defaults/mintpy.yaml`. It is recommended to copy over to the ~/.config/dask/ directory before running MintPy in parallel on your HPC system. One can then modify this file according to their own vailable computing environment, i.e. naming the configuration by 1) the job scheduler type (lsf, pbs, slurm) or 2) HPC name (stampede, comet, pegasus, etc.). The latter is useful if you are using MintPy on multiple HPC systems or in different resource schemes.

```bash
cp $MINTPY_HOME/mintpy/defaults/mintpy.yaml ~/.config/dask/mintpy.yaml
```

Note on `DASK_ROOT_CONFIG`: if you would like to not use the ~/.config/dask/ directory to store your configuration files for some reason, you can also set the `DASK_ROOT_CONFIG` environment variable to a custom directory. Any YAML files in this directory will be searched and added to the dask configuration object for the dask_jobqueue cluster object. However, files in `DASK_ROOT_CONFIG` directory has lower priority than the `~/.config/dask/` directory. **It is generally not recommended to set a custom DASK_ROOT_CONFIG variable.**

### Configuration parameters in `~/.config/dask/mintpy.yaml` ###

**name** - Name of the worker job as it will appear to the job scheduler. Any values are perfectly fine.

**cores** - Total number of cores you would like your parallel code to use. Note that this is not the number of cores per job or per worker, but the TOTAL number of cores you code can possibly use. It is recommended to check with your HPC cluster admin, or your HPC cluster documentation to find the number of cores that are available per compute node on your system. The value of this parameter is not recommended to exceed the total cores per compute node. **Note that this is a required parameter for dask_jobqueue cluster objects and must have a non-null value.**

**memory** - Total amount of memory you would like your code to use. Like the `cores` parameter, this is not memory per job or per worker, but total memory you code could use, and it is recommended to check what the maximum memory available on a compute node on you system is before setting this value. For `MintPy` is recommended to set this to a reasonably high value so as to ensure the dask workers do not run our memory while performing their calculations, which will result in the processing routine failing. A good metric found by testing on several `SLURM` clusters is to allocate at least 2GB of memory per core for smaller datasets (such as the `MintPy` test data) and as much as 5GB per core for much larger datasets. **Note that this is a required parameter for dask_jobqueue cluster objects and must have a non-null value. For systems that do not use the —mem directive in their job specifications, please see the section on the ‘header-skip’.**

**processes** - Number of python processes you would like your code to use PER JOB. For compute nodes with many cores, each core can often handle multiple running processes at a time which allow for faster computations and lower running times, without effecting slowing down the queuing time by requesting additional physical resources (cores and memory). If left as `null`, this value will automatically be set to `~sqrt(cores)` so that the number of processes and threads per process are approximately equal.  Processes can be thought of as the number of workers a given compute core can handle simultaneously, so, by setting `processes: 1` as in the default `mintpy.yaml` file, each core handles a single worker at a time, and then begins working on the next one after completing the previous one. If `processes: 5`, each core works on 5 workers at a time, and begins handling another worker as soon as any one of the previous 5 have completed. **It is generally recommended by dask to keep the number of processes low and instead to request additional physical resources, as specifying too many processes on some systems has been known to lead to failed jobs and workers.**

**walltime** - Maximum amount of time your job is permitted to run on a compute node before it is killed by the node. Walltimes are required so as to ensure users and jobs cannot hog nodes for extended periods of time without doing any substantial computations. Most queues have maximum allowed compute times, after which still running jobs will be terminated, but this is often on the order of 48+ hours. It is recommended to use relatively short walltimes for `MintPy` (the test dataset uses `wall time: 00:30:00` (30 minutes)) to have a higher priority in the job queue, but not too short, or your job may never connect to the dask scheduler and finish its computation (when in doubt, start small and work towards larger wall times as needed). Note that the walltime format changes slightly depending on the job scheduler that is in place on your HPC system. For `LSF` schedulers, walltime is accepted and expected to be in `HH:MM` format, while `SLURM` and `PBS` schedulers only accept walltime in `HH:MM:SS` format.

**queue** - Scheduler queue to submit your jobs to. Most systems have several different queues with different amounts and types of resources you can request. For this reason, it is important to check with your HPC cluster admin or documentation to determine the most appropriate queue to `MintPy` on and what resources are available on that queue. Requesting for more resources (core, memory, nodes, etc.) on a queue than are available often leads to your jobs being rejected from queue, and can, on some systems, lead to fines or account suspensions.

**project** - Project allocation associated with your account. Most HPC systems bill projects based on the amount of resources requested, or allocate projects certain allowances based on grants or other funding, so it is important to specify the proper project for your jobs to be associated with.

**interface** - Network interface being used by your HPC cluster. For all tested clusters, this has been the `ib0` interface. Available network interfaces can be fond by running `ipconfig` on your login node.

**job-extra** - this allows you to specify additional scheduler parameters that are not support out of the box by dask. For example, some `SLURM` clusters require the presence of the `--nodes` SBATCH directive in order to be accepted by the job scheduler, but dask does not use this directive directly. It can instead be specified in the `job-extra` array, and will be written with the other directive at submit time.

*There are a multitude of other configuration options that can be specified to dask to further customize the way in which jobs are run and executed on your HPC cluster, but the above are the most commonly used ones for MintPy. For further details and info on possible config options, see the [dask_jobqueue documentation](https://jobqueue.dask.org/en/latest/configuration.html)*

### Control parameters in `smallbaselineApp.cfg` ###

There are 5 options related to the dask features that can be controlled in MintPy via the `smallbaselineApp.cfg` file:

```cfg
mintpy.networkInversion.parallel  = auto #[yes / no], turn ON or OFF the parallel processing with dask, auto for no.
mintpy.networkInversion.cluster   = auto #[lsf / pbs / slurm / local], job scheduler in your HPC system auto for local.
mintpy.networkInversion.config    = auto #[name / no], name of the configuration section in YAML file, auto for no (to use the same name as the cluster type specified above)
mintpy.networkInversion.numWorker = auto #[int > 1], number of worker to submit and run, auto for 4 (local) or 40 (non-local), set to "all" to use all available cores.
mintpy.networkInversion.walltime  = auto #[HH:MM], walltime to be used for each dask job, auto for 00:40.
```
