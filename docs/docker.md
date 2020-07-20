## Running as Docker container

Thanks to Andre Theron for putting together an [MintPy container on DockerHub](https://hub.docker.com/r/andretheronsa/mintpy). [Docker](https://docs.docker.com/get-started/) allows you to run MintPy in a dedicated container (essentially an efficient virtual machine). [Here](https://docs.docker.com/install/) is the instruction to install docker.

To pull the MintPy container from Dockerhub to your local machine: 

```
docker pull andretheronsa/mintpy:latest
```

To start an interactive shell session in the container from the terminal, with bash for example: 

```
docker run -it andretheronsa/mintpy:latest bash
```

To map data on the host (local) machine to the container use [volumes](https://docs.docker.com/storage/volumes/):

```
docker run -it -v /path/to/data/dir:/home/work/ andretheronsa/mintpy:latest bash
```

Background processing is possible using something like:  

```
docker run -it -v /path/to/data/dir:/home/work/ andretheronsa/mintpy:latest python /home/python/MintPy/mintpy/smallbaselineApp.py /home/work/smallbaselineApp.cfg
```

### Notes ###

+ The container may have strong permissions for directories you map to it.   

+ Container was built on `insarlab/main` - should be updated with new releases.  

+ Needs further testing and improvement - can be made smaller (use Alpine instead of Debian...)  
