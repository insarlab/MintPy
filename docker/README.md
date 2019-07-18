## Running as Docker container

[Docker](https://docs.docker.com/get-started/) allows you to run MintPy in a dedicated container (essentially an efficient virtual machine). Follow the link to install [installed](https://docs.docker.com/install/).

Pull the container from Dockerhub to your local machine using: 

```
docker pull andretheronsa/mintpy:latest
```

Run an interactive shell session in the container with: 

```
docker run -it andretheronsa/mintpy:latest bash
```

+ `-it` Opens an interactive session and links to your terminal

+ `andretheronsa/mintpy:latest` Specifies the container

+ `bash` Specefies which command to run - bash in this case opens a bash session - but any command can be used.

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

+ Container was built on insarlab/master - should be updated with new releases.

+ Needs further testing and improvement - can be made smaller (use Alpine instead of Debian...)
