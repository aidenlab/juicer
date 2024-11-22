# Dockerized versions of Juicer 2

## Available Container Images

### 1. Custom Kasm Workspace Image (aidenlab/juicer-kasm)
This version provides Juicer2 and Juicebox in a full Ubuntu Desktop environment using the [Kasm framework](https://www.kasmweb.com/docs/latest/how_to/building_images.html).

#### Running via Browser Interface
Access the container through a browser with a fully functional Ubuntu Desktop environment:
```bash
docker run --shm-size=512m -p 6901:6901 -e VNC_PW=password aidenlab/juicer-kasm:latest
```
Then navigate to: `https://<IP>:6901`
- Username: kasm_user
- Password: password

#### Interactive Terminal Access
```bash
docker run --rm -it --shm-size=512m --entrypoint /bin/bash aidenlab/juicer-kasm:latest
```


### 2. Standard Docker Image (aidenlab/juicer:v2.0.1)
This is a lightweight version without the desktop environment, ideal for command-line usage and pipeline integration.

```bash
# Pull the image
docker pull aidenlab/juicer:v2.0.1

# Run Juicer
docker run --rm --gpus 1 -v /path/to/data:/data aidenlab/juicer:v2.0.1 [options]

# Run interactively
docker run --rm --gpus 1 -it --entrypoint /bin/bash -v /path/to/data:/data aidenlab/juicer:v2.0.1
juicer.sh [options]
```

### 3. Singularity/Apptainer Support
Juicer can also be run using Singularity/Apptainer, which is commonly available on HPC systems.

```bash
# Pull the Singularity image from sylabs.io repository
singularity pull --arch amd64 library://aidenlab/juicer/juicer:2.0.1

# or Convert Docker image to Singularity
singularity pull juicer.sif docker://aidenlab/juicer:v2.0.1

# Run Juicer 
./juicer_2.0.1.sif [options]

# or 
singularity run --nv juicer_2.0.1.sif [options]

# Run interactively from the container
singularity shell --nv juicer_2.0.1.sif
juicer.sh [options]

```

## Running Juicer with Demo Data

### Using any container version:
```bash
# Execute the following command in an interactive terminal within the container:
/aidenlab/download-and-run-demo.sh

# or
singularity exec juicer_2.0.1.sif /aidenlab/download-and-run-demo.sh
```

## Resource Requirements
- Minimum recommended memory: 16GB RAM
- Recommended storage: 50GB free space for demo data
- For Kasm version: Additional 512MB shared memory (--shm-size=512m)

## More Documentation
- [Juicer Wiki](https://github.com/aidenlab/juicer/wiki)
- [Juicer GitHub Repository](https://github.com/aidenlab/juicer)
- [Kasm Documentation](https://www.kasmweb.com/docs/latest/how_to/building_images.html)

## Support and Issues
Please report any issues on our [GitHub Issues page](https://github.com/aidenlab/juicer/issues).

