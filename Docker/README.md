# Dockerized version of Juicer 2

## Custom Kasm framework Workspace image
### Dockerfile.kasm
This is the Dockerized version of Juicer2 and Juicebox using the [Kasm framework](https://www.kasmweb.com/docs/latest/how_to/building_images.html) base docker image. 

### Instructions of running Juicer 2 with demo data standalone:

#### Accessing the container via a browser using fully functional Ubuntu Desktop environment : https://\<IP\>:6901
``` bash
docker run --shm-size=512m -p 6901:6901 -e VNC_PW=password aidenlab/juicer-kasm:latest
```
The container is now accessible via a browser : https://\<IP\>:6901

User : kasm_user <BR>
Password: password
<BR>

#### Opening an iteractive terminal within the container:
``` bash
docker run --rm -it --shm-size=512m --entrypoint bash aidenlab/juicer-kasm:latest
```

#### Downloading demo data and running Juicer:
```bash
# In a terminal window execute:
/aidenlab/scripts/download-and-run-demo.sh
```

## More Documentation
Please see [the wiki](https://github.com/aidenlab/juicer/wiki) for extensive documentation or 
https://github.com/aidenlab/juicer

