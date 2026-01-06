# Dockerfile for Edoxaban Model

## Build image
To build the latest image use:
```bash
docker build -f Dockerfile -t matthiaskoenig/edoxaban:0.7.1 -t matthiaskoenig/edoxaban:latest .
```

## Push images
The image is pushed to dockerhub: [Docker Hub – edoxaban](https://hub.docker.com/repository/docker/matthiaskoenig/edoxaban/general)

```bash
docker login
docker push --all-tags matthiaskoenig/edoxaban
```
