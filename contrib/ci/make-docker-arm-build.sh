#!/bin/bash

docker system prune -f

# build arm image:
docker buildx build --no-cache --platform linux/arm64 --output type=registry \
  -t dealii/dealii:master-jammy-arm \
                --build-arg IMG=jammy \
                --build-arg VER=master \
                --build-arg NJOBS=4 \
                github

# combine images:
docker pull dealii/dealii:master-jammy --platform amd64
docker tag dealii/dealii:master-jammy dealii/dealii:master-jammy-amd64
docker push dealii/dealii:master-jammy-amd64
docker buildx imagetools create -t dealii/dealii:master-jammy \
       dealii/dealii:master-jammy-arm \
       dealii/dealii:master-jammy-amd64

# cleanup:
docker system prune -f
