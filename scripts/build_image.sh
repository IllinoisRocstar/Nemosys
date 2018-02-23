#!/bin/bash

IMAGE_NAME=nemosys_build
CONTAINER_NAME=nemosys_test

docker stop $CONTAINER_NAME

docker rm $CONTAINER_NAME

docker build -t $IMAGE_NAME Nemosys/

docker run -itd -v /Projects/IR/Users/misenholt/dev/Nemosys/Nemosys/python:/Nemosys/python \
                -v /Projects/IR/Users/misenholt/dev/Nemosys/Nemosys/testing:/Nemosys/testing \
                --name $CONTAINER_NAME $IMAGE_NAME

docker exec -it $CONTAINER_NAME bash -c "./build.sh /Nemosys /Nemosys/contrib/nemosys_tpls.tar.gz"
docker exec -it $CONTAINER_NAME bash -c "cd build && make test VERBOSE=1"
