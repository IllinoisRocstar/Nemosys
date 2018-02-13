FROM ubuntu:16.04

RUN apt-get update

RUN apt-get -y install build-essential cmake libvtk6-dev libproj-dev libcgns-dev libmetis-dev libhdf5-dev libfltk1.3-dev liblapack-dev libgmp-dev libjpeg-dev libsm-dev libice-dev gfortran

RUN apt-get -y install python3.5 python3-pip

RUN pip3 install --upgrade pip

COPY . /Nemosys

WORKDIR /Nemosys/

RUN pip3 install -r requirements.txt

#RUN ./build.sh $PWD $PWD/contrib/nemosys_tpls.tar.gz

ENTRYPOINT /bin/bash
