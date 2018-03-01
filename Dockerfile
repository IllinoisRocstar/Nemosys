FROM ubuntu:16.04

RUN apt-get update

RUN apt-get -y install build-essential cmake libvtk6-dev libproj-dev libcgns-dev libmetis-dev libhdf5-dev libfltk1.3-dev liblapack-dev libgmp-dev libjpeg-dev libsm-dev libice-dev gfortran

RUN apt-get -y install python3.5-dev python3-pip

RUN apt-get -y install swig

RUN apt-get -y install vim

RUN pip3 install --upgrade pip

COPY . /Nemosys

WORKDIR /Nemosys/

#RUN pip3 install -r requirements.txt

ENV PYTHONPATH=/Nemosys/python LD_LIBRARY_PATH=/Nemosys/build/lib

#RUN ./build.sh $PWD $PWD/contrib/nemosys_tpls.tar.gz

ENTRYPOINT /bin/bash
