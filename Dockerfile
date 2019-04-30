FROM ubuntu:16.04

RUN apt update

RUN apt -y install build-essential cmake gfortran
RUN apt -y install zlib1g-dev libfreetype6-dev "libfltk1.3-dev" libxmu-dev libxi-dev
RUN apt -y install libhdf5-dev liblapack-dev libjpeg-dev

RUN apt -y install libcgns-dev libmetis-dev libexodusii-dev
RUN apt -y install "python3.5-dev" python3-pip "python2.7-dev" python-pip swig
#RUN apt -y install libvtk6-dev libproj-dev libgmp-dev libsm-dev libice-dev
#RUN apt -y install vim

#RUN pip3 install --upgrade pip

COPY ./contrib/nemosys_tpls.tar.gz /Nemosys/
COPY ./scripts/build.sh /Nemosys/

WORKDIR /Nemosys/

#RUN pip3 install -r requirements.txt

#ENV PYTHONPATH=/Nemosys/python LD_LIBRARY_PATH=/Nemosys/build/lib

RUN ./build.sh $PWD/nemosys_tpls.tar.gz /Nemosys-Deps

ENTRYPOINT /bin/bash
