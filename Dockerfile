FROM ubuntu:16.04

RUN apt update

# Compilers and build tools.
RUN apt -y install build-essential cmake gfortran libopenmpi-dev

# System TPLs needed by NEMoSys dependencies.
RUN apt -y install zlib1g-dev libfreetype6-dev "libfltk1.3-dev" libxmu-dev libxi-dev
RUN apt -y install libhdf5-dev liblapack-dev libjpeg-dev

# System TPLs needed by NEMoSys.
RUN apt -y install libcgns-dev libmetis-dev libexodusii-dev
RUN apt -y install "python3.5-dev" python3-pip "python2.7-dev" python-pip swig

#RUN apt -y install libvtk6-dev libproj-dev libgmp-dev libsm-dev libice-dev
#RUN apt -y install vim

# Run setuptools upgrade twice in case of "distribute" switch
RUN pip2 install --upgrade pip
RUN pip2 install --upgrade setuptools
RUN pip2 install --upgrade setuptools
RUN pip3 install --upgrade pip
RUN pip3 install --upgrade setuptools
RUN pip3 install --upgrade setuptools

# Copy the TPLs and the build script.
COPY ./contrib/nemosys_tpls.tar.gz /Nemosys/
COPY ./scripts/build.sh /Nemosys/

WORKDIR /Nemosys/

# Install NEMoSys dependencies
RUN ./build.sh $PWD/nemosys_tpls.tar.gz /Nemosys-Deps

ENTRYPOINT /bin/bash
