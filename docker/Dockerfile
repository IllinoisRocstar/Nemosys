FROM centos:7

RUN yum update -y
RUN yum upgrade -y

RUN yum install https://dl.fedoraproject.org/pub/epel/epel-release-latest-7.noarch.rpm -y
RUN yum groupinstall -y "Development Tools"

# Compilers and build tools.
RUN yum install -y cmake3 gcc-c++.x86_64 make freetype-devel.x86_64 fltk-devel.x86_64
RUN yum install -y cgnslib-devel hdf5-devel hdf-devel.x86_64 boost169-devel.x86_64
RUN yum install metis-devel.x86_64 lapack-devel exodusii-devel.x86_64 -y
RUN yum install libXmu-devel libXi-devel.x86_64 -y

# need these three for building opencascade - Make sure to build static
# RUN yum install tcl-devel.x86_64 -y
# RUN yum install tk-devel.x86_64 -y
# RUN yum install mesa-libGL-devel.x86_64 mesa-libGL.x86_64 -y
RUN yum install centos-release-scl scl-utils -y
RUN yum install devtoolset-7 -y
# RUN echo "scl enable devtoolset-7 bash" >> /etc/skel/.bashrc
RUN echo "source /opt/rh/devtoolset-7/enable" >> /etc/skel/.bashrc

RUN yum install openmpi-devel -y
RUN echo "source /etc/profile.d/modules.sh" >> /etc/skel/.bashrc
RUN echo "module load mpi/openmpi-x86_64" >> /etc/skel/.bashrc

COPY . /home/TPLs/

RUN yum install /home/TPLs/gmsh-4.5.1-Linux64.rpm -y
RUN yum install /home/TPLs/openfoam7-7.20190902-Linux.rpm -y
RUN yum install /home/TPLs/VTK-7.1.0-Linux.rpm -y
RUN yum install /home/TPLs/netgen-6.1-dev-Linux.rpm -y

RUN echo "source /opt/OpenFOAM-7/etc/bashrc" >> /etc/skel/.bashrc

RUN useradd -U --shell /bin/bash -m -k /etc/skel -d /home/user -u 35553 -c "Build Account" build
USER build

ENV BASH_ENV "/etc/skel/.bashrc"

WORKDIR /home/user/
ENTRYPOINT /bin/bash
