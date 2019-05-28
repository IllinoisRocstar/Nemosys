#!/bin/bash

# This script will build NEMoSys and all of the necessary dependencies
# given a tarball of projects with gmsh, occt, hdf5, cgns, netgen and vtk.
# Usage:
#    ./build.sh PATH_TO_NEMOSYS_TPLS_TARBALL PATH_TO_INSTALL_DIR


set -x

NEMOSYS_DEPS_BUILD_DIR=/tmp/nemosys_build
NEMOSYS_TARBALL_PATH=$1
NEMOSYS_DEPS_INSTALL_PATH=$2

# get number of threads for building
num_threads=$(nproc)

# clean
rm -rf ${NEMOSYS_DEPS_BUILD_DIR}
mkdir ${NEMOSYS_DEPS_BUILD_DIR}
rm -rf ${NEMOSYS_DEPS_INSTALL_PATH}

# extract
cd ${NEMOSYS_DEPS_BUILD_DIR}
cp ${NEMOSYS_TARBALL_PATH} .
tar xzf nemosys_tpls.tar.gz

# build netgen
cd ${NEMOSYS_DEPS_BUILD_DIR}/nemosys_tpls
tar xzf netgen-mesher-git.tar.gz
cd netgen-mesher-git
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${NEMOSYS_DEPS_INSTALL_PATH}/netgen \
      -DUSE_GUI=OFF ..
make -j${num_threads}
make install

# build occt
cd ${NEMOSYS_DEPS_BUILD_DIR}/nemosys_tpls
tar xzf opencascade-7.3.0.tgz
cd opencascade-7.3.0
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${NEMOSYS_DEPS_INSTALL_PATH}/opencascade \
      -DBUILD_MODULE_Draw=OFF \
      -DBUILD_MODULE_Visualization=OFF \
      -DBUILD_MODULE_ApplicationFramework=OFF ..
make -j${num_threads}
make install

# build gmsh
cd ${NEMOSYS_DEPS_BUILD_DIR}/nemosys_tpls
tar xzf gmsh-gmsh_4_2_3.tar.gz
cd gmsh-gmsh_4_2_3
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${NEMOSYS_DEPS_INSTALL_PATH}/gmsh \
      -DCMAKE_PREFIX_PATH=${NEMOSYS_DEPS_INSTALL_PATH}/opencascade \
      -DENABLE_BUILD_LIB=ON -DENABLE_BUILD_SHARED=ON -DENABLE_PRIVATE_API=ON \
      -DDEFAULT=ON -DENABLE_CGNS=OFF -DENABLE_NETGEN=OFF -DENABLE_HXT=OFF ..
make lib shared -j${num_threads}
make install

# build vtk
cd ${NEMOSYS_DEPS_BUILD_DIR}/nemosys_tpls
tar xzf vtk-7.1.0.tar.gz
cd VTK-7.1.0
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=${NEMOSYS_DEPS_INSTALL_PATH}/vtk ..
make -j${num_threads}
make install

# build boost
cd ${NEMOSYS_DEPS_BUILD_DIR}/nemosys_tpls
tar xzf boost_1_68_0.tar.gz
cd boost_1_68_0
./bootstrap.sh --prefix=${NEMOSYS_DEPS_INSTALL_PATH}/boost
./b2 install --prefix=${NEMOSYS_DEPS_INSTALL_PATH}/boost --with=all
