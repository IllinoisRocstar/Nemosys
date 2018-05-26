#!/bin/bash

# This script will build Nemosys and all of the necessary dependencies
# given a tarball of projects with gmsh, madlib, and cgns.
# Usage: 
#	- without simmetrix support
#		./build.sh PATH_TO_NEMOSYS PATH_TO_NEMOSYS_TPLS_TARBALL PATH_TO_INSTALL_DIR
#	- with simmetrix support
#		./ build.sh PATH_TO_NEMOSYS PATH_TO_NEMOSYS_TPLS_TARBALL PATH_TO_INSTALL_DIR PATH_TO_HDF5


set -x

NEMOSYS_DEPS_BUILD_DIR=/tmp/nemosys_build
NEMOSYS_PROJECT_PATH=$1
#NEMOSYS_DEPS_INSTALL_PATH=$NEMOSYS_PROJECT_PATH/install
NEMOSYS_DEPS_INSTALL_PATH=$3
NEMOSYS_TARBALL_PATH=$2

# clean
rm -rf $NEMOSYS_DEPS_BUILD_DIR
mkdir $NEMOSYS_DEPS_BUILD_DIR
rm -rf $NEMOSYS_DEPS_INSTALL_PATH

# extract
cd $NEMOSYS_DEPS_BUILD_DIR
cp $NEMOSYS_TARBALL_PATH .
tar zxf nemosys_tpls.tar.gz

# build netgen
#cd $NEMOSYS_DEPS_BUILD_DIR/nemosys_tpls
#tar xzf netgen-5.3.1.tar.gz
#cd netgen-5.3.1
#./configure --prefix=$NEMOSYS_DEPS_INSTALL_PATH/netgen/opt/netgen --exec-prefix=$NEMOSYS_DEPS_INSTALL_PATH/netgen/opt/netgen --with-tcl=/usr/lib/tcl8.5/ --with-tk=/usr/lib/tk8.5/ --with-togl=/usr/lib --enable-shared --enable-nglib
#make
#make install

# build netgen
cd $NEMOSYS_DEPS_BUILD_DIR/nemosys_tpls
tar xzf netgen-mesher-git.tar.gz
cd netgen-mesher-git
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$NEMOSYS_DEPS_INSTALL_PATH/netgen -DUSE_GUI=OFF ..
make
make install

# build gmsh
cd $NEMOSYS_DEPS_BUILD_DIR/nemosys_tpls
tar zxf gmsh-2.15.0-source.tgz
cd gmsh-2.15.0-source
mkdir lib
cd lib
cmake -DCMAKE_INSTALL_PREFIX=$NEMOSYS_DEPS_INSTALL_PATH/gmsh -DDEFAULT=0 -DENABLE_BUILD_LIB=1 -DENABLE_BUILD_SHARED=1 ..
make lib shared install/fast -j8
cp ../Mesh/meshPartitionObjects.h $NEMOSYS_DEPS_INSTALL_PATH/gmsh/include/gmsh/
cp ../Mesh/Generator.h $NEMOSYS_DEPS_INSTALL_PATH/gmsh/include/gmsh/

# build madlib
cd $NEMOSYS_DEPS_BUILD_DIR/nemosys_tpls
tar zxf madlib-1.3.0.tar.gz
cd madlib-1.3.0/
./configure --prefix=$NEMOSYS_DEPS_INSTALL_PATH/madlib --enable-moveIt --enable-benchmarks --enable-ann --enable-gmsh --with-gmsh-prefix=$NEMOSYS_DEPS_INSTALL_PATH/gmsh
make -j8
make install
ls Mesh/MeshDataBase.h Mesh/MeshDataBaseInterface.h Mesh/MeshDataBaseIterators.h Mesh/MeshDataBaseAttachable.h Mesh/MeshDataBaseMiniMesh.h Mesh/MshTags.h Mesh/MeshDataBaseIO.h Mesh/CheckOrientation.h Common/MAdMessage.h Common/MAdSingleton.h Geo/GmshEntities.h Geo/Physical.h Adapt/utils/NodalDataManager.h | xargs -I {} cp {} $NEMOSYS_DEPS_INSTALL_PATH/madlib/include/MAdLib

HDF5_DIR=$4

# build cgns
cd $NEMOSYS_DEPS_BUILD_DIR/nemosys_tpls
tar zxf cgns.tar.gz
cd CGNS
mkdir build
cd build
CMAKE_PREFIX_PATH=$HDF5_DIR cmake -DCGNS_ENABLE_HDF5=ON -DCMAKE_INSTALL_PREFIX=$NEMOSYS_DEPS_INSTALL_PATH/cgns ..
make -j8
make install
