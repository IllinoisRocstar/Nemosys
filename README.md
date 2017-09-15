NEMoSys
----------
This is the NEMoSys.

### Build Instructions ###

Need to `apt install` at least the following dependencies:

* build-essential
* cmake
* libvtk6-dev
* libproj-dev
* libcgns-dev
* libmetis-dev
* libhdf5-dev
* libfltk1.3-dev
* liblapack-dev
* libgmp-dev
* libjpeg-dev
* libsm-dev
* libice-dev
* gfortran

You can use `build.sh` found in the root project directory to automatically
build the appropriate dependencies given the path to this directory and the
tarball of dependencies, or you can build each dependency manually. You can
find all the dependencies at:

/Projects/IR/Users/msafdari/share/nemosys_tpls.tar.gz

Assuming your building everything manually and $NEMOSYS_PROJECT_PATH is the
path to the local install directory, you start by extracting the whole archive:

```
$ cd $HOME/appropriate/project/path
$ cp /Projects/IR/Users/msafdari/share/nemosys_tpls.tar.gz .
$ tar zxf nemosys_tpls.tar.gz
$ cd nemosys_tpls
```

#### Building Gmsh ####

Unpack Gmsh from the `neomsys_tpls` directory:

```
$ tar zxf gmsh-2.15.0-source.tgz
$ cd gmsh-2.15.0-source
```

Build Gmsh by running the following commands:

```
$ mkdir lib
$ cd lib
$ cmake -DCMAKE_INSTALL_PREFIX=$NEMOSYS_PROJECT_PATH/gmsh -DDEFAULT=0
        -DENABLE_BUILD_LIB=1 -DENABLE_BUILD_SHARED=1 ..
$ make lib shared install/fast -j8
$ cp ./Mesh/meshPartitionObjects.h $NEMOSYS_PROJECT_PATH/gmsh/include/gmsh/
```

#### Building madlib ####

Unpack madlib from the `neomsys_tpls` directory:

```
$ tar zxf madlib-1.3.0.tar.gz
$ cd madlib-1.3.0/
```

Build madlib:

```
$ ./configure --prefix=$NEMOSYS_PROJECT_PATH/madlib --enable-moveIt
              --enable-benchmarks --enable-ann
              --enable-gmsh --with-gmsh-prefix=$NEMOSYS_PROJECT_PATH/gmsh
$ make -j8
$ make install
```

Afterwards, a number of header files from the madlib source directory will
need to be manually copied into `$NEMOSYS_PROJECT_PATH/madlib/include/MAdLib`:

* Mesh/MeshDataBase.h
* Mesh/MeshDataBaseInterface.h
* Mesh/MeshDataBaseIterators.h
* Mesh/MeshDataBaseAttachable.h
* Mesh/MeshDataBaseMiniMesh.h
* Mesh/MshTags.h
* Mesh/MeshDataBaseIO.h
* Mesh/CheckOrientation.h
* Common/MAdMessage.h
* Common/MAdSingleton.h
* Geo/GmshEntities.h
* Geo/Physical.h
* Adapt/utils/NodalDataManager.h

#### Building CGNS ####

Unpack Gmsh from the `neomsys_tpls` directory:

```
$ tar zxf cgns.tar.gz
$ cd CGNS
```

Build CGNS:

```
$ mkdir build && cd build
$ cmake -DCMAKE_INSTALL_PREFIX=$NEMOSYS_PROJECT_PATH/cgns ..
$ make -j8
$ make install
```

#### Building Nemosys ####

With the local dependencies in place, you should be able to build Nemosys
using the following commands:

```
$ mkdir build
$ cd build
$ CMAKE_PREFIX_PATH=../install/madlib:../install/gmsh:../install/cgns cmake ..
$ make -j8
```

There may be some cmake and build warnings. You can get rid of the cmake
warnings by building all of the dependencies from /Projects, installing them
locally, and configuring with cmake, but this is unnecessary and some projects
will take an extremely long time to build.
