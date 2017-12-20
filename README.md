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

Once these dependencies are installed, the easiest way to build the rest of the 
package is with the script 'build.sh'. Assume $NEMOSYS_PROJECT_PATH is the path to
the local install directory. Execute the following:

```
$ cd $NEMOSYS_PROJECT_PATH
$ ./build.sh $PWD $PWD/contrib/nemosys_tpls.tar.gz

```

If this fails, you can try building the tpls independently
Extract the whole archive as such:

```
$ cd $NEMOSYS_PROJECT_PATH/
$ tar zxf contrib/nemosys_tpls.tar.gz 
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
$ cmake -DCMAKE_INSTALL_PREFIX=$NEMOSYS_PROJECT_PATH/install/gmsh -DDEFAULT=0
        -DENABLE_BUILD_LIB=1 -DENABLE_BUILD_SHARED=1 ..
$ make lib shared install/fast -j8
$ cp ./Mesh/meshPartitionObjects.h $NEMOSYS_PROJECT_PATH/install/gmsh/include/gmsh/
```

#### Building madlib ####

Unpack madlib from the `neomsys_tpls` directory:

```
$ tar zxf madlib-1.3.0.tar.gz
$ cd madlib-1.3.0/
```

Build madlib:

```
$ ./configure --prefix=$NEMOSYS_PROJECT_PATH/install/madlib --enable-moveIt
              --enable-benchmarks --enable-ann
              --enable-gmsh --with-gmsh-prefix=$NEMOSYS_PROJECT_PATH/install/gmsh
$ make -j8
$ make install
```

Afterwards, a number of header files from the madlib source directory will
need to be manually copied into `$NEMOSYS_PROJECT_PATH/install/madlib/include/MAdLib`:

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
$ cmake -DCMAKE_INSTALL_PREFIX=$NEMOSYS_PROJECT_PATH/install/cgns ..
$ make -j8
$ make install
```

#### Building Netgen ####
```
Unpack Netgen from the `nemosys_tpls` directory:

```
$ tar xzf netgen-meshter-git.tar.gz
$ cd netgen-mesher-git
```

Build Netgen:

$ mkdir build && cd build
$ cmake -DCMAKE_INSTALL_PREFIX=$NEMOSYS_PROJECT_PATH/install/netgen -DUSE_GUI=OFF ..
$ make
$ make install
```

#### Building Nemosys ####

With the local dependencies in place, you should be able to build Nemosys and all utilities
using the following commands:

```
$ mkdir build
$ cd build
$ CMAKE_PREFIX_PATH=../install/madlib:../install/gmsh:../install/cgns cmake -DBUILD_UTILS=ON ..
$ make -j8
```
