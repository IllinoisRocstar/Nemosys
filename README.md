NEMoSys
----------
The **N**uclear **E**nergy **Mo**deling **Sys**tem is a modular, extensible resource 
designed to be used in typical application development systems as well as in distributed
web-services environments. The focus of the project is on providing a framework for robust,
automated mesh generation, mesh quality analysis, adaptive mesh refinement and data transfer
between arbitrary meshes. Python bindings to the Nemosys library can also be enabled.

## Getting Started ##
To acquire NEMosys, you can download it from Illinois Rocstar's GitHub
or clone it with the following command:
```
$ git clone git@git.illinois.rocstar:Nemosys/Nemosys.git
```
## Build Instructions ##
### Build Dependencies ###
You will need to `apt install` at least the following dependencies:

* build-essential
* cmake
* libproj-dev
* libmetis-dev
* libfltk1.3-dev
* liblapack-dev
* libgmp-dev
* libjpeg-dev
* libsm-dev
* libice-dev
* gfortran
* libxt-dev
* zlib1g-dev
* tcl-dev
* tk-dev
* libxmu-dev
* python-dev
* swig (if you want python bindings)
* an MPI compiler

Once these dependencies are installed, the easiest way to build the required third party
libraries is with the script `build.sh`. Assume $NEMOSYS_PROJECT_PATH is the path to Nemosys, 
and $NEMOSYS_INSTALL_PATH is the desired installation location. Make sure to use absolute paths 
and execute the following:
```
$ NEMOSYS_PROJECT_PATH=/full/path/to/Nemosys
$ NEMOSYS_INSTALL_PATH=/full/path/to/install
$ $NEMOSYS_PROJECT_PATH/scripts/build.sh $NEMOSYS_PROJECT_PATH $NEMOSYS_PROJECT_PATH/contrib/nemosys_tpls.tar.gz $NEMOSYS_INSTALL_PATH
```

### Build Nemosys ###
Now, we can compile the Nemosys library, create its python bindings and other utilities: 
```
$ cd $NEMOSYS_PROJECT_PATH
$ mkdir build && cd build
$ export CC=mpicc
$ export CXX=mpicxx
$ CMAKE_PREFIX_PATH=$NEMOSYS_INSTALL_PATH/madlib:$NEMOSYS_INSTALL_PATH/gmsh:$NEMOSYS_INSTALL_PATH/netgen cmake -DCMAKE_INSTALL_PREFIX=$NEMOSYS_INSTALL_PATH -DENABLE_PYTHON_BINDINGS=ON -DENABLE_BUILD_UTILS=ON -DENABLE_TESTING=ON -DENABLE_MPI=ON -DBUILD_SHARED_LIBS=ON .. 
$ make -j6 (or however many threads you'd like to use)
$ make install (sudo if install location requires it)
$ export LD_LIBRARY_PATH=$NEMOSYS_INSTALL_PATH/Nemosys/lib:$LD_LIBRARY_PATH
$ export PYTHONPATH=$NEMOSYS_INSTALL_PATH/Nemosys/python/lib/python2.7/site-packages:$PYTHONPATH
```
Executing the commands above will build all libraries, executables and bindings. The libraries are
installed in `$NEMOSYS_INSTALL_PATH/Nemosys/lib`. Executables are installed in 
`$NEMOSYS_INSTALL_PATH/Nemosys/bin`. If python bindings are enabled, the `pyNemosys` module files are
installed in `$NEMOSYS_INSTALL_PATH/Nemosys/python/lib/python2.7/site-packages`.
The `pyNemosys` module can be imported in python as `import pyNemosys`. The build configuration 
can modified through the CMake curses interface (ccmake) or by passing the command line options to cmake.

## Testing Nemosys ##
From the build directory, execute the following command to test the installation:
```
$ make test
```
This will execute several tests found in `$NEMOSYS_PROJECT_PATH/testing`.

### Manually Build Third Party Libraries ###
If execution of `build.sh` fails, or you have already installed some of the dependecies,
you can try building the remaining tpls independently
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
$ cmake -DCMAKE_INSTALL_PREFIX=$NEMOSYS_INSTALL_PATH/gmsh -DDEFAULT=0
        -DENABLE_BUILD_LIB=1 -DENABLE_BUILD_SHARED=1 ..
$ make lib shared install/fast -j8
$ cp ./Mesh/meshPartitionObjects.h $NEMOSYS_INSTALL_PATH/gmsh/include/gmsh/
```

#### Building madlib ####
Unpack madlib from the `neomsys_tpls` directory:
```
$ tar zxf madlib-1.3.0.tar.gz
$ cd madlib-1.3.0/
```
Build madlib:
```
$ ./configure --prefix=$NEMOSYS_INSTALL_PATH/madlib --enable-moveIt
              --enable-benchmarks --enable-ann
              --enable-gmsh --with-gmsh-prefix=$NEMOSYS_INSTAL_PATH/gmsh
$ make -j8
$ make install
```
Afterwards, a number of header files from the madlib source directory will
need to be manually copied into `$NEMOSYS_INSTALL_PATH/madlib/include/MAdLib`:

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


#### Building Netgen ####
Unpack Netgen from the `nemosys_tpls` directory:
```
$ tar xzf netgen-meshter-git.tar.gz
$ cd netgen-mesher-git
```
Build Netgen:
```
$ mkdir build && cd build
$ cmake -DCMAKE_INSTALL_PREFIX=$NEMOSYS_PROJECT_PATH/install/netgen -DUSE_GUI=OFF ..
$ make -j8
$ make install
```

See the building Nemosys section to proceed from this point and complete the build.

## License ##
This project is licensed under the University of Illinois/NCSA Open Source License- see LICENSE for
more details.
