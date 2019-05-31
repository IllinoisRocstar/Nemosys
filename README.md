NEMoSys
----------
The **N**uclear **E**nergy **Mo**deling **Sys**tem is a modular, extensible
resource designed for use in typical application development systems as well as
distributed web-services environments. The project focus is providing a
framework for robust, automated mesh generation, mesh quality analysis, adaptive
mesh refinement, and data transfer between arbitrary meshes. Python bindings to
the NEMoSys library can also be enabled.

## Version ##
Version 0.29.0

NEMoSys follows semantic versioning. The versions will be major.minor.patch.
We will:
* Increase the patch version for bug fixes, security fixes, and code
documentation. Backwards compatible; no breaking changes.
* Increase the minor version for new features and additions to the library’s
interface. Backwards compatible; no breaking changes.
* Increase the major version for breaking changes to the library’s interface or
breaking changes to behavior.

## Getting Started ##
To acquire NEMoSys, you can download it from Illinois Rocstar's GitHub or clone
it with the following command:
```
$ git clone git@git.illinois.rocstar:Nemosys/Nemosys.git
```

## Build Instructions ##
### Build Dependencies ###
You will need to `apt install` at least the following dependencies:

* build-essential
* cmake 
* gfortran
* libopenmpi-dev (or other MPI compiler)
* zlib1g-dev
* libfreetype6-dev
* libfltk1.3-dev
* libxmu-dev
* libxi-dev
* libhdf5-dev
* liblapack-dev
* libjpeg-dev
* libcgns-dev
* libmetis-dev

Optional dependencies for additional functionality:

* libexodusii-dev
* python3.5-dev
* python3-pip
* python2.7-dev
* python-pip
* swig

Once these dependencies are installed, the easiest way to build the required
third-party libraries is with the `build.sh` script. Assume
`$NEMOSYS_PROJECT_PATH` is the path to NEMoSys, and `$NEMOSYS_INSTALL_PATH` is
the desired installation location. Make sure to use absolute paths and execute
the following:
```
$ NEMOSYS_PROJECT_PATH=/full/path/to/Nemosys/source
$ NEMOSYS_DEPS_INSTALL_PATH=/full/path/to/dependency/install
$ ${NEMOSYS_PROJECT_PATH}/scripts/build.sh \
      ${NEMOSYS_PROJECT_PATH}/contrib/nemosys_tpls.tar.gz \
      ${NEMOSYS_DEPS_INSTALL_PATH}
```

### Build NEMoSys ###
Now, we can compile the NEMoSys library, create its Python bindings, and build
other utilities: 
```
$ NEMOSYS_INSTALL_PATH=/full/path/to/Nemosys/install
$ cd ${NEMOSYS_PROJECT_PATH}
$ mkdir build && cd build
$ cmake .. \
        -DCMAKE_PREFIX_PATH="\
${NEMOSYS_DEPS_INSTALL_PATH}/opencascade;\
${NEMOSYS_DEPS_INSTALL_PATH}/gmsh;\
${NEMOSYS_DEPS_INSTALL_PATH}/vtk;\
${NEMOSYS_DEPS_INSTALL_PATH}/netgen" \
        -DCMAKE_INSTALL_PREFIX=${NEMOSYS_INSTALL_PATH} \
        -DENABLE_BUILD_UTILS=ON \
        -DENABLE_TESTING=ON \
        -DBUILD_SHARED_LIBS=ON \
        -DENABLE_EXODUS=ON \
        -DENABLE_PYTHON_BINDINGS=ON
$ make -j$(nproc) (or however many threads you'd like to use)
$ make install (sudo if install location requires it)
```
Executing the commands above will build all libraries, executables, and
bindings. The libraries are installed in `$NEMOSYS_INSTALL_PATH/Nemosys/lib`.
Executables are installed in `$NEMOSYS_INSTALL_PATH/Nemosys/bin`. If Python
bindings are enabled, the `pyNemosys` module is installed for the user. The
`pyNemosys` module can be imported in Python as `import pyNemosys`. The build
configuration can be modified through the CMake Curses interface, `ccmake`, or
by passing the command line options to `cmake`.

### Enabling cfMesh ###
**cfMesh** is an open-source meshing engine implemented on top of **OpenFOAM**. NEMoSys comes with a fully integrated cfMesh-based meshing module. 
To enable the capability, the NEMOSys should be compiled with `ENABLE_CFMSH=ON` and also pointed to the location of the installation of the cfMesh's meshLibrary headers (`${PATH_TO_CFMSH_INCPATH}`). Following OpenFOAM conventions, headers are usually kept in `meshLibrary\lnInclude` folder.
Note that cfMesh dependens on the OpenFOAM, so before starting the compile process make sure to load OpenFOAM environment variables. Depending on the version, OpenFOAM can be loaded by sourcing the bashrc, or cshrc scripts provided in the `OpenFoam-x.y/etc/`. Refer to the OpenFOAM documentation for furhter instructions. After OpenFOAM environment is loaded, enable cfMesh build by adding these lines to the cmake command:
```
-DENABLE_CFMSH=ON -DCFMSH_INCPATH=${PATH_TO_CFMSH_INCPATH}
```

## Testing NEMoSys ##
From the build directory, execute the following command to test the
installation:
```
$ make test
```
This will execute several tests found in `$NEMOSYS_PROJECT_PATH/testing`.

### Manually Build Third Party Libraries ###
If execution of `build.sh` fails, or you have already installed some of the
dependencies, you can try building the remaining TPLs independently. Extract the
whole archive as such:
```
$ cd ${NEMOSYS_PROJECT_PATH}
$ tar zxf contrib/nemosys_tpls.tar.gz 
$ cd nemosys_tpls
```

#### Building OpenCASCADE ####
Unpack OpenCASCADE from the `nemosys_tpls` directory:
```
$ tar xzf opencascade-7.3.0.tgz
$ cd opencascade-7.3.0
```
Build and install OpenCASCADE:
```
$ mkdir build
$ cd build
$ cmake .. \
        -DCMAKE_INSTALL_PREFIX=${NEMOSYS_DEPS_INSTALL_PATH}/opencascade \
        -DBUILD_MODULE_Draw=OFF \
        -DBUILD_MODULE_Visualization=OFF \
        -DBUILD_MODULE_ApplicationFramework=OFF
$ make -j$(nproc)
$ make install
```

#### Building Gmsh ####
Gmsh depends on OpenCASCADE.
Once installed, unpack Gmsh from the `neomsys_tpls` directory:
```
$ tar xzf gmsh-gmsh_4_2_3.tar.gz
$ cd gmsh-gmsh_4_2_3
```
Build and install Gmsh by running the following commands:
```
$ mkdir build
$ cd build
$ cmake .. \
        -DCMAKE_INSTALL_PREFIX=${NEMOSYS_DEPS_INSTALL_PATH}/gmsh \
        -DCMAKE_PREFIX_PATH=${NEMOSYS_DEPS_INSTALL_PATH}/opencascade \
        -DENABLE_BUILD_LIB=ON -DENABLE_BUILD_SHARED=ON -DENABLE_PRIVATE_API=ON \
        -DDEFAULT=ON -DENABLE_CGNS=OFF -DENABLE_NETGEN=OFF -DENABLE_HXT=OFF
$ make lib shared -j$(nproc)
$ make install -j$(nproc)
```

#### Building VTK ####
Unpack VTK from the `neomsys_tpls` directory:
```
$ tar xzf vtk-7.1.0.tar.gz
$ cd vtk-7.1.0
```
Build and install VTK by running the following commands:
```
$ mkdir build
$ cd build
$ cmake .. \
        -DCMAKE_INSTALL_PREFIX=${NEMOSYS_DEPS_INSTALL_PATH}/vtk
$ make -j$(nproc)
$ make install
```

#### Building Netgen ####
Unpack Netgen from the `nemosys_tpls` directory:
```
$ tar xzf netgen-meshter-git.tar.gz
$ cd netgen-mesher-git
```
Build and install Netgen:
```
$ mkdir build && cd build
$ cmake -DCMAKE_INSTALL_PREFIX=${NEMOSYS_DEPS_INSTALL_PATH}/netgen \
        -DUSE_GUI=OFF ..
$ make -j$(nproc)
$ make install
```

See the building NEMoSys section to proceed from this point and to complete the
build.
