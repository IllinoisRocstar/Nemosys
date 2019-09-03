NEMoSys
----------
The **N**uclear **E**nergy **Mo**deling **Sys**tem is a modular, extensible
resource designed for use in typical application development systems as well as
distributed web-services environments. The project focus is providing a
framework for robust, automated mesh generation, mesh quality analysis, adaptive
mesh refinement, and data transfer between arbitrary meshes. Python bindings to
the NEMoSys library can also be enabled.

## Version ##
Version 0.42.0

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

## Build Options ##
The following table contains all available CMake options to
configure NEMoSys functionality. The necessary third-party library is listed
in the notes section.

| Option name            | Option description              | Default | Notes                            |
|------------------------|---------------------------------|---------|----------------------------------|
| ENABLE_TESTING         | Enable testing                  | ON      |                                  |
| ENABLE_METIS           | Enable Metis interface          | ON      | Requires METIS                   |
| ENABLE_NETGEN          | Enable Netgen interface         | ON      | Requires Netgen                  |
| ENABLE_BUILD_UTILS     | Build utilities                 | OFF     |                                  |
| ENABLE_MPI             | Enable MPI support              | OFF     | Requires MPI compiler            |
| ENABLE_EXODUS          | Enable EXODUS II extensions     | OFF     | Requires Exodus II               |
| ENABLE_EPIC            | Enable EPIC preprocessor        | OFF     | Requires ENABLE_EXODUS           |
| ENABLE_PYTHON_BINDINGS | Enable Python bindings          | OFF     | Requires Python and SWIG         |
| ENABLE_SIMMETRIX       | Enable Simmetrix Meshing engine | OFF     | Requires Simmetrix (UNSUPPORTED) |
| ENABLE_CFMSH           | Enable cfMesh Meshing engine    | OFF     | Requires OpenFOAM     |
| ENABLE_DTK             | Enable DTK extensions           | OFF     | UNSUPPORTED                      |

### Enabling cfMesh ###
**cfMesh** is an open-source meshing engine implemented on top of **OpenFOAM**.
NEMoSys comes with a fully integrated cfMesh-based meshing module. To enable
the capability, the NEMoSys should be compiled with `ENABLE_CFMSH=ON`.

**Note:** cfMesh depends on
OpenFOAM, so before starting the compile process make sure to load OpenFOAM
environment variables. Depending on the version, OpenFOAM can be loaded by
sourcing the bashrc, or cshrc scripts provided in the `OpenFoam-x.y/etc/`.
Refer to the OpenFOAM documentation for further instructions. After OpenFOAM
environment is loaded, enable cfMesh build by adding this line to the cmake
command:
```
-DENABLE_CFMSH=ON
```

### Enabling Simmetrix (UNSUPPORTED) ###
**Simmetrix** is a commercial meshing engine developed by Simmetrix Inc.
(http://www.simmetrix.com/). To enable NEMoSys interface to the Simmetrix,
compile with `ENABLE_SIMMETRIX=ON` and set the location of the Simmetrix
libraries (`${PATH_TO_SIMMETRIX_LIBS}`). Ensure the folder structure set
by Simmetrix Inc. is maintained as it is used to find the headers
automatically. On the command line, the additional CMake options are:
```
-DENABLE_SIMMETRIX=ON -DSIMMETRIX_LIB_DIR=${PATH_TO_SIMMETRIX_LIBS}
```

## Unix Building Instructions ##
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
bindings. The libraries are installed in `$NEMOSYS_INSTALL_PATH/lib`.
Executables are installed in `$NEMOSYS_INSTALL_PATH/bin`. If Python
bindings are enabled, the `pyNemosys` module is installed for the user. The
`pyNemosys` module can be imported in Python as `import pyNemosys`. The build
configuration can be modified through the CMake Curses interface, `ccmake`, or
by passing the command line options to `cmake`.

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
        -DBUILD_DOC_Overview=OFF \
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

## Windows Build Instructions ##
The dependencies are similar to a UNIX build of NEMoSys with the addition of 
boost. An archive of pre-built Windows dependencies is available with some
exceptions:

Note: When downloading pre-built libraries, ensure to select the version
matching the bitness (32 or 64) and the MSVC compiler version
(14.1, 14.0, etc.) used.

**boost**: https://sourceforge.net/projects/boost/files/boost-binaries/

**netCDF** (only if Exodus/Epic support is enabled): https://www.unidata.ucar.edu/downloads/netcdf/index.jsp

**SWIG** (only if Python bindings are enabled): http://swig.org/download.html

### boost Dependency ###
boost will install by default to `C:\local\boost_#_##_#` where the `#` are
the boost version. The boost location must be specified to CMake with the
`BOOST_ROOT` variable.

### netCDF Dependency ###
The netCDF installation is not detected automatically. Pass the location to the
build system through the `CMAKE_PREFIX_PATH` variable.

### Archived Dependencies ###
Extract the archive to a custom location (`%TPL_DIR%). The dependencies must
be given to CMake through the `CMAKE_PREFIX_PATH` variable.
```
-DCMAKE_PREFIX_PATH="%TPL_DIR%\cgns;%TPL_DIR%\exodusii;%TPL_DIR%\gmsh;%TPL_DIR%\hdf5\cmake\hdf5;%TPL_DIR%\metis;%TPL_DIR%\netgen;%TPL_DIR%\vtk;%TPL_DIR%\zlib"
```
Exclude any you wish to use a custom installation.

### Python on Windows ###
Python bindings are generated using SWIG. A SWIG installer is needed on the
system and must be available to CMake and available to the Python environment.

This is accomplished in two parts: (1) add the `SWIG_EXECUTABLE` variable
pointing to the SWIG executable (`swig.exe`) for CMake and (2) add the folder
containing the SWIG executable to the system `PATH` for Python.
A pre-built Windows executable can be downloaded from the
SWIG website: http://swig.org/download.html

Note: On a system with multiple Python installations, the Cmake variable
`PYTHON_EXECUTABLE` can be set to the `python.exe` executable to force the
specific environment.


### Build NEMoSys ###
With the dependencies specified above installed, we can compile NEMoSys with
Python bindings and build tools with the following command from a MSVC command
prompt:

```
> set NETCDF_LOCATION=C:\full\path\to\netCDF\install
> set BOOST_LOCATION=C:\full\path\to\boost\install
> set SWIG_EXE=C:\full\path\to\SWIG\executable\swig.exe
> set NEMOSYS_INSTALL_PATH=C:\full\path\to\Nemosys\install
> cd %NEMOSYS_PROJECT_PATH%
> md build && cd build
> cmake .. ^
        -G "Ninja" ^
        -DBOOST_ROOT="%BOOST_LOCATION%" ^
        -DSWIG_EXECUTABLE="%SWIG_EXE_LOCATION%" ^
        -DCMAKE_INSTALL_PREFIX="%NEMOSYS_INSTALL_PATH%" ^
        -DBUILD_SHARED_LIBS=ON ^
        -DENABLE_TESTING=ON ^
        -DENABLE_BUILD_UTILS=ON ^
        -DENABLE_PYTHON_BINDINGS=ON ^
        -DCMAKE_PREFIX_PATH="%TPL_DIR%\cgns;%TPL_DIR%\exodusii;%TPL_DIR%\gmsh;%TPL_DIR%\hdf5\cmake\hdf5;%TPL_DIR%\metis;%TPL_DIR%\netgen;%TPL_DIR%\vtk;%TPL_DIR%\zlib;%NETCDF_LOCATION%"
> ninja -j %NUMBER_OF_PROCESSORS% (or however many threads you'd like to use)
> ninja install
```
Executing the commands above will build all libraries, executables, and
bindings. The libraries are installed in `%NEMOSYS_INSTALL_PATH%\lib`.
Executables are installed in `%NEMOSYS_INSTALL_PATH%\bin`.

## Testing NEMoSys ##
From the build directory, execute the following command to test the
installation:
```
UNIX:
$ make test
WINDOWS:
> ninja test
```

This will execute several tests found in `$NEMOSYS_PROJECT_PATH/testing`.
