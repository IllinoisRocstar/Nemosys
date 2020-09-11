NEMoSys
----------
The **N**uclear **E**nergy **Mo**deling **Sys**tem is a modular, extensible
resource designed for use in typical application development systems as well as
distributed web-services environments. The project focus is providing a
framework for robust, automated mesh generation, mesh quality analysis, adaptive
mesh refinement, and data transfer between arbitrary meshes. Python bindings to
the NEMoSys library can also be enabled.

## Version ##
Version 0.51.7

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
| ENABLE_MPI             | Enable MPI support              | OFF     | Requires MPI compiler            |
| ENABLE_TESTING         | Enable testing                  | ON      |                                  |
| ENABLE_BUILD_UTILS     | Build utilities                 | OFF     |                                  |
| ENABLE_PYTHON_BINDINGS | Enable Python bindings          | OFF     | Requires Python and SWIG         |
| ENABLE_CFMSH           | Enable cfMesh Meshing engine    | OFF     | Requires OpenFOAM                |
| ENABLE_CGNS            | Enable CGNS extensions          | OFF     | Requires CGNS                    |
| ENABLE_CONSRV_SURFACE_TRANSFER | Enable conservative surface transfer | OFF | Requires IMPACT         |
| ENABLE_CONSRV_VOLUME_TRANSFER | Enable conservative volume transfer | OFF | Requires MPI              |
| ENABLE_DTK             | Enable DTK extensions           | OFF     | UNSUPPORTED                      |
| ENABLE_EXODUS          | Enable EXODUS II extensions     | OFF     | Requires Exodus II               |
| ENABLE_EPIC            | Enable EPIC preprocessor        | OFF     | Requires ENABLE_EXODUS           |
| ENABLE_HDF5            | Enable HDF5 extensions          | OFF     | Requires HDF5                    |
| ENABLE_METIS           | Enable Metis partitioner        | ON      | Requires METIS                   |
| ENABLE_NETGEN          | Enable Netgen meshing engine    | ON      | Requires Netgen                  |
| ENABLE_OMEGAH          | Enable Omega_h mesh adaptivity  | ON      |                                  |
| ENABLE_OMEGAH_CUDA     | Enable GPU for Omega_h          | OFF     | Requires Kokkos                  |
| ENABLE_OPENCASCADE     | Enable OpenCASACADE support     | ON      | Requires OpenCASCADE (OCCT)      |
| ENABLE_SIMMETRIX       | Enable Simmetrix Meshing engine | OFF     | Requires Simmetrix (UNSUPPORTED) |
| ENABLE_TEMPLATE_MESH   | Enable meshing templates        | ON      |                                  |
| ENABLE_MLAMR           | Enable machine learning based AMR | OFF | Requires Frugally-deep library     |

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

### Enabling machine learning based AMR ###
NEMoSys adaptive mesh refinement module for CFD is now equipped with machine
learning support. This module allows users to use trained machine learning
models for adaptive mesh refinement. To enable this capability, the NEMoSys 
should be compiled with `ENABLE_MLAMR=ON`.

User will also need a header only library 
[frugally-deep](https://github.com/Dobiasd/frugally-deep), which loads python 
trained ML models in C++. User will need to provide installation path for this
library using `-DCMAKE_PREFIX_PATH` flag along with other dependecies of 
NEMoSys. frugally-deep library also requires Tensorflow 2.1.0 installed.
```
-DENABLE_MLAMR=ON
-DCMAKE_PREFIX_PATH="\
${NEMOSYS_DEPS_INSTALL_PATH}/opencascade;\
${NEMOSYS_DEPS_INSTALL_PATH}/gmsh;\
${NEMOSYS_DEPS_INSTALL_PATH}/vtk;\
${NEMOSYS_DEPS_INSTALL_PATH}/netgen; \
/Install/path/to/frugally-deep" \
```

### Enabling CUDA for Omega_h ###
Omega_h can be built with CUDA support using the Kokkos backend, assuming Kokkos
is built with CUDA support. Currently, only Kokkos version 2 is supported.
To enable this, make sure that the following flag is set:
```
-DENABLE_OMEGAH_CUDA=ON
```
and that `CMAKE_PREFIX_PATH` (or `$PATH`) contains
`${NEMOSYS_DEPS_INSTALL_PATH}/kokkos/lib/CMake`. Note that both Kokkos and
NEMoSys must be built as shared libraries (`-DBUILD_SHARED_LIBS=ON`), and that
any code that requires the `oshGeoMesh.H` header must then also be compiled for
CUDA.

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

#### Note
We no longer maintain the build script. The file and the archive of TPLs in `contribs` folder 
are now deprecated. Instead, we publish the latest tested versions of TPLs used in the 
project. The list includes following items:

* Gmsh 4.5.1
* Netgen v6.2-dev (commit hash a2f434ebbf)
* OpenCASCADE 7.3.0
* VTK 7.1.0
* Boost 1.68.0
* OpenFoam version 4, 5, 6, or 7
* Kokkos 2

Some of the TPLs need to be compiled in specific configurations. Directions for
such TPLs are provided in the **Manually Build Third Party Libraries** section.

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
        -DENABLE_PYTHON_BINDINGS=ON \
        -DCMAKE_BUILD_TYPE=Release 
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
Throughout this section, we assume `NEMOSYS_DEPS_INSTALL_PATH` is
pointing to the location of the installation of the TPLs.

#### Building OpenCASCADE ####
Checkout proper version of OpenCASCADE, build and install the project by 
running the following commands:
```
$ mkdir build
$ cd build
$ cmake .. \
        -DCMAKE_INSTALL_PREFIX=${NEMOSYS_DEPS_INSTALL_PATH}/opencascade \
        -DBUILD_DOC_Overview=OFF \
        -DBUILD_MODULE_Draw=OFF \
        -DBUILD_MODULE_Visualization=OFF \
        -DBUILD_MODULE_ApplicationFramework=OFF
        -DBUILD_LIBRARY_TYPE=STATIC
$ make -j$(nproc)
$ make install
```

#### Building Gmsh ####
Gmsh depends on OpenCASCADE. Checkout proper version of Gmsh, 
build and install by running the following commands:
```
$ mkdir build
$ cd build
$ cmake .. \
        -DCMAKE_INSTALL_PREFIX=${NEMOSYS_DEPS_INSTALL_PATH}/gmsh \
        -DCMAKE_PREFIX_PATH=${NEMOSYS_DEPS_INSTALL_PATH}/opencascade \
        -DENABLE_BUILD_LIB=OFF -DENABLE_BUILD_SHARED=ON -DENABLE_PRIVATE_API=ON \
        -DDEFAULT=ON -DENABLE_CGNS=OFF -DENABLE_NETGEN=OFF -DENABLE_HXT=ON \
        -DENABLE_FLTK=ON -DENABLE_OCC_STATIC=ON -DENABLE_BUILD_DYNAMIC=ON \
        -DENABLE_OPENMP=ON
$ make lib shared -j$(nproc)
$ make install -j$(nproc)
```

#### Building VTK ####
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
Build and install Netgen by running the following commands:
```
$ mkdir build && cd build
$ cmake -DCMAKE_INSTALL_PREFIX=${NEMOSYS_DEPS_INSTALL_PATH}/netgen \
        -DUSE_GUI=OFF ..
$ make -j$(nproc)
$ make install
```

#### Building Kokkos ####
Build and install Kokkos version 2 with CUDA backend using the following:
```
$ cd build
$ cmake ${KOKKOS_SRC} \
$       -DCMAKE_CXX_COMPILER=${KOKKOS_SRC}/bin/nvcc_wrapper \
$       -DCMAKE_INSTALL_PREFIX=${NEMOSYS_DEPS_INSTALL_PATH}/kokkos \
$       -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
$       -DCKOKKOS_ARCH=${CUDA_ARCH_CC} \
$       -DKOKKOS_ENABLE_CUDA=ON \
$       -DKOKKOS_ENABLE_CUDA_LAMBDA=ON \
$       -DBUILD_SHARED_LIBS=ON
$ make install
```
where `${KOKKOS_SRC}` is the source directory and `${CUDA_ARCH_CC}` refers to
the architecture and compute capability of your GPU, from the list
`Kepler30 Kepler32 Kepler35 Kepler37 Maxwell50 Maxwell52 Maxwell53 Pascal60
Pascal61 Volta70 Volta72`.

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

## Repository Installation ##

*Note* : This is for internal IR user only.
You can install experimental built of the latest NEMoSys for Ubuntu 18 and CentOS7 using following directions.

### Ubunte 18.04 LTS ###

1- Install following packages/dependencies

```
sudo apt update
sudo apt install apt-transport-https ca-certificates curl software-properties-common wget
```
2- Add repository
```
curl -fsSL http://nemosys-repository.illinois.rocstar/nemosys-repository-pub.gpg | sudo apt-key add -
sudo add-apt-repository "deb http://nemosys-repository.illinois.rocstar/ bionic main"
```
3- Install OpenFoam
```
sudo sh -c "wget -O - https://dl.openfoam.org/gpg.key | apt-key add -"
sudo add-apt-repository http://dl.openfoam.org/ubuntu
```
4- Install Nemosys
```
sudo apt-get install nemosys
```
5- Test installation
```
source /opt/openfoam7/etc/bachrc
nemosysRun
```
you should see:
```
Usage: nemosysRun input.json
```

### CentOS 7 ###

1- Add NEMoSys Repository
```
echo "[NEMoSys]
name=NEMoSys
baseurl=http://nemosys-rpm-repository.illinois.rocstar
enabled=1
gpgcheck=0
" >> /etc/yum.repos.d/nemosys.repo
```

2- Install Nemosys
```
sudo yum install nemosys
```
3- Test installation
```
source /opt/openfoam7/etc/bachrc
nemosysRun
```
you should see:
```
Usage: nemosysRun input.json
```
