rm -rf build
mkdir build && cd build
CMAKE_PREFIX_PATH=../install/madlib:../install/gmsh:../install/cgns cmake -DRUN_SWIG=ON -DPYTHON_BINDINGS=ON ..
make -j8 -Wall
