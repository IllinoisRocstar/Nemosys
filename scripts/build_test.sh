export PYTHONPATH=/Nemosys/python
export LD_LIBRARY_PATH=/Nemosys/build/lib

cd build
CMAKE_PREFIX_PATH=../install/madlib:../install/gmsh:../install/cgns cmake -DRUN_SWIG=ON -DPYTHON_BINDINGS=ON ..
make -j8 -Wall
cd ../testing
python3 test_scripts/test_pyNemosys.py
