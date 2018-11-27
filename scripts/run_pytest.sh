export PYTHONPATH=/Nemosys/python
export LD_LIBRARY_PATH=/Nemosys/build/lib

cd /Nemosys/testing/test_scripts

python3 -m unittest test_pyNemosys.TestPyNemosys

