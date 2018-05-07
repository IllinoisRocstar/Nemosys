rm /Nemosys/python/pyNemosys.py
rm /Nemosys/python/pyNemosys_wrap.cpp
rm -rf /Nemosys/build/python
ls /Nemosys/python
cd /Nemosys/build
cmake ..
make
