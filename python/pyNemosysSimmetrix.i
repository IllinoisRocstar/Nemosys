%module pyNemosys
%include "pyNemosys.i"
%{
//#include "RocRestartDriver.H"
#include "simmetrixGen.H"
#include "simmetrixParams.H"
#include "cgnsAnalyzer.H"
#include "meshStitcher.H"
#include "meshPartitioner.H"
%}

%include "simmetrixParams.H"

%include "simmetrixGen.H"

%include "meshStitcher.H"

%include "meshPartitioner.H"
