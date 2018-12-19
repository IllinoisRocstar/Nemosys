#ifndef ep16Prep_H
#define ep16Prep_H

#ifdef HAVE_EPIC

// Nemosys
#include "inputGen.H"
#include "exoMesh.H"

// other
#include <iostream>
#include "jsoncons/json.hpp"

using jsoncons::json;


// A preprocessor class derived from generic
// high-fidelity input generation class inputGen
// This derived class will implement methods to generate
// input files specific to EPIC 2016 code.
class ep16Prep : public inputGen
{
public:
    ep16Prep(json inputjson):_shortForm(true),_mdb(NULL),inputGen(inputjson) {};
    ~ep16Prep();

    static ep16Prep* readJSON(std::string jsonFName);
    static ep16Prep* readJSON(json inputjson);
    
    // basic methods
public:
    void readJSON();
    void process();

    // feature specifications
public:
    void addNdeSet() {};
    void addElmSet() {};

    // attribute assignement
public:
    void addMat() {};
    void addBC() {};
    void addIC() {};
    bool addMisc() {};

    // internal processing
public:
    void wrtPre(std::string _tsk, std::string __tsk);
    void wrtCmnt(std::string cmnt);
    void wrtMsh(std::string _tsk, std::string __tsk);
    void wrtNdeCrds() {};
    void wrtElmCon() {};
    void wrtPost() {};
    void wrtNdeSet() {};
    void wrtElmSet() {};
    void wrtMisc(std::string _tsk, std::string __tsk);

private:
    bool _shortForm;
    EXOMesh::exoMesh* _mdb;
    std::map<std::string,int> _mat;
    std::map<std::string,std::string> _bcs;

};

#endif

#endif
