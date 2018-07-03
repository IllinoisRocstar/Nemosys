#ifndef ep16Prep_H
#define ep16Prep_H

#ifdef HAVE_EPIC

// Nemosys
#include "inputGen.H"

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
    ep16Prep(json inputjson):_shortForm(true),inputGen(inputjson) {};
    ~ep16Prep(){};

    static ep16Prep* readJSON(std::string jsonFName);
    static ep16Prep* readJSON(json inputjson);
    
    // basic methods
public:
    void readJSON();
    void process() {};
    void write() {};

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
private:
    void wrtPre() {};
    void wrtCmnt() {};
    void wrtMsh() {};
    void wrtNdeCrds() {};
    void wrtElmCon() {};
    void wrtPost() {};
    void wrtNdeSet() {};
    void wrtElmSet() {};
    void wrtMisc() {};

private:
    bool _shortForm;
};

#endif

#endif
