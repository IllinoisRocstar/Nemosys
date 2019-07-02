#ifndef ep16Prep_H
#define ep16Prep_H

#ifdef HAVE_EPIC

// Nemosys headers
#include "nemosys_export.h"
#include "inputGen.H"
#include "exoMesh.H"

// other
#include <iostream>


// A preprocessor class derived from generic
// high-fidelity input generation class inputGen
// This derived class will implement methods to generate
// input files specific to EPIC 2016 code.
class NEMOSYS_EXPORT ep16Prep : public inputGen {
  public:
    explicit ep16Prep(const jsoncons::json &inputjson)
        : _shortForm(false), _mdb(nullptr), inputGen(inputjson) {};

    ~ep16Prep() override;

    static ep16Prep *readJSON(const std::string &jsonFName);
    static ep16Prep *readJSON(const jsoncons::json &inputjson);

  // basic methods
  public:
    void readJSON();
    void process() override;
    void read(const std::string &fname) override;
    void close(const std::string &fname) override;

  // feature specifications
  public:
    void addNdeSet() override {};
    void addElmSet() override {};

  // attribute assignment
  public:
    void addMat() override {};
    void addBC() override {};
    void addIC() override {};
    bool addMisc() override { return false; };

  // internal processing
  public:
    void wrtPre(const std::string &_tsk, const std::string &__tsk) override;
    void wrtCmnt(const std::string &cmnt) override;
    void wrtMsh(const std::string &_tsk, const std::string &__tsk) override;
    void wrtNdeCrds() override {};
    void wrtElmCon() override {};
    void wrtPost() override {};
    void wrtNdeSet() override {};
    void wrtElmSet() override {};
    void wrtMisc(const std::string &_tsk, const std::string &__tsk) override;
    void edit(const std::string &_tsk, const std::string &__tsk) override;

  private:
    std::string _fname;
    bool _shortForm;
    EXOMesh::exoMesh *_mdb;
    std::map<std::string, int> _mat;
    std::map<std::string, std::string> _bcs;
    std::stringstream _buffer;
};

#endif

#endif
