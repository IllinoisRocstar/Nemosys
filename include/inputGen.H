#ifndef inputGen_H
#define inputGen_H

// headers
#include <sstream>
#include <map>
#include <jsoncons/json.hpp>

using jsoncons::json;

class meshBase;

enum inpFileType {INPGEN_TXT, INPGEN_CGNS, INPGEN_VTK};

// class declerations
// inputGen is an abstract class that provides a general-purpose,
// tempalate structure for generating input needed for the 
// high fidelity physics codes.
class inputGen {

    // constructors and destructors
public:
    inputGen(json jsnStrm):_jstrm(jsnStrm), _key("generic"){};
    ~inputGen() {};

    // basic methods
public:
    void _set_key(std::string key) {_key = key;}
    virtual void setNameType(std::string fname, inpFileType ftyp, std::string key=std::string());
    virtual void setOrder(std::vector<std::string>& ord, std::string key=std::string());
    virtual std::vector<std::string> getOrder(std::string key=std::string());
    virtual void pushOrder(std::string ord, std::string key=std::string());
    virtual void setMsh(meshBase* mb, std::string key=std::string());
    virtual void process() = 0;
    virtual void write() = 0;

    // feature specifications
public:
    virtual void addNdeSet() = 0;
    virtual void addElmSet() = 0;

    // attribute assignement
public:
    virtual void addMat() = 0;
    virtual void addBC() = 0;
    virtual void addIC() = 0;
    virtual bool addMisc() = 0;

    // internal processing
protected:
    virtual void wrtPre() = 0;
    virtual void wrtCmnt() = 0;
    virtual void wrtMsh() = 0;
    virtual void wrtNdeCrds() = 0;
    virtual void wrtElmCon() = 0;
    virtual void wrtPost() = 0;
    virtual void wrtNdeSet() = 0;
    virtual void wrtElmSet() = 0;
    virtual void wrtMisc() = 0;

protected:
    json _jstrm;
    std::stringstream _inp;
    std::string _key;
    std::map<std::string, std::vector<meshBase*> > _mb;
    std::map<std::string, std::vector<std::string> > _ord;
    std::map<std::string, std::vector<std::string> > _fn;
    std::map<std::string, std::vector<inpFileType> > _tpe;
};

#endif

