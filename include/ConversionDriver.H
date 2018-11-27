#ifndef CONVERSIONDRIVER_H
#define CONVERSIONDRIVER_H

#include <NemDriver.H>

#ifdef HAVE_EXODUSII
#include "exoMesh.H"
#endif

class ConversionDriver : public NemDriver
{
  public:
    
    ConversionDriver(std::string srcmsh, std::string trgmsh, std::string method, std::string ofname, json inputjson);
    ~ConversionDriver();

    static ConversionDriver* readJSON(json inputjson);
    static ConversionDriver* readJSON(std::string ifname);

#ifdef HAVE_EXODUSII
    void procExo(json ppJson, std::string fname, EXOMesh::exoMesh* em);
#endif


  private:
    meshBase* source;
};


#endif
