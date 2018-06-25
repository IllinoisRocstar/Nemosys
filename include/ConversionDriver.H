#ifndef CONVERSIONDRIVER_H
#define CONVERSIONDRIVER_H

#include <NemDriver.H>

class ConversionDriver : public NemDriver
{
  public:
    
    ConversionDriver(std::string srcmsh, std::string trgmsh, std::string method, std::string ofname, json inputjson);

    static ConversionDriver* readJSON(json inputjson);
    static ConversionDriver* readJSON(std::string ifname);

    ~ConversionDriver();

  private:
    meshBase* source;
};


#endif
