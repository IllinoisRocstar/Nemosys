#ifndef TRANSFERDRIVER_H
#define TRANSFERDRIVER_H

#include <NemDriver.H>

class TransferDriver : public NemDriver
{
  public:
    TransferDriver(std::string srcmsh, std::string trgmsh, std::string method,
                   std::string ofname, bool checkQuality);

    TransferDriver(std::string srcmsh, std::string trgmsh, std::string method,
                   std::vector<std::string> arrayNames, std::string ofname,
                   bool checkQuality);

    static TransferDriver* readJSON(json inputjson);
    static TransferDriver* readJSON(std::string ifname);

    ~TransferDriver();

  private:
    meshBase* source;
    meshBase* target; 
};


#endif
