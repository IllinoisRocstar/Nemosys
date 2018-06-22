#ifndef ROCPREPDRIVER_H
#define ROCPREPDRIVER_H

#include <NemDriver.H>

class RocPrepDriver : public NemDriver
{
  public:
    RocPrepDriver(std::string& fname);
    ~RocPrepDriver(){}
    static RocPrepDriver* readJSON(json inputjson);
};
#endif
