#ifndef ROCPREPDRIVER_H
#define ROCPREPDRIVER_H

#include <NemDriver.H>

class RocPrepDriver : public NemDriver
{
  public:
    RocPrepDriver(std::string& fname, int numPartitions);
    // for testing recreate
    RocPrepDriver(const std::vector<std::string>& fnames);
    ~RocPrepDriver(){}
    static RocPrepDriver* readJSON(json inputjson);
};
#endif
