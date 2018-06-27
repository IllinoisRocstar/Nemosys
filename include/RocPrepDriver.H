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

std::vector<std::string> getCgFNames(const std::string& case_dir, 
                                     const std::string& prefix,
                                     const std::string& base_t);


#endif
