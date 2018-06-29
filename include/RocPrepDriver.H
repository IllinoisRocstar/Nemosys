#ifndef ROCPREPDRIVER_H
#define ROCPREPDRIVER_H

#include <NemDriver.H>

class RocPrepDriver : public NemDriver
{
  public:
    RocPrepDriver(std::string& fname, int numPartitions);
    ~RocPrepDriver();
    static RocPrepDriver* readJSON(json inputjson);

  private:
    meshBase* mesh;
    std::vector<meshBase*> partitions;

};

std::vector<std::string> getCgFNames(const std::string& case_dir, 
                                     const std::string& prefix,
                                     const std::string& base_t);


#endif
