#ifndef REFINEDRIVER_H
#define REFINEDRIVER_H

#include <NemDriver.H>

class RefineDriver : public NemDriver
{
  public:
    RefineDriver(std::string _mesh, std::string method, std::string arrayName,
                 double dev_mult, bool maxIsmin, double edgescale, 
                 std::string ofname, bool transferData);
    RefineDriver(std::string _mesh, std::string method, 
                 double edgescale, std::string ofname, bool transferData);
  
    RefineDriver(std::string _mesh, std::string method, std::string arrayName, int order,
                 std::string ofname, bool transferData);
    

    static RefineDriver* readJSON(json inputjson);
    static RefineDriver* readJSON(std::string ifname);

    ~RefineDriver();

    meshBase* getMesh() { return mesh; }

  private:
    meshBase* mesh;

};

#endif
