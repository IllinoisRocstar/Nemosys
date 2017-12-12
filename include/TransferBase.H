#ifndef TRANSFERBASE_H
#define TRANSFERBASE_H

#include <meshBase.H>

class TransferBase
{
  // constructors and destructors
  public:
    TransferBase():source(NULL),target(NULL)
    { 
      std::cout << "TransferBase constructed" << std::endl;
    }
    virtual ~TransferBase()
    { 
      std::cout << "TransferBase destroyed" << std::endl;
    }
 
  // 'virtual' constructor
  public:
    static TransferBase* Create(std::string method, meshBase* source, meshBase* target);
 
  // transfer methods
  public:
    // transfer point data with given id from source to target
    virtual int runPD(vtkPointData* pd, int arrayID) = 0;
    virtual int runPD(const std::vector<int>& arrayIDs) = 0;
    virtual int runPD(int arrayID) = 0;
    // transfer cell data with given id from source to target
    virtual int runCD(vtkCellData* cd, int arrayID,
                      std::vector<std::vector<double>>& targetCenters) = 0;
    virtual int runCD(const std::vector<int>& arrayIDs) = 0;
    // transfer all data from source to target
    virtual int run() = 0; 

  protected:
    meshBase* source;
    meshBase* target;   
};



#endif
