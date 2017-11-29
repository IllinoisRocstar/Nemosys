#ifndef TRANSFER_H
#define TRANSFER_H

#include <meshBase.H>
#include <ANN.h>
#include <chrono>

class Transfer
{
  // constructors and destructors
  public:
    Transfer():source(NULL),target(NULL)
    { 
      std::cout << "Transfer constructed" << std::endl;
    }
    virtual ~Transfer()
    { 
      std::cout << "Transfer destroyed" << std::endl;
    }
 
  // 'virtual' constructor
  public:
    static Transfer* Create(std::string method, meshBase* source, meshBase* target);
 
  // transfer methods
  public:
    // transfer point data with given id from source to target
    virtual int runPD(vtkPointData* pd, int arrayID) = 0;
    // transfer cell data with given id from source to target
    virtual int runCD(vtkCellData* cd, int arrayID) = 0;
    // transfer all data from source to target
    virtual int run() = 0; 

  protected:
    meshBase* source;
    meshBase* target;   
};

class FETransfer : public Transfer
{
  public:
    FETransfer(meshBase* _source, meshBase* _target)
    {
      source = _source;
      target = _target;
      std::cout << "FETransfer constructed" << std::endl;
    }

    ~FETransfer()
    {
      std::cout << "FETransfer destroyed" << std::endl;
    }
    int runPD(vtkPointData* pd, int arrayID);
    int runCD(vtkCellData* cd, int arrayID){}; // TODO: not implemented
    int run();
};

#endif
