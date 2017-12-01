#ifndef FETRANSFER_H
#define FETRANSFER_H

#include <TransferBase.H>

class FETransfer : public TransferBase
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
    int runPD(int arrayID);
    int runCD(vtkCellData* cd, int arrayID){}; // TODO: not implemented
    int run();
};

#endif
