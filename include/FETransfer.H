#ifndef FETRANSFER_H
#define FETRANSFER_H

#include <TransferBase.H>

using namespace nemAux;

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

    // transfer point data from source to target
    int runPD(vtkPointData* pd, int arrayID);
    // transfer point data from source to provided target points
    // this method is used when transferring converted cellToPoint data
    // from source to cell centers in target 
    std::vector<std::vector<double>> 
      runPD(std::vector<std::vector<double>>& sourceData, 
            std::vector<std::vector<double>>& targetPnts);
    // transfer point data with given ids from source to target
    int runPD(const std::vector<int>& arrayIDs);
    // transfer point data with given id from source to target
    int runPD(int arrayID);
    // transfer cell data with given id from source to target
    int runCD(vtkCellData* cd, int arrayID, 
              std::vector<std::vector<double>>& targetCenters);
    // transfer cell data with given ids from source to target
    int runCD(const std::vector<int>& arrayIDs);
    // transfer all cell and point data from source to target
    int run();
};

#endif
