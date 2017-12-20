#ifndef FETRANSFER_H
#define FETRANSFER_H

#include <TransferBase.H>

using namespace nemAux;

// This class is used for data transfer between meshes based on the element transfer method

class FETransfer : public TransferBase
{
  public:
    FETransfer(meshBase* _source, meshBase* _target);

    ~FETransfer()
    {
      std::cout << "FETransfer destroyed" << std::endl;
    }

  // point data transfer
  public:
    /* Transfers point data with arrayID from source mesh to target
       The algorithm is as follows;
        1) For each point in the target mesh, find the cell of the source
           mesh in which it exists.
            - using a cell locator
            - if cell locator fails, find the nearest neighbor in the source mesh
              and all cells sharing this neighbor point. Check if the target point is 
              in any of these neighboring cells
        2) When the cell is identified, evaluate the weights for interpolation of the 
           solution to the target point and perform the interpolation. */
    int runPD(vtkPointData* pd, int arrayID);
    // transfer point data from source to provided target points
    // this method is used when transferring converted cellToPoint data
    // from source to cell centers in target 
    std::vector<std::vector<double>> 
      runPD(vtkSmartPointer<vtkCellLocator> cellLocator,
            std::vector<std::vector<double>>& sourceData, 
            std::vector<std::vector<double>>& targetPnts);
    // transfer point data with given ids from source to target
    int runPD(const std::vector<int>& arrayIDs);

  // cell data transfer
  public:
    /* Transfer cell data from source mesh to target
       The algorithm is as follows:
        1)  Convert the cell data on the source mesh by inverse-distance 
            weighted averaging of data at cells sharing given point
              - cell data is assumed to be perscribed at cell centers
        2)  Compute the centers of cell in the target mesh
        3)  Transfer the converted cell-point data from the source mesh
            to the cell centers of the target mesh using the runPD methods*/
    int runCD(vtkCellData* cd, int arrayID, 
              std::vector<std::vector<double>>& targetCenters);
    // transfer cell data with given ids from source to target
    int runCD(const std::vector<int>& arrayIDs);
    // transfer all cell and point data from source to target
    int run();
};

#endif
