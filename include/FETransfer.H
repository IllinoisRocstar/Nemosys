#ifndef FETRANSFER_H
#define FETRANSFER_H

#include "nemosys_export.h"
#include "TransferBase.H"

#include <vtkDoubleArray.h>
#include <vtkGenericCell.h>

class NEMOSYS_EXPORT FETransfer : public TransferBase {
 public:
  FETransfer(meshBase *_source, meshBase *_target);

  ~FETransfer() override { std::cout << "FETransfer destroyed" << std::endl; }

 public:
  static FETransfer *Create(meshBase *source, meshBase *target) {
    return new FETransfer(source, target);
  }

  static std::shared_ptr<FETransfer> CreateShared(meshBase *source,
                                                  meshBase *target) {
    return std::shared_ptr<FETransfer>(FETransfer::Create(source, target));
  }

 public:
  /**
   * @brief Transfers point data with arrayID from source mesh to target
   * The algorithm is as follows;
   * 1) For each point in the target mesh, find the cell of the source
   *      mesh in which it exists.
   *      - using a cell locator
   *      - if cell locator fails, find the nearest neighbor in the source mesh
   *        and all cells sharing this neighbor point. Check if the target point
   * is in any of these neighboring cells 2) When the cell is identified,
   * evaluate the weights for interpolation of the solution to the target point
   * and perform the interpolation.
   * @param arrayIDs array of array ids for transfer
   * @param newnames optional array of names for transferred data:
   * newnames[i] specifies the transferred array name in the target at
   * arrayId[i]
   */
  int transferPointData(const std::vector<int> &arrayIDs,
                        const std::vector<std::string> &newnames =
                            std::vector<std::string>()) override;

  int transferPointData(int pointId,
                        vtkSmartPointer<vtkGenericCell> containingCell,
                        std::vector<vtkSmartPointer<vtkDoubleArray>> &dasSource,
                        std::vector<vtkSmartPointer<vtkDoubleArray>> &dasTarget,
                        bool flip);

  // cell data transfer
 public:
  /**
   * @brief Transfer cell data from source mesh to target
   * The algorithm is as follows:
   *  1)  Convert the cell data on the source mesh by inverse-distance
   *      weighted averaging of data at cells sharing given point
   *        - cell data is assumed to be prescribed at cell centers
   *  2)  Compute the centers of cell in the target mesh
   *  3)  Transfer the converted cell-point data from the source mesh
   *      to the cell centers of the target mesh using the runPD methods
   * @param arrayIDs array of array ids for transfer
   * @param newnames optional array of names for transferred data:
   * newnames[i] specifies the transferred array name in the target at
   * arrayId[i]
   */
  int transferCellData(const std::vector<int> &arrayIDs,
                       const std::vector<std::string> &newnames =
                           std::vector<std::string>()) override;

  int transferCellData(
      int i, vtkSmartPointer<vtkGenericCell> genCell,
      std::vector<vtkSmartPointer<vtkDoubleArray>> &dasSourceToPoint,
      std::vector<vtkSmartPointer<vtkDoubleArray>> &dasTarget);

  /**
   * @brief Transfer all cell and point data from source to target.
   */
  int run(const std::vector<std::string> &newnames =
              std::vector<std::string>()) override;

 private:
  /**
   * @brief Thread safe implementation using static point locator.
   */
  void getClosestSourceCell(double x[3], bool flip, vtkIdType &id,
                            vtkSmartPointer<vtkGenericCell> closestCell,
                            double closestPoint[3], int &subId,
                            double &minDist2, double *&weights);
};

#endif
