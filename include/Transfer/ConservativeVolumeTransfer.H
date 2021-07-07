#ifndef NEMOSYS_CONSERVATIVEVOLUMETRANSFER_H_
#define NEMOSYS_CONSERVATIVEVOLUMETRANSFER_H_

#include "Transfer/ConservativeSurfaceTransfer.H"
#include "Transfer/TransferBase.H"
#include "nemosys_export.h"

#include <vtkQuadratureSchemeDefinition.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

#include <cstring>

class NEMOSYS_EXPORT ConservativeVolumeTransfer : public TransferBase {

public:
  ConservativeVolumeTransfer(meshBase *_source, meshBase *_target);

  ~ConservativeVolumeTransfer() override {
    std::cout << "Conservative Volume Transfer destroyed" << std::endl;
  }

public:
  static ConservativeVolumeTransfer *Create(meshBase *_source,
                                            meshBase *_target) {
    return new ConservativeVolumeTransfer(_source, _target);
  }

  static std::shared_ptr<ConservativeVolumeTransfer>
  CreateShared(meshBase *_source, meshBase *_target) {
    return std::shared_ptr<ConservativeVolumeTransfer>(
        ConservativeVolumeTransfer::Create(_source, _target));
  }

public:
  int transferPointData(
      const std::vector<int> &arrayIDs,
      const std::vector<std::string> &newnames = std::vector<std::string>());

  int transferPointData(int i, vtkSmartPointer<vtkGenericCell> genCell,
                        std::vector<vtkSmartPointer<vtkDoubleArray>> &dasSource,
                        std::vector<vtkSmartPointer<vtkDoubleArray>> &dasTarget,
                        bool flip) {
    return 0;
  }

public:
  int transferCellData(
      const std::vector<int> &arrayIDs,
      const std::vector<std::string> &newnames = std::vector<std::string>()) {
    return 0;
  }

  int transferCellData(
      int i, vtkSmartPointer<vtkGenericCell> genCell,
      std::vector<vtkSmartPointer<vtkDoubleArray>> &dasSourceToPoint,
      std::vector<vtkSmartPointer<vtkDoubleArray>> &dasTarget) {
    return 0;
  }

  int run(const std::vector<std::string> &newnames =
              std::vector<std::string>()) override;

  int initialize();

  int transferSurfacePointData(const std::vector<int> &arrayIDs);

public:
  /**
   * @brief Enable surface transfer. This is necessary for meshes with a shared
   * boundary.
   */
  void enableSurfaceTransfer() { surfaceTransferEnabled = true; }

  /**
   * @brief Disable surface transfer.
   */
  void disableSurfaceTransfer() { surfaceTransferEnabled = false; }

  /**
   * @brief Construct supermesh.
   */
  int constructSupermesh();
  /**
   * @brief Construct mass matrix corresponding to target.
   */
  int constructMassMatrix();
  /**
   * @brief Construct mixed mass matrix corresponding to the source and target.
   * Requires a constructed supermesh.
   */
  int constructMixedMassMatrix();
  /**
   * @brief Construct load vector corresponding to a given array and component.
   */
  int constructLoadVector(int arrayId, int componentId);
  /**
   * @brief Transfer array given its id.
   */
  int transfer(int arrayId);

  /**
   * @brief Convertes supermesh to an unstructured grid.
   * @see getSupermeshGrid
   */
  int convertSupermeshToUnstructuredGrid();

  /**
   * @brief Get the source grid in the transfer.
   */
  vtkSmartPointer<vtkUnstructuredGrid> getSourceGrid() { return sourceGrid; }
  /**
   * @brief Get the target grid in the transfer.
   */
  vtkSmartPointer<vtkUnstructuredGrid> getTargetGrid() { return targetGrid; }

  /**
   * @brief Get the supermesh grid in the transfer. Call
   * convertSupermeshToUnstructuredGrid before this query.
   */
  vtkSmartPointer<vtkUnstructuredGrid> getSupermeshGrid() {
    return supermeshGrid;
  }

private:
  double TET8[32] = {
      0.3281633025163817,  0.3281633025163817,  0.3281633025163817,
      0.01551009245085488, 0.1080472498984286,  0.1080472498984286,
      0.1080472498984286,  0.6758582503047142,  0.3281633025163817,
      0.3281633025163817,  0.01551009245085488, 0.3281633025163817,
      0.1080472498984286,  0.1080472498984286,  0.6758582503047142,
      0.1080472498984286,  0.3281633025163817,  0.01551009245085488,
      0.3281633025163817,  0.3281633025163817,  0.1080472498984286,
      0.6758582503047142,  0.1080472498984286,  0.1080472498984286,
      0.01551009245085488, 0.3281633025163817,  0.3281633025163817,
      0.3281633025163817,  0.6758582503047142,  0.1080472498984286,
      0.1080472498984286,  0.1080472498984286};

  double TET8W[8] = {0.1362178425370874, 0.1137821574629126, 0.1362178425370874,
                     0.1137821574629126, 0.1362178425370874, 0.1137821574629126,
                     0.1362178425370874, 0.1137821574629126};

  double TET4[16] = {.5854101966249684, .1381966011250105, .1381966011250105,
                     .1381966011250105, .1381966011250105, .5854101966249684,
                     .1381966011250105, .1381966011250105, .1381966011250105,
                     .1381966011250105, .5854101966249684, .1381966011250105,
                     .1381966011250105, .1381966011250105, .1381966011250105,
                     .5854101966249684};

  double TET4W[4] = {0.250000000000000, 0.250000000000000, 0.250000000000000,
                     0.250000000000000};

  double TET10[40] = {
      .7784952948213300, .0738349017262234, .0738349017262234,
      .0738349017262234, .0738349017262234, .7784952948213300,
      .0738349017262234, .0738349017262234, .0738349017262234,
      .0738349017262234, .7784952948213300, .0738349017262234,
      .0738349017262234, .0738349017262234, .0738349017262234,
      .7784952948213300,

      .4062443438840510, .4062443438840510, .0937556561159491,
      .0937556561159491, .4062443438840510, .0937556561159491,
      .4062443438840510, .0937556561159491, .4062443438840510,
      .0937556561159491, .0937556561159491, .4062443438840510,

      .0937556561159491, .4062443438840510, .4062443438840510,
      .0937556561159491, .0937556561159491, .4062443438840510,
      .0937556561159491, .4062443438840510, .0937556561159491,
      .0937556561159491, .4062443438840510, .4062443438840510,
  };

  double TET10W[10] = {
      .0476331348432089, .0476331348432089, .0476331348432089,
      .0476331348432089,

      .1349112434378610, .1349112434378610, .1349112434378610,
      .1349112434378610, .1349112434378610, .1349112434378610,
  };

  long numSourceNodes;
  long numTargetNodes;

  long initSourceTetId;
  long initTargetTetId;

  vtkSmartPointer<vtkUnstructuredGrid> sourceGrid;
  vtkSmartPointer<vtkUnstructuredGrid> targetGrid;
  vtkSmartPointer<vtkUnstructuredGrid> supermeshGrid;

  vtkSmartPointer<vtkQuadratureSchemeDefinition> quadrature;

  // Eigen requires indexing via a signed type (can't use size_t)
  Eigen::SparseMatrix<double, Eigen::ColMajor, long> massMatrix;
  Eigen::SparseMatrix<double, Eigen::ColMajor, long> mixedMassMatrix;
  Eigen::VectorXd loadVector;

  std::vector<std::vector<double>> children_a;
  std::vector<std::vector<double>> children_b;

  std::vector<std::vector<double>> tets_c;
  std::vector<long> parents_a;
  std::vector<long> parents_b;

  void getLibSupermeshData(vtkDataSet *data, long &nnodes, int &dim,
                           long &nelements, int &loc, double *&positions,
                           long *&enlist, long &initTetId);

  // Given cell array id, interpolate to points in source mesh and return
  // corresponding point array id in source.
  int interpolateCellDataToPoints(int cellArrayId);

  vtkSmartPointer<vtkPolyData> extractSurface(vtkUnstructuredGrid *grid);

  void volumeCheck();

  bool surfaceTransferEnabled = false;
};

#endif  // NEMOSYS_CONSERVATIVEVOLUMETRANSFER_H_
