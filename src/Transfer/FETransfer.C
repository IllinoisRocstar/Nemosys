#include "AuxiliaryFunctions.H"
#include "Transfer/FETransfer.H"

#include <vtkCellData.h>
#include <vtkCellTypes.h>
#include <vtkPointData.h>
#include <vtkStaticPointLocator.h>

#ifdef HAVE_OPENMP
#  include <omp.h>
#endif

using nemAux::operator-;  // for vector subtraction.
using nemAux::operator*;  // for vector multiplication.

FETransfer::FETransfer(meshBase *_source, meshBase *_target) {
  source = _source;
  target = _target;

  // since we are deprecating meshBase
  // srcPointLocator = source->buildStaticPointLocator();
  srcPointLocator = vtkSmartPointer<vtkStaticPointLocator>::New();
  srcPointLocator->SetDataSet(source->getDataSet());
  srcPointLocator->BuildLocator();
  // trgPointLocator = target->buildStaticPointLocator();
  trgPointLocator = vtkSmartPointer<vtkStaticPointLocator>::New();
  trgPointLocator->SetDataSet(target->getDataSet());
  trgPointLocator->BuildLocator();

  // srcCellLocator = source->buildStaticCellLocator();
  srcCellLocator = vtkSmartPointer<vtkStaticCellLocator>::New();
  srcCellLocator->SetDataSet(source->getDataSet());
  srcCellLocator->BuildLocator();
  // trgCellLocator = target->buildStaticCellLocator();
  trgCellLocator = vtkSmartPointer<vtkStaticCellLocator>::New();
  trgCellLocator->SetDataSet(target->getDataSet());
  trgCellLocator->BuildLocator();

  auto sourceCellTypes = vtkSmartPointer<vtkCellTypes>::New();
  auto targetCellTypes = vtkSmartPointer<vtkCellTypes>::New();

  std::cout << "FETransfer constructed" << std::endl;
}

/* transfers point data with arrayID from source mesh to target
   The algorithm is as follows:
    1) For each point in the target mesh, find the cell of the source
       mesh in which it resides.
        - using a cell locator
        - if cell locator fails, find the nearest neighbor in the source mesh
          and all cells sharing this neighbor point. Check if the target point
   is in any of these neighboring cells 2) When the cell is identified, evaluate
   the weights for interpolation of the solution to the target point and perform
   the interpolation.
*/
int FETransfer::transferPointData(const std::vector<int> &arrayIDs,
                                  const std::vector<std::string> &newnames) {
  if (arrayIDs.empty()) {
    std::cerr << "No arrays are selected for interpolation" << std::endl;
    throw;
  }

  vtkSmartPointer<vtkPointData> pd = source->getDataSet()->GetPointData();
  // clean target data of duplicate names if no newnames specified
  if (newnames.empty()) {
    int numArr = pd->GetNumberOfArrays();
    for (int arrayID : arrayIDs) {
      if (arrayID >= numArr) {
        std::cerr << "Selected arrayID is out of bounds\n";
        std::cerr << "There are " << numArr << " point data arrays"
                  << std::endl;
        throw;
      }
      target->unsetPointDataArray(pd->GetArrayName(arrayID));
    }
  }

  std::vector<vtkSmartPointer<vtkDoubleArray>> dasSource(arrayIDs.size());
  std::vector<vtkSmartPointer<vtkDoubleArray>> dasTarget(arrayIDs.size());

  // initializing arrays storing interpolated data
  for (int id = 0; id < arrayIDs.size(); ++id) {
    // get desired point data array from source to be transferred to target
    vtkSmartPointer<vtkDoubleArray> daSource =
        vtkDoubleArray::SafeDownCast(pd->GetArray(arrayIDs[id]));
    // get tuple length of given data`
    int numComponent = daSource->GetNumberOfComponents();
    // declare data array to be populated with values at target points
    vtkSmartPointer<vtkDoubleArray> daTarget =
        vtkSmartPointer<vtkDoubleArray>::New();
    // names and sizing
    if (newnames.empty()) {
      daTarget->SetName(pd->GetArrayName(arrayIDs[id]));
    } else {
      daTarget->SetName(newnames[id].c_str());
    }
    daTarget->SetNumberOfComponents(numComponent);
    daTarget->SetNumberOfTuples(target->getNumberOfPoints());

    std::cout << "Number of components : " << numComponent << std::endl;
    std::cout << "Number of points     : " << target->getNumberOfPoints()
              << std::endl
              << std::endl;
    dasSource[id] = daSource;
    dasTarget[id] = daTarget;
  }
  nemAux::Timer timer;
  timer.start();
#ifdef HAVE_OPENMP
  int numTargetPoints = target->getNumberOfPoints();
  std::cout << "Max number of threads available in transfer "
            << omp_get_max_threads() << std::endl;
// transfer point data
#  pragma omp parallel default(none) \
      shared(dasSource, dasTarget, numTargetPoints, cout)
  {
#  pragma omp single
    {
      std::cout << "Number of threads used in transfer "
                << omp_get_num_threads() << std::endl;
      // the following methods are thread safe IF called from a single thread
      // first
      vtkSmartPointer<vtkGenericCell> tempCell =
          vtkSmartPointer<vtkGenericCell>::New();
      source->getDataSet()->GetCell(0, tempCell);
      target->getDataSet()->GetCell(0, tempCell);
      double x[3];
      source->getDataSet()->GetPoint(0, x);
      target->getDataSet()->GetPoint(0, x);
      auto tempCellIds = vtkSmartPointer<vtkIdList>::New();
      source->getDataSet()->GetPointCells(0, tempCellIds);
      target->getDataSet()->GetPointCells(0, tempCellIds);
    }
    // generic cell buffer used by locator
    vtkSmartPointer<vtkGenericCell> containingCell =
        vtkSmartPointer<vtkGenericCell>::New();
#  pragma omp for schedule(static)
    for (int iPnt = 0; iPnt < numTargetPoints; ++iPnt) {
      transferPointData(iPnt, containingCell, dasSource, dasTarget, false);
    }
  }
#else
  vtkSmartPointer<vtkGenericCell> containingCell =
      vtkSmartPointer<vtkGenericCell>::New();
  for (int iPnt = 0; iPnt < target->getNumberOfPoints(); ++iPnt) {
    transferPointData(iPnt, containingCell, dasSource, dasTarget, false);
  }
#endif
  timer.stop();
  std::cout << "TRANSFER TIME : " << timer.elapsed() << std::endl;
  for (int id = 0; id < arrayIDs.size(); ++id) {
    target->getDataSet()->GetPointData()->AddArray(dasTarget[id]);
  }
  if (checkQual) {
    std::vector<vtkSmartPointer<vtkDoubleArray>> newDasSource(arrayIDs.size());
    for (int id = 0; id < arrayIDs.size(); ++id) {
      int numComponent = dasTarget[id]->GetNumberOfComponents();
      // declare data array to be populated with values at target points
      vtkSmartPointer<vtkDoubleArray> newDaSource =
          vtkSmartPointer<vtkDoubleArray>::New();
      newDaSource->SetNumberOfComponents(numComponent);
      newDaSource->SetNumberOfTuples(target->getNumberOfPoints());
      newDasSource[id] = newDaSource;
    }

    vtkSmartPointer<vtkGenericCell> genCell =
        vtkSmartPointer<vtkGenericCell>::New();
    for (int i = 0; i < source->getNumberOfPoints(); ++i) {
      transferPointData(i, genCell, dasTarget, newDasSource, true);
    }

    for (int id = 0; id < arrayIDs.size(); ++id) {
      int numComponent = newDasSource[id]->GetNumberOfComponents();
      double diffsum = 0.0;
      for (int i = 0; i < source->getNumberOfPoints(); ++i) {
        auto *comps_old = new double[numComponent];
        auto *comps_new = new double[numComponent];
        dasSource[id]->GetTuple(i, comps_old);
        newDasSource[id]->GetTuple(i, comps_new);
        for (int j = 0; j < numComponent; ++j) {
          double diff = std::fabs((comps_new[j] - comps_old[j]) / comps_old[j]);
          diffsum += std::isnan(diff) ? 0.0 : diff * diff;
        }
        delete[] comps_old;
        delete[] comps_new;
      }
      double rmse =
          std::sqrt(diffsum / (numComponent * source->getNumberOfPoints()));
      std::cout << "RMS Error in Nodal Transfer: "
                << (!(std::isnan(rmse) || std::isinf(rmse)) ? rmse : 0)
                << std::endl;
    }
  }
  return 0;
}

int FETransfer::transferPointData(
    int pointId, vtkSmartPointer<vtkGenericCell> containingCell,
    std::vector<vtkSmartPointer<vtkDoubleArray>> &dasSource,
    std::vector<vtkSmartPointer<vtkDoubleArray>> &dasTarget, const bool flip) {
  // id of the cell containing source/target mesh point

  double x[3];             // target point
  double closestPoint[3];  // nearest source point
  vtkIdType cellId;
  int subId;
  double minDist2;
  double *weights;

  if (!flip) {
    target->getDataSet()->GetPoint(pointId, x);
    getClosestSourceCell(x, flip, cellId, containingCell, closestPoint, subId,
                         minDist2, weights);
  } else {
    source->getDataSet()->GetPoint(pointId, x);
    getClosestSourceCell(x, flip, cellId, containingCell, closestPoint, subId,
                         minDist2, weights);
  }
  for (int i = 0; i < dasSource.size(); ++i) {
    int numComponent = dasSource[i]->GetNumberOfComponents();
    auto *comps = new double[numComponent];
    std::vector<double> interps(numComponent, 0.0);
    for (int m = 0; m < containingCell->GetNumberOfPoints(); ++m) {
      int pntId = containingCell->GetPointId(m);
      dasSource[i]->GetTuple(pntId, comps);
      for (int h = 0; h < numComponent; ++h) {
        interps[h] += comps[h] * weights[m];
      }
    }
    delete[] comps;
    // adding interpolated value to data of cell
    dasTarget[i]->SetTuple(pointId, interps.data());
  }
  return 0;
}

/* Transfer cell data from source mesh to target
   The algorithm is as follows:
    1)  Convert the cell data on the source mesh by inverse-distance
        weighted averaging of data at cells sharing given point
          - cell data is assumed to be prescribed at cell centers
    2)  Compute the centers of cell in the target mesh
    3)  Transfer the converted cell-point data from the source mesh
        to the cell centers of the target mesh using the runPD methods
*/
int FETransfer::transferCellData(const std::vector<int> &arrayIDs,
                                 const std::vector<std::string> &newnames) {
  if (arrayIDs.empty()) {
    std::cerr << "no arrays selected for interpolation" << std::endl;
    exit(1);
  }

  vtkSmartPointer<vtkCellData> cd = source->getDataSet()->GetCellData();
  // clean target data of duplicate names if no newnames specified
  if (newnames.empty()) {
    int numArr = cd->GetNumberOfArrays();
    for (int arrayID : arrayIDs) {
      if (arrayID >= numArr) {
        std::cerr << "ERROR: arrayID is out of bounds\n";
        std::cerr << "There are " << numArr << " cell data arrays" << std::endl;
        exit(1);
      }
      target->unsetCellDataArray(cd->GetArrayName(arrayID));
    }
  }
  std::vector<vtkSmartPointer<vtkDataArray>> dasSource(arrayIDs.size());
  std::vector<vtkSmartPointer<vtkDoubleArray>> dasTarget(arrayIDs.size());
  for (int id = 0; id < arrayIDs.size(); ++id) {
    // get desired cell data array from source to be transferred to target
    vtkSmartPointer<vtkDataArray> daSource = cd->GetArray(arrayIDs[id]);
    // get tuple length of given data
    int numComponent = daSource->GetNumberOfComponents();
    // declare data array to be populated with values at target points
    vtkSmartPointer<vtkDoubleArray> daTarget =
        vtkSmartPointer<vtkDoubleArray>::New();
    // names and sizing
    if (newnames.empty())
      daTarget->SetName(cd->GetArrayName(arrayIDs[id]));
    else
      daTarget->SetName(newnames[id].c_str());
    daTarget->SetNumberOfComponents(numComponent);
    daTarget->SetNumberOfTuples(target->getNumberOfCells());
    dasSource[id] = daSource;
    dasTarget[id] = daTarget;
  }
  // straight forward transfer without weighted averaging by locating target
  // cell in source mesh and assigning cell data
  if (!continuous) {
    std::cout << "Non-continuous cell data transfer invoked" << std::endl;
    vtkSmartPointer<vtkGenericCell> genCell =
        vtkSmartPointer<vtkGenericCell>::New();
    for (int i = 0; i < target->getNumberOfCells(); ++i) {
      std::vector<double> targetCenter = target->getCellCenter(i);
      // id of the cell containing source mesh point
      vtkIdType id;
      int subId;
      double minDist2;
      // find closest point and closest cell to x
      double closestPoint[3];
      double *x = targetCenter.data();
      double *weights;
      getClosestSourceCell(x, false, id, genCell, closestPoint, subId, minDist2,
                           weights);
      /*
      srcCellLocator->FindClosestPoint(x, closestPoint, genCell, id, subId,
                                       minDist2);
                                       */
      if (id >= 0) {
        if (minDist2 > c2cTrnsDistTol)
          std::cout << "Warning: For cell at " << x[0] << " " << x[1] << " "
                    << x[2] << " closest cell point found is at "
                    << closestPoint[0] << " " << closestPoint[1] << " "
                    << closestPoint[2] << " with distance " << minDist2
                    << ", Cell IDs: source " << id << " target " << i
                    << std::endl;
        for (int j = 0; j < dasSource.size(); ++j) {
          int numComponent = dasSource[j]->GetNumberOfComponents();
          // double comps[numComponent];
          std::vector<double> comps;
          comps.resize(numComponent, 0.);
          dasSource[j]->GetTuple(id, &comps[0]);
          dasTarget[j]->SetTuple(i, &comps[0]);
        }
      } else {
        std::cerr << "Could not locate target cell " << i
                  << " from in the source mesh! Check the source mesh."
                  << std::endl;
        exit(1);
      }
    }
  } else  // transfer with weighted averaging
  {
    std::cout << "Continuous cell data transfer invoked" << std::endl;
    // ---------------------- Convert source cell data to point data -------- //
    std::vector<vtkSmartPointer<vtkDoubleArray>> dasSourceToPoint(
        arrayIDs.size());

    // cellId container for cells sharing a point
    vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();

    // initializing arrays storing interpolated data
    for (int id = 0; id < arrayIDs.size(); ++id) {
      // get desired cell data array from source to be transferred to target
      vtkSmartPointer<vtkDataArray> daSource = cd->GetArray(arrayIDs[id]);
      // get tuple length of given data
      int numComponent = daSource->GetNumberOfComponents();
      // initialize cellToPoint source array
      vtkSmartPointer<vtkDoubleArray> daSourceToPoint =
          vtkSmartPointer<vtkDoubleArray>::New();
      daSourceToPoint->SetNumberOfComponents(numComponent);
      daSourceToPoint->SetNumberOfTuples(source->getNumberOfPoints());

      for (int i = 0; i < source->getNumberOfPoints(); ++i) {
        // find cells sharing point i
        source->getDataSet()->GetPointCells(i, cellIds);
        int numSharedCells = cellIds->GetNumberOfIds();
        double totW = 0;
        double W;
        // contains averaged/weighted data
        std::vector<double> interps(numComponent, 0.0);
        for (int j = 0; j < numSharedCells; ++j) {
          int cellId = cellIds->GetId(j);
          // double comps[numComponent];
          std::vector<double> comps;
          comps.resize(numComponent, 0.);
          daSource->GetTuple(cellId, &comps[0]);
          // compute distance from point to cell center
          W = 1. / nemAux::l2_Norm(source->getCellCenter(cellId) -
                                   source->getPoint(i));
          // average over shared cells, weighted by inverse distance to center
          for (int k = 0; k < numComponent; ++k) {
            interps[k] += W * comps[k];
          }
          totW += W;
        }
        interps = (1.0 / totW) * interps;
        daSourceToPoint->SetTuple(i, interps.data());
      }
      dasSourceToPoint[id] = daSourceToPoint;
    }
    // genCell used by locator
    vtkSmartPointer<vtkGenericCell> genCell =
        vtkSmartPointer<vtkGenericCell>::New();
    for (int i = 0; i < target->getNumberOfCells(); ++i) {
      transferCellData(i, genCell, dasSourceToPoint, dasTarget);
    }
  }
  for (int id = 0; id < arrayIDs.size(); ++id) {
    target->getDataSet()->GetCellData()->AddArray(dasTarget[id]);
  }
  return 0;
}

int FETransfer::transferCellData(
    int i, vtkSmartPointer<vtkGenericCell> genCell,
    std::vector<vtkSmartPointer<vtkDoubleArray>> &dasSourceToPoint,
    std::vector<vtkSmartPointer<vtkDoubleArray>> &dasTarget) {
  // getting point from target and setting as query
  std::vector<double> targetCenter = target->getCellCenter(i);
  // id of the cell containing source mesh point
  vtkIdType id;
  int subId;
  double minDist2;
  // find closest point and closest cell to x
  double closestPoint[3];
  double *x = targetCenter.data();
  double *weights;
  getClosestSourceCell(x, false, id, genCell, closestPoint, subId, minDist2,
                       weights);
  /*
  srcCellLocator->FindClosestPoint(x, closestPoint, genCell, id, subId,
                                   minDist2);
                                   */
  if (id >= 0) {
    // passed to evaluate position if called
    double pcoords[3];
    // parameters for interpolation
    // int pntId;
    int result;
    auto *weights = new double[genCell->GetNumberOfPoints()];
    result = genCell->EvaluatePosition(x, nullptr, subId, pcoords, minDist2,
                                       weights);
    if (result > 0) {
      for (int id = 0; id < dasSourceToPoint.size(); ++id) {
        int numComponent = dasSourceToPoint[id]->GetNumberOfComponents();
        auto *comps = new double[numComponent];
        std::vector<double> interps(numComponent, 0.0);
        for (int m = 0; m < genCell->GetNumberOfPoints(); ++m) {
          int pntId = genCell->GetPointId(m);
          dasSourceToPoint[id]->GetTuple(pntId, comps);
          for (int h = 0; h < numComponent; ++h) {
            interps[h] += comps[h] * weights[m];
          }
        }
        delete[] comps;
        // adding interpolated value to data of cell
        dasTarget[id]->SetTuple(i, interps.data());
      }
    } else if (result == 0) {
      delete[] weights;
      std::cerr
          << "Could not locate point from target mesh in any cells sharing"
          << " its nearest neighbor in the source mesh" << std::endl;
      exit(1);
    } else {
      delete[] weights;
      std::cerr
          << "problem encountered evaluating position of point from target"
          << " mesh with respect to cell in source mesh" << std::endl;
      exit(1);
    }
  } else {
    std::cerr << "Could not locate center of cell " << i
              << " from target in source mesh" << std::endl;
    exit(1);
  }
  return 0;
}

int FETransfer::run(const std::vector<std::string> &newnames) {
  if (!(source && target)) {
    std::cerr << "source and target meshes must be initialized" << std::endl;
    exit(1);
  }

  // transferring point data
  int numArr = source->getDataSet()->GetPointData()->GetNumberOfArrays();
  if (numArr > 0) {
    std::vector<int> arrayIDs(numArr);
    std::cout << "Transferring point arrays: " << std::endl;
    for (int i = 0; i < numArr; ++i) {
      arrayIDs[i] = i;
      std::cout << "\t" << source->getDataSet()->GetPointData()->GetArrayName(i)
                << "\n";
    }
    transferPointData(arrayIDs, newnames);
  } else {
    std::cout << "no point data found" << std::endl;
  }

  // transferring cell data
  numArr = source->getDataSet()->GetCellData()->GetNumberOfArrays();
  if (numArr > 0) {
    std::vector<int> arrayIDs(numArr);
    std::cout << "Transferring cell arrays: " << std::endl;
    for (int i = 0; i < numArr; ++i) {
      arrayIDs[i] = i;
      std::cout << "\t" << source->getDataSet()->GetCellData()->GetArrayName(i)
                << "\n";
    }
    transferCellData(arrayIDs, newnames);
  } else {
    std::cout << "no cell data found" << std::endl;
  }

  return 0;
}

void FETransfer::getClosestSourceCell(
    double x[3], bool flip, vtkIdType &id,
    vtkSmartPointer<vtkGenericCell> closestCell, double closestPoint[3],
    int &subId, double &minDist2, double *&weights) {
  auto locator = flip ? trgCellLocator : srcCellLocator;
  auto localSource = flip ? target : source;

  // using fast builtin locator findcell method
  double pcoords[3];
  weights = new double[20];
  id = locator->FindCell(x, 1e-10, closestCell, pcoords, weights);
  weights = new double[closestCell->GetNumberOfPoints()];
  closestCell->EvaluatePosition(x, closestPoint, subId, pcoords, minDist2,
                                weights);

  // if failed try searching
  // TODO: find a better way to compute initial bound
  double bound = 1e-3;
  if (id < 0) {
    while (bound < 1e+3) {
      double bbox[6] = {x[0] - bound, x[0] + bound, x[1] - bound,
                        x[1] + bound, x[2] - bound, x[2] + bound};
      int subId;
      double pcoords[3];
      auto cellIds = vtkSmartPointer<vtkIdList>::New();
      locator->FindCellsWithinBounds(bbox, cellIds);
      for (int i = 0; i < cellIds->GetNumberOfIds(); ++i) {
        id = cellIds->GetId(i);
        localSource->getDataSet()->GetCell(id, closestCell);
        weights = new double[closestCell->GetNumberOfPoints()];
        closestCell->EvaluatePosition(x, closestPoint, subId, pcoords, minDist2,
                                      weights);
        if (minDist2 < 1e-9) {
          return;
        }
      }
      bound *= 10;
    }
  }

  if (id < 0) {
    std::cerr << "Could not find containing cell for point : " << x[0] << " "
              << x[1] << " " << x[2] << std::endl;
    closestCell = nullptr;
    minDist2 = vtkMath::Inf();
    throw;
  }
}
