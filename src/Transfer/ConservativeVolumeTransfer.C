#include "Transfer/ConservativeVolumeTransfer.H"

#include "AuxiliaryFunctions.H"
#include "Mesh/vtkMesh.H"

#include <libsupermesh-c.h>

#include <vtkCell.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkDoubleArray.h>
#include <vtkMergePoints.h>
#include <vtkMeshQuality.h>
#include <vtkPointData.h>
#include <vtkTetra.h>

ConservativeVolumeTransfer::ConservativeVolumeTransfer(meshBase *_source,
                                                       meshBase *_target) {
  source = _source;
  target = _target;

  // Set grids, allocate vectors/matrices, quadrature scheme and construct
  // transfer objects
  sourceGrid = vtkUnstructuredGrid::SafeDownCast(source->getDataSet());
  targetGrid = vtkUnstructuredGrid::SafeDownCast(target->getDataSet());

  numSourceNodes = sourceGrid->GetNumberOfPoints();
  numTargetNodes = targetGrid->GetNumberOfPoints();

  massMatrix = Eigen::SparseMatrix<double, Eigen::ColMajor, long>(
      numTargetNodes, numTargetNodes);
  mixedMassMatrix = Eigen::SparseMatrix<double, Eigen::ColMajor, long>(
      numTargetNodes, numSourceNodes);

  // set quadrature scheme
  quadrature = vtkSmartPointer<vtkQuadratureSchemeDefinition>::New();
  quadrature->Initialize(VTK_TETRA, 4, 4, TET4, TET4W);
  // quadrature->Initialize(VTK_TETRA, 4, 8, TET8, TET8W);
}

int ConservativeVolumeTransfer::run(const std::vector<std::string> &newnames) {
  constructSupermesh();
  constructMassMatrix();
  constructMixedMassMatrix();

  // transfer all point data
  vtkSmartPointer<vtkPointData> sourcePointData = sourceGrid->GetPointData();
  for (int arrayIdx = 0; arrayIdx < sourcePointData->GetNumberOfArrays();
       ++arrayIdx) {
    transfer(arrayIdx);
  }

  /*
  vtkSmartPointer<vtkCellData> sourceCellData = sourceGrid->GetCellData();
  for(int arrayIdx = 0; arrayIdx < sourceCellData->GetNumberOfArrays();
  ++arrayIdx)
  {
    // TODO implement transferCellData
  }
  */
  return 0;
}

int ConservativeVolumeTransfer::initialize() {
  constructSupermesh();
  constructMassMatrix();
  constructMixedMassMatrix();
}

int ConservativeVolumeTransfer::transferPointData(
    const std::vector<int> &arrayIDs,
    const std::vector<std::string> &newnames) {
  /*
  if(surfaceTransferEnabled)
  {
    // extract source surface
    vtkSmartPointer<vtkPolyData> sourceSurface = extractSurface(sourceGrid);
    vtkSmartPointer<vtkPolyData> targetSurface = extractSurface(targetGrid);

    vtkMesh* sourceMesh = new vtkMesh(sourceSurface, "sourceSurface.vtp");
    vtkMesh* targetMesh = new vtkMesh(targetSurface, "targetSurface.vtp");

    sourceMesh->write("source_surface_test.vtp");
    targetMesh->write("target_surface_test.vtp");

    ConservativeSurfaceTransfer* surfaceTransfer = new
  ConservativeSurfaceTransfer(sourceMesh, targetMesh);

    surfaceTransfer->transferPointData(arrayIDs);

    targetMesh->write("target_surface_test.vtp");
  }
  */

  for (const int &arrayId : arrayIDs) {
    transfer(arrayId);
  }
  return 0;
}

int ConservativeVolumeTransfer::constructSupermesh() {
  nemAux::Timer supermeshTimer;
  std::cout << "Constructing supermesh..." << std::endl;
  supermeshTimer.start();

  /*
   * A - source
   * B - target
   */

  long nnodes_a;
  int dim_a;
  long nelements_a;
  int loc_a;
  double *positions_a;
  long *enlist_a;

  nemAux::Timer sourceConversionTimer;
  sourceConversionTimer.start();

  getLibSupermeshData(source->getDataSet(), nnodes_a, dim_a, nelements_a, loc_a,
                      positions_a, enlist_a, initSourceTetId);

  sourceConversionTimer.stop();
  std::cout << "source converted to libsupermesh format in "
            << sourceConversionTimer.elapsed() << " milliseconds." << std::endl;

  long nnodes_b;
  int dim_b;
  long nelements_b;
  int loc_b;
  double *positions_b;
  long *enlist_b;

  nemAux::Timer targetConversionTimer;
  targetConversionTimer.start();

  getLibSupermeshData(target->getDataSet(), nnodes_b, dim_b, nelements_b, loc_b,
                      positions_b, enlist_b, initTargetTetId);

  targetConversionTimer.stop();
  std::cout << "target converted to libsupermesh format in "
            << targetConversionTimer.elapsed() << " milliseconds." << std::endl;

  // volumeCheck();

  nemAux::Timer intersectionSearchTimer;
  intersectionSearchTimer.start();

  libsupermesh_tree_intersection_finder_set_input(
      &nnodes_a, &dim_a, &nelements_a, &loc_a, &nnodes_b, &dim_b, &nelements_b,
      &loc_b, positions_a, enlist_a, positions_b, enlist_b);

  intersectionSearchTimer.stop();
  std::cout << "intersection search completed in "
            << intersectionSearchTimer.elapsed() << std::endl;

  nemAux::Timer intersectionCalcTimer;
  intersectionCalcTimer.start();

  long nindices;
  libsupermesh_tree_intersection_finder_query_output(&nindices);

  long nelements = nelements_a;
  long *indices = new long[nindices];
  long *ind_ptr = new long[nelements_a + 1];
  libsupermesh_tree_intersection_finder_get_output(&nelements, &nindices,
                                                   indices, ind_ptr);

  auto getTet = [&](long i, long *enlist, double *positions) -> double * {
    double *tet = new double[3 * 4];
    for (int j = 0; j < 4; ++j) {
      // i*4 + j : index of node j in tet i
      // multiply this by 3 to account for offset in positions array
      long ptIndex = enlist[i * 4 + j] * 3;
      double pt[3] = {positions[ptIndex], positions[ptIndex + 1],
                      positions[ptIndex + 2]};
      int offset = j * 3;
      tet[offset] = pt[0];
      tet[offset + 1] = pt[1];
      tet[offset + 2] = pt[2];
    }
    return tet;
  };

  double *tets_c_buf = new double[1000];

  auto pushData = [&](std::vector<double> tet_c, long parent_a, long parent_b) {
    tets_c.push_back(tet_c);
    parents_a.push_back(parent_a);
    parents_b.push_back(parent_b);
  };

  for (long ai = 0; ai < nelements_a; ++ai) // cells in A
  {
    double *tet_a = getTet(ai, enlist_a, positions_a);

    for (long bptr = ind_ptr[ai]; bptr < ind_ptr[ai + 1];
         ++bptr) // indices of intersection candidates in B
    {
      long bj = indices[bptr]; // get index in b
      double *tet_b = getTet(bj, enlist_b, positions_b);

      int n_tets_c;
      // TODO if n_tets_c > buffer size -> reallocate
      libsupermesh_intersect_tets_real(tet_a, tet_b, tets_c_buf, &n_tets_c);
      if (!n_tets_c)
        continue;

      for (int ck = 0; ck < n_tets_c; ++ck) {
        int offs = ck * 12;

        std::vector<double> tet_c;
        tet_c.insert(tet_c.end(), &tets_c_buf[offs], &tets_c_buf[offs + 12]);

        // push tet and its parent indices in A and B (source and target)
        pushData(tet_c, ai, bj);
      }
    }
  }

  intersectionCalcTimer.stop();
  std::cout << "intersection calculation completed in "
            << intersectionCalcTimer.elapsed() << " milliseconds." << std::endl;

  double totalVol = 0.;
  for (auto tet : tets_c) {
    double vol;
    libsupermesh_tetrahedron_volume(&tet[0], &vol);
    if (vol < 0)
      std::cout << "NEGATIVE VOL" << std::endl;
    totalVol += vol;
  }
  std::cout << std::setprecision(15) << "supermesh volume : " << totalVol
            << std::endl;
  std::cout << "number of elems  : " << tets_c.size() << std::endl;

  // free
  delete[] indices;
  delete[] ind_ptr;
  delete[] tets_c_buf;

  supermeshTimer.stop();
  std::cout << "Supermesh constructed in " << supermeshTimer.elapsed()
            << " milliseconds." << std::endl;
  return 0;
}

int ConservativeVolumeTransfer::constructMassMatrix() {
  nemAux::Timer timer;
  std::cout << "Constructing mass matrix..." << std::endl;
  timer.start();

  const double *weights = quadrature->GetQuadratureWeights();

  typedef Eigen::Triplet<double> Triplet;
  std::vector<Triplet> tripletList;

  // TODO start from initTargetTetId
  for (long cellId = initTargetTetId; cellId < targetGrid->GetNumberOfCells();
       ++cellId) {
    if (targetGrid->GetCellType(cellId) != VTK_TETRA)
      continue;
    vtkTetra *tet = vtkTetra::SafeDownCast(targetGrid->GetCell(cellId));

    double detJ = vtkMeshQuality::TetVolume(tet);

    for (int qi = 0; qi < quadrature->GetNumberOfQuadraturePoints(); ++qi) {
      const double *psi = quadrature->GetShapeFunctionWeights(qi);
      double wi = weights[qi];

      for (int r = 0; r < 4; ++r) {
        for (int s = 0; s < 4; ++s) {
          // TODO check that point ids are zero indexed
          long rGlob = tet->GetPointId(r);
          long sGlob = tet->GetPointId(s);
          tripletList.push_back(
              Triplet(rGlob, sGlob, wi * psi[r] * psi[s] * detJ));
        }
      }
    }
  }

  massMatrix.setFromTriplets(tripletList.begin(), tripletList.end());

  timer.stop();
  std::cout << "Mass matrix constructed in " << timer.elapsed()
            << " milliseconds." << std::endl;

  return 0;
}

int ConservativeVolumeTransfer::constructMixedMassMatrix() {
  nemAux::Timer timer;
  timer.start();
  std::cout << "Constructing mixed mass matrix ..." << std::endl;

  const double *weights = quadrature->GetQuadratureWeights();

  typedef Eigen::Triplet<double> Triplet;
  std::vector<Triplet> tripletList(tets_c.size());

  for (long subCellIdx = 0; subCellIdx < tets_c.size(); ++subCellIdx) {
    std::vector<double> childTet = tets_c[subCellIdx];

    long sourceCellId = parents_a[subCellIdx];
    long targetCellId = parents_b[subCellIdx];

    // TODO check indexing for initSource/initTargetTetIds
    vtkTetra *sourceParentTet = vtkTetra::SafeDownCast(
        sourceGrid->GetCell(initSourceTetId + sourceCellId));
    vtkTetra *targetParentTet = vtkTetra::SafeDownCast(
        targetGrid->GetCell(initTargetTetId + targetCellId));

    double childVertices[4][3];
    for (int i = 0; i < 4; ++i) {
      int offs = i * 3;
      childVertices[i][0] = childTet[offs];
      childVertices[i][1] = childTet[offs + 1];
      childVertices[i][2] = childTet[offs + 2];
    }

    double sourceParentVertices[4][3];
    double targetParentVertices[4][3];
    for (int i = 0; i < 4; ++i) {
      double x[3];

      sourceParentTet->GetPoints()->GetPoint(i, x);
      for (int j = 0; j < 3; ++j)
        sourceParentVertices[i][j] = x[j];

      targetParentTet->GetPoints()->GetPoint(i, x);
      for (int j = 0; j < 3; ++j)
        targetParentVertices[i][j] = x[j];
    }

    // TODO check vol = detJ
    double detJ;
    libsupermesh_tetrahedron_volume(&childTet[0], &detJ);
    // if(detJ < 1e-20) continue;

    for (int qi = 0; qi < quadrature->GetNumberOfQuadraturePoints(); ++qi) {
      const double *psi = quadrature->GetShapeFunctionWeights(qi);
      double wi = weights[qi];

      double x[3] = {0., 0., 0.};

      // get geometric coordinate
      // TODO make sure x is contained in subtet
      for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 3; ++j)
          x[j] += psi[i] * childVertices[i][j];
      }

      // TODO check source/target psis are positive and sum to 1
      // source shape function val
      double srcPsi[4] = {};
      vtkTetra::BarycentricCoords(
          x, sourceParentVertices[0], sourceParentVertices[1],
          sourceParentVertices[2], sourceParentVertices[3], srcPsi);

      // target shape function val
      double tgtPsi[4] = {};
      vtkTetra::BarycentricCoords(
          x, targetParentVertices[0], targetParentVertices[1],
          targetParentVertices[2], targetParentVertices[3], tgtPsi);

      for (int r = 0; r < 4; ++r) {
        for (int s = 0; s < 4; ++s) {
          // TODO check that point ids are zero indexed
          long rGlobal = targetParentTet->GetPointId(r);
          long sGlobal = sourceParentTet->GetPointId(s);
          tripletList.push_back(
              Triplet(rGlobal, sGlobal, wi * tgtPsi[r] * srcPsi[s] * detJ));
        }
      }
    }
  }

  mixedMassMatrix.setFromTriplets(tripletList.begin(), tripletList.end());

  timer.stop();

  std::cout << "Mixed mass matrix constructed in " << timer.elapsed()
            << " milliseconds." << std::endl;

  return 0;
}

int ConservativeVolumeTransfer::constructLoadVector(int arrayId,
                                                    int componentId) {
  nemAux::Timer timer;

  std::cout << "Constructing load vector ..." << std::endl;
  timer.start();

  vtkSmartPointer<vtkDoubleArray> sourceArray = vtkDoubleArray::SafeDownCast(
      sourceGrid->GetPointData()->GetArray(arrayId));

  // allocate and populate source values vector
  Eigen::VectorXd sourceValuesVector = Eigen::VectorXd(numSourceNodes);
  auto sourcePointValues = vtkDoubleArray::SafeDownCast(
      sourceGrid->GetPointData()->GetArray(arrayId));
  for (int i = 0; i < numSourceNodes; ++i) {
    sourceValuesVector(i) = sourcePointValues->GetComponent(i, componentId);
  }

  loadVector = mixedMassMatrix * sourceValuesVector;
  timer.stop();

  std::cout << "Load vector constructed in " << timer.elapsed()
            << " milliseconds." << std::endl;

  return 0;
}

int ConservativeVolumeTransfer::transfer(int arrayId) {
  nemAux::Timer timer;
  timer.start();
  std::cout << "Transferring ..." << std::endl;

  vtkSmartPointer<vtkDataArray> sourceArray =
      sourceGrid->GetPointData()->GetArray(arrayId);

  auto targetPointValues = vtkSmartPointer<vtkDoubleArray>::New();
  targetPointValues->SetName(sourceArray->GetName());
  targetPointValues->SetNumberOfComponents(
      sourceArray->GetNumberOfComponents());

  Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>,
                  Eigen::COLAMDOrdering<int>>
      solver;
  solver.analyzePattern(massMatrix);
  solver.factorize(massMatrix);

  for (int componentId = 0; componentId < sourceArray->GetNumberOfComponents();
       ++componentId) {
    constructLoadVector(arrayId, componentId);

    Eigen::VectorXd targetValuesVector = solver.solve(loadVector);

    for (int targetNodeId = 0; targetNodeId < numTargetNodes; ++targetNodeId) {
      targetPointValues->InsertComponent(targetNodeId, componentId,
                                         targetValuesVector(targetNodeId));
    }
  }

  targetGrid->GetPointData()->AddArray(targetPointValues);

  timer.stop();
  std::cout << "Transfer completed in " << timer.elapsed() << " milliseconds."
            << std::endl;

  return 0;
}

int ConservativeVolumeTransfer::interpolateCellDataToPoints(int cellArrayId) {
  vtkSmartPointer<vtkDataArray> sourceCellValues =
      sourceGrid->GetCellData()->GetArray(cellArrayId);

  auto sourcePointValues = vtkSmartPointer<vtkDoubleArray>::New();
  sourcePointValues->SetName(sourceCellValues->GetName());
  sourcePointValues->SetNumberOfComponents(
      sourceCellValues->GetNumberOfComponents());

  for (int pointId = 0; pointId < sourceGrid->GetNumberOfPoints(); ++pointId) {
    auto incidentCellIds = vtkSmartPointer<vtkIdList>::New();
    sourceGrid->GetPointCells(pointId, incidentCellIds);

    double *tuple = new double[sourceCellValues->GetNumberOfComponents()]();

    double volSum = 0.;

    // get contribution from each incident cell
    for (int i = 0; i < incidentCellIds->GetNumberOfIds(); ++i) {
      int cellId = incidentCellIds->GetId(i);
      vtkTetra *tet = vtkTetra::SafeDownCast(sourceGrid->GetCell(cellId));

      double vol = vtkMeshQuality::TetVolume(tet);
      volSum += vol;
      // ... at each component
      for (int compIdx = 0; compIdx < sourceCellValues->GetNumberOfComponents();
           ++compIdx) {
        tuple[compIdx] += vol * sourceCellValues->GetComponent(cellId, compIdx);
      }
    }
    double volAvg = volSum / incidentCellIds->GetNumberOfIds();

    // divide by volume average
    for (int compIdx = 0; compIdx < sourceCellValues->GetNumberOfComponents();
         ++compIdx) {
      tuple[compIdx] = tuple[compIdx] / volAvg;
    }
  }

  return 0;
}

void ConservativeVolumeTransfer::getLibSupermeshData(
    vtkDataSet *data, long &nnodes, int &dim, long &nelements, int &loc,
    double *&positions, long *&enlist, long &initTetId) {
  initTetId = 0;
  for (int cellId = 0; cellId < data->GetNumberOfCells(); ++cellId) {
    if (data->GetCellType(cellId) == VTK_TETRA) {
      initTetId = cellId;
      break;
    }
  }

  long numNodes = data->GetNumberOfPoints();
  // assume remaining cells are tets
  long numTets = data->GetNumberOfCells() - initTetId;

  nnodes = numNodes;
  dim = 3;
  nelements = numTets;
  loc = 4;

  // coordinates of nodes
  positions = new double[dim * nnodes];
  // element-node graph
  enlist = new long[loc * nelements];

  // get coordinates
  for (int ptIdx = 0; ptIdx < numNodes; ++ptIdx) {
    double pt[3];
    data->GetPoint(ptIdx, pt);

    int offset = ptIdx * 3;
    positions[offset] = pt[0];
    positions[offset + 1] = pt[1];
    positions[offset + 2] = pt[2];
  }

  // get element-node graph
  for (int cellId = initTetId; cellId < data->GetNumberOfCells(); ++cellId) {
    vtkCell *tet = data->GetCell(cellId);

    int offset = (cellId - initTetId) * 4;
    enlist[offset] = tet->GetPointId(0);
    enlist[offset + 1] = tet->GetPointId(1);
    enlist[offset + 2] = tet->GetPointId(2);
    enlist[offset + 3] = tet->GetPointId(3);
  }
}

int ConservativeVolumeTransfer::convertSupermeshToUnstructuredGrid() {
  supermeshGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

  auto mergePoints = vtkSmartPointer<vtkMergePoints>::New();
  auto supermeshGridPoints = vtkSmartPointer<vtkPoints>::New();
  double bounds[6];
  sourceGrid->GetBounds(bounds);

  mergePoints->InitPointInsertion(supermeshGridPoints, bounds, tets_c.size());

  // insert points
  for (auto tet : tets_c) {
    for (int ptIdx = 0; ptIdx < 4; ++ptIdx) {
      int offs = ptIdx * 3;

      double x[3] = {tet[offs], tet[offs + 1], tet[offs + 2]};
      vtkIdType ptId;
      int inserted = mergePoints->InsertUniquePoint(x, ptId);
    }
  }

  supermeshGrid->SetPoints(supermeshGridPoints);

  auto tetArray = vtkSmartPointer<vtkCellArray>::New();
  // add tets
  for (auto tet : tets_c) {
    auto tetra = vtkSmartPointer<vtkTetra>::New();

    Eigen::Vector3d faceVertices[3];
    Eigen::Vector3d topVertex;

    vtkIdType ptIds[4];

    // get points of face and "top" vertex
    for (int ptIdx = 0; ptIdx < 4; ++ptIdx) {
      int offs = ptIdx * 3;

      double x[3] = {tet[offs], tet[offs + 1], tet[offs + 2]};

      if (ptIdx < 3)
        faceVertices[ptIdx] = Eigen::Vector3d(x[0], x[1], x[2]);
      else
        topVertex = Eigen::Vector3d(x[0], x[1], x[2]);

      ptIds[ptIdx] = mergePoints->IsInsertedPoint(x);

      if (ptIds[ptIdx] < 0)
        std::cout << "INVALID POINT" << std::endl;
    }

    Eigen::Vector3d cross = (faceVertices[2] - faceVertices[0])
                                .cross(faceVertices[1] - faceVertices[0]);
    cross = cross / cross.norm();

    Eigen::Vector3d topDirection = topVertex - faceVertices[0];
    topDirection = topDirection / topDirection.norm();

    // if negative, reverse first 3 point ids
    if (cross.dot(topDirection) > 0)
      std::swap(ptIds[0], ptIds[2]);

    tetra->GetPointIds()->SetId(0, ptIds[0]);
    tetra->GetPointIds()->SetId(1, ptIds[1]);
    tetra->GetPointIds()->SetId(2, ptIds[2]);
    tetra->GetPointIds()->SetId(3, ptIds[3]);

    // take first 3 - orient so that corresponding normal points to 4th
    tetArray->InsertNextCell(tetra);
  }

  supermeshGrid->SetCells(VTK_TETRA, tetArray);

  std::cout << "Converted grid has : " << supermeshGrid->GetNumberOfCells()
            << " cells." << std::endl;

  if (supermeshGrid->GetNumberOfCells() != tets_c.size()) {
    std::cout << "Reference grid has " << tets_c.size() << " cells."
              << std::endl;
  }

  return 0;
}

void ConservativeVolumeTransfer::volumeCheck() {
  for (int i = initSourceTetId; i < sourceGrid->GetNumberOfCells(); ++i) {
    vtkTetra *tet = vtkTetra::SafeDownCast(sourceGrid->GetCell(i));
    if (vtkMeshQuality::TetVolume(tet) < 0.) {
      std::cout << "source cell " << i << " has negative volume" << std::endl;
    }
  }
  for (int i = initTargetTetId; i < targetGrid->GetNumberOfCells(); ++i) {
    vtkTetra *tet = vtkTetra::SafeDownCast(targetGrid->GetCell(i));
    if (vtkMeshQuality::TetVolume(tet) < 0.) {
      std::cout << "target cell " << i << " has negative volume" << std::endl;
    }
  }
}

vtkSmartPointer<vtkPolyData>
ConservativeVolumeTransfer::extractSurface(vtkUnstructuredGrid *grid) {
  // extract surface
  auto surfaceFilter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
  surfaceFilter->SetPassThroughPointIds(true);
  surfaceFilter->SetOriginalPointIdsName("original_point_ids");

  surfaceFilter->SetInputData(grid);
  surfaceFilter->Update();

  vtkPolyData *surface = surfaceFilter->GetOutput();

  // polydata appears to be freed once dataset filter is freed, make deep copy
  auto surfaceCopy = vtkSmartPointer<vtkPolyData>::New();
  surfaceCopy->DeepCopy(surface);

  return surfaceCopy;
}
