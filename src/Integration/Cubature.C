#include "Integration/Cubature.H"

#include <vtkQuadraturePointsGenerator.h>
#include <vtkMeshQuality.h>
#include <vtkInformationQuadratureSchemeDefinitionVectorKey.h>
#include <vtkCellData.h>
#include <vtkCellTypes.h>
#include <vtkInformation.h>
#include <vtkXMLPolyDataWriter.h>
#include "Mesh/vtkMesh.H" // for writeVTFile

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "AuxiliaryFunctions.H"
#include <iostream>

// Table 10.4 Quadrature for unit tetrahedra in http://me.rice.edu/~akin/Elsevier/Chap_10.pdf
// OR
// Table 6.3.1 and 6.3.1 in http://www.cs.rpi.edu/~flaherje/pdf/fea6.pdf
// also
// https://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.AppI.d/AFEM.AppI.pdf

// The following arrays are shape functions evaluated at quadrature points supporting second
// order integration for the indicated cell type as well as the quadrature point weights.

// TODO: Add a check to preserve cubature information on mesh for future use
// TODO: Add L2-norm function to simplify implementation in other classes where this is used
double TRI3[] =
    {
        .666666666666667, .166666666666667, .166666666666667,
        .166666666666667, .666666666666667, .166666666666667,
        .166666666666667, .166666666666667, .666666666666667
    };
//double TRI3W[] = {0.166666666666666, 0.166666666666666, 0.166666666666666};
double TRI3W[] = {0.333333333333333, 0.333333333333333, 0.333333333333333};

// second order
double TET4[] =
    {
        .585410196624969, .138196601125011, .138196601125011, .138196601125011,
        .138196601125011, .585410196624969, .138196601125011, .138196601125011,
        .138196601125011, .138196601125011, .585410196624969, .138196601125011,
        .138196601125011, .138196601125011, .138196601125011, .585410196624969
    };
//double TET4W[] = {0.041666666666667, 0.041666666666667,
//                  0.041666666666667, 0.041666666666667};
double TET4W[] = {0.25, 0.25, 0.25, 0.25};

double HEX8[] =
    {
        0.490562612162344, 0.131445855765802, 0.0352208109008645,
        0.131445855765802, 0.131445855765802, 0.0352208109008645,
        0.00943738783765592, 0.0352208109008645, 0.131445855765802,
        0.490562612162344, 0.131445855765802, 0.0352208109008645,
        0.0352208109008645, 0.131445855765802, 0.0352208109008645,
        0.00943738783765592, 0.131445855765802, 0.0352208109008645,
        0.131445855765802, 0.490562612162344, 0.0352208109008645,
        0.00943738783765592, 0.0352208109008645, 0.131445855765802,
        0.0352208109008645, 0.131445855765802, 0.490562612162344,
        0.131445855765802, 0.00943738783765592, 0.0352208109008645,
        0.131445855765802, 0.0352208109008645, 0.131445855765802,
        0.0352208109008645, 0.00943738783765592, 0.0352208109008645,
        0.490562612162344, 0.131445855765802, 0.0352208109008645,
        0.131445855765802, 0.0352208109008645, 0.131445855765802,
        0.0352208109008645, 0.00943738783765592, 0.131445855765802,
        0.490562612162344, 0.131445855765802, 0.0352208109008645,
        0.0352208109008645, 0.00943738783765592, 0.0352208109008645,
        0.131445855765802, 0.131445855765802, 0.0352208109008645,
        0.131445855765802, 0.490562612162344, 0.00943738783765592,
        0.0352208109008645, 0.131445855765802, 0.0352208109008645,
        0.0352208109008645, 0.131445855765802, 0.490562612162344,
        0.131445855765802
    };
double HEX8W[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};


GaussCubature::GaussCubature(vtkDataSet *_dataSet)
    : dataSet(_dataSet), numVolCells(0), totalComponents(0)
{
  dataSet->GetCellData()->RemoveArray("QuadratureOffSet");
  constructGaussMesh();
}

GaussCubature::GaussCubature(vtkDataSet *_dataSet,
                             const std::vector<int> &_arrayIDs)
    : dataSet(_dataSet), numVolCells(0), arrayIDs(_arrayIDs),
      totalComponents(0)
{
  dataSet->GetCellData()->RemoveArray("QuadratureOffSet");
  constructGaussMesh();
  interpolateToGaussPoints();
}

GaussCubature *GaussCubature::Create(vtkDataSet *_dataSet)
{
  return new GaussCubature(_dataSet);
}

GaussCubature *GaussCubature::Create(vtkDataSet *_dataSet,
                                     const std::vector<int> &arrayIDs)
{
  return new GaussCubature(_dataSet, arrayIDs);
}

std::unique_ptr<GaussCubature> GaussCubature::CreateUnique(vtkDataSet *_dataSet)
{
  return std::unique_ptr<GaussCubature>(GaussCubature::Create(_dataSet));
}

std::unique_ptr<GaussCubature>
GaussCubature::CreateUnique(vtkDataSet *_dataSet,
                            const std::vector<int> &arrayIDs)
{
  return std::unique_ptr<GaussCubature>(
      GaussCubature::Create(_dataSet, arrayIDs));
}

std::shared_ptr<GaussCubature>
GaussCubature::CreateShared(vtkDataSet *_dataSet)
{
  std::shared_ptr<GaussCubature> cuby;
  cuby.reset(GaussCubature::Create(_dataSet));
  return cuby;
}

std::shared_ptr<GaussCubature>
GaussCubature::CreateShared(vtkDataSet *_dataSet,
                            const std::vector<int> &arrayIDs)
{
  std::shared_ptr<GaussCubature> cuby;
  cuby.reset(GaussCubature::Create(_dataSet, arrayIDs));
  return cuby;
}


void GaussCubature::constructGaussMesh()
{
  // check whether arrayIDs exist in mesh
  {
    vtkSmartPointer<vtkPointData> pd = dataSet->GetPointData();
    int numArr = pd->GetNumberOfArrays();
    for (int arrayID : arrayIDs)
    {
      if (arrayID >= numArr)
      {
        std::cerr << "Array ID " << arrayID << " not found in dataset"
                  << std::endl;
        exit(1);
      }
    }
  }

  // Get the dictionary key     
  auto key = vtkQuadratureSchemeDefinition::DICTIONARY();

  // Get the cell types used by the data set
  auto cellTypes = vtkSmartPointer<vtkCellTypes>::New();
  dataSet->GetCellTypes(cellTypes);
  int nCellTypes = cellTypes->GetNumberOfTypes();

  // create offset array and store the dictionary within
  auto offsets = vtkSmartPointer<vtkIdTypeArray>::New();
  std::string basename = "QuadratureOffset";
  offsets->SetName(basename.c_str());
  auto info = vtkSmartPointer<vtkInformation>::New();
  info = offsets->GetInformation();

  for (int typeId = 0; typeId < nCellTypes; ++typeId)
  {
    int cellType = cellTypes->GetCellType(typeId);
    // Initialize quadrature scheme definition for given cell type
    auto def = vtkSmartPointer<vtkQuadratureSchemeDefinition>::New();
    switch (cellType)
    {
      case VTK_TRIANGLE:
        def->Initialize(VTK_TRIANGLE, 3, 3, TRI3, TRI3W);
        break;
      case VTK_TETRA:
        def->Initialize(VTK_TETRA, 4, 4, TET4, TET4W);
        break;
      case VTK_HEXAHEDRON:
        def->Initialize(VTK_HEXAHEDRON, 8, 8, HEX8, HEX8W);
        break;
      default:
        std::cerr << "Error: Cell type: " << cellType << " found "
                  << "with no quadrature definition provided" << std::endl;
        exit(1);
    }
    // the definition must appear in the dictionary associated with
    // the offset array
    key->Set(info, def, cellType);
  }
  // get dictionary size 
  int dictSize = key->Size(info);
  dict = new vtkQuadratureSchemeDefinition *[dictSize];
  key->GetRange(info, dict, 0, 0, dictSize);
  vtkIdType numCells = dataSet->GetNumberOfCells();

  offsets->SetNumberOfTuples(numCells);

#ifdef HAVE_OPENMP
  std::vector<vtkIdType> chunkOffsets(omp_get_max_threads() + 1, 0);
  #pragma omp parallel default(none) shared(numCells, chunkOffsets, offsets)
  {
    vtkIdType subOffset = 0;
    auto genCell = vtkSmartPointer<vtkGenericCell>::New();

    #pragma omp single
    {
      // GetCell method is thread safe when first called from single thread
      dataSet->GetCell(0, genCell);
    }

    #pragma omp for schedule(static)
    for (vtkIdType cellId = 0; cellId < numCells; ++cellId) {
      offsets->SetValue(cellId, subOffset);

      dataSet->GetCell(cellId, genCell);
      int cellType = genCell->GetCellType();
      if (cellType >= VTK_TETRA) numVolCells += 1;

      vtkQuadratureSchemeDefinition *celldef = dict[cellType];
      subOffset += celldef->GetNumberOfQuadraturePoints();
    }

    int tid = omp_get_thread_num();
    chunkOffsets[tid + 1] = subOffset;
    #pragma omp barrier
    #pragma omp single
    {
      for (int i = 1; i < chunkOffsets.size(); ++i) {
        chunkOffsets[i] = chunkOffsets[i] + chunkOffsets[i - 1];
      }
    }

    #pragma omp for schedule(static)
    for (vtkIdType cellId = 0; cellId < numCells; ++cellId) {
      vtkIdType offset = chunkOffsets[tid] + offsets->GetValue(cellId);
      offsets->SetValue(cellId, offset);
    }
  }
#else
  vtkIdType offset = 0;

  for (int cellid = 0; cellid < numCells; ++cellid) {
    offsets->SetValue(cellid, offset);
    int cellType = dataSet->GetCell(cellid)->GetCellType();
    if (cellType >= VTK_TETRA) {
      numVolCells += 1;
    }
    vtkQuadratureSchemeDefinition *celldef = dict[cellType];
    offset += celldef->GetNumberOfQuadraturePoints();
  }
#endif

  dataSet->GetCellData()->AddArray(offsets);

  auto pointGen = vtkSmartPointer<vtkQuadraturePointsGenerator>::New();

  pointGen->SetInputArrayToProcess
      (0, 0, 0,
       vtkDataObject::FIELD_ASSOCIATION_CELLS,
       "QuadratureOffset");
  pointGen->SetInputData(dataSet);
  gaussMesh = vtkSmartPointer<vtkPolyData>::New();
  gaussMesh = vtkPolyData::SafeDownCast(pointGen->GetOutput());
  pointGen->Update();
}


double GaussCubature::computeCellVolume(vtkSmartPointer<vtkGenericCell> genCell,
                                        int cellType) const
{
  switch (cellType)
  {
    case VTK_TRIANGLE:
    {
      return vtkMeshQuality::TriangleArea(genCell);
    }
    case VTK_QUAD:
    {
      return vtkMeshQuality::QuadArea(genCell);
    }
    case VTK_TETRA:
    {
      return vtkMeshQuality::TetVolume(genCell);
    }
    case VTK_HEXAHEDRON:
    {
      return vtkMeshQuality::HexVolume(genCell);
    }
    default:
    {
      std::cerr << "Error: Cell type: " << cellType << "found "
                << "with no quadrature definition provided" << std::endl;
      exit(1);
    }
  }
}


// FIXME: Can remove switch if you just use the general method implemented in tetra's case
double GaussCubature::computeJacobian(vtkSmartPointer<vtkGenericCell> genCell,
                                      int cellType) const
{
  switch (cellType)
  {
    case VTK_TRIANGLE:
    {
      return vtkMeshQuality::TriangleArea(genCell);
    }
    case VTK_QUAD:
    {
      return vtkMeshQuality::QuadArea(genCell) / 4.0;
    }
    case VTK_TETRA:
    {
      return vtkMeshQuality::TetVolume(genCell);
    }
    case VTK_HEXAHEDRON:
    {
      return vtkMeshQuality::HexVolume(genCell) / 8.0;
    }
    default:
    {
      std::cerr << "Error: Cell type: " << cellType << "found "
                << "with no quadrature definition provided" << std::endl;
      exit(1);
    }
  }
}


int GaussCubature::getOffset(int cellID) const
{
  vtkIdType offsets[1];
  vtkIdTypeArray::FastDownCast(
      dataSet->GetCellData()->GetArray("QuadratureOffset"))
      ->GetTypedTuple(cellID, offsets);
  return offsets[0];
}


pntDataPairVec GaussCubature::getGaussPointsAndDataAtCell(int cellID)
{
  if (arrayIDs.empty())
  {
    std::cerr << "no array have been selected for interpolation" << std::endl;
    exit(1);
  }

  int numDataArr = gaussMesh->GetPointData()->GetNumberOfArrays();

  if (numDataArr == 0)
  {
    interpolateToGaussPoints();
  }

  // get number of gauss points in cell from dictionary
  int numGaussPoints = dict[dataSet->GetCell(
      cellID)->GetCellType()]
      ->GetNumberOfQuadraturePoints();
  // get offset from nodeMesh for lookup of gauss points in polyData
  int offset = getOffset(cellID);

  //pntDataPairVec container;
  pntDataPairVec container(numGaussPoints);

  vtkSmartPointer<vtkPointData> pd = gaussMesh->GetPointData();
  for (int i = 0; i < numGaussPoints; ++i)
  {
    double x_tmp[3];
    gaussMesh->GetPoint(offset + i, x_tmp);
    std::vector<double> gaussPnt(x_tmp, x_tmp + 3);
    std::vector<double> data(totalComponents);
    int currcomp = 0;
    for (int j = 0; j < numComponents.size(); ++j)
    {
      double *comps = new double[numComponents[j]];
      pd->GetArray(j)->GetTuple(offset + i, comps);
      for (int k = 0; k < numComponents[j]; ++k)
      {
        data[currcomp] = comps[k];
        ++currcomp;
      }
      delete[] comps;
    }
    container[i] = std::move(std::make_pair(gaussPnt, data));
  }
  return container;
}


int GaussCubature::interpolateToGaussPointsAtCell
    (const int cellID,
     vtkSmartPointer<vtkGenericCell> genCell,
     const std::vector<vtkSmartPointer<vtkDataArray>> &das,
     std::vector<vtkSmartPointer<vtkDoubleArray>> &daGausses) const
{
  // putting current cell into genCell
  dataSet->GetCell(cellID, genCell);
  // getting cellType information for lookup in map
  int cellType = dataSet->GetCellType(cellID);
  // get quadrature weights for this cell type
  const double *shapeFunctionWeights = dict[cellType]->GetShapeFunctionWeights();
  // number of gauss points in this cell
  int numGaussPoints = dict[cellType]->GetNumberOfQuadraturePoints();
  // get offset from nodeMesh for lookup of gauss points in polyData
  int offset = getOffset(cellID);
  // interpolation loop
  for (int j = 0; j < numGaussPoints; ++j)
  {
    for (int id = 0; id < das.size(); ++id)
    {
      int numComponent = das[id]->GetNumberOfComponents();
      auto *comps = new double[numComponent];
      std::vector<double> interps(numComponent, 0.0);
      for (int m = 0; m < genCell->GetNumberOfPoints(); ++m)
      {
        int pntId = genCell->GetPointId(m);
        das[id]->GetTuple(pntId, comps);
        for (int h = 0; h < numComponent; ++h)
        {
          interps[h] += comps[h] *
                        shapeFunctionWeights[j * genCell->GetNumberOfPoints() + m];
        }
      }
      delete[] comps;
      // adding interpolated value to data of cell
      daGausses[id]->SetTuple(j + offset, interps.data());
    }
  }
  return numGaussPoints;
}


void GaussCubature::interpolateToGaussPoints()
{
  if (arrayIDs.empty())
  {
    std::cerr << "no arrays selected for interpolation" << std::endl;
    exit(1);
  }

  std::vector<vtkSmartPointer<vtkDoubleArray>> daGausses(arrayIDs.size());
  std::vector<vtkSmartPointer<vtkDataArray>> das(arrayIDs.size());
  numComponents.resize(arrayIDs.size());
  // initializing arrays storing interpolated data
  for (int id = 0; id < arrayIDs.size(); ++id)
  {
    // get desired point data array to be interpolated to gauss points
    vtkSmartPointer<vtkDataArray> da = dataSet->GetPointData()->GetArray(
        arrayIDs[id]);
    // get tuple length of given data
    int numComponent = da->GetNumberOfComponents();
    // declare data array to be populated with values at gauss points
    vtkSmartPointer<vtkDoubleArray> daGauss = vtkSmartPointer<vtkDoubleArray>::New();
    // names and sizing
    daGauss->SetName(
        dataSet->GetPointData()->GetArrayName(arrayIDs[id]));
    daGauss->SetNumberOfComponents(numComponent);
    daGauss->SetNumberOfTuples(gaussMesh->GetNumberOfPoints());
    das[id] = da;
    daGausses[id] = daGauss;
    numComponents[id] = numComponent;
    totalComponents += numComponent;
  }
  int numCells = dataSet->GetNumberOfCells();
#ifdef HAVE_OPENMP
  #pragma omp parallel default(none) shared(numCells, das, daGausses)
  {
    vtkSmartPointer<vtkGenericCell> genCell =
        vtkSmartPointer<vtkGenericCell>::New();
    #pragma omp single
    {
      // GetCell method is thread safe when first called from single thread
      dataSet->GetCell(0, genCell);
    }
    #pragma omp for schedule(static)
    for (int i = 0; i < numCells; ++i) {
      interpolateToGaussPointsAtCell(i, genCell, das, daGausses);
    }
  }
#else
  // generic cell to store given cell in dataSet
  vtkSmartPointer<vtkGenericCell> genCell = vtkSmartPointer<vtkGenericCell>::New();
  for (int i = 0; i < dataSet->GetNumberOfCells(); ++i)
  {
    interpolateToGaussPointsAtCell(i, genCell, das, daGausses);
  }
#endif
  for (int id = 0; id < arrayIDs.size(); ++id)
  {
    gaussMesh->GetPointData()->AddArray(daGausses[id]);
  }
}


void GaussCubature::interpolateToGaussPoints(
    const std::vector<std::string> &newArrayNames) {
  if (newArrayNames.empty()) {
    std::cerr << "no arrays selected for interpolation" << std::endl;
    exit(1);
  }

  std::vector<vtkSmartPointer<vtkDoubleArray>> daGausses(newArrayNames.size());
  std::vector<vtkSmartPointer<vtkDataArray>> das(newArrayNames.size());
  // initializing arrays storing interpolated data
  for (int id = 0; id < newArrayNames.size(); ++id) {
    // get desired point data array to be interpolated to gauss points
    vtkSmartPointer<vtkDataArray> da =
        dataSet->GetPointData()->GetArray(&(newArrayNames[id])[0u]);
    // get tuple length of given data
    int numComponent = da->GetNumberOfComponents();
    // declare data array to be populated with values at gauss points
    vtkSmartPointer<vtkDoubleArray> daGauss = vtkSmartPointer<vtkDoubleArray>::New();
    // names and sizing
    daGauss->SetName(&(newArrayNames[id])[0u]);
    daGauss->SetNumberOfComponents(numComponent);
    daGauss->SetNumberOfTuples(gaussMesh->GetNumberOfPoints());
    das[id] = da;
    daGausses[id] = daGauss;
  }
  int numCells = dataSet->GetNumberOfCells();
#ifdef HAVE_OPENMP
  #pragma omp parallel default(none) shared(numCells, das, daGausses)
  {
    // generic cell to store given cell in dataSet
    vtkSmartPointer<vtkGenericCell> genCell =
        vtkSmartPointer<vtkGenericCell>::New();
    #pragma omp single
    {
      // GetCell method is thread safe when first called from single thread
      dataSet->GetCell(0, genCell);
    }
    #pragma omp for schedule(static)
    for (int i = 0; i < numCells; ++i) {
      interpolateToGaussPointsAtCell(i, genCell, das, daGausses);
    }
  }
#else
  // generic cell to store given cell in dataSet
  vtkSmartPointer<vtkGenericCell> genCell =
      vtkSmartPointer<vtkGenericCell>::New();
  for (int i = 0; i < dataSet->GetNumberOfCells(); ++i) {
    interpolateToGaussPointsAtCell(i, genCell, das, daGausses);
  }
#endif
  for (int id = 0; id < arrayIDs.size(); ++id)
  {
    gaussMesh->GetPointData()->AddArray(daGausses[id]);
  }
}

void GaussCubature::integrateOverCell
    (int cellID,
     vtkSmartPointer<vtkGenericCell> genCell,
     vtkSmartPointer<vtkPointData> pd,
     std::vector<vtkSmartPointer<vtkDoubleArray>> &integralData,
     std::vector<std::vector<double>> &totalIntegralData) const
{
  // putting cell from nodeMesh into genCell
  dataSet->GetCell(cellID, genCell);
  // getting cellType for looking up numGaussPoints in dictionary
  // as well as computing scaled Jacobian
  int cellType = genCell->GetCellType();

  // get number of gauss points in cell from dictionary
  int numGaussPoints;
  const double* quadWeights;
  numGaussPoints = dict[cellType]->GetNumberOfQuadraturePoints();
  quadWeights = dict[cellType]->GetQuadratureWeights();

  // computing Jacobian for integration
  double jacobian = computeJacobian(genCell, cellType);
  // get quadrature weights for this cell type
  // get offset from nodeMesh for lookup of gauss points in polyData
  int offset = getOffset(cellID);
  // holds integrated data for each array
  std::vector<std::vector<double>> data(arrayIDs.size());
  // integration loop
  for (int j = 0; j < integralData.size(); ++j) {
    int numComponent = integralData[j]->GetNumberOfComponents();
    data[j].resize(numComponent, 0.0);
    auto comps = new double[numComponent];
    for (int i = 0; i < numGaussPoints; ++i) {
      pd->GetArray(j)->GetTuple(offset + i, comps);
      for (int k = 0; k < numComponent; ++k) {
        // TODO: generalize to support surface integration
        if (genCell->GetCellDimension() == 3) {
          data[j][k] += comps[k] * quadWeights[i];//*jacobian;
        } else
          data[j][k] += 0.0;
      }
    }
    delete[] comps;
    for (int k = 0; k < numComponent; ++k) {
      data[j][k] *= jacobian;
      totalIntegralData[j][k] += data[j][k];
    }
    // adding integrated value to data of cell
    integralData[j]->SetTuple(cellID, data[j].data());
  }
}


void GaussCubature::integrateOverCell
    (int cellID,
     vtkSmartPointer<vtkGenericCell> genCell,
     vtkSmartPointer<vtkPointData> pd,
     std::vector<vtkSmartPointer<vtkDoubleArray>> &integralData,
     std::vector<std::vector<double>> &totalIntegralData,
     const std::vector<std::string> &newArrayNames,
     bool computeRMSE) const
{
  // putting cell from nodeMesh into genCell
  dataSet->GetCell(cellID, genCell);
  // getting cellType for looking up numGaussPoints in dictionary
  int cellType = dataSet->GetCell(cellID)->GetCellType();
  // get number of gauss points in cell from dictionary
  int numGaussPoints = dict[cellType]->GetNumberOfQuadraturePoints();
  // computing Jacobian for integration
  double jacobian = computeJacobian(genCell, cellType);
  // computing volume for RMSE
  double volume = computeCellVolume(genCell, cellType);
  // get quadrature weights for this cell type
  const double *quadWeights = dict[cellType]->GetQuadratureWeights();
  // get offset from nodeMesh for lookup of gauss points in polyData
  int offset = getOffset(cellID);
  // holds integrated data for each array
  std::vector<std::vector<double>> data(arrayIDs.size());
  // integration loop
  for (int j = 0; j < integralData.size(); ++j)
  {
    int numComponent = integralData[j]->GetNumberOfComponents();
    data[j].resize(numComponent);
    auto *comps = new double[numComponent];
    for (int i = 0; i < numGaussPoints; ++i)
    {
      pd->GetArray(&(newArrayNames[j])[0u])->GetTuple(offset + i, comps);
      for (int k = 0; k < numComponent; ++k)
      {
        // TODO: generalize to support surface integration
        if (genCell->GetCellDimension() == 3)
        {
          data[j][k] += comps[k] * quadWeights[i];
        }
        else
          data[j][k] += 0.0;
      }
    }
    delete[] comps;
    // taking sqrt of integrated value (RMSE)
    for (int k = 0; k < numComponent; ++k)
    {
      data[j][k] *= jacobian;
      data[j][k] = (computeRMSE ? std::sqrt(data[j][k] / volume) : data[j][k]);
      //if (computeRMSE)
      //  data[j][k] = std::sqrt(data[j][k] / volume);
      totalIntegralData[j][k] += data[j][k];
    }
    // adding integrated value to data of cell
    integralData[j]->SetTuple(cellID, data[j].data());
  }
}

std::vector<std::vector<double>> GaussCubature::integrateOverAllCells()
{
  int numCells = dataSet->GetNumberOfCells();

  if (gaussMesh->GetPointData()->GetNumberOfArrays() == 0)
  {
    interpolateToGaussPoints();
  }

  vtkSmartPointer<vtkPointData> pd = gaussMesh->GetPointData();
  std::vector<vtkSmartPointer<vtkDoubleArray>> cellIntegralData(arrayIDs.size());
  std::vector<std::vector<double>> totalIntegralData(arrayIDs.size());
  for (int id = 0; id < arrayIDs.size(); ++id)
  {
    std::string arrName(
        dataSet->GetPointData()->GetArrayName(arrayIDs[id]));
    arrName.append("Integral");
    vtkSmartPointer<vtkDoubleArray> integralDatum = vtkSmartPointer<vtkDoubleArray>::New();
    integralDatum->SetName(&arrName[0u]);
    integralDatum->SetNumberOfComponents(numComponents[id]);
    integralDatum->SetNumberOfTuples(dataSet->GetNumberOfCells());
    cellIntegralData[id] = integralDatum;
    totalIntegralData[id].resize(numComponents[id], 0);
  }

#ifdef HAVE_OPENMP
  #pragma omp parallel default(none) shared(numCells, pd, cellIntegralData, totalIntegralData, cerr)
  {
    auto genCell = vtkSmartPointer<vtkGenericCell>::New();

    #pragma omp single
    {
      // GetCell method is thread safe when first called from single thread
      dataSet->GetCell(0, genCell);
    }
    auto partialIntegralData = totalIntegralData;

    #pragma omp for schedule(static)
    for (int i = 0; i < numCells; ++i) {
      integrateOverCell(i, genCell, pd, cellIntegralData, partialIntegralData);
    }

    #pragma omp critical
    {
      // total integral data has same dimensions as cell integral data ...
      // sum partial integrals
      for (int id = 0; id < partialIntegralData.size(); ++id) {
        for (int k = 0; k < partialIntegralData[id].size(); ++k) {
          totalIntegralData[id][k] += partialIntegralData[id][k];
        }
      }
    }
  }
#else
  auto genCell = vtkSmartPointer<vtkGenericCell>::New();
  for (int i = 0; i < numCells; ++i) {
    integrateOverCell(i, genCell, pd, cellIntegralData, totalIntegralData);
  }
#endif

  for (int id = 0; id < arrayIDs.size(); ++id) {
    dataSet->GetCellData()->AddArray(cellIntegralData[id]);
  }
  return totalIntegralData;
}


std::vector<std::vector<double>>
GaussCubature::integrateOverAllCells(
    const std::vector<std::string> &newArrayNames,
    bool computeRMSE)
{
  interpolateToGaussPoints(newArrayNames);

  vtkSmartPointer<vtkPointData> pd = gaussMesh->GetPointData();
  std::vector<vtkSmartPointer<vtkDoubleArray>> cellIntegralData(newArrayNames.size());
  std::vector<std::vector<double>> totalIntegralData(newArrayNames.size());
  for (int id = 0; id < newArrayNames.size(); ++id)
  {
    std::string name = newArrayNames[id] + "Integral";
    vtkSmartPointer<vtkDoubleArray> integralDatum = vtkSmartPointer<vtkDoubleArray>::New();
    integralDatum->SetName(&name[0u]);
    int numComponent = pd->GetArray(
        &(newArrayNames[id])[0u])->GetNumberOfComponents();
    integralDatum->SetNumberOfComponents(numComponent);
    integralDatum->SetNumberOfTuples(dataSet->GetNumberOfCells());
    cellIntegralData[id] = integralDatum;
    totalIntegralData[id].resize(numComponent, 0);
  }

  int numCells = dataSet->GetNumberOfCells();

#ifdef HAVE_OPENMP
  #pragma omp parallel default(none) \
  shared(numCells, pd, cellIntegralData, totalIntegralData, newArrayNames, computeRMSE)
  {
    auto genCell = vtkSmartPointer<vtkGenericCell>::New();
    #pragma omp single
    {
      // GetCell method is thread safe when first called from single thread
      dataSet->GetCell(0, genCell);
    }
    auto partialIntegralData = totalIntegralData;
    #pragma omp for schedule(static)
    for (int i = 0; i < numCells; ++i) {
      integrateOverCell(i, genCell, pd, cellIntegralData, partialIntegralData,
                        newArrayNames, computeRMSE);
    }
    #pragma omp critical
    {
      // total integral data has same dimensions as cell integral data ...
      // sum partial integrals
      for(int id = 0; id < partialIntegralData.size(); ++id) {
        for(int k = 0; k < partialIntegralData[id].size(); ++k) {
          totalIntegralData[id][k] += partialIntegralData[id][k];
        }
      }
    }
  }
#else
  vtkSmartPointer<vtkGenericCell> genCell = vtkSmartPointer<vtkGenericCell>::New();
  for (int i = 0; i < numCells; ++i) {
    integrateOverCell(i, genCell, pd, cellIntegralData, totalIntegralData,
                      newArrayNames, computeRMSE);
  }
#endif

  for (int id = 0; id < newArrayNames.size(); ++id)
  {
    dataSet->GetCellData()->AddArray(cellIntegralData[id]);
  }
  return totalIntegralData;
}


void GaussCubature::writeGaussMesh(const char *name) const
{
  if (gaussMesh)
    writeVTFile<vtkXMLPolyDataWriter>(name, gaussMesh);
  else
  {
    std::cerr << "Gauss point mesh has not been constructed" << std::endl;
    exit(1);
  }
}
