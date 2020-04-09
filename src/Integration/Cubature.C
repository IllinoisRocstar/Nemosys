#include "Cubature.H"

#include <vtkQuadraturePointsGenerator.h>
#include <vtkMeshQuality.h>
#include <vtkInformationQuadratureSchemeDefinitionVectorKey.h>
#include <vtkCellData.h>
#include <vtkInformation.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkMesh.H> // for writeVTFile


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


GaussCubature::GaussCubature(meshBase *_nodeMesh)
    : nodeMesh(_nodeMesh), numVolCells(0), totalComponents(0)
{
  nodeMesh->unsetCellDataArray("QuadratureOffSet");
  constructGaussMesh();
}

GaussCubature::GaussCubature(meshBase *_nodeMesh,
                             const std::vector<int> &_arrayIDs)
    : nodeMesh(_nodeMesh), numVolCells(0), arrayIDs(_arrayIDs),
      totalComponents(0)
{
  nodeMesh->unsetCellDataArray("QuadratureOffSet");
  constructGaussMesh();
  interpolateToGaussPoints();
}

GaussCubature *GaussCubature::Create(meshBase *nodeMesh)
{
  return new GaussCubature(nodeMesh);
}

GaussCubature *GaussCubature::Create(meshBase *nodeMesh,
                                     const std::vector<int> &arrayIDs)
{
  return new GaussCubature(nodeMesh, arrayIDs);
}

std::unique_ptr<GaussCubature> GaussCubature::CreateUnique(meshBase *nodeMesh)
{
  return std::unique_ptr<GaussCubature>(GaussCubature::Create(nodeMesh));
}

std::unique_ptr<GaussCubature>
GaussCubature::CreateUnique(meshBase *nodeMesh,
                            const std::vector<int> &arrayIDs)
{
  return std::unique_ptr<GaussCubature>(
      GaussCubature::Create(nodeMesh, arrayIDs));
}

std::shared_ptr<GaussCubature>
GaussCubature::CreateShared(meshBase *nodeMesh)
{
  std::shared_ptr<GaussCubature> cuby;
  cuby.reset(GaussCubature::Create(nodeMesh));
  return cuby;
}

std::shared_ptr<GaussCubature>
GaussCubature::CreateShared(meshBase *nodeMesh,
                            const std::vector<int> &arrayIDs)
{
  std::shared_ptr<GaussCubature> cuby;
  cuby.reset(GaussCubature::Create(nodeMesh, arrayIDs));
  return cuby;
}


void GaussCubature::constructGaussMesh()
{
  // check whether arrayIDs exist in mesh
  {
    vtkSmartPointer<vtkPointData> pd = nodeMesh->getDataSet()->GetPointData();
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
  vtkInformationQuadratureSchemeDefinitionVectorKey *key =
      vtkQuadratureSchemeDefinition::DICTIONARY();

  // Get the cell types used by the data set
  vtkSmartPointer<vtkCellTypes> cellTypes
      = vtkSmartPointer<vtkCellTypes>::New();
  nodeMesh->getDataSet()->GetCellTypes(cellTypes);
  int nCellTypes = cellTypes->GetNumberOfTypes();

  // create offset array and store the dictionary within
  vtkSmartPointer<vtkIdTypeArray> offsets
      = vtkSmartPointer<vtkIdTypeArray>::New();
  std::string basename = "QuadratureOffset";
  offsets->SetName(basename.c_str());
  vtkSmartPointer<vtkInformation> info = vtkSmartPointer<vtkInformation>::New();
  info = offsets->GetInformation();

  for (int typeId = 0; typeId < nCellTypes; ++typeId)
  {
    int cellType = cellTypes->GetCellType(typeId);
    // Initialize quadrature scheme definition for given cell type
    vtkSmartPointer<vtkQuadratureSchemeDefinition> def
        = vtkSmartPointer<vtkQuadratureSchemeDefinition>::New();
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
  offsets->SetNumberOfTuples(nodeMesh->getDataSet()->GetNumberOfCells());
  vtkIdType offset = 0;

  for (int cellid = 0;
       cellid < nodeMesh->getDataSet()->GetNumberOfCells(); ++cellid)
  {
    offsets->SetValue(cellid, offset);
    int cellType = nodeMesh->getDataSet()->GetCell(cellid)->GetCellType();
    if (cellType >= VTK_TETRA)
      numVolCells += 1;
    vtkQuadratureSchemeDefinition *celldef = dict[cellType];
    offset += celldef->GetNumberOfQuadraturePoints();
  }

  nodeMesh->getDataSet()->GetCellData()->AddArray(offsets);

  vtkSmartPointer<vtkQuadraturePointsGenerator> pointGen =
      vtkSmartPointer<vtkQuadraturePointsGenerator>::New();

  pointGen->SetInputArrayToProcess
      (0, 0, 0,
       vtkDataObject::FIELD_ASSOCIATION_CELLS,
       "QuadratureOffset");
  pointGen->SetInputData(nodeMesh->getDataSet());
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
      nodeMesh->getDataSet()->GetCellData()->GetArray("QuadratureOffset"))
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
  int numGaussPoints = dict[nodeMesh->getDataSet()->GetCell(
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
  nodeMesh->getDataSet()->GetCell(cellID, genCell);
  // getting cellType information for lookup in map
  int cellType = nodeMesh->getDataSet()->GetCellType(cellID);
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
    vtkSmartPointer<vtkDataArray> da = nodeMesh->getDataSet()->GetPointData()->GetArray(
        arrayIDs[id]);
    // get tuple length of given data
    int numComponent = da->GetNumberOfComponents();
    // declare data array to be populated with values at gauss points
    vtkSmartPointer<vtkDoubleArray> daGauss = vtkSmartPointer<vtkDoubleArray>::New();
    // names and sizing
    daGauss->SetName(
        nodeMesh->getDataSet()->GetPointData()->GetArrayName(arrayIDs[id]));
    daGauss->SetNumberOfComponents(numComponent);
    daGauss->SetNumberOfTuples(gaussMesh->GetNumberOfPoints());
    das[id] = da;
    daGausses[id] = daGauss;
    numComponents[id] = numComponent;
    totalComponents += numComponent;
  }
  // generic cell to store given cell in nodeMesh->getDataSet()
  vtkSmartPointer<vtkGenericCell> genCell = vtkSmartPointer<vtkGenericCell>::New();
  for (int i = 0; i < nodeMesh->getNumberOfCells(); ++i)
  {
    interpolateToGaussPointsAtCell(i, genCell, das, daGausses);
  }
  for (int id = 0; id < arrayIDs.size(); ++id)
  {
    gaussMesh->GetPointData()->AddArray(daGausses[id]);
  }
}


void GaussCubature::interpolateToGaussPoints(
    const std::vector<std::string> &newArrayNames)
{
  if (newArrayNames.empty())
  {
    std::cerr << "no arrays selected for interpolation" << std::endl;
    exit(1);
  }

  std::vector<vtkSmartPointer<vtkDoubleArray>> daGausses(newArrayNames.size());
  std::vector<vtkSmartPointer<vtkDataArray>> das(newArrayNames.size());
  // initializing arrays storing interpolated data
  for (int id = 0; id < newArrayNames.size(); ++id)
  {
    // get desired point data array to be interpolated to gauss points
    vtkSmartPointer<vtkDataArray> da
        = nodeMesh->getDataSet()->GetPointData()->GetArray(
            &(newArrayNames[id])[0u]);
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
  // generic cell to store given cell in nodeMesh->getDataSet()
  vtkSmartPointer<vtkGenericCell> genCell = vtkSmartPointer<vtkGenericCell>::New();
  for (int i = 0; i < nodeMesh->getNumberOfCells(); ++i)
  {
    interpolateToGaussPointsAtCell(i, genCell, das, daGausses);
  }
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
  nodeMesh->getDataSet()->GetCell(cellID, genCell);
  // getting cellType for looking up numGaussPoints in dictionary
  // as well as computing scaled Jacobian
  int cellType = nodeMesh->getDataSet()->GetCell(cellID)->GetCellType();
  // get number of gauss points in cell from dictionary
  int numGaussPoints = dict[cellType]->GetNumberOfQuadraturePoints();
  // computing Jacobian for integration
  double jacobian = computeJacobian(genCell, cellType);
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
    data[j].resize(numComponent, 0.0);
    auto *comps = new double[numComponent];
    for (int i = 0; i < numGaussPoints; ++i)
    {
      pd->GetArray(j)->GetTuple(offset + i, comps);
      for (int k = 0; k < numComponent; ++k)
      {
        // TODO: generalize to support surface integration
        if (genCell->GetCellDimension() == 3)
        {
          data[j][k] += comps[k] * quadWeights[i];//*jacobian;
        }
        else
          data[j][k] += 0.0;
      }
    }
    delete[] comps;
    for (int k = 0; k < numComponent; ++k)
    {
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
  nodeMesh->getDataSet()->GetCell(cellID, genCell);
  // getting cellType for looking up numGaussPoints in dictionary
  int cellType = nodeMesh->getDataSet()->GetCell(cellID)->GetCellType();
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
  if (gaussMesh->GetPointData()->GetNumberOfArrays() == 0)
  {
    interpolateToGaussPoints();
  }

  vtkSmartPointer<vtkPointData> pd = gaussMesh->GetPointData();
  std::vector<vtkSmartPointer<vtkDoubleArray>> integralData(arrayIDs.size());
  std::vector<std::vector<double>> totalIntegralData(arrayIDs.size());
  for (int id = 0; id < arrayIDs.size(); ++id)
  {
    std::string arrName(
        nodeMesh->getDataSet()->GetPointData()->GetArrayName(arrayIDs[id]));
    arrName.append("Integral");
//    std::cout << arrName << std::endl;
    vtkSmartPointer<vtkDoubleArray> integralDatum = vtkSmartPointer<vtkDoubleArray>::New();
    integralDatum->SetName(&arrName[0u]);
    integralDatum->SetNumberOfComponents(numComponents[id]);
    integralDatum->SetNumberOfTuples(nodeMesh->getNumberOfCells());
    integralData[id] = integralDatum;
    totalIntegralData[id].resize(numComponents[id], 0);
  }
  vtkSmartPointer<vtkGenericCell> genCell = vtkSmartPointer<vtkGenericCell>::New();
  for (int i = 0; i < nodeMesh->getNumberOfCells(); ++i)
  {
    integrateOverCell(i, genCell, pd, integralData, totalIntegralData);
  }

  for (int id = 0; id < arrayIDs.size(); ++id)
  {
    nodeMesh->getDataSet()->GetCellData()->AddArray(integralData[id]);
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
  std::vector<vtkSmartPointer<vtkDoubleArray>> integralData(
      newArrayNames.size());
  std::vector<std::vector<double>> totalIntegralData(newArrayNames.size());
  for (int id = 0; id < newArrayNames.size(); ++id)
  {
    std::string name = newArrayNames[id] + "Integral";
    vtkSmartPointer<vtkDoubleArray> integralDatum = vtkSmartPointer<vtkDoubleArray>::New();
    integralDatum->SetName(&name[0u]);
    int numComponent = pd->GetArray(
        &(newArrayNames[id])[0u])->GetNumberOfComponents();
    integralDatum->SetNumberOfComponents(numComponent);
    integralDatum->SetNumberOfTuples(nodeMesh->getNumberOfCells());
    integralData[id] = integralDatum;
    totalIntegralData[id].resize(numComponent, 0);
  }
  vtkSmartPointer<vtkGenericCell> genCell = vtkSmartPointer<vtkGenericCell>::New();
  for (int i = 0; i < nodeMesh->getNumberOfCells(); ++i)
  {
    integrateOverCell(i, genCell, pd, integralData, totalIntegralData,
                      newArrayNames, computeRMSE);
  }

  for (int id = 0; id < newArrayNames.size(); ++id)
  {
    nodeMesh->getDataSet()->GetCellData()->AddArray(integralData[id]);
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

//void GaussCubature::constructGaussMesh(const std::vector<int>& arrayIDs)
//{
//  // building poly data
//  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
//  vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();
//  for (int i = 0; i < nodeMesh->getDataSet()->GetNumberOfCells(); ++i)
//  {
//    std::vector<std::vector<double>> gaussPoints = getGaussPointsAtCell(i);
//    for (int j = 0; j < gaussPoints.size(); ++j)
//    {
//      vtkIdType id[1];
//      id[0] = points->InsertNextPoint(gaussPoints[j].data());
//      vertices->InsertNextCell(1,id);
//    }  
//  }
//  vtkSmartPointer<vtkPolyData> gaussMesh = vtkSmartPointer<vtkPolyData>::New();
//  gaussMesh->SetPoints(points);
//  gaussMesh->SetVerts(vertices);
//  interpolateToGaussPoints(gaussMesh,arrayIDs);
//  writeVTFile<vtkXMLPolyDataWriter> ("gaussTest.vtp",gaussMesh);   
//}


//  // allocate storate for polygons
//  //  gaussMesh->Allocate();
//  for (int i = 0; i < nodeMesh->getNumberOfCells(); ++i)
//  {
//    //gaussMesh->InsertNextCell(VTK_POLYGON,polyCellIds); 
//    polyPnt += interpolateToGaussPointsAtCell(i,genCell,das,daGausses,numComponents,polyPnt);
//    //polyPnt += numGaussPoints;
//  }
//
//int GaussCubature::getNumGaussPointsForCellType(int cellType)
//{
//  int numGaussPoints;
//  switch(cellType)
//  {
//    case VTK_TRIANGLE:
//    {
//      numGaussPoints = 3;
//      break;
//    }
//    case VTK_TETRA:
//    { 
//      numGaussPoints = 4;
//      break;
//    }
//    default:
//    {
//      std::cerr << "Error: Cell type: " << cellType << "found "
//                << "with no quadrature definition provided" << std::endl;
//      exit(1);
//    }
//  }
//  return numGaussPoints; 
//}

//void GaussCubature::buildMap()
//{
//  // building quadrature scheme map
//  vtkSmartPointer<vtkCellTypes> cellTypes 
//    = vtkSmartPointer<vtkCellTypes>::New();
//  nodeMesh->getDataSet()->GetCellTypes(cellTypes);
//  int nCellTypes = cellTypes->GetNumberOfTypes(); 
//  for (int i = 0; i < nCellTypes; ++i)
//  {
//    int cellType = cellTypes->GetCellType(i);
//    nGaussForCellTMap[cellType] = getNumGaussPointsForCellType(cellType);
//  }
//}


//std::vector<std::vector<double>> GaussCubature::getGaussPointsAtCell(int cellID)
//{
//  // get cell type for quadrature scheme definition
//  int cellType = nodeMesh->getDataSet()->GetCellType(cellID);
//  // this vector holds the shape function evaluated at parametric coord of gauss point
//  std::vector<double> shapeFuncAtGauss;
//  std::vector<double>::iterator beg = shapeFuncAtGauss.begin();
//  // dispatch number of gauss points by cell type
//  int numGaussPoints;
//  switch(cellType)
//  {
//    case VTK_TRIANGLE:
//    {
//      shapeFuncAtGauss.insert(beg, TRI3, TRI3+9);
//      numGaussPoints = 3;
//      break;
//    }
//    case VTK_TETRA:
//    { 
//      shapeFuncAtGauss.insert(beg, TET4, TET4+16); 
//      numGaussPoints = 4;
//      break;
//    }
//    default:
//    {
//      std::cerr << "Error: Cell type: " << cellType << "found "
//                << "with no quadrature definition provided" << std::endl;
//      exit(1);
//    }
//  }
//
//  vtkSmartPointer<vtkGenericCell> genCell = vtkSmartPointer<vtkGenericCell>::New();
//  nodeMesh->getDataSet()->GetCell(cellID,genCell);
//  int numPointsInCell = genCell->GetNumberOfPoints();
//  std::vector<std::vector<double>> gaussPoints;
//  gaussPoints.resize(numGaussPoints);
//  for (int j = 0; j < gaussPoints.size(); ++j)
//  {
//    gaussPoints[j].resize(3,0);
//    for (int k = 0; k < numPointsInCell; ++k)
//    {
//      int pntID = genCell->GetPointId(k);
//      double x[3];
//      nodeMesh->getDataSet()->GetPoint(pntID,x); 
//      for (int i = 0; i < 3; ++i)
//      {
//        gaussPoints[j][i] += shapeFuncAtGauss[j*numGaussPoints + k]*x[i]; 
//      }
//    }
//  }
//  return gaussPoints; 
//}


