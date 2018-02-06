#include <Cubature.H>



double TRI3 [] = 
{
.666666666666666, .166666666666666, .166666666666666,
.166666666666666, .666666666666666, .166666666666666,
.166666666666666, .166666666666666, .666666666666666
};
double TRI3W [] = {0.333333333333333,0.333333333333333,0.333333333333333};

double TET4 [] = 
{
.585410196624968, .138196601125010, .138196601125010, .138196601125010,
.138196601125010, .585410196624968, .138196601125010, .138196601125010,
.138196601125010, .138196601125010, .585410196624968, .138196601125010,
.138196601125010, .138196601125010, .138196601125010, .585410196624968 
};
double TET4W [] = {0.25,0.25,0.25,0.25};

GaussCubature::GaussCubature(meshBase* _nodeMesh):nodeMesh(_nodeMesh),totalComponents(0)
{ 
  nodeMesh->unsetCellDataArray("QuadratureOffSet");
  constructGaussMesh();
} 

GaussCubature::GaussCubature(meshBase* _nodeMesh, const std::vector<int>& _arrayIDs)
  : nodeMesh(_nodeMesh), arrayIDs(_arrayIDs),totalComponents(0)
{
  nodeMesh->unsetCellDataArray("QuadratureOffSet");
  constructGaussMesh();
  interpolateToGaussPoints(); 
} 

GaussCubature* GaussCubature::Create(meshBase* nodeMesh)
{
  return new GaussCubature(nodeMesh);
}

GaussCubature* GaussCubature::Create(meshBase* nodeMesh, const std::vector<int>& arrayIDs)
{
  return new GaussCubature(nodeMesh, arrayIDs);
}

std::unique_ptr<GaussCubature> GaussCubature::CreateUnique(meshBase* nodeMesh)
{
  return std::unique_ptr<GaussCubature>(GaussCubature::Create(nodeMesh));
}

std::unique_ptr<GaussCubature>
GaussCubature::CreateUnique(meshBase* nodeMesh, const std::vector<int>& arrayIDs)
{
  return std::unique_ptr<GaussCubature>(GaussCubature::Create(nodeMesh, arrayIDs));
}

std::shared_ptr<GaussCubature> 
GaussCubature::CreateShared(meshBase* nodeMesh)
{
  std::shared_ptr<GaussCubature> cuby;
  cuby.reset(GaussCubature::Create(nodeMesh));
  return cuby; 
}

std::shared_ptr<GaussCubature>
GaussCubature::CreateShared(meshBase* nodeMesh, const std::vector<int>& arrayIDs)
{
  std::shared_ptr<GaussCubature> cuby;
  cuby.reset(GaussCubature::Create(nodeMesh,arrayIDs));
  return cuby;
}



double GaussCubature::computeScaledJacobian(vtkSmartPointer<vtkGenericCell> genCell, int cellType)
{
  double jacobi = 0;
  switch(cellType)
  {
    case VTK_TRIANGLE:
    {
      double x1[3],x2[3],x3[3];
      genCell->Points->GetPoint(0,x1);
      genCell->Points->GetPoint(1,x2);
      genCell->Points->GetPoint(2,x3);
      jacobi = vtkTriangle::TriangleArea(x1,x2,x3);
      break;
    }
    case VTK_TETRA:
    { 
      double pcoords[3];
      double derivs[12];
      genCell->InterpolateDerivs(pcoords,derivs);
      Matrix3d jacobianMat = Matrix3d::Zero();
      double x[3];
      //double *m[3]; 
      //double m0[3],m1[3],m2[3];
      //m[0] = m0; m[1] = m1; m[2] = m2;
      //for (int i = 0; i < 3; ++i)
      //  m0[i] = m1[i] = m2[i] = 0.0;

      for (int j = 0; j < 4; ++j)
      {
        genCell->Points->GetPoint(j,x);
        for (int i = 0; i < 3; ++i)
        {
          jacobianMat(i,0) += x[i] * derivs[j];
          jacobianMat(i,1) += x[i] * derivs[j+4];
          jacobianMat(i,2) += x[i] * derivs[j+8];
          //m0[i] += x[i] * derivs[j];
          //m1[i] += x[i] * derivs[j+4];
          //m2[i] += x[i] * derivs[j+8];
        }
      }
      //vtkTetra* tetcell = vtkTetra::SafeDownCast(genCell);//->SetCellTypeToTetra();
      //std::vector<std::vector<double>> x(4);
      //double x1[3], x2[3], x3[3], x4[3];
      //genCell->Points->GetPoint(0, x1);
      //genCell->Points->GetPoint(1, x2);
      //genCell->Points->GetPoint(2, x3);
      //genCell->Points->GetPoint(3, x4);

      //jacobi = tetcell->ComputeVolume(x1,x2,x3,x4);
      jacobi = jacobianMat.determinant()/6.0;
     //jacobi = vtkMath::Determinant3x3(m[0],m[1],m[2]);
      break;
    }
    default:
    {
      std::cerr << "Error: Cell type: " << cellType << "found "
                << "with no quadrature definition provided" << std::endl;
      exit(1);
    }
  }
  return jacobi;
}


// get offset from nodeMesh for lookup of gauss points in polydata
int GaussCubature::getOffset(int cellID)
{
  int offset;
  {
    double offsets[1];
    nodeMesh->getDataSet()->GetCellData()->GetArray("QuadratureOffset")->GetTuple(cellID,offsets);
    offset = offsets[0];
  }
  return offset;
}

pntDataPairVec GaussCubature::getGaussPointsAndDataAtCell(int cellID)
{
  if (arrayIDs.size() == 0)
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
  int numGaussPoints = dict[nodeMesh->getDataSet()->GetCell(cellID)->GetCellType()]
                        ->GetNumberOfQuadraturePoints();
  // get offset from nodeMesh for lookup of gauss points in polydata
  int offset = getOffset(cellID);
  //pntDataPairVec container;
  pntDataPairVec container(numGaussPoints);

  vtkSmartPointer<vtkPointData> pd = gaussMesh->GetPointData();
  for (int i = 0; i < numGaussPoints; ++i)
  {
    double x_tmp[3];
    gaussMesh->GetPoint(offset+i,x_tmp);
    std::vector<double> gaussPnt(x_tmp,x_tmp + 3); 
    //std::vector<std::vector<double>> data(numComponents.size());
    std::vector<double> data(totalComponents);
    //std::vector<VectorXd> data(numComponents.size());
		int currcomp = 0;
    for (int j = 0; j < numComponents.size(); ++j)
    {
      //data[j].resize(numComponents[j]);
      double comps[numComponents[j]];
      pd->GetArray(j)->GetTuple(offset+i,comps);
      //std::vector<double> comps(comps_tmp,comps_tmp+numComponents[j]);
     // data[j] = std::move(comps);
    	for (int k = 0; k < numComponents[j]; ++k)
			{
        data[currcomp] = comps[k];
				//data[j](k) = comps[k];
			  ++currcomp;
      }
		}
    container[i] = std::move(std::make_pair(gaussPnt,data));
  }
  
  return container;
}

void GaussCubature::integrateOverCell
                      (int cellID,
                       vtkSmartPointer<vtkGenericCell> genCell,
                       vtkSmartPointer<vtkPointData> pd,
                       std::vector<vtkSmartPointer<vtkDoubleArray>>& integralData,
                       std::vector<std::vector<double>>& totalIntegralData)
{
  // putting cell from nodemesh into gencell
  nodeMesh->getDataSet()->GetCell(cellID,genCell);
  // getting celltype for looking up numGaussPoints in dictionary 
  // as well as computing scaled jacobian
  int cellType = nodeMesh->getDataSet()->GetCell(cellID)->GetCellType();
  // get number of gauss points in cell from dictionary
  int numGaussPoints = dict[cellType]
                        ->GetNumberOfQuadraturePoints();
  // computing scaled jacobian for proper integration
  double scaledJacobian = computeScaledJacobian(genCell,cellType);
  // get quadrature weights for this cell type
  const double* quadWeights = dict[cellType]->GetQuadratureWeights();
  // get offset from nodeMesh for lookup of gauss points in polydata
  int offset = getOffset(cellID);
  // holds integrated data for each array
  std::vector<std::vector<double>> data(arrayIDs.size());
  // integration loop
  for (int j = 0; j < arrayIDs.size(); ++j)
  {
    data[j].resize(numComponents[j]);
    double comps[numComponents[j]];
    for (int i = 0; i < numGaussPoints; ++i)
    {
      pd->GetArray(j)->GetTuple(offset+i,comps);
      for (int k = 0; k < numComponents[j]; ++k)
      {
        data[j][k] += comps[k]*quadWeights[i]*scaledJacobian;
        totalIntegralData[j][k] += data[j][k]; 
      }
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
    std::string arrName(nodeMesh->getDataSet()->GetPointData()->GetArrayName(arrayIDs[id]));
    arrName.append("Integral");
    std::cout << arrName << std::endl;
    vtkSmartPointer<vtkDoubleArray> integralDatum = vtkSmartPointer<vtkDoubleArray>::New();
    integralDatum->SetName(&arrName[0u]);
    integralDatum->SetNumberOfComponents(numComponents[id]);
    integralDatum->SetNumberOfTuples(nodeMesh->getNumberOfCells());
    integralData[id] = integralDatum;
    totalIntegralData[id].resize(numComponents[id],0);
  }
  vtkSmartPointer<vtkGenericCell> genCell = vtkSmartPointer<vtkGenericCell>::New();
  for (int i = 0; i < nodeMesh->getNumberOfCells(); ++i)
  {
    integrateOverCell(i,genCell,pd,integralData,totalIntegralData);
  }
  
	for (int id = 0; id < arrayIDs.size(); ++id)
  {	
    nodeMesh->getDataSet()->GetCellData()->AddArray(integralData[id]);
  }
  return totalIntegralData;
}


int GaussCubature::interpolateToGaussPointsAtCell
                    (const int cellID,
                     vtkSmartPointer<vtkGenericCell> genCell,
                     const std::vector<vtkSmartPointer<vtkDataArray>>& das,
                     std::vector<vtkSmartPointer<vtkDoubleArray>>& daGausses)
{
  // putting current cell into genCell
  nodeMesh->getDataSet()->GetCell(cellID, genCell);
  // getting cellType information for lookup in map
  int cellType = nodeMesh->getDataSet()->GetCellType(cellID);
  // number of gauss points in this cell
  //int numGaussPoints = dict[cellType]->GetNumberOfQuadraturePoints();
  int numGaussPoints = dict[cellType]->GetNumberOfQuadraturePoints();
  // get offset from nodeMesh for lookup of gauss points in polydata
  int offset = getOffset(cellID);
  // interpolation loop
  for (int j = 0; j < numGaussPoints; ++j)
  {
    double x[3];
    gaussMesh->GetPoint(j+offset,x); 
    // parameters for interpolation
    int subId; // not used
    double minDist2; // not used
    double pcoords[3];
    double weights[genCell->GetNumberOfPoints()];
    genCell->EvaluatePosition(x,NULL,subId,pcoords,minDist2,weights);
    for (int id = 0; id < das.size(); ++id)
		{
      double comps[numComponents[id]];
			std::vector<double> interps(numComponents[id],0.0);
			for (int m = 0; m < genCell->GetNumberOfPoints(); ++m)
			{
				int pntId = genCell->GetPointId(m);
				das[id]->GetTuple(pntId, comps);
				for (int h = 0; h < numComponents[id]; ++h)
				{
					interps[h] += comps[h]*weights[m]; 
				}
			}
      // adding interpolated value to data of cell
			daGausses[id]->SetTuple(j+offset,interps.data());
    }
	}
  return numGaussPoints;
}

void GaussCubature::interpolateToGaussPoints()
{
  if (arrayIDs.size() == 0)
  {
    std::cerr << "no arrays selected for interpolation" << std::endl;
    exit(1);
  }
  std::vector<vtkSmartPointer<vtkDoubleArray>> daGausses(arrayIDs.size());
	std::vector<vtkSmartPointer<vtkDataArray>> das(arrayIDs.size());
	numComponents.resize(arrayIDs.size());
	for (int id = 0; id < arrayIDs.size(); ++id)
	{
    // get desired point data array to be interpolated to gauss points
    vtkSmartPointer<vtkDataArray> da = nodeMesh->getDataSet()->GetPointData()->GetArray(arrayIDs[id]);
    // get tuple length of given data
    int numComponent = da->GetNumberOfComponents();
    // declare data array to be populated with values at gauss points
		vtkSmartPointer<vtkDoubleArray> daGauss = vtkSmartPointer<vtkDoubleArray>::New();
    // names and sizing
    daGauss->SetName(nodeMesh->getDataSet()->GetPointData()->GetArrayName(arrayIDs[id]));
    daGauss->SetNumberOfComponents(numComponent);
		daGauss->SetNumberOfTuples(gaussMesh->GetNumberOfPoints());
    das[id] = da;
		daGausses[id] = daGauss;
		numComponents[id] = numComponent;	
	  totalComponents += numComponent;
  }
  // generic cell to store given cell in nodeMesh->getDataSet()
  vtkSmartPointer<vtkGenericCell> genCell = vtkSmartPointer<vtkGenericCell>::New();       
  // number of points in polydata to which data has been interpolated
//  int polyPnt = 0;
  for (int i = 0; i < nodeMesh->getNumberOfCells(); ++i)
  {
    interpolateToGaussPointsAtCell(i,genCell,das,daGausses);
  }
	for (int id = 0; id < arrayIDs.size(); ++id)
  {	
    gaussMesh->GetPointData()->AddArray(daGausses[id]);
  }
}

void GaussCubature::constructGaussMesh()
{
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
    switch(cellType)
    {
      case VTK_TRIANGLE:
        def->Initialize(VTK_TRIANGLE,3,3,TRI3,TRI3W);
        break;
      case VTK_TETRA:
        def->Initialize(VTK_TETRA, 4, 4, TET4,TET4W);
        break;
      default:
        std::cerr << "Error: Cell type: " << cellType << "found "
                  << "with no quadrature definition provided" << std::endl;
        exit(1);
    }
    // the definition must apear in the dictionary associated with 
    // the offset array
    key->Set(info, def, cellType); 
  }
  // get dictionary size 
  int dictSize = key->Size(info);
  dict = new vtkQuadratureSchemeDefinition * [dictSize];
  key->GetRange(info, dict, 0, 0, dictSize);
  offsets->SetNumberOfTuples(nodeMesh->getDataSet()->GetNumberOfCells());
  vtkIdType offset = 0;
  
  for (int cellid = 0; cellid < nodeMesh->getDataSet()->GetNumberOfCells() ; ++cellid)
  {
    offsets->SetValue(cellid, offset);
    vtkCell* cell = nodeMesh->getDataSet()->GetCell(cellid);
    int cellType = cell->GetCellType();
    vtkQuadratureSchemeDefinition* celldef = dict[cellType];
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




void GaussCubature::writeGaussMesh(const char* name)
{
  if (gaussMesh)
    writeVTFile<vtkXMLPolyDataWriter> (name,gaussMesh);   
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


