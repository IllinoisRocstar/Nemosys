#include <Cubature.H>


GaussCubature::GaussCubature(meshBase* _nodeMesh)
{
  nodeMesh = _nodeMesh;
  nodeMesh->unsetCellDataArray("QuadratureOffset");
	constructGaussMesh(); 
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
  vtkSmartPointer<vtkDataArray> data 
    = nodeMesh->getDataSet()->GetCellData()->GetArray(basename.c_str());

  offsets->SetName(basename.c_str());
  nodeMesh->getDataSet()->GetCellData()->AddArray(offsets);
  
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
        def->Initialize(VTK_TRIANGLE,3,3,TRI3);
        break;
      case VTK_TETRA:
        def->Initialize(VTK_TETRA, 4, 4, TET4);
        break;
      default:
        std::cout << "Error: Cell type: " << cellType << "found "
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

void GaussCubature::interpolateToGaussPoints(const std::vector<int>& arrayIDs)
{

  vtkSmartPointer<vtkGenericCell> genCell = vtkSmartPointer<vtkGenericCell>::New();       
  for (int id = 0; id < arrayIDs.size(); ++id)
  {
    // get desired point data array to be interpolated to gauss points
    vtkSmartPointer<vtkDataArray> da = nodeMesh->getDataSet()->GetPointData()->GetArray(arrayIDs[id]);
    // get tuple length of given data
    int numComponent = da->GetNumberOfComponents();
    // declare data array to be populated for polydata
	  vtkSmartPointer<vtkDoubleArray> daGauss = vtkSmartPointer<vtkDoubleArray>::New();
    // names and sizing
    daGauss->SetName(nodeMesh->getDataSet()->GetPointData()->GetArrayName(arrayIDs[id]));
    daGauss->SetNumberOfComponents(numComponent);
    // generic cell to store given cell in nodeMesh->getDataSet()

    // number of points in polydata to which data has been interpolated
    int polyPnt = 0;
    // allocate storate for polygons
    gaussMesh->Allocate();
    for (int i = 0; i < nodeMesh->getNumberOfCells(); ++i)
    {
      // putting current cell into genCell
      nodeMesh->getDataSet()->GetCell(i, genCell);
      // getting cellType information for lookup in dictionary
      int cellType = nodeMesh->getDataSet()->GetCellType(i);
      // number of gauss points in this cell
      int numGaussPoints = dict[cellType]->GetNumberOfQuadraturePoints();
      // id list for 'gauss' cells in polydata 
      vtkSmartPointer<vtkIdList> polyCellIds = vtkSmartPointer<vtkIdList>::New();
      for (int j = 0; j < numGaussPoints; ++j)
      {
        polyCellIds->InsertNextId(j+polyPnt); 
        double x[3];
        gaussMesh->GetPoint(j+polyPnt,x); 
        // parameters for interpolation
        int subId; // not used
        double minDist2; // not used
        double pcoords[3];
        double comps[numComponent];
        double weights[genCell->GetNumberOfPoints()];
        genCell->EvaluatePosition(x,NULL,subId,pcoords,minDist2,weights); 
        double interps[numComponent];
        for (int m = 0; m < genCell->GetNumberOfPoints(); ++m)
        {
          int pntId = genCell->GetPointId(m);
          da->GetTuple(pntId, comps);
          for (int h = 0; h < numComponent; ++h)
          {
            interps[h] += comps[h]*weights[m]; 
          }
        }
        daGauss->InsertNextTuple(interps);
      }
      gaussMesh->InsertNextCell(VTK_POLYGON,polyCellIds); 
      polyPnt += numGaussPoints;
    }
    std::cout << "NUM POLY POINTS: " << polyPnt << std::endl; 
    gaussMesh->GetPointData()
      ->SetActiveScalars(nodeMesh->getDataSet()->GetPointData()->GetArrayName(arrayIDs[id]));
    gaussMesh->GetPointData()->SetScalars(daGauss);    
  }      
}
