#include <FETransfer.H>


int FETransfer::runPD(vtkPointData* pd, int arrayID)
{
  // removing point data in target if it exists with same name
  {
    std::string srcArrName = pd->GetArrayName(arrayID);
    std::string trgArrName = 
      (
        target->getDataSet()->GetPointData()->GetArrayName(arrayID) ?
        target->getDataSet()->GetPointData()->GetArrayName(arrayID) : "NULL"
      );
    if (!srcArrName.compare(trgArrName))
      target->unsetPointDataArray(arrayID);
  }

  // extracting data array 
  std::vector<std::vector<double>> transferData;
  vtkDataArray* da = pd->GetArray(arrayID);
  // finding dim of data in array 
  int numComponent = da->GetNumberOfComponents();
  transferData.resize(target->getNumberOfPoints());

  for (int j = 0; j < target->getNumberOfPoints(); ++j)
  {
    transferData[j].resize(numComponent,0.0);
    // getting point from target and setting as query
    double x[3]; 
    target->getDataSet()->GetPoint(j,x);
    // searching for nearest neighbor of query in source mesh
    int nnID = source->getDataSet()->FindPoint(x);
    if (nnID < 0 )
    {
      std::cout << "Error finding neighbor point in source mesh for " 
                << " point " << j << "in target mesh" << std::endl;
      exit(1);
    }
    else
    {
      // searching for cells in source that contain nn of query
      vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
      source->getDataSet()->GetPointCells(nnID, cellIds);
      // getting first cell from list of cells containing nn of query
      int it = cellIds->GetId(0);
      vtkCell* cell = source->getDataSet()->GetCell(it); 
      int cellType = cell->GetCellType();
      int numPointsInCell = cell->GetNumberOfPoints();
      // declaring variables for use in interpolation function
      int subId;
      double minDist2; // not used
      double pcoords[3];
      double weights[numPointsInCell]; 
      int result; 
      switch(cellType)
      {
        /*case VTK_TRIANGLE:
        { 
          vtkTriangle* triCell = vtkTriangle::SafeDownCast(cell);
          vtkIdList* ids = triCell->GetPointIds();
          double x1[3],x2[3],x3[3];
          double tol2 = 0.000001;
          int in = triCell->PointInTriangle(x,x1,x2,x3,tol); 
          if (in)
          {
            // getting interpolation weights at query point in cell
            result = EvaluatePosition(x,NULL,subId,pcoords,minDist2,weights);
            if ( result == -1)  
            {
              std::cout << "problem encountered evaluating position of point from target"
                        << " mesh with respect to cell in source mesh" << std::endl;
              exit(1);
            }
            break;
          }
        }*/
        
        case VTK_TETRA:
        {
          vtkTetra* tetCell = vtkTetra::SafeDownCast(cell);
          vtkIdList* ids = tetCell->GetPointIds();
          double x1[3],x2[3],x3[3],x4[4];
          source->getDataSet()->GetPoint(ids->GetId(0),x1);
          source->getDataSet()->GetPoint(ids->GetId(1),x2);
          source->getDataSet()->GetPoint(ids->GetId(2),x3);
          source->getDataSet()->GetPoint(ids->GetId(3),x4);
          // getting interpolation weights at query point in cell
          result = tetCell->BarycentricCoords(x,x1,x2,x3,x4,weights);
          if ( result == 0)  
          {
            std::cout << "problem encountered evaluating barycentric coordinates. "
                      << "Tetrahedron is degenerate " << std::endl;
            exit(1);
          }
          break;
        }
        
        default:
        {
          // getting interpolation weights at query point in cell
          result = cell->EvaluatePosition(x,NULL,subId,pcoords, minDist2,weights);
          if ( result == -1)  
          {
            std::cout << "problem encountered evaluating position of point from target"
                      << " mesh with respect to cell in source mesh" << std::endl;
            exit(1);
          }
        }
      }
      if (result == 1)
      {
        for (int m = 0; m < cell->GetNumberOfPoints(); ++m)
        {
          int pntId = cell->GetPointId(m);
          double comps[numComponent];
          da->GetTuple(pntId, comps);
          
          for (int h = 0; h < numComponent; ++h)
          {
            transferData[j][h] += comps[h]*weights[m]; 
          }
        }

      }
    }
  }
  target->setPointDataArray(pd->GetArrayName(arrayID),transferData);
  return 0;
}  

int FETransfer::runPD(int arrayID)
{
  if (!(source && target))
  {
    std::cout << "source and target meshes must be initialized" << std::endl;
    exit(1);
  }

    vtkSmartPointer<vtkPointData> pd = vtkSmartPointer<vtkPointData>::New();
    pd = source->getDataSet()->GetPointData();
    if (pd)
    {
      int numArr = pd->GetNumberOfArrays();
      if (arrayID >= numArr)
      {
        std::cout << "ERROR: arrayID is out of bounds" << std::endl;
        std::cout << "There are " << numArr << " point data arrays" << std::endl;
        exit(1);
      } 
      runPD(pd, arrayID);
    }
    else
    {
      std::cout << "no point data found" << std::endl;
      return 1;
    }
  return 0;
}

int FETransfer::run()
{
  if (!(source && target))
  {
    std::cout << "source and target meshes must be initialized" << std::endl;
    exit(1);
  }

  // transferring point data
  {
    vtkSmartPointer<vtkPointData> pd = vtkSmartPointer<vtkPointData>::New();
    pd = source->getDataSet()->GetPointData();
    if (pd)
    {
      // finding number of data arrays
      int numArr = pd->GetNumberOfArrays();
    
      for (int i = 0; i < numArr; ++i)
      { 
        runPD(pd, i);
      }
    }
    else
    {
      std::cout << "no point data found" << std::endl;
      return 1;
    }
  }

  // transferring cell data
  /*{
    vtkSmartPointer<vtkCellData> cd = vtkSmartPointer<vtkCellData>::New();
    cd = source->getDataSet()->GetCellData();
    if (cd)
    {
      // finding number of data arrays
      int numArr = cd->GetNumberOfArrays();
      for (int i = 0; i < numArr; ++i)
      {
        runCD(cd, i);
      }
    }
  }*/

  return 0;
}
