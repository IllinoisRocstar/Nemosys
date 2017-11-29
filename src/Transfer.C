#include <Transfer.H>

Transfer* Transfer::Create(std::string method, meshBase* _source, meshBase* _target)
{
  if (!method.compare("FE"))
  {
    FETransfer* transobj = new FETransfer( _source , _target);
    return transobj; 
  }  
}

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
    std::vector<double> queryPt = target->getPoint(j); 
    double x[3]; 
    x[0] = queryPt[0]; x[1] = queryPt[1]; x[2] = queryPt[2];
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
      std::vector<int> commonCells = source->getCellsWithPoint(nnID);  
      auto it = commonCells.begin(); 
      int subId;
      double minDist2; // not used
      double pcoords[3];
      vtkCell* cell = source->getDataSet()->GetCell(*it); 
      int cellType = cell->GetCellType();
      int numPointsInCell = cell->GetNumberOfPoints();
      bool finished = 0;
      switch(cellType)
      {
        case VTK_TETRA:
        {
          vtkTetra* TetCell = vtkTetra::SafeDownCast(cell);
          vtkIdList* ids = TetCell->GetPointIds();
          double x1[3],x2[3],x3[3],x4[4];
          source->getDataSet()->GetPoint(ids->GetId(0),x1);
          source->getDataSet()->GetPoint(ids->GetId(1),x2);
          source->getDataSet()->GetPoint(ids->GetId(2),x3);
          source->getDataSet()->GetPoint(ids->GetId(3),x4);
          double baryCoord[4]; 
          TetCell->BarycentricCoords(x,x1,x2,x3,x4,baryCoord);
          for (int h = 0; h < numComponent; ++h)
          {
            for (int m = 0; m < 4; ++m)
            {
              int pntId = ids->GetId(m);//cell->GetPointId(m);
              double comps[numComponent];
              da->GetTuple(pntId, comps);
              transferData[j][h] += comps[h]*baryCoord[m];
            }
          }
          finished = 1;
          break;
        }
        default:
        {
          double weights[numPointsInCell];
          // getting interpolation weights at query point in cell
          int result = cell->
                          EvaluatePosition(x,NULL,subId,pcoords,minDist2,weights);
          if (result == 1)
          {
            for (int h = 0; h < numComponent; ++h)
            {
              for (int m = 0; m < cell->GetNumberOfPoints(); ++m)
              {
                int pntId = cell->GetPointId(m);
                double comps[numComponent];
                da->GetTuple(pntId, comps);
                transferData[j][h] += comps[h]*weights[m];
              }
            }
            finished = 1;
          }
          else if ( result == -1)  
          {
            std::cout << "problem encountered evaluating position of point from target"
                      << " mesh with respect to cell in source mesh" << std::endl;
            exit(1);
          }

        }

      }
    }
  }
  target->setPointDataArray(pd->GetArrayName(arrayID),transferData);
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
