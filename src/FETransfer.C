#include <FETransfer.H>

FETransfer::FETransfer(meshBase* _source, meshBase* _target)
{
  source = _source;
  srcCellLocator = source->buildLocator();
  target = _target;
  trgCellLocator = target->buildLocator();
  std::cout << "FETransfer constructed" << std::endl;
}

/* transfers point data with arrayID from source mesh to target
   The algorithm is as follows;
    1) For each point in the target mesh, find the cell of the source
       mesh in which it exists.
        - using a cell locator
        - if cell locator fails, find the nearest neighbor in the source mesh
          and all cells sharing this neighbor point. Check if the target point is 
          in any of these neighboring cells
    2) When the cell is identified, evaluate the weights for interpolation of the 
       solution to the target point and perform the interpolation.
*/

int FETransfer::runPD(vtkPointData* pd, int arrayID)
{
  std::cout << "Transferring point array " << pd->GetArrayName(arrayID) << std::endl;
  // extracting data array 
  std::vector<std::vector<double>> transferData;
  vtkDataArray* da = pd->GetArray(arrayID);
  // finding dim of data in array 
  int numComponent = da->GetNumberOfComponents();
  transferData.resize(target->getNumberOfPoints());
  // genCell used by locator
  vtkSmartPointer<vtkGenericCell> genCell = vtkSmartPointer<vtkGenericCell>::New();       
  // TODO REMOVE
  //std::ofstream outputStream1("pointInterpTest");

  for (int j = 0; j < target->getNumberOfPoints(); ++j)
  {
    transferData[j].resize(numComponent,0.0);
    // getting point from target and setting as query
    double x[3];
    target->getDataSet()->GetPoint(j,x);
    // id of the cell containing source mesh point
    vtkIdType id;
    int subId;
    double minDist2; 
    int result;
    // find closest point and closest cell to x
    double closestPoint[3];
    srcCellLocator->FindClosestPoint(x, closestPoint, genCell,id,subId,minDist2);
    if (id >= 0)
    {
      // passed evaulate position if called
      double pcoords[3];
      // parameters for interpolation
      double comps[numComponent];
      int pntId;
      double weights[genCell->GetNumberOfPoints()];
      result = genCell->EvaluatePosition(x,NULL,subId,pcoords,minDist2,weights); 

    //////////////////////////////////////////////////////
   /*   {
      std::vector<double> y(3);
      std::vector<double> xVec(x,x+3);
      for (int f = 0; f < genCell->GetNumberOfPoints(); ++f)
      { 
        pntId = genCell->GetPointId(f);
        std::vector<double> coords;
        coords = source->getPoint(pntId);
        for (int g = 0; g < 3; ++g)
        {
          y[g] += coords[g]*weights[f];
        }
      }
      std::vector<double> diff = y-xVec;
      outputStream1 << L2_Norm(diff) << " " << diff[0] << " " << diff[1] << " " << diff[2] << std::endl;
      }*/
    ///////////////////////////////////////////////////////
      if (result > 0)
      {
        for (int m = 0; m < genCell->GetNumberOfPoints(); ++m)
        {
          pntId = genCell->GetPointId(m);
          da->GetTuple(pntId, comps);
          for (int h = 0; h < numComponent; ++h)
          {
            transferData[j][h] += comps[h]*weights[m]; 
          }
        }
      }
      else if (result == 0)
      {
        std::cout << "Could not locate point from target mesh in any cells sharing"
                  << " its nearest neighbor in the source mesh" << std::endl;
        exit(1);
      }
      else   
      {
        std::cout << "problem encountered evaluating position of point from target"
                  << " mesh with respect to cell in source mesh" << std::endl;
        exit(1);
      }
    }
    else
    {
      std::cout << "Could not locate point from target in source mesh" << std::endl;
      exit(1);
    }
  }
  target->setPointDataArray(pd->GetArrayName(arrayID),transferData);
  if (checkQual)
  {
    std::vector<std::vector<double>> sourcePnts(source->getNumberOfPoints());
    std::vector<double> oldData(source->getNumberOfPoints()*numComponent);
    for (int i = 0; i < source->getNumberOfPoints(); ++i)
    {
      double comps[numComponent];
      sourcePnts[i] = source->getPoint(i);
      da->GetTuple(i, comps);
      for (int j = 0; j < numComponent; ++j)
      {
        oldData[i*numComponent + j] = comps[j];
      }
    }
    //std::string name = pd->GetArrayName(arrayID);
    //name += "backInterp";
    //std::vector<std::vector<double>> newData = runPD(trgCellLocator, transferData, sourcePnts);
    //source->setPointDataArray(&name[0u],newData);
    //std::vector<double> scaleSourcePnts = flatten(sourcePnts);
    //scaleSourcePnts = hadamard(scaleSourcePnts,scaleSourcePnts);
    //scaleSourcePnts = 298374.0*scaleSourcePnts;
    //std::vector<std::vector<double>> scaleSourcePntsFold = fold(scaleSourcePnts,3);
    //source->setPointDataArray("pointCoords", scaleSourcePntsFold);
    //std::vector<double> newData = flatten(runPD(target, transferData, sourcePnts));
    std::vector<double> newData1 = flatten(runPD(trgCellLocator, transferData, sourcePnts));//newData);
    double a,b,diff;
    //std::ofstream outputstream("transferQualityCheck.txt");
    //outputstream << "points with shit interpolation for " << pd->GetArrayName(arrayID) << std::endl
    //             << "Point" << " oldVal " << "newVal " << "diff " << std::endl;
    double sum = 0;
    for (int i = 0; i < newData1.size(); ++i)
    {
      a = newData1[i];
      b = oldData[i];
      diff = std::fabs((a-b)/b);
      sum += diff;
      //if (diff > 1e-4)
      //{
      //  outputstream << i << " " << b << " " << a << " " << diff << std::endl;
      //} 
    }
    std::cout << "Average Error in Nodal Transfer: " << sum/newData1.size() << std::endl;
  }
  return 0;
}  

int FETransfer::runPD(const std::vector<int>& arrayIDs)
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
    for (int i = 0; i < arrayIDs.size(); ++i)
    {
      if (arrayIDs[i] >= numArr)
      {
        std::cout << "ERROR: arrayID is out of bounds" << std::endl;
        std::cout << "There are " << numArr << " point data arrays" << std::endl;
        exit(1);
      }

      target->unsetPointDataArray(pd->GetArrayName(arrayIDs[i]));
       
      runPD(pd, arrayIDs[i]);
    }
  }
  else
  {
    std::cout << "no point data found" << std::endl;
    return 1;
  }
  return 0;
}

std::vector<std::vector<double>> 
  FETransfer::runPD(vtkSmartPointer<vtkCellLocator> cellLocator, 
                    std::vector<std::vector<double>>& sourceData, 
                    std::vector<std::vector<double>>& targetPnts)
{
  // extracting data array 
  std::vector<std::vector<double>> transferData;
  // finding dim of data in array 
  int numComponent = sourceData[0].size();
  transferData.resize(targetPnts.size());
  // genCell used by locator
  vtkSmartPointer<vtkGenericCell> genCell = vtkSmartPointer<vtkGenericCell>::New();       

  for (int j = 0; j < targetPnts.size(); ++j)
  {
    transferData[j].resize(numComponent,0.0);
    // getting point from target and setting as query
    // id of the cell containing source mesh point
    vtkIdType id;
    int subId;
    double minDist2;
    // find closest point and closest cell to x
    double closestPoint[3];
    double* x = targetPnts[j].data();
    cellLocator->FindClosestPoint(x, closestPoint, genCell,id,subId,minDist2);
    if (id >= 0)
    {
      // passed to evaulate position if called
      double pcoords[3];
      // parameters for interpolation
      double comps[numComponent];
      int pntId;
      int result;
      double weights[genCell->GetNumberOfPoints()];
      result = genCell->EvaluatePosition(x,NULL,subId,pcoords,minDist2,weights); 
      if (result > 0)
      {
        for (int m = 0; m < genCell->GetNumberOfPoints(); ++m)
        {
          pntId = genCell->GetPointId(m);
          for (int h = 0; h < numComponent; ++h)
          {
            transferData[j][h] += sourceData[pntId][h]*weights[m]; 
          }
        }
      }
      else if (result == 0)
      {
        std::cout << "Could not locate point from target mesh in any cells sharing"
                  << " its nearest neighbor in the source mesh" << std::endl;
        exit(1);
      }
      else   
      {
        std::cout << "problem encountered evaluating position of point from target"
                  << " mesh with respect to cell in source mesh" << std::endl;
        exit(1);
      }
    }
    else
    {
      std::cout << "Could not locate point " 
                << j << " from target in source mesh" << std::endl;
      exit(1);
    }
  }
  return transferData;
}

/* Transfer cell data from source mesh to target
   The algorithm is as follows:
    1)  Convert the cell data on the source mesh by inverse-distance 
        weighted averaging of data at cells sharing given point
          - cell data is assumed to be perscribed at cell centers
    2)  Compute the centers of cell in the target mesh
    3)  Transfer the converted cell-point data from the source mesh
        to the cell centers of the target mesh using the runPD methods
*/
int FETransfer::runCD(vtkCellData* cd, int arrayID, 
                      std::vector<std::vector<double>>& targetCenters)
{
 
  // ---------------------- Convert source cell data to point data -------- // 
  // extract cell data 
  vtkDataArray* da = cd->GetArray(arrayID);
  // finding dim of data in array 
  int numComponent = da->GetNumberOfComponents();

  // cellToPointData holds the cell data converted to points
  std::vector<std::vector<double>> cellToPointData;
  cellToPointData.resize(source->getNumberOfPoints());
  // cellId container for cells sharing a point
  vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New(); 
  
  for (int i = 0; i < source->getNumberOfPoints(); ++i)
  {
    cellToPointData[i].resize(numComponent,0);
    // find cells sharing point i
    source->getDataSet()->GetPointCells(i, cellIds);
    int numSharedCells = cellIds->GetNumberOfIds(); 
    double totW = 0;
    double W;
    for (int j = 0; j < numSharedCells; ++j)
    {
      int cellId = cellIds->GetId(j);
      double comps[numComponent];
      da->GetTuple(cellId, comps);
      // compute distance from point to cell center
      W = 1./L2_Norm(source->getCellCenter(cellId)
                     - source->getPoint(i));
      // average over shared cells, weighted by inverse distance to center
      for (int k = 0; k < numComponent; ++k)
      {
        cellToPointData[i][k] += W*comps[k];
      }
      totW += W; 
    } 
    cellToPointData[i] = (1.0/totW)*cellToPointData[i];
  }
  std::string arrName = cd->GetArrayName(arrayID);
  // ------------------- Transfer point data to cell centers of target ------// 

  std::vector<std::vector<double>> transferData = runPD(srcCellLocator, cellToPointData, targetCenters);
  target->setCellDataArray(&arrName[0u], transferData);
  return 0;
}

int FETransfer::runCD(const std::vector<int>& arrayIDs)
{

  if (!(source && target))
  {
    std::cout << "source and target meshes must be initialized" << std::endl;
    exit(1);
  }
  
  vtkSmartPointer<vtkCellData> cd = vtkSmartPointer<vtkCellData>::New();
  cd = source->getDataSet()->GetCellData();
  if (cd)
  {
    std::vector<std::vector<double>> targetCenters(target->getNumberOfCells());
    for (int j = 0; j < target->getNumberOfCells(); ++j)
    {
      targetCenters[j] = target->getCellCenter(j);
    }
    int numArr = cd->GetNumberOfArrays();
    for (int i = 0; i < arrayIDs.size(); ++i)
    {
      if (arrayIDs[i] >= numArr)
      {
        std::cout << "ERROR: arrayID is out of bounds" << std::endl;
        std::cout << "There are " << numArr << " cell data arrays" << std::endl;
        exit(1);
      }
      target->unsetCellDataArray(cd->GetArrayName(arrayIDs[i]));
      runCD(cd, arrayIDs[i], targetCenters);
    }
  }
  else
  {
    std::cout << "no cell data found" << std::endl;
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
  int numArr = source->getDataSet()->GetPointData()->GetNumberOfArrays();
  if (numArr > 0)
  {
    std::vector<int> arrayIDs(numArr);
    for (int i = 0; i < numArr; ++i)
    {
      arrayIDs[i] = i;
    }
    runPD(arrayIDs);
  }
  else
  {
    std::cout << "no point data found" << std::endl;
  }
  

  // transferring cell data
  numArr = source->getDataSet()->GetCellData()->GetNumberOfArrays();
  if (numArr > 0)
  {
    std::vector<int> arrayIDs(numArr);
    for (int i = 0; i < numArr; ++i)
    {
      arrayIDs[i] = i;
    }
    runCD(arrayIDs);
  }
  else
  {
    std::cout << "no cell data found" << std::endl;
  }

  return 0;
}

