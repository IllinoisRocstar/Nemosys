#include <FETransfer.H>


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
  std::cout << "Transferring point array " << arrayID << std::endl;
  // extracting data array 
  std::vector<std::vector<double>> transferData;
  vtkDataArray* da = pd->GetArray(arrayID);
  // finding dim of data in array 
  int numComponent = da->GetNumberOfComponents();
  transferData.resize(target->getNumberOfPoints());
  // constructing cell locator for source mesh
  vtkSmartPointer<vtkCellLocator> cellLocator = 
    vtkSmartPointer<vtkCellLocator>::New();
  cellLocator->SetDataSet(source->getDataSet()); 
  cellLocator->BuildLocator();
	// gencell used by locator
	vtkSmartPointer<vtkGenericCell> gencell = vtkSmartPointer<vtkGenericCell>::New();				
	// passed to locator and evaulate position if called
	double pcoords[3];
	// creating buffer for weights array because we don't know how many weights to 
	// allocate before the cell is located
	double weights[10];
	// cell used when cell locator fails
	vtkCell* cell;
	// id of the cell containing source mesh point
  int id;
	// cellIds for when cell locator fails
	vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
  // parameters for interpolation
	double comps[numComponent];
	int pntId;
	int cellType;
	int numPointsInCell;
  int subId;
  double minDist2; // not used
	int result = 0;
  double x[3];
  std::vector<double> Weights; 
	for (int j = 0; j < target->getNumberOfPoints(); ++j)
  {
    transferData[j].resize(numComponent,0.0);
    // getting point from target and setting as query
    target->getDataSet()->GetPoint(j,x);
		// locating cell in source containing point in target
		id = cellLocator->FindCell(x, .5, gencell,pcoords,weights); 
		if (id >= 0)
		{
      for (int m = 0; m < gencell->GetNumberOfPoints(); ++m)
      {
        pntId = gencell->GetPointId(m);
        da->GetTuple(pntId, comps);
        for (int h = 0; h < numComponent; ++h)
        {
          transferData[j][h] += comps[h]*weights[m]; 
        }
      }
		
		}
		// if cell locator could not find cell containing target point
		else 
		{
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
    	  source->getDataSet()->GetPointCells(nnID, cellIds);
    	  // checking if query is inside of any of these neighbor cells
				// and calculating weights if so
        for (int i = 0; i < cellIds->GetNumberOfIds(); ++i)
				{
					id = cellIds->GetId(i);
    	  	cell = source->getDataSet()->GetCell(id); 
					cellType = cell->GetCellType();
    	  	numPointsInCell = cell->GetNumberOfPoints();
    	  	double weights[numPointsInCell];
          Weights.resize(numPointsInCell); 
    	    result = cell->EvaluatePosition(x,NULL,subId,pcoords, minDist2,weights);
					if (result == 1) 
          {
            for (int k = 0; k < numPointsInCell; ++k)
            {
              Weights[k] = weights[k];
            }
            break;
				  }
        }
			
				if (result == 0)
				{
					std::cout << "Could not locate point from target mesh in any cells sharing"
										<< " its nearest neighbor in the source mesh" << std::endl;
					exit(1);
				}
	
    	  else if (result == 1)
    	  {
    	    for (int m = 0; m < numPointsInCell; ++m)
    	    {
    	      pntId = cell->GetPointId(m);
    	      da->GetTuple(pntId, comps);
    	      
    	      for (int h = 0; h < numComponent; ++h)
    	      {
    	        transferData[j][h] += comps[h]*Weights[m]; 
    	      }
    	    }
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
  FETransfer::runPD(std::vector<std::vector<double>>& sourceData, 
                    std::vector<std::vector<double>>& targetPnts)
{
  // extracting data array 
  std::vector<std::vector<double>> transferData;
  // finding dim of data in array 
  int numComponent = sourceData[0].size();
  transferData.resize(targetPnts.size());
  // constructing cell locator for source mesh
  vtkSmartPointer<vtkCellLocator> cellLocator = 
    vtkSmartPointer<vtkCellLocator>::New();
  cellLocator->SetDataSet(source->getDataSet()); 
  cellLocator->BuildLocator();
	// gencell used by locator
	vtkSmartPointer<vtkGenericCell> gencell = vtkSmartPointer<vtkGenericCell>::New();				
	// passed to locator and evaulate position if called
	double pcoords[3];
	// creating buffer for weights array because we don't know how many weights to 
	// allocate before the cell is located
	double weights[10];
	// cell used when cell locator fails
	vtkCell* cell;
	// id of the cell containing source mesh point
  int id;
	// cellIds for when cell locator fails
	vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
  // parameters for interpolation
	double comps[numComponent];
	int pntId;
	int cellType;
	int numPointsInCell;
  int subId;
  double minDist2; // not used
	int result = 0;
  std::vector<double> Weights;
	for (int j = 0; j < targetPnts.size(); ++j)
  {
    transferData[j].resize(numComponent,0.0);
    // getting point from target and setting as query
    double* x = targetPnts[j].data(); 
		// locating cell in source containing point in target
		id = cellLocator->FindCell(x, .5, gencell,pcoords,weights); 
		if (id >= 0)
		{
      for (int m = 0; m < gencell->GetNumberOfPoints(); ++m)
      {
        pntId = gencell->GetPointId(m);
        for (int h = 0; h < numComponent; ++h)
        {
          transferData[j][h] += sourceData[pntId][h]*weights[m]; 
        }
      }
		
		}
		// if cell locator could not find cell containing target point
		else 
		{
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
    	  source->getDataSet()->GetPointCells(nnID, cellIds);
    	  // checking if query is inside of any of these neighbor cells
				// and calculating weights if so
				for (int i = 0; i < cellIds->GetNumberOfIds(); ++i)
				{
					id = cellIds->GetId(i);
    	  	cell = source->getDataSet()->GetCell(id); 
					cellType = cell->GetCellType();
    	  	numPointsInCell = cell->GetNumberOfPoints();
    	  	double weights[numPointsInCell]; 
          Weights.resize(numPointsInCell);
    	    result = cell->EvaluatePosition(x,NULL,subId,pcoords, minDist2,weights);
					if (result == 1)
          {
            for (int k = 0; k < numPointsInCell; ++k)
            {
              Weights[k] = weights[k];
            }
				    break;
          }
        }
			
				if (result == 0)
				{
					std::cout << "Could not locate point from target mesh in any cells sharing"
										<< " its nearest neighbor in the source mesh" << std::endl;
					exit(1);
				}
	
    	  else if (result == 1)
    	  {
    	    for (int m = 0; m < cell->GetNumberOfPoints(); ++m)
    	    {
    	      pntId = cell->GetPointId(m);
    	      for (int h = 0; h < numComponent; ++h)
    	      {
    	        transferData[j][h] += sourceData[pntId][h]*Weights[m]; 
    	      }
    	    }
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
  double comps[numComponent];
  // cellToPointData holds the cell data converted to points
  std::vector<std::vector<double>> cellToPointData;
  cellToPointData.resize(source->getNumberOfPoints());
  // cellId container for cells sharing a point
  vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New(); 
  int numSharedCells;
  double totW, W;
  for (int i = 0; i < source->getNumberOfPoints(); ++i)
  {
    cellToPointData[i].resize(numComponent,0);
    // find cells sharing point i
    source->getDataSet()->GetPointCells(i, cellIds);
    numSharedCells = cellIds->GetNumberOfIds(); 
    totW = 0;
    for (int j = 0; j < numSharedCells; ++j)
    {
      int cellId = cellIds->GetId(j);
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

  std::vector<std::vector<double>> transferData = runPD(cellToPointData, targetCenters);
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
