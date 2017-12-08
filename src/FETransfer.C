#include <FETransfer.H>
// TODO: Figure out why all data is being transferred regardless

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
  // constructing cell locator for source mesh
  vtkSmartPointer<vtkCellLocator> cellLocator = 
    vtkSmartPointer<vtkCellLocator>::New();
  cellLocator->SetDataSet(source->getDataSet()); 
  cellLocator->BuildLocator();
  double scale_len;
  {
    std::vector<double> meanStdev = getMeanStdev(source->getCellLengths());
    scale_len = meanStdev[0]*.5;
  }
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
	int result;
	for (int j = 0; j < target->getNumberOfPoints(); ++j)
  {
    transferData[j].resize(numComponent,0.0);
    // getting point from target and setting as query
    double x[3]; 
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
    	    result = cell->EvaluatePosition(x,NULL,subId,pcoords, minDist2,weights);
					if (result == 1) break;
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
    	      da->GetTuple(pntId, comps);
    	      
    	      for (int h = 0; h < numComponent; ++h)
    	      {
    	        transferData[j][h] += comps[h]*weights[m]; 
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
