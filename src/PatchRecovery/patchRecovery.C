#include <patchRecovery.H>


PatchRecovery::PatchRecovery(meshBase* nodeMesh, int _order, const std::vector<int>& arrayIDs)
	: order(_order)
{
  cubature = GaussCubature::CreateUnique(nodeMesh, arrayIDs);
    
}

//TODO: stop being an idiot and change the type of pntDataPairVec to hold eigen vec type for data

// extract coordinates and data from cell by pushing into coords and data
// vectors accross patch
void PatchRecovery::extractAxesAndData(pntDataPairVec& pntsAndData, 
																       std::vector<std::vector<double>>& coords,
                                       std::vector<VectorXd>& data,
                                       const std::vector<int>& numComponents,
                                       int& pntNum)
{

	// number of components for each data being considered
	// {x,y,z} -> {{d1_1, ... ,d1_m} , ... , {dk_1, ... , dk_m}} ; i_m != k_m in general

	for (int i = 0; i < pntsAndData.size(); ++i)
	{
		//coords.push_back(std::move(pntsAndData[i].first));
    coords[pntNum] = std::move(pntsAndData[i].first);
    printVec(coords[pntNum]);
    int currcomp = 0;
    for (int j = 0; j < numComponents.size(); ++j)
    {
      for (int k = 0; k < numComponents[j]; ++k)
      {
        data[currcomp](pntNum) = pntsAndData[i].second[currcomp];
        ++currcomp;
      }
    }
	  pntNum += 1;
  }
}
													


// pair type for coordinate and data there
//typedef std::pair<std::vector<double>, std::vector<double>> pntDataPair;

// holds gauss points and data at these points as pairs
//typedef std::vector<pntDataPair> pntDataPairVec;

void PatchRecovery::recoverNodalSolution()
{
  std::cout << "WARNING: mesh is assumed to be properly numbered" << std::endl;
  meshBase* nodeMesh = cubature->getNodeMesh();
  int numPoints = nodeMesh->getNumberOfPoints();
  vtkSmartPointer<vtkIdList> patchCellIDs = vtkSmartPointer<vtkIdList>::New();
  vtkQuadratureSchemeDefinition** dict = cubature->getDict();
	
	std::vector<int> numComponents = cubature->getNumComponents();
  int totalComponents = cubature->getTotalComponents();

	for (int i = 0; i < numPoints; ++i)
  {
   
    double tmp[1];
    nodeMesh->getDataSet()->GetPointData()->GetArray(1)->GetTuple(i,tmp);
    // get ids of cells in patch of node
    nodeMesh->getDataSet()->GetPointCells(i,patchCellIDs);
    // get total number of gauss points in patch 
    int numPatchPoints = 0;
    for (int k = 0; k < patchCellIDs->GetNumberOfIds(); ++k)
		{
      int cellType = nodeMesh->getDataSet()->GetCell(patchCellIDs->GetId(k))->GetCellType();
      numPatchPoints += dict[cellType]->GetNumberOfQuadraturePoints(); 
    }
    // coordinates of each gauss point in patch
    std::vector<std::vector<double>> coords(numPatchPoints);
		// vector of components of data at all gauss points in patch
    std::vector<VectorXd> data(totalComponents);
    for (int k = 0; k < totalComponents; ++k)
    {
      data[k].resize(numPatchPoints);
    }

    int pntNum = 0;
    for (int j = 0; j < patchCellIDs->GetNumberOfIds(); ++j)
    {
      pntDataPairVec pntsAndData = cubature->getGaussPointsAndDataAtCell(patchCellIDs->GetId(j));
			extractAxesAndData(pntsAndData, coords, data, numComponents, pntNum);

    }

		// construct orthopoly from coords
		std::unique_ptr<orthoPoly3D> patchPolyApprox 
      = orthoPoly3D::CreateUnique(order,coords);//(new orthoPoly3D(order, coords));
		// get coordinate of node generating patch
    std::vector<double> genNodeCoord = nodeMesh->getPoint(i);
		for (int k = 0; k < totalComponents; ++k)
		{
      // get polynomial approximant of component over patch
      patchPolyApprox->computeA(data[k]);
      std::cout << patchPolyApprox->eval(genNodeCoord) << std::endl;
    }

  } 
}

		//int numGaussPoints = 0;
		//for (int j = 0; j < patchCellIDs->GetNumberOfIds(); ++j)
		//{
		//	int cellType = nodeMesh->getDataSet()->GetCell(patchCellIDs->GetId(j))->GetCellType();
		//	numGaussPoints += dict[cellType]->GetNUmberOfQuadraturePoints();
		//}		

