#include <patchRecovery.H>


PatchRecovery::PatchRecovery(meshBase* nodeMesh, int _order, const std::vector<int>& arrayIDs)
	: order(_order)
{
  cubature = GaussCubature::CreateUnique(nodeMesh, arrayIDs);
}



// extract coordinates and data from cell by pushing into coords and data
// vectors accross patch
void PatchRecovery::extractAxesAndData(pntDataPairVec& pntsAndData, 
																       std::vector<std::vector<double>>& coords,
                                       std::vector<VectorXd>& data,
                                       const std::vector<int>& numComponents,
                                       int& pntNum)
{
	for (int i = 0; i < pntsAndData.size(); ++i)
	{
    coords[pntNum] = std::move(pntsAndData[i].first);
   // printVec(coords[pntNum]);
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
  // getting node mesh from cubature
  meshBase* nodeMesh = cubature->getNodeMesh();
  int numPoints = nodeMesh->getNumberOfPoints();
  // initializing id list for patch cells
  vtkSmartPointer<vtkIdList> patchCellIDs = vtkSmartPointer<vtkIdList>::New();
  // getting cubature scheme dictionary for indexing
  vtkQuadratureSchemeDefinition** dict = cubature->getDict();
	std::vector<int> numComponents = cubature->getNumComponents();
  // getting number of doubles representing all data at point
  int totalComponents = cubature->getTotalComponents();
  // storage for new point data
  std::vector<vtkSmartPointer<vtkDoubleArray>> newPntData(numComponents.size());
 
  // initializing double arrays for new data 
  for (int i = 0; i < numComponents.size(); ++i)
  {
    std::string name = nodeMesh->getDataSet()->GetPointData()
                               ->GetArrayName(cubature->getArrayIDs()[i]);
    name += "New";
    newPntData[i] = vtkSmartPointer<vtkDoubleArray>::New();
    newPntData[i]->SetName(&name[0u]);
    newPntData[i]->SetNumberOfComponents(numComponents[i]);  
    newPntData[i]->SetNumberOfTuples(numPoints); 
  }

  // looping over all points, looping over patches per point
	for (int i = 0; i < numPoints; ++i)
  {
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

		// construct polyapprox from coords
		std::unique_ptr<polyApprox> patchPolyApprox 
      = polyApprox::CreateUnique(order,std::move(coords));//(new orthoPoly3D(order, coords));
		// get coordinate of node that generates patch
    std::vector<double> genNodeCoord = nodeMesh->getPoint(i);
    int currComp = 0;
    for (int k = 0; k < numComponents.size(); ++k)
    {
      double comps[numComponents[k]];
      for (int l = 0; l < numComponents[k]; ++l)
      {
        patchPolyApprox->computeCoeff(data[currComp++]);
        comps[l] = patchPolyApprox->eval(genNodeCoord);
      }
      newPntData[k]->InsertTuple(i,comps);
    }
  }
  for (int k = 0; k < numComponents.size(); ++k)
  {
    nodeMesh->getDataSet()->GetPointData()->AddArray(newPntData[k]);
  }
}

		//for (int k = 0; k < totalComponents; ++k)
		//{

    //  // get polynomial approximant of component over patch
    //  patchPolyApprox->computeCoeff(data[k]);
    //  for (; currArr < numComponents[currArr]; ++currArr)
    //  newPntData[currArr]->InsertNextTuple(patchPolyApprox->eval(genNodeCoord));
    //    
    //  //std::cout << patchPolyApprox->eval(genNodeCoord) << std::endl;
    //}
		
    //int numGaussPoints = 0;
		//for (int j = 0; j < patchCellIDs->GetNumberOfIds(); ++j)
		//{
		//	int cellType = nodeMesh->getDataSet()->GetCell(patchCellIDs->GetId(j))->GetCellType();
		//	numGaussPoints += dict[cellType]->GetNUmberOfQuadraturePoints();
		//}		

