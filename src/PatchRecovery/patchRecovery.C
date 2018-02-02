#include <patchRecovery.H>


PatchRecovery::PatchRecovery(meshBase* nodeMesh, int _order, const std::vector<int>& arrayIDs)
	: order(_order)
{
  cubature = GaussCubature::CreateUnique(nodeMesh, arrayIDs);
    
}

//TODO: stop being an idiot and change the type of pntDataPairVec to hold eigen vec type for data

// extract coordinates and data from cell by pushing into coords and data
void PatchRecovery::extractAxes(pntDataPairVec& pntsAndData, 
																std::vector<std::vector<double>>& coords)
{

	// number of components for each data being considered
	// {x,y,z} -> {{d1_1, ... ,d1_m} , ... , {dk_1, ... , dk_m}} ; i_m != k_m in general
	std::vector<int> numComponents = cubature->getNumComponents();			

	for (int i = 0; i < pntsAndData.size(); ++i)
	{
		coords.push_back(std::move(pntsAndData[i].first));
	}
}
													
// pair type for coordinate and data there
//typedef std::pair<std::vector<double>,std::vector<std::vector<double>>> pntDataPair;

// holds gauss points and data at these points as pairs
//typedef std::vector<pntDataPair> pntDataPairVec;


//void PatchRecovery::recoverNodalSolution(const std::vector<int>& arrayIDs)
//{
//  std::cout << "WARNING: mesh is assumed to be properly numbered" << std::endl;
//  meshBase* nodeMesh = cubature->getNodeMesh();
//  int numPoints = nodeMesh->getNumberOfPoints();
//  vtkSmartPointer<vtkIdList> patchCellIDs = vtkSmartPointer<vtkIdList>::New();
//  vtkQuadratureSchemeDefinition** dict = cubature->getDict();
//	
//	std::vector<int> numComponents = cubature->getNumComponents();
//
//	for (int i = 0; i < numPoints; ++i)
//  {
//    
//		std::vector<std::vector<double>> coords;
//		std::vector<std::vector<VectorXd>> data;
//
//    // get ids of cells in patch of node
//    nodeMesh->getDataSet()->GetPointCells(i,patchCellIDs);
//		int numGaussPoints = 0;
//		for (int j = 0; j < patchCellIDs->GetNumberOfIds(); ++j)
//		{
//			int cellType = nodeMesh->getDataSet()->GetCell(patchCellIDs->GetId(j))->GetCellType();
//			numGaussPoints += dict[cellType]->GetNUmberOfQuadraturePoints();
//		}		
//
//    for (int j = 0; j < patchCellIDs->GetNumberOfIds(); ++j)
//    {
//		
//      pntDataPairVec pntsAndData = cubature->getGaussPointsAndDataAtCell(patchCellIDs->GetId(j));
//			extractAxesAndData(pntsAndData, coords, data);
//    }
//
//		//TODO: SOMETHING LIKE THIS
//		// construct orthopoly from coords
//		std::unique_ptr<orthoPoly3D> patchPolyApprox(new orthoPoly3D(order, coords));
//		for (int k = 0; k < numComponents.size(); ++k)
//		{	
//			for (int l = 0; l < numComponents[k].size(); ++k)
//			{
//				patchPolyApprox->computeA(data[k][l]);	
//			
//			}
//	
//		}
//
//  } 
//
//}
