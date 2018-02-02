#include <patchRecovery.H>

PatchRecovery::PatchRecovery(meshBase* nodeMesh, const std::vector<int>& arrayIDs)
{
  cubature = GaussCubature::CreateUnique(nodeMesh, arrayIDs);
    
}

void PatchRecovery::recoverNodalSolution(const std::vector<int>& arrayIDs)
{
  std::cout << "WARNING: mesh is assumed to be properly numbered" << std::endl;
  meshBase* nodeMesh = cubature->getNodeMesh();
  int numPoints = nodeMesh->getNumberOfPoints();
  vtkSmartPointer<vtkIdList> patchCellIDs = vtkSmartPointer<vtkIdList>::New();
  for (int i = 0; i < numPoints; ++i)
  {
    
    // get ids of cells in patch of node
    nodeMesh->getDataSet()->GetPointCells(i,patchCellIDs);
    for (int j = 0; j < patchCellIDs->GetNumberOfIds(); ++j)
    {
      pntDataPairVec pntsAndData = cubature->getGaussPointsAndDataAtCell(patchCellIDs->GetId(j));
      // TODO: efficient way to extract x,y,z in std::vec and data into VectorXd
      //       note move semantics can be used because we don't care about 
      //       pntDataPairVec once we've got the info   

    } 


  } 

}
