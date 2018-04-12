#include <patchRecovery.H>

//TODO: To define orthogonal polynomials over patches of a structured grid that has
//      deformed from a rectilinear grid, a conformal mapping must be applied to transform
//      the deformed mesh back to a rectilinear grid. Only then can we satisfy the 
//      requirements for a tensor-product construction of multivariate orthogonal polynomials.

PatchRecovery::PatchRecovery(meshBase* nodeMesh, int _order, const std::vector<int>& arrayIDs)
  : order(_order)
{
  cubature = GaussCubature::CreateUnique(nodeMesh, arrayIDs);
}


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
                          
std::vector<double> 
PatchRecovery::getMinMaxCoords(const std::vector<std::vector<double>>& coords)
{
  std::vector<double> minMaxXYZ(6);
  for (int i = 0; i < 3; ++i)
  {
    minMaxXYZ[2*i] = coords[0][i];
    minMaxXYZ[2*i + 1] = minMaxXYZ[2*i];
  }

  for (int i = 0; i < coords.size(); ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      // min
      if (minMaxXYZ[2*j] > coords[i][j])
      {
        minMaxXYZ[2*j] = coords[i][j];
      }
      // max
      else if (minMaxXYZ[2*j+1] < coords[i][j])
      {
        minMaxXYZ[2*j+1] = coords[i][j];
      }
    }  
  }  
  return minMaxXYZ;
}

void PatchRecovery::regularizeCoords(std::vector<std::vector<double>>& coords,
                                     std::vector<double>& genNodeCoord)
{ 
  std::vector<double> minMaxXYZ = getMinMaxCoords(coords);
  for (int j = 0; j < 3; ++j)
  {
    double coord_min = minMaxXYZ[2*j];
    double bound = minMaxXYZ[2*j+1] - coord_min;
    for (int i = 0; i < coords.size(); ++i)
    {
      coords[i][j] = -1 + 2*(coords[i][j] - coord_min)/bound;
    }
    genNodeCoord[j] = -1 + 2*(genNodeCoord[j] - coord_min)/bound;
  }
}

void PatchRecovery::recoverNodalSolution(bool ortho)
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

  int totPatchPoints = 0; 
  // looping over all points, looping over patches per point
  for (int i = 0; i < numPoints; ++i) //FIXME
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

    if (patchCellIDs->GetNumberOfIds() < 2)
    {
      std::cerr << "Only " << patchCellIDs->GetNumberOfIds() 
                << " cell in patch of point " << i << std::endl;
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

    // get coordinate of node that generates patch
    std::vector<double> genNodeCoord = nodeMesh->getPoint(i);
    // regularizing coordinates for preconditioning of basis matrix
    regularizeCoords(coords, genNodeCoord);
    if (!ortho)
    {
      // construct polyapprox from coords
      std::unique_ptr<polyApprox> patchPolyApprox 
        = polyApprox::CreateUnique(order,std::move(coords));//(new orthoPoly3D(order, coords));
      int currComp = 0;
      for (int k = 0; k < numComponents.size(); ++k)
      {
        double comps[numComponents[k]];
        for (int l = 0; l < numComponents[k]; ++l)
        {
          patchPolyApprox->computeCoeff(data[currComp]);
          comps[l] = patchPolyApprox->eval(genNodeCoord);
          patchPolyApprox->resetCoeff();
          ++currComp;
        }
        newPntData[k]->InsertTuple(i,comps);
      }
    }
    else
    {
      for (int k = 0; k < patchCellIDs->GetNumberOfIds(); ++k)
      {
        std::cout << "point " << i << " patch cell: " << patchCellIDs->GetId(k) << std::endl;
      }

      std::unique_ptr<orthoPoly3D> patchPolyApprox
        = orthoPoly3D::CreateUnique(order,std::move(coords)); 
      int currComp = 0;
      for (int k = 0; k < numComponents.size(); ++k)
      {
        double comps[numComponents[k]];
        for (int l = 0; l < numComponents[k]; ++l)
        {
          std::cout << "computing coefficients of ortho polys" << std::endl;
          patchPolyApprox->computeA(data[currComp]);
          comps[l] = patchPolyApprox->eval(genNodeCoord);
          patchPolyApprox->resetA();
          ++currComp;
        }
        newPntData[k]->InsertTuple(i,comps);
      }
    }
  }
  for (int k = 0; k < numComponents.size(); ++k)
  {
    nodeMesh->getDataSet()->GetPointData()->AddArray(newPntData[k]);
  }
}

// TODO: check for whether recovered solution exists
std::vector<std::vector<double>> PatchRecovery::computeNodalError()
{
  std::cout << "WARNING: mesh is assumed to be properly numbered" << std::endl;
  // getting node mesh from cubature
  meshBase* nodeMesh = cubature->getNodeMesh();
  std::vector<int> arrayIDs = cubature->getArrayIDs();
  int numPoints = nodeMesh->getNumberOfPoints();
  // initializing id list for patch cells
  vtkSmartPointer<vtkIdList> patchCellIDs = vtkSmartPointer<vtkIdList>::New();
  // getting cubature scheme dictionary for indexing
  vtkQuadratureSchemeDefinition** dict = cubature->getDict();
  std::vector<int> numComponents = cubature->getNumComponents();
  // getting number of doubles representing all data at a given point
  int totalComponents = cubature->getTotalComponents();
  // storage for new point data
  std::vector<vtkSmartPointer<vtkDoubleArray>> newPntData(numComponents.size());
  // storage for error data
  std::vector<vtkSmartPointer<vtkDoubleArray>> errorPntData(numComponents.size());
  // reference point data
  vtkSmartPointer<vtkPointData> pd = nodeMesh->getDataSet()->GetPointData(); 

  std::vector<std::string> errorNames(numComponents.size());
 
  // initializing double arrays for new data
  for (int i = 0; i < numComponents.size(); ++i)
  {
    std::string name = pd->GetArrayName(arrayIDs[i]);
    std::string newName = name + "New";
    newPntData[i] = vtkSmartPointer<vtkDoubleArray>::New();
    newPntData[i]->SetName(&newName[0u]);
    newPntData[i]->SetNumberOfComponents(numComponents[i]);  
    newPntData[i]->SetNumberOfTuples(numPoints); 
  
    std::string errName = name + "Error";
    errorPntData[i] = vtkSmartPointer<vtkDoubleArray>::New();
    errorPntData[i]->SetName(&errName[0u]);
    errorPntData[i]->SetNumberOfComponents(numComponents[i]);  
    errorPntData[i]->SetNumberOfTuples(numPoints);
    errorNames[i] = errName; 
  }
 
  // storage for nodal average element sizes
  vtkSmartPointer<vtkDoubleArray> nodeSizes = vtkSmartPointer<vtkDoubleArray>::New();
  nodeSizes->SetName("nodeSizes");
  nodeSizes->SetNumberOfComponents(1);
  nodeSizes->SetNumberOfTuples(numPoints);

  // storage for generic cell if needed
  vtkSmartPointer<vtkGenericCell> genCell = vtkSmartPointer<vtkGenericCell>::New(); 

  // looping over all points, looping over patches per point
  for (int i = 0; i < numPoints; ++i)
  {
    // get ids of cells in patch of node
    nodeMesh->getDataSet()->GetPointCells(i,patchCellIDs);
    // get total number of gauss points in patch and assign element size to 
    // patch generating node 
    int numPatchPoints = 0;
    // also get average size of elements in patch
    double nodeSize = 0;
    for (int k = 0; k < patchCellIDs->GetNumberOfIds(); ++k)
    {
      // put current patch cell into gencell
      nodeMesh->getDataSet()->GetCell(patchCellIDs->GetId(k),genCell);
      int cellType = nodeMesh->getDataSet()->GetCell(patchCellIDs->GetId(k))->GetCellType();
      numPatchPoints += dict[cellType]->GetNumberOfQuadraturePoints();
      //nodeSize += cbrt(2.356194490192344*cubature->computeCellVolume(genCell, cellType));   
      nodeSize += std::sqrt(genCell->GetLength2());
    }
    // patch-averaged node size
    nodeSize /= patchCellIDs->GetNumberOfIds();
    nodeSizes->InsertTuple(i, &nodeSize);
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

    // get coordinate of node that generates patch
    std::vector<double> genNodeCoord = nodeMesh->getPoint(i);
    // regularizing coordinates for preconditioning of basis matrix
    regularizeCoords(coords, genNodeCoord);
    // construct polyapprox from coords
    std::unique_ptr<polyApprox> patchPolyApprox 
      = polyApprox::CreateUnique(order,std::move(coords));//(new orthoPoly3D(order, coords));
    int currComp = 0;
    for (int k = 0; k < numComponents.size(); ++k)
    {
      double comps[numComponents[k]];
      double refComps[numComponents[k]];
      double errorComps[numComponents[k]];
      pd->GetArray(arrayIDs[k])->GetTuple(i,refComps);
      for (int l = 0; l < numComponents[k]; ++l)
      {
        patchPolyApprox->computeCoeff(data[currComp]);
        comps[l] = patchPolyApprox->eval(genNodeCoord);
        errorComps[l] = std::pow(comps[l] - refComps[l],2);  
        patchPolyApprox->resetCoeff();
        ++currComp;
      }
      newPntData[k]->InsertTuple(i,comps);
      errorPntData[k]->InsertTuple(i,errorComps);
    }
  }
  for (int k = 0; k < numComponents.size(); ++k)
  {
//    pd->AddArray(newPntData[k]);
    pd->AddArray(errorPntData[k]);
  }
 
  pd->AddArray(nodeSizes); 
  return cubature->integrateOverAllCells(errorNames, 1);  
}
