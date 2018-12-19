#include <Z2ErrorSizeField.H>
#include <patchRecovery.H>
#include <vtkCellData.h>

Z2ErrorSizeField::Z2ErrorSizeField(meshBase* _mesh, int arrayID)
{
  initialize(_mesh, arrayID, 0,0, "Z2ErrorSF");
  order = mesh->getOrder();
  std::cout << "Z2ErrorSizeField constructed" << std::endl; 
}

double Z2ErrorSizeField::computeNodalError(int arrayID)
{
  std::vector<int> arrayIDs(1);
  arrayIDs[0] = arrayID;
  std::unique_ptr<PatchRecovery> recoverObj
    = std::unique_ptr<PatchRecovery> (new PatchRecovery(mesh,order,arrayIDs));
  return recoverObj->computeNodalError()[0][0];
}

void Z2ErrorSizeField::computeSizeField(int arrayID)
{
  double aveError = computeNodalError(arrayID)/mesh->getNumberOfCells();

  std::string errorName = mesh->getDataSet()->GetPointData()->GetArrayName(arrayID);
  errorName += "ErrorIntegral";
  
  vtkSmartPointer<vtkDataArray> errorIntegrals 
    = mesh->getDataSet()->GetCellData()->GetArray(&errorName[0u]);
  vtkSmartPointer<vtkDataArray> nodeSizes
    = mesh->getDataSet()->GetPointData()->GetArray("nodeSizes");

  vtkSmartPointer<vtkGenericCell> genCell = vtkSmartPointer<vtkGenericCell>::New();

  std::vector<double> sizeField(mesh->getNumberOfCells());

  for (int i = 0; i < mesh->getNumberOfCells(); ++i)
  {
    double interpSize = 0; 
    mesh->getDataSet()->GetCell(i,genCell);
    double center[3];
    double weights[genCell->GetNumberOfPoints()];
    int subId;    // dummy
    double x[3];  // dummy
    genCell->GetParametricCenter(center);
    genCell->EvaluateLocation(subId,center,x,weights);
    for (int j = 0; j < genCell->GetNumberOfPoints(); ++j)
    {
      int pntID = genCell->GetPointId(j);
      double size[1];
      nodeSizes->GetTuple(pntID,size);
      interpSize += size[0]*weights[j];
    }
    double error[1];
    errorIntegrals->GetTuple(i,error);
    sizeField[i] = interpSize/(std::pow(error[0]/aveError,1.0/order))*sizeFactor;
  }
  mesh->setCellDataArray(&sfname[0u], sizeField);
  mesh->setSFBool(1);
}
