#include <meshBase.H>
#include <stdlib.h>

int main(int argc, char* argv[])
{
  
  if (argc < 4) 
  {
    std::cout << "Usage: " << argv[0] << " file1.vtu file2.vtu TOL" << std::endl;
    exit(1);
  }

  meshBase* mesh1 = meshBase::Create(argv[1]);
  meshBase* mesh2 = meshBase::Create(argv[2]);
  double tol = atof(argv[3]);


  if (mesh1->getNumberOfPoints() != mesh2->getNumberOfPoints() ||
      mesh1->getNumberOfCells() != mesh2->getNumberOfCells())
  {
    std::cout << "Meshes don't have the same number of points or cells" << std::endl;
    exit(1);
  }

  for (int i = 0; i < mesh1->getNumberOfPoints(); ++i)
  {
    std::vector<double> coord1 = mesh1->getPoint(i);
    std::vector<double> coord2 = mesh2->getPoint(i);
    for (int j = 0; j < 3; ++j)
    {
      if (std::fabs(coord1[j]-coord2[j]) > tol)
      {
        std::cout << "Meshes differ in point coordinates" << std::endl;
        exit(1);
      }
    } 
  }

  for (int i = 0; i < mesh1->getNumberOfCells(); ++i)
  {
    std::vector<std::vector<double>> cell1 = mesh1->getCellVec(i);
    std::vector<std::vector<double>> cell2 = mesh2->getCellVec(i);
    if (cell1.size() != cell2.size())
    {
      std::cout << "Meshes differ in cells" << std::endl;
      exit(1); 
    }
    for (int j = 0; j < cell1.size(); ++j)
    {
      for (int k = 0; k < 3; ++k)
      {
        if (std::fabs(cell1[j][k] - cell2[j][k]) > tol)
        {
          std::cout << "Meshes differ in cells" << std::endl;
          exit(1); 
        }
      }
    }   
  }

  vtkSmartPointer<vtkPointData> pd1 = vtkSmartPointer<vtkPointData>::New();
  pd1 = mesh1->getDataSet()->GetPointData();
  vtkSmartPointer<vtkPointData> pd2 = vtkSmartPointer<vtkPointData>::New();
  pd2 = mesh2->getDataSet()->GetPointData(); 
  int numArr1 = pd1->GetNumberOfArrays(); 
  int numArr2 = pd2->GetNumberOfArrays(); 
 
  if (numArr1 != numArr2)
  {
    std::cout << "Meshes have different numbers of point data" << std::endl;
    exit(1);
  }

  for (int i = 0; i < numArr1; ++i)
  {
    vtkDataArray* da1 = pd1->GetArray(i);
    vtkDataArray* da2 = pd2->GetArray(i);
    int numComponent = da1->GetNumberOfComponents();
    for (int j = 0; j < mesh1->getNumberOfPoints(); ++j)
    {
      double comps1[numComponent];
      double comps2[numComponent];
      da1->GetTuple(j, comps1);
      da2->GetTuple(j, comps2);
      for (int k = 0; k < numComponent; ++k)
      {
        if (std::fabs(comps1[k] - comps2[k]) > tol)
        {
          std::cout << "Meshes differ in point data values" << std::endl;
          exit(1);
        }
      }
    }
  }
  if (mesh1) delete mesh1;
  if (mesh2) delete mesh2;
  std::cout << "Meshes are the same" << std::endl;
  return EXIT_SUCCESS;
}
