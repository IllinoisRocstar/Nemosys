// standard headers
#include <cstring>
#include <string.h>
#include <iostream>
#include <memory>

// Nemosys headers
#include <meshBase.H>
#include <vtkThreshold.h>
#include <vtkGeometryFilter.h>
#include <vtkCellData.h>

int main(int argc, char* argv[])
{
  if (argc==1 || (argc==2 && !std::strcmp(argv[1], "-h")) ) {
    std::cout << "Usage: " << argv[0] 
              << " gmshFilename" << std::endl
              << "Eg) gmsh2Cobalt myMesh.msh" << std::endl;
    return 0;
  }

  std::string fname(argv[1]);
  std::string trimmedName;
  if (fname.find(".msh") != -1) 
  {
    // get trimmed filename
    size_t lastindex = fname.find_last_of("."); 
    trimmedName = fname.substr(0, lastindex);
  }
  else
  {
    std::cout << "File must be in .msh format." << std::endl;
    exit(1);
  }

  
  // create meshbase object
  std::shared_ptr<meshBase> myMesh = meshBase::CreateShared(fname);

  // declare array to store element type
  vtkSmartPointer<vtkDataArray> elementTypeArray = vtkSmartPointer<vtkIdTypeArray>::New();;
  elementTypeArray->SetNumberOfComponents(1);
  elementTypeArray->SetName("elemType");

  // get element types
  double tmp[1];  // storing inserted tuples
  for (int iCell = 0; iCell < myMesh->getDataSet()->GetNumberOfCells(); iCell++)
  {
    if (myMesh->getDataSet()->GetCellType(iCell) != VTK_TETRA &&
        myMesh->getDataSet()->GetCellType(iCell) != VTK_TRIANGLE)
    {
      std::cout << "Only triangle and tetrahedral elements supported." << std::endl;
      exit(1);
    }
    tmp[0] = myMesh->getDataSet()->GetCellType(iCell);
    elementTypeArray->InsertNextTuple(tmp);
  }
  myMesh->getDataSet()->GetCellData()->AddArray(elementTypeArray);

  // threshold to get only volume cells
  vtkSmartPointer<vtkThreshold> volThreshold = vtkSmartPointer<vtkThreshold>::New();
  volThreshold->SetInputData(myMesh->getDataSet());
  volThreshold->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "elemType");
  volThreshold->ThresholdBetween(VTK_TETRA, VTK_TETRA);
  volThreshold->Update();
  // store output in unstructured grid
  vtkSmartPointer<vtkUnstructuredGrid> volUG = volThreshold->GetOutput();
  // create meshbase object
  std::shared_ptr<meshBase> volUG_mb = meshBase::CreateShared(volUG, "extractedVolume.vtu");
  volUG_mb->write();  // write to file

  // threshold to get only surface cells
  vtkSmartPointer<vtkThreshold> surfThreshold = vtkSmartPointer<vtkThreshold>::New();
  surfThreshold->SetInputData(myMesh->getDataSet());
  surfThreshold->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "elemType");
  surfThreshold->ThresholdBetween(VTK_TRIANGLE, VTK_TRIANGLE);
  surfThreshold->Update();
  // store output in unstructured grid
  vtkSmartPointer<vtkUnstructuredGrid> surfUG = surfThreshold->GetOutput();
  // geometry filter to polydata object
  vtkSmartPointer<vtkGeometryFilter> geomFilter = vtkGeometryFilter::New();
  geomFilter->SetInputData(surfUG);
  geomFilter->Update();
  vtkSmartPointer<vtkPolyData> surfPD = geomFilter->GetOutput();
  // create meshbase object
  std::shared_ptr<meshBase> surfPD_mb = meshBase::CreateShared(surfPD, "extractedSurface.vtp");
  surfPD_mb->write();  // write to file

  // create cobalt output names
  std::string mapName = trimmedName + ".cgi";
  std::string grdName = trimmedName + ".cgr";

  // write to cobalt format
  volUG_mb->writeCobalt(surfPD_mb.get(), mapName, grdName);

  std::cout << "Application ended successfully!\n";
  return 0;
}

