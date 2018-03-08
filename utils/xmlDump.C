//
// DumpXMLFile - report on the contents of an XML or legacy vtk file
//  Usage: DumpXMLFile XMLFile1 XMLFile2 ...
//         where
//         XMLFile is a vtk XML file of type .vtu, .vtp, .vts, .vtr,
//         .vti, .vto
//


#include <meshBase.H> 

int main (int argc, char *argv[])
{
  if (argc < 2)
  {
    std::cerr<< " VTK XML file processing utility v1.0" << std::endl
             << std::endl
             << " Looks into input vtk file(s) and reports on its contents. "<< std::endl
             << std::endl
             << " Usage: " << argv[0] << " XMLFile1 XMLFile2 ..." << std::endl;
  }

  // Process each file on the command line
  int f = 1;
  while (f < argc)
  {
    std::string fname(argv[f]);
    std::unique_ptr<meshBase> myVTK = meshBase::CreateUnique(fname);

    // report statistics
    myVTK->report();

    vtkSmartPointer<vtkPointData> arr = myVTK->getDataSet()->GetPointData();

    if (arr->GetNumberOfArrays())
    { 
       // access one of point data
      std::cout << "Total number of data = " 
                << arr->GetNumberOfArrays();
      std::cout << " Number of tuple = " << arr->GetArray(0)->GetNumberOfTuples() 
                << " number of component = " << arr->GetArray(0)->GetNumberOfComponents() << std::endl;
    }
    vtkSmartPointer<vtkCellData> arr1 = myVTK->getDataSet()->GetCellData();
    if (arr1->GetNumberOfArrays())
    {
      // access one of cell data
      std::cout << "Total number of data = " 
                << arr1->GetNumberOfArrays();
      std::cout << " Number of tuple = " << arr->GetArray(0)->GetNumberOfTuples() 
                << " number of component = " << arr->GetArray(0)->GetNumberOfComponents() << std::endl;
    }
    f++;
  }
  return EXIT_SUCCESS;
}
