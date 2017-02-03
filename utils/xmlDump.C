//
// DumpXMLFile - report on the contents of an XML or legacy vtk file
//  Usage: DumpXMLFile XMLFile1 XMLFile2 ...
//         where
//         XMLFile is a vtk XML file of type .vtu, .vtp, .vts, .vtr,
//         .vti, .vto
//


#include <vtkAnalyzer.H> 

int main (int argc, char *argv[])
{
  if (argc < 2)
  {
    std::cerr<< " VTK XML file processing utility v1.0" << std::endl
             << std::endl
             << " Looks into input vtk file(s) and reports ont contents. "<< std::endl
             << std::endl
             << " Usage: " << argv[0] << " XMLFile1 XMLFile2 ..." << std::endl;
  }

  // Process each file on the command line
  int f = 1;
  while (f < argc)
  {
    vtkAnalyzer* myVTK;
    myVTK = new vtkAnalyzer(argv[f]);
    myVTK->read();
     
    // report statistics
    myVTK->report();

    // access one of point data
    std::vector<std::vector<double> > arr;
    int numTuple, numComponent;
    std::cout << "Total number of data = " 
              << myVTK->getPointDataArray(0, arr, numTuple, numComponent);
    std::cout << " Number of tuple = " << numTuple << " number of component = " << numComponent << std::endl;

    // access one of cell data
    numTuple = 0;
    numComponent = 0;
    std::vector<std::vector<double> > arr2;
    std::cout << "Total number of data = " 
              << myVTK->getCellDataArray(1, arr2, numTuple, numComponent);
    std::cout << " Number of tuple = " << numTuple << " number of component = " << numComponent << std::endl;
    f++;
  }
  return EXIT_SUCCESS;
}
