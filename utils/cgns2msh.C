#include "cgnsAnalyzer.H"
#include <cstring>
#include <string.h>
#include <iostream>
#include <memory>
// MAdLib headers 
#include "ModelInterface.h"

typedef std::shared_ptr<cgnsAnalyzer> cgPtr;
 
int main(int argc, char* argv[])
{
  std::string cgFileName;

  // check input
  if (argc>3 || (argc==2 && !std::strcmp(argv[1], "-h")) ){
    std::cout << "Usage: " << argv[0] << " cgnsFileName mshFileName" << std::endl;
    return 0;
  }

  // setting to default if needed
  if (argc == 1)
    cgFileName = "test.cgns";
  else
    cgFileName = argv[1];

  // writing the sample cgns file
  cgPtr cgObj(new cgnsAnalyzer(cgFileName));
  if (argc == 1)
  {
    //cgObj->writeSampleStructured();
    cgObj->writeSampleUnstructured();
    std::cout << "Successfuly wrote grid to file test.cgns\n";
  }
  // reading the sample file
  cgObj->loadGrid();

  // get the list of solution data
  std::vector<std::string> list;
  cgObj->getSolutionDataNames(list);
  std::cout << "Number of solution field = " << list.size() << std::endl;

  // access one of data
  std::vector<double> slnData;
  int dataType;
  dataType = cgObj->getSolutionData(list[0], slnData);  // returns: 3 cell data, 2 for vertex
  std::cout << "Accessing the first dataset.\n";
  std::cout << list[0] << " has " << slnData.size()
            << " data in it.\n"
            << "data[0] = " << slnData[0] 
            << "\ndataType = " << (dataType==3 ? "Elemental":"Nodal")
            << std::endl;

  // exporting mesh to the MAdLib
  MAd::pGModel model = NULL;
  MAd::GM_create(&model,"");
  MAd::pMesh mesh = M_new(model);
  cgObj->exportToMAdMesh(mesh);

  // investigating sanity of the mesh
  MAd::M_info(mesh);

  // converting to msh format
  std::string mshFileName;
  mshFileName = (argc==3) ? argv[2]:"output.msh";
  M_writeMsh(mesh, mshFileName.c_str(), 2, NULL);

  std::cout << "The file was converted to Gmesh format successfully!" << std::endl;
  return 0;
}


