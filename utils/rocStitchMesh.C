// standard headers
#include <cstring>
#include <string.h>
#include <iostream>
#include <memory>

// Nemosys headers
#include <meshBase.H>
#include <meshPartitioner.H>
#include <meshStitcher.H>
#include <cgnsWriter.H>
#include <vtkAppendFilter.h>

/*   Main Function */ 
int main(int argc, char* argv[])
{
  // check input
  int nInCgFile = 0;
  
  std::vector<std::string> cgFileName;
  if (argc==1 || (argc==2 && !std::strcmp(argv[1], "-h")) ) {
    std::cout << "Usage: " << argv[0] 
              << " nCgFile CgFileName0" << std::endl;
    return 0;
  }
  std::string::size_type sz;   // alias of size_t
  nInCgFile = std::stoi(argv[1],&sz);
  meshStitcher* stitcher = new meshStitcher(nInCgFile, argv[2]);
  meshBase* trgVTK = stitcher->getStitchedMB();
 
  if (stitcher){ delete stitcher; trgVTK = 0;}
  
  //std::unique_ptr<meshBase> mesh1 = meshBase::CreateUnique(argv[1]);
  //std::unique_ptr<meshBase> mesh2 = meshBase::CreateUnique(argv[2]);
  //mesh1->unsetCellDataArray("mdot_old");
  //vtkSmartPointer<vtkAppendFilter> appender
  //  = vtkSmartPointer<vtkAppendFilter>::New();
  //appender->AddInputData(mesh1->getDataSet());
  //appender->AddInputData(mesh2->getDataSet());
  //appender->Update();
  //std::unique_ptr<meshBase> mesh3
  //  = std::unique_ptr<meshBase>(meshBase::Create(appender->GetOutput(), "stitched_ifluid_b_ni.vtu"));
  //mesh3->report();
  //mesh3->write();
  std::cout << "Application ended successfully!\n";
  return 0;
}
