// standard headers
#include <cstring>
#include <string.h>
#include <iostream>
#include <memory>
#include <ctime>

// Nemosys headers
#include "gridTransfer.H"

// MAdLib headers 
#include "MAdLib.h"

/* auxiliary functions */
void getRegCenters(MAd::pMesh msh, std::vector<double>& regCntCrds);
std::vector<double> getCOG(std::vector<double>& regCntCrds);
void helpExit();

void helpExit()
{
  std::cout << "NEMoSys Webservice Backend  (version 1.0) \n"
	    << "A grid-to-grid transfer utility \n"
	    << "\n"
	    << "Usage: grid2gridTransfer --command  sourceGridFile [targetGridFile]\n" 
	    << "where command can be :\n"
	    << "help or h :  provides current help.\n"
	    << "transCGNS : to transfer quantities between CGNS grids.\n"
	    << " statCGNS : gives statistics about the grid.\n"
	    << "  chkCGNS : runs MAdLib checks on the grid.\n"
	    << " lslnCGNS : lists soltion names exisiting on a CGNS grid.\n"
	    << std::endl;
  exit(0);
}

/*   Main Function */ 
int main(int argc, char* argv[])
{
  // check input
  if ((argc<=2) || (argc==2 && !strcmp(argv[1], "--h"))
                || (argc==2 && !strcmp(argv[1], "--help")))
    helpExit();
  
  std::vector<std::string> cmd;
  cmd.push_back(argv[1]);
  // switch processing

  // --tCGNS
  if (!strcmp(cmd[0].c_str(),"--transCGNS"))
  {
    // input processing
    if (argc < 4)
      helpExit();
    std::vector<std::string> cgFileName;
    cgFileName.push_back(argv[2]);
    cgFileName.push_back(argv[3]);
    std::cout << "Transfering between the grids #################################################\n";
    std::cout << "Transfering from " << cgFileName[0] << " -> " << cgFileName[1] << std::endl;
    // reading source CGNS file
    gridTransfer* transObj = new gridTransfer(cgFileName[0], cgFileName[1]);
    transObj->loadSrcCgSeries(1); 
    transObj->exportMeshToMAdLib("src");
    transObj->convertToVTK("src", true);
    //transObj->exportNodalDataToMAdLib();
    transObj->exportSrcToGModel();
    // reading target grid
    transObj->loadTrgCg(); 
    transObj->exportMeshToMAdLib("trg");
    transObj->convertToMsh("trg");
    // transfer physical values between the meshes 
    transObj->transfer();
    // write target with solutions to vtk 
    transObj->convertToVTK("trg", true);
    // write target to CGNSL 
    transObj->writeTrgCg("target.cgns");
    std::cout << "Transfer finished successfully.\n";
  }

  // --statCGNS
  else if (!strcmp(cmd[0].c_str(),"--statCGNS"))
  {
    // input processing
    if (argc < 3)
      helpExit();
    std::vector<std::string> cgFileName;
    cgFileName.push_back(argv[2]);
    std::cout << "Acquiring grid statistics ####################################################\n";
    // reading source CGNS file
    gridTransfer* transObj = new gridTransfer(cgFileName[0], "dummy");
    transObj->loadSrcCgSeries(1); 
    transObj->exportMeshToMAdLib("src");
    transObj->gridStats();
    std::cout << "Statistics generated finished successfully.\n";
  }

  // --chkCGNS
  else if (!strcmp(cmd[0].c_str(),"--chkCGNS"))
  {
    // input processing
    if (argc < 3)
      helpExit();
    std::vector<std::string> cgFileName;
    cgFileName.push_back(argv[2]);
    std::cout << "Checking the grid ###########################################################\n";
    // reading source CGNS file
    gridTransfer* transObj = new gridTransfer(cgFileName[0], "dummy");
    transObj->loadSrcCgSeries(1); 
    transObj->exportMeshToMAdLib("src");
    transObj->gridCheck();
    std::cout << "Finished checking grid successfully.\n";
  }

  // --lslnCGNS
  else if (!strcmp(cmd[0].c_str(),"--lslnCGNS"))
  {

  }
  else
  {
    std::cout << "Error: " << cmd[0] << " is not a valid command.\n";
    std::cout << "For more help type " << argv[0] << " -h \n";
    throw;
  }

  // ending application
  std::cout << "Request processed successfully.\n";
  return 0;
}



///////////////////////////////////////////////////////////////////////
//                       AUX FUNCTIONS                               // 
///////////////////////////////////////////////////////////////////////

/* get element/cell/region center coordinates */
void getRegCenters(MAd::pMesh msh, std::vector<double>& regCntCrds)
{
   MAd::RIter ri = M_regionIter(msh);
   int rCnt = 0;
   while (MAd::pRegion pr = RIter_next(ri)) 
   {
     double xc[3];
     MAd::R_center(pr, xc);
     regCntCrds.push_back(xc[0]);
     regCntCrds.push_back(xc[1]);
     regCntCrds.push_back(xc[2]);
     /* 
     std::cout << "Region " 
               << rCnt++ 
               << " center coordinate = "
               << xc[0] << " "
               << xc[1] << " "
               << xc[2] << " "
               << std::endl;
     */
   }
} 

