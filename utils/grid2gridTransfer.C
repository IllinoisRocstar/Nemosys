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
	    << "Usage: grid2gridTransfer --command  sourceGridFile [targetGridFile] [extraOptions]\n" 
	    << "where command can be :\n"
	    << "help or h :  provides current help.\n"
	    << "transCGNS : to transfer quantities between CGNS grids.\n"
	    << " statCGNS : gives statistics about the grid.\n"
	    << "  chkCGNS : runs MAdLib checks on the grid.\n"
	    << " cgns2stl : skins a CGNS grid and writes it into a stl file.\n"
	    << " cgns2vtk : converts cgns grid into vtk format.\n"
	    << "  listSln : lists current existing solution names defined on the grid and their type.\n"
            << " plotHist : generats mesh quality histogram and writes it to hist.json\n"
	    << std::endl;
  exit(0);
}

/*   Main Function */ 
int main(int argc, char* argv[])
{
  // handshake message
  std::cout << "NEMoSys Webserver (v1.0)\n";
  std::cout << "Processing request....\n";

  // check input
  if ((argc<=2) || (argc==2 && !strcmp(argv[1], "--h"))
                || (argc==2 && !strcmp(argv[1], "--help")))
    helpExit();
  
  std::vector<std::string> cmd;
  cmd.push_back(argv[1]);
  // switch processing

  // --transCGNS
  if (!strcmp(cmd[0].c_str(),"--transCGNS"))
  {
    // input processing
    bool calcErr = false;
    if (argc < 5)
    {
      helpExit();
      std::cout << "Example: " << argv[0]
                << "--transCGNS source.cgns target.cgns target_with_solution.cgns [--withErr]"
                << std::endl;
    }
    std::vector<std::string> cgFileName;
    cgFileName.push_back(argv[2]);
    cgFileName.push_back(argv[3]);
    cgFileName.push_back(argv[4]);
    if (argc==6 && !strcmp(argv[5],"--withErr"))
      calcErr = true;
    std::cout << "Transfering between the grids ##########\n";
    std::cout << "Transfering from " << cgFileName[0] << " -> " << cgFileName[1] << std::endl;
    // reading source CGNS file
    gridTransfer* transObj = new gridTransfer(cgFileName[0], cgFileName[1]);
    transObj->loadSrcCg(); 
    transObj->exportMeshToMAdLib("src");
    transObj->convertToVTK("src", true);
    //transObj->exportNodalDataToMAdLib();
    transObj->exportToGModel("src");
    // reading target grid
    transObj->loadTrgCg(); 
    transObj->exportMeshToMAdLib("trg");
    transObj->convertToMsh("trg");
    // transfer physical values between the meshes 
    transObj->transfer();
    // write target with solutions to vtk 
    transObj->convertToVTK("trg", true);
    // write target to CGNSL 
    transObj->writeTrgCg(cgFileName[2]);
    // report on solution transfer quality
    if (calcErr) {
      std::vector<std::string> slnList;
      transObj->getSolutionDataNames(slnList);
      int slnIdx = 1;
      std::cout << setw(3) << " # " 
		<< setw(12) << "   Name   " 
		<< setw(12) << "    Type    "
		<< setw(14) << "  RMS Error \n";
      std::cout << setw(3) << "---" 
		<< setw(12) << "----------" 
		<< setw(12) << "----------" 
		<< setw(14) << "------------\n";
      for (auto it=slnList.begin(); it!=slnList.end(); it++){
	solution_type_t slnType;
	std::vector<double> slnData;
	slnType = transObj->getSolutionData(*it, slnData);
	double accuracy = transObj->calcTransAcc(*it);
	std::cout << setw(3) << slnIdx++ 
		  << setw(12) << *it 
		  << setw(12) << ((slnType==0)? "Nodal":"Elemental") 
		  << setw(13) << accuracy
		  << std::endl;
      }
    }
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
    std::cout << "Acquiring grid statistics ############\n";
    // reading source CGNS file
    gridTransfer* transObj = new gridTransfer(cgFileName[0], "dummy");
    transObj->loadSrcCg(); 
    transObj->exportMeshToMAdLib("src");
    transObj->gridStats();
    std::cout << "Statistics generated successfully.\n";
  }

  // --plotHist : added by Woo for generating CGNS
  else if (!strcmp(cmd[0].c_str(),"--plotHist"))
  {
    // input processing
    if (argc < 3)
      helpExit();
    std::vector<std::string> cgFileName;
    cgFileName.push_back(argv[2]);
    std::cout << "Acquiring histogram stats ############\n";
    // reading source CGNS file
    gridTransfer* transObj = new gridTransfer(cgFileName[0], "dummy");
    transObj->loadSrcCg(); 
    transObj->exportMeshToMAdLib("src");
    transObj->gridHist();
    std::cout << "Statistics generated successfully.\n";
  }


  // --chkCGNS
  else if (!strcmp(cmd[0].c_str(),"--chkCGNS"))
  {
    // input processing
    if (argc < 3)
      helpExit();
    std::vector<std::string> cgFileName;
    cgFileName.push_back(argv[2]);
    std::cout << "Checking the grid #########\n";
    // reading source CGNS file
    gridTransfer* transObj = new gridTransfer(cgFileName[0], "dummy");
    transObj->loadSrcCg(); 
    transObj->exportMeshToMAdLib("src");
    transObj->gridCheck();
    std::cout << "Finished checking grid successfully.\n";
  }

  // --cgns2stl
  else if (!strcmp(cmd[0].c_str(),"--cgns2stl"))
  {
    // input processing
    if (argc < 5) {
      std::cout << "Example : " << argv[0] 
                << " --cgns2stl srouce.cgns output_path [src or trg]"
                << std::endl;
      helpExit();
    }
    std::vector<std::string> cgFileName;
    cgFileName.push_back(argv[2]);
    std::string dist = argv[3];
    std::string gridName = argv[4];
    std::string newname = argv[5];
    if (!strcmp(gridName.c_str(), "src") && !strcmp(gridName.c_str(), "trg"))
    {
      std::cout << "The forth switch should be src or trg " << std::endl;
      helpExit();
    }
    std::cout << "Converting to stl ########\n";
    // reading the CGNS file
    gridTransfer* transObj = new gridTransfer(cgFileName[0], cgFileName[0]);
    if (!strcmp(gridName.c_str(), "src")) 
    {
      transObj->loadSrcCg(); 
    } else {
      transObj->loadTrgCg();
    }
    transObj->exportMeshToMAdLib(gridName);
    transObj->convertToSTL(gridName, dist,newname);
    std::cout << "Finished conversion of " << gridName << " to STL successfully.\n";
  }

  // --cgns2json
  else if (!strcmp(cmd[0].c_str(),"--cgns2json"))
  {
    // input processing
    if (argc < 5) {
      std::cout << "Example : " << argv[0] 
                << " --cgns2json srouce.cgns output_path [src or trg]"
                << std::endl;
      helpExit();
    }
    std::vector<std::string> cgFileName;
    cgFileName.push_back(argv[2]);
    std::string dist = argv[3];
    std::string gridName = argv[4];
    if (!strcmp(gridName.c_str(), "src") && !strcmp(gridName.c_str(), "trg"))
    {
      std::cout << "The forth switch should be src or trg " << std::endl;
      helpExit();
    }
    std::cout << "Converting to json ########\n";
    // reading the CGNS file
    gridTransfer* transObj = new gridTransfer(cgFileName[0], cgFileName[0]);
    if (!strcmp(gridName.c_str(), "src")) 
    {
      transObj->loadSrcCg(); 
    } else {
      transObj->loadTrgCg();
    }
    transObj->exportMeshToMAdLib(gridName);
    transObj->convertToJSON(gridName, dist, false);
    std::cout << "Finished conversion of " << gridName << " to JSON successfully.\n";
  }

  // --cgns2vtk
  else if (!strcmp(cmd[0].c_str(),"--cgns2vtk"))
  {
    // input processing
    if (argc < 5) 
    {
      std::cout << "Example " << argv[0]
                << "--cgns2vtk grid.cgns path_to_distenation_vtk [src or trg] "
                << std::endl;
      helpExit();
    }
    std::vector<std::string> cgFileName;
    cgFileName.push_back(argv[2]);
    std::string dist = argv[3];
    std::string gridName = argv[4];
    if (!strcmp(gridName.c_str(), "src") && !strcmp(gridName.c_str(), "trg"))
    {
      std::cout << "The forth switch should be src or trg " << std::endl;
      helpExit();
    }
    std::cout << "Converting to vtk ########\n";
    // reading the CGNS file
    gridTransfer* transObj = new gridTransfer(cgFileName[0], cgFileName[0]);
    if (!strcmp(gridName.c_str(), "src"))
      transObj->loadSrcCg(); 
    else
      transObj->loadTrgCg();
    transObj->exportMeshToMAdLib(gridName);
    transObj->convertToVTK(gridName, true, dist);
    std::cout << "Finished conversion of " << gridName << " to VTK successfully.\n";
  }

  // --listSln
  else if (!strcmp(cmd[0].c_str(),"--listSln"))
  {
    // input processing
    if (argc < 3)
      helpExit();
    std::vector<std::string> cgFileName;
    cgFileName.push_back(argv[2]);
    std::cout << "Reading solution names ########\n";
    // reading the CGNS file
    gridTransfer* transObj = new gridTransfer(cgFileName[0], "dummy");
    transObj->loadSrcCg(); 
    std::vector<std::string> slnList;
    transObj->getSolutionDataNames(slnList);
    int slnIdx = 1;
    std::cout << setw(3) << " # " << setw(12) << "   Name   " << setw(16) << "    Type    \n";
    std::cout << setw(3) << "---" << setw(12) << "----------" << setw(16) << "-----------\n";
    for (auto it=slnList.begin(); it!=slnList.end(); it++){
      solution_type_t slnType;
      std::vector<double> slnData;
      slnType = transObj->getSolutionData(*it, slnData);
      std::cout << setw(3) << slnIdx++ << setw(12) << *it << setw(15) << ((slnType==0)? "NODAL":"Elemental") << std::endl;
    }
    std::cout << "Finished reading solution names successfully.\n";
  }

  // wrong switch
  else
  {
    std::cout << "Error: " << cmd[0] << " is not a valid command.\n";
    std::cout << "For more help type " << argv[0] << " -h \n";
    throw;
  }

  // ending application
  std::cout << "Request processed successfully #######\n";
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

