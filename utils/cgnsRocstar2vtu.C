/*   Program write_grid_str.c    */
/*
Creates simple 3-D structured grid and writes it to a
CGNS file.

Example compilation for this program is (change paths!):

cc -I ../CGNS_CVS/cgnslib -c write_grid_str.c
cc -o write_grid_str_c write_grid_str.o -L ../CGNS_CVS/cgnslib/LINUX -lcgns

(../CGNS_CVS/cgnslib/LINUX/ is the location where the compiled
library libcgns.a is located)
*/

// standard headers
#include <cstring>
#include <string.h>
#include <iostream>
#include <memory>
#include <algorithm>
#include <sstream>
#include <vector>
#include <dirent.h>
#include <sys/types.h>

// Nemosys headers
#include "cgnsAnalyzer.H"
#include "vtkAnalyzer.H"
#include "baseInterp.H"
#include "meshPartitioner.H"
#include "cgnsWriter.H"

// MAdLib headers 
#include "ModelInterface.h"
#include "MAdLib.h"
#include "NodalDataManager.h"

// typedefs
//typedef std::shared_ptr<cgnsAnalyzer> cgPtr;


/* auxiliary functions */
void getRegCenters(MAd::pMesh msh, std::vector<double>& regCntCrds);

std::vector <std::string> read_directory(  std::string path );

void get_file_names(std::vector<std::string> &names, std::string distinct);

int get_core_count(std::vector<std::string> names);

int stoi(std::string str);


/*   Main Function */ 
int main(int argc, char* argv[])
{
  // check input
  int nInCgFile = 0;
  int nOutCgFile = 0;
  
  
  if (argc==1 || (argc==2 && !std::strcmp(argv[1], "-h")) ) {
    std::cout << "Usage: " << argv[0] 
              << " /path/to/directory begginingOfUniqueFileName" << std::endl;
    return 0;
  }

  std::string directory = argv[1];

  std::vector<std::string> names;

  names = read_directory(directory);

  get_file_names(names, argv[2]);


  std::string::size_type sz;   // alias of size_t
  nInCgFile = get_core_count(names);
  nOutCgFile = nInCgFile; // setting the number of output files the same as input
  std::vector<std::string> cgFileName(nInCgFile);


//Start for loop of each processor
for(int outer = 0; outer < names.size(); outer += nInCgFile)
{
	std::cout << "Outer: " << outer << std::endl;
  for (int iCg=0; iCg < nInCgFile; iCg++)
     cgFileName[iCg] = names.at(iCg + outer);



  std::cout << "Reading input files #################################################\n";
  // reading cgns file1
  cgnsAnalyzer* cgObj1 = new cgnsAnalyzer(cgFileName[0]);
  cgObj1->loadGrid();
  cgObj1->checkVertex();
  // attaching partion id to the mesh
  std::vector<double> slnData(cgObj1->getNElement(),0);
  cgObj1->appendSolutionData("partitionOld", slnData, ELEMENTAL, cgObj1->getNElement(), 1);
  // adding new cgns files
  for (int iCg=1; iCg < nInCgFile; iCg++) {
    cgnsAnalyzer* cgObj2 = new cgnsAnalyzer(cgFileName[iCg]);
    cgObj2->loadGrid();
    // appending data
    std::vector<double> slnData(cgObj2->getNElement(),iCg);
    cgObj2->appendSolutionData("partitionOld", slnData, ELEMENTAL, cgObj2->getNElement(), 1);
    // stitching meshes
    cgObj1->stitchMesh(cgObj2, true);
  }
  if (nInCgFile > 1)
     std::cout << "Meshes stitched successfully!\n";


  std::cout << "Exporting mesh to MAdLib format #####################################\n";
  // exporting mesh to the MAdLib
  MAd::pGModel model = NULL;
  MAd::GM_create(&model,"");
  MAd::pMesh mesh = M_new(model);
  cgObj1->exportToMAdMesh(mesh);
  cgObj1->classifyMAdMeshOpt(mesh);
  
  // writing the mesh to gmsh and convert to vtk
  std::cout <<"Writing to gmsh format.\n";
  M_writeMsh(mesh, "stitched.msh", 2, NULL);
  std::cout << "Converting from gmsh to vtk format.\n";
  GModel* trgGModel;
  trgGModel = new GModel("stitched"); 
  trgGModel->readMSH("stitched.msh");
  trgGModel->writeVTK("stitched.vtk", false, true);
  std::string output = "stitched" + std::to_string(outer/nInCgFile) + ".vtk";
  trgGModel->writeVTK(output.c_str(), false, true);

  // write physical quantities to vtk file
  std::cout << "Writing physical quantities to vtk file.\n";
  vtkAnalyzer* trgVTK;
  trgVTK = new vtkAnalyzer((char*)"stitched.vtk");
  trgVTK->read();
  trgVTK->report();
  // figure out what is existing on the stitched grid
  int outNData, outNDim;
  std::vector<std::string> slnNameList;
  std::vector<std::string> appSlnNameList;
  cgObj1->getSolutionDataNames(slnNameList);  
  cgObj1->getAppendedSolutionDataName(appSlnNameList);
  slnNameList.insert(slnNameList.end(),
                     appSlnNameList.begin(), appSlnNameList.end());
  // write all data into vtk file
  for (auto is=slnNameList.begin(); is<slnNameList.end(); is++)
  {
    std::vector<double> physData;
    cgObj1->getSolutionDataStitched(*is, physData, outNData, outNDim);
    solution_type_t dt = cgObj1->getSolutionDataObj(*is)->getDataType();
    if (dt == NODAL)      
      trgVTK->setPointDataArray((*is).c_str(), 1, physData);
    else
      trgVTK->setCellDataArray((*is).c_str(), 1, physData);
  }
  trgVTK->report();
  trgVTK->write("stitchedPhys.vtu");
  output = "stitchedPhys" + std::to_string(outer/nInCgFile) + ".vtu";
  std::vector<char> char_array(output.begin(), output.end());
  char_array.push_back(0);
  trgVTK->write(&char_array[0]);
  
}
  std::cout << "Application ended successfully!\n";
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


int stoi(std::string str)
{
   std::stringstream ss(str);
   int n;
   ss << str;
   ss >> n;
   return n;
}

int get_core_count(std::vector<std::string> names)
{
  int previous = 0;
  int current = 0;


  for(int i = 0; i < names.size() && previous <= current; i++)
  {
    previous = current;
    current = std::stoi(names.at(i).substr(names.at(i).length() - 9, 4));

  }

  return previous + 1;
}

void get_file_names(std::vector<std::string> &names, std::string distinct)
{
  bool already_erased = false;

  for(int i = 0; i < names.size();i++)
  {
    already_erased = false;
    if(distinct.compare(names.at(i).substr(0,distinct.size())) != 0)
    {
      names.erase(names.begin() + i);
      i--;
      already_erased = true;
    }

    if( !already_erased && names.at(i).substr(names.at(i).length() - 4 ) != "cgns")
    {
      names.erase(names.begin() + i);
      i--;
    }

  }
}





// read_directory()
//   Return an ASCII-sorted vector of filename entries in a given directory.
//   If no path is specified, the current working directory is used.
//
//   Always check the value of the global 'errno' variable after using this
//   function to see if anything went wrong. (It will be zero if all is well.)
//
std::vector <std::string> read_directory( std::string path )
  {
  std::vector <std::string> result;
  dirent* de;
  DIR* dp;
  errno = 0;
  dp = opendir( path.empty() ? "." : path.c_str() );
  if (dp)
    {
    while (true)
      {
      errno = 0;
      de = readdir( dp );
      if (de == NULL) break;
      result.push_back( std::string( de->d_name ) );
      }
    closedir( dp );
    std::sort( result.begin(), result.end() );
    }
  return result;
  }


