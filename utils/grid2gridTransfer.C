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

// typedefs
//typedef std::shared_ptr<cgnsAnalyzer> cgPtr;


/* auxiliary functions */
void getRegCenters(MAd::pMesh msh, std::vector<double>& regCntCrds);
std::vector<double> getCOG(std::vector<double>& regCntCrds);


/*   Main Function */ 
int main(int argc, char* argv[])
{
  // check input
  if (argc!=3 || (argc==2 && !strcmp(argv[1], "-h")) ) {
    std::cout << "Usage: " << argv[0] 
              << " srcCGNSFile disCGNSFile" << std::endl;
    return 0;
  }

  std::vector<std::string> cgFileName;
  cgFileName.push_back(argv[1]);
  cgFileName.push_back(argv[2]);
  std::cout << "Transfering from " << cgFileName[0] << " -> " << cgFileName[1] << std::endl;

  std::cout << "Reading input files #################################################\n";
  // reading source CGNS file
  gridTransfer* transObj = new gridTransfer(cgFileName[0], cgFileName[1]);
  transObj->loadSrcCgSeries(4); 
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
  transObj->writeTrgCg("test.cgns");

  /*
  // write CGNS files for the new grid
  for (int iCg=0; iCg<nOutCgFile; iCg++)
  {
     //std::ostringstream exp;
     //exp << "test_" << iCg << ".cgns";
     int1->clearCache();
     std::string fCgName;
     fCgName =cgFileName[iCg];
     std::size_t pos = fCgName.find_last_of("/");
     fCgName = fCgName.substr(pos+1);
     std::cout << "Writing remeshed " << fCgName << std::endl;
     // define elementary information
     cgnsWriter* cgWrtObj = new cgnsWriter(fCgName, cgObj1->getBaseName(), 3, 3);
     cgWrtObj->setUnits(cgObj1->getMassUnit(), cgObj1->getLengthUnit(),
			cgObj1->getTimeUnit(), cgObj1->getTemperatureUnit(),
			cgObj1->getAngleUnit());
     cgWrtObj->setBaseItrData(cgObj1->getBaseItrName(), cgObj1->getNTStep(), cgObj1->getTimeStep());
     cgWrtObj->setZoneItrData(cgObj1->getZoneItrName(), cgObj1->getGridCrdPntr(), cgObj1->getSolutionPntr());
     cgWrtObj->setZone(cgObj1->getZoneName(iCg), cgObj1->getZoneType());
     cgWrtObj->setNVrtx(mPart->getNNdePart(iCg));
     cgWrtObj->setNCell(mPart->getNElmPart(iCg));
     // define coordinates
     cgWrtObj->setGridXYZ(mPart->getCrds(iCg, MAd::M_getVrtXCrds(mesh)), 
			  mPart->getCrds(iCg, MAd::M_getVrtYCrds(mesh)), 
			  mPart->getCrds(iCg, MAd::M_getVrtZCrds(mesh)));
     // define connctivity
     cgWrtObj->setSection(cgObj1->getSectionName(), 
			  (ElementType_t) cgObj1->getElementType(), 
			  mPart->getConns(iCg));
     // define vertex and cell data 
     std::map<std::string, GridLocation_t> slnNLMap = cgObj1->getSolutionNameLocMap();
     for (auto is=slnNLMap.begin(); is!=slnNLMap.end(); is++)
       cgWrtObj->setSolutionNode(is->first, is->second);
     // write skelleton of the file
     cgWrtObj->writeGridToFile();
     // write individual data fields
     std::map<int,std::pair<int,keyValueList> > slnMap = cgObj1->getSolutionMap();
     std::vector<GridLocation_t> gLoc = cgObj1->getSolutionGridLocations();
     std::vector<std::string> slnName = cgObj1->getSolutionNodeNames();
     std::vector<double> regCntCrdsPart = mPart->getElmSlnVec(iCg, regCntCrdsNew, 3);
     //cog = getCOG(regCntCrdsPart);
     int iSol = -1;
     for (auto is=slnMap.begin(); is!=slnMap.end(); is++)
     {
       std::pair<int,keyValueList> slnPair = is->second;
       int slnIdx = slnPair.first;
       keyValueList fldLst = slnPair.second;
       for (auto ifl=fldLst.begin(); ifl!=fldLst.end(); ifl++)
       {
	 iSol++;
	 std::vector<double> stitPhysData;
	 std::vector<double> partPhysData;
	 int nData;
	 if (gLoc[iSol] == Vertex)
	 {
	   nData = mPart->getNNdePart(iCg);
	   ma->getMeshData(ifl->second, &stitPhysData);
	   partPhysData = mPart->getNdeSlnScalar(iCg, stitPhysData);
	 } else {
	   nData = mPart->getNElmPart(iCg);
	   std::vector<double> oldPhysData;
	   int nDataT, nDimT;
	   cgObj1->getSolutionDataStitched(ifl->second, oldPhysData, nDataT, nDimT);
	   int1->interpolate(mPart->getNElmPart(iCg), regCntCrdsPart, oldPhysData, partPhysData);      
           //std::cout << "Minimum element = " 
           //          << *std::min_element(partPhysData.begin(), partPhysData.end())
           //          << "\n Maximum element = " 
           //          << *std::max_element(partPhysData.begin(), partPhysData.end())
           //          << std::endl;
	 }
	 std::cout << "Writing "
		   << nData 
		   << " to "
		   << ifl->second
		   << " located in "
		   << slnName[iSol]
		   << std::endl;
	 // write to file
	 cgWrtObj->writeSolutionField(ifl->second, slnName[iSol], RealDouble, &partPhysData[0]);
       }
     }
     delete cgWrtObj;
  }
  */

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

/* returns the cartesian coordinates for the geometric center */
std::vector<double> getCOG(std::vector<double>& regCntCrds)
{
  int nNde = regCntCrds.size()/3;
  double x,y,z;
  x=0.0;
  y=0.0;
  z=0.0;
  for (int iNde=0; iNde<nNde; iNde++)
  {
    x += regCntCrds[iNde*3];
    y += regCntCrds[iNde*3 + 1];
    z += regCntCrds[iNde*3 + 2];
  }
  std::vector<double> cog;
  cog.push_back(x/nNde);
  cog.push_back(y/nNde);
  cog.push_back(z/nNde);
  std::cout << "Region center of geometry (" << cog[0]
            << ", " << cog[1] << " , " << cog[2] << ")\n";
  return(cog);
}
