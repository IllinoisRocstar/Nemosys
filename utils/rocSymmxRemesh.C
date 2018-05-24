// standard headers
#include <cstring>
#include <string.h>
#include <iostream>
#include <memory>

// Nemosys headers
#include <cgnsAnalyzer.H>
#include <meshBase.H>
#include <symmxGen.H>
#include <symmxParams.H>
#include <meshPartitioner.H>
#include <meshStitcher.H>
#include <cgnsWriter.H>

/*   Main Function */ 
int main(int argc, char* argv[])
{
  // check input
  int nInCgFile = 0;
  int nOutCgFile = 0;
  
  std::vector<std::string> cgFileName;
  if (argc==1 || (argc==2 && !std::strcmp(argv[1], "-h")) ) {
    std::cout << "Usage: " << argv[0] 
              << " nCGNSFileIn nCGNSFileOut inCgFileName1 inCgFileName2 ..." << std::endl;
    return 0;
  }
  std::string::size_type sz;   // alias of size_t
  nInCgFile = std::stoi(argv[1],&sz);
  nOutCgFile = nInCgFile; // setting the number of output files the same as input  
  //nOutCgFile = std::stoi(argv[2],&sz); 
  for (int iCg=0; iCg < nInCgFile; iCg++)
     cgFileName.push_back(argv[2+iCg]);

  meshStitcher* stitcher = new meshStitcher(cgFileName);
  cgnsAnalyzer* cgObj1 = stitcher->getStitchedCGNS();

  ///////////////////// REMESH WITH SIMMETRIX //////////////////////////////////////////

  meshBase* trgVTK = stitcher->getStitchedMB();
  std::cout << "Extracting surface mesh #############################################\n";
  meshBase* trgVTKsurf = meshBase::Create(trgVTK->extractSurface(), "skinMeshPart.vtp"); 
 
  // remeshing with simmetrix
  trgVTKsurf->write("skinMeshPart.stl");
  symmxParams* params = new symmxParams(); 
  params->logFName = "symmxGen.log";
  params->licFName = "/home/users/snatesh/snatesh-storage/SYMMETRIX/LICENSE/simmodsuite.lic";
  params->features = "geomsim_core,meshsim_surface,meshsim_volume,geomsim_discrete";
  meshBase* newMesh = meshBase::generateMesh("skinMeshPart.stl","simmetrix",params);
  trgVTK->setContBool(0);

  ///////////////////// PARTITION THE MESH WITH METIS /////////////////////////////////

  std::cout << " Partitioning the mesh with METIS.\n";
  meshPartitioner* mPart = new meshPartitioner(newMesh);
  mPart->partition(nOutCgFile);
  
  // write CGNS files for the new grid
  for (int iCg=0; iCg<nOutCgFile; iCg++)
  {
    std::string fCgName;
    fCgName =cgFileName[iCg];
    std::size_t pos = fCgName.find_last_of("/");
    fCgName = fCgName.substr(pos+1);
    fCgName = trim_fname(fCgName,"New.cgns");
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
    std::vector<std::vector<double>> comp_crds(newMesh->getVertCrds()); 
    cgWrtObj->setGridXYZ(mPart->getCrds(iCg, comp_crds[0]), 
		                     mPart->getCrds(iCg, comp_crds[1]), 
		                     mPart->getCrds(iCg, comp_crds[2]));
    // define connctivity
    cgWrtObj->setSection(cgObj1->getSectionName(), 
		   (ElementType_t) cgObj1->getElementType(), 
		   mPart->getConns(iCg));

    // write partitioned vtk files with data transfered from stitched mesh
    std::vector<int> vtkConn(mPart->getConns(iCg));
    for (auto it = vtkConn.begin(); it != vtkConn.end(); ++it)
    {
      *it -= 1;
    }
    std::stringstream vtkname;
    vtkname << "tmp" << iCg << ".vtu";
    meshBase* tmpvtk = meshBase::Create(mPart->getCrds(iCg, comp_crds[0]),
                                        mPart->getCrds(iCg, comp_crds[1]),
                                        mPart->getCrds(iCg, comp_crds[2]),
                                        vtkConn, VTK_TETRA, vtkname.str());
    trgVTK->transfer(tmpvtk, "Consistent Interpolation");
    tmpvtk->write();


    // define vertex and cell data 
    std::map<std::string, GridLocation_t> slnNLMap = cgObj1->getSolutionNameLocMap();
    for (auto is=slnNLMap.begin(); is!=slnNLMap.end(); is++)
      cgWrtObj->setSolutionNode(is->first, is->second);
    // write skeleton of the file
    cgWrtObj->writeGridToFile();

    /////////////////////// WRITE SOLUTION DATA TO CGNS ///////////////////////////////

    // write individual data fields
    std::map<int,std::pair<int,keyValueList> > slnMap = cgObj1->getSolutionMap();
    std::vector<GridLocation_t> gLoc = cgObj1->getSolutionGridLocations();
    std::vector<std::string> slnName = cgObj1->getSolutionNodeNames();

    int iSol = -1;
    for (auto is=slnMap.begin(); is!=slnMap.end(); is++)
    {
      std::pair<int,keyValueList> slnPair = is->second;
      int slnIdx = slnPair.first;
      keyValueList fldLst = slnPair.second;
      for (auto ifl=fldLst.begin(); ifl!=fldLst.end(); ifl++)
      {
	      iSol++;
	      std::vector<double> partPhysData;
	      int nData;
	      if (gLoc[iSol] == Vertex)
	      {
	       tmpvtk->getPointDataArray(ifl->second, partPhysData);              
	      } 
        else 
        {
          tmpvtk->getCellDataArray(ifl->second, partPhysData);
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
    delete tmpvtk;       
  }
 
  if (trgVTKsurf) delete trgVTKsurf; 
  if (newMesh) delete newMesh;
  if (stitcher){ delete stitcher; cgObj1 = 0; trgVTK = 0;}
  if (params) delete params;
  if (mPart) delete mPart;
  
  std::cout << "Application ended successfully!\n";

  return 0;
}
