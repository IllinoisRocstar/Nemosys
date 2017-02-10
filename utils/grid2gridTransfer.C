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

  std::cout << "Reading input files #################################################\n";
  // reading source CGNS file
  gridTransfer* transObj = new gridTransfer(cgFileName[0], cgFileName[1]);
  transObj->loadSrcCgSeries(1); 

  // reading input CGNS file
  transObj->exportMeshToMAdLib("src");
  transObj->convertToVtk("src", true);
  transObj->exportNodalDataToMAdLib();
  transObj->exportSrcToGModel();

  // reading target grid
  transObj->loadTrgCgSeries(1); 
  transObj->dummy();
  

  /*
  cgObj1->loadGrid();
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

  // checking quality of stitching
  //cgObj1->checkVertex();
  //cgObj1->checkElmConn(3);

  // experiment with mesh partitioner
  //meshPartitioner* mPart1 = new meshPartitioner(cgObj1);
  //mPart1->partition(4);

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

  std::cout << "Optimizing the mesh #############################################\n";
  cgObj1->classifyMAdMeshBnd(mesh); // registering boundaries for the optimization step

  // perparing for element quantity interpolation
  std::vector<double> regCntCrdsOld;
  getRegCenters(mesh, regCntCrdsOld);
  //std::vector<double> tmp = getCOG(regCntCrdsOld);
  basicInterpolant* int1 = new basicInterpolant(3, M_numRegions(mesh), 1, regCntCrdsOld);

  // prepare the mesh for the optimization
  std::cout <<"Checking mesh sanity.\n";
  MAd::checkMesh(mesh);
  MAd::M_info(mesh);

  // defining addaptive refinement parameters
  MAd::PWLSField * sizeField = new MAd::PWLSField(mesh);
  sizeField->setCurrentSize();
  sizeField->scale(1.5);
  MAd::MeshAdapter* ma = new MAd::MeshAdapter(mesh,sizeField);
  //ma->uglyTheMesh(0.5,20);
  ma->printParameters();

  // attach all nodal data to the mesh optimizer
  int nData, nDim;
  std::vector<std::string> cgSlnNameList;
  std::vector<std::string> cgAppSlnNameList;
  cgObj1->getSolutionDataNames(cgSlnNameList);  
  cgObj1->getAppendedSolutionDataName(cgAppSlnNameList);
  cgSlnNameList.insert(cgSlnNameList.end(),
                     cgAppSlnNameList.begin(), cgAppSlnNameList.end());
  for (auto is=cgSlnNameList.begin(); is<cgSlnNameList.end(); is++)
  {
    std::vector<double> physData;
    cgObj1->getSolutionDataStitched(*is, physData, nData, nDim);
    solution_type_t dt = cgObj1->getSolutionDataObj(*is)->getDataType();
    if (dt == NODAL) {
      std::cout << "Reading " << *is << std::endl;
      ma->registerData(*is, physData); 
    }
  }

  // optimize the mesh
  std::cout << "Statistics before optimization: \n";
  ma->printStatistics(std::cout);
  std::cout << "Optimizing the mesh ...\n";
  ma->run();
  std::cout << "Statistics after optimization: \n";
  ma->printStatistics(std::cout);
  // write bulk mesh
  cgObj1->unclassifyMAdMeshBnd(mesh); // remove registered boundaries for proper output
  MAd::M_writeMsh (mesh, "optimizedMesh.msh", 2, NULL);
  std::cout << "Volume Mesh "; 
  MAd::M_info(mesh);
  // convert optimized Gmsh file to vtk
  std::cout << "Converting from optimized gmsh to vtk format.\n";
  trgGModel = new GModel("stitchedOpt"); 
  trgGModel->readMSH("optimizedMesh.msh");
  trgGModel->writeVTK("optimizedMesh.vtk", false, true);

  // skinning the mesh
  std::cout << "Computing surface mesh.\n"; 
  std::vector<int> skinElmIds;
  MAd::pGModel tmpMdl = NULL;
  MAd::GM_create(&tmpMdl,"");
  MAd::pMesh skinMesh = M_new(tmpMdl);
  mesh->skin_me(skinMesh, skinElmIds);
  std::cout << "Num of elements in the mesh = " << MAd::M_numRegions(skinMesh) << std::endl;
  std::cout << "Num of triangles in the mesh = " << MAd::M_numTriangles(skinMesh) << std::endl;
  std::cout << "Minimum surface element index = " << *std::min_element(skinElmIds.begin(),skinElmIds.end()) << std::endl;
  std::cout << "Maximum surface element index = " << *std::max_element(skinElmIds.begin(),skinElmIds.end()) << std::endl;
  skinMesh->classify_unclassified_entities();
  skinMesh->destroyStandAloneEntities();
  //MAd::M_writeMsh (skinMesh, "skinMesh.msh", 2, NULL);
  MAd::SaveGmshMesh (skinMesh, "skinMesh.msh", 2, false, NULL);

  MAd::M_info(skinMesh);

  std::cout << "Transfering solution values #############################################\n";
  // write physical quantities to vtk file
  std::cout << "Writing transfered physical quantities to vtk file.\n";
  // get nodal data after refinement and write them
  vtkAnalyzer* trgVTK2;
  trgVTK2 = new vtkAnalyzer((char*)"optimizedMesh.vtk");
  trgVTK2->read();
  std::cout << "Writing cell and nodal data....\n";
  std::vector<double> regCntCrdsNew;
  getRegCenters(mesh, regCntCrdsNew);
  //std::vector<double> cog = getCOG(regCntCrdsNew);
  int nNewElm = M_numRegions(mesh);
  for (auto is=cgSlnNameList.begin(); is<cgSlnNameList.end(); is++)
  {
    solution_type_t dt = cgObj1->getSolutionDataObj(*is)->getDataType();
    if (dt == NODAL) {
      std::cout << "Writing nodal " << *is << std::endl;
      std::vector<double> physData;
      ma->getMeshData((*is), &physData);
      trgVTK2->setPointDataArray((*is).c_str(), 1, physData);
      //MAd::NodalDataManagerSgl::instance().writeData((*is),((*is)+".pos").c_str());
    } else {
      //gs field is wiered in irocstar files we dont write it back
      if (!(*is).compare("gs")){
        continue;
      }
      std::cout << "Writing cell-based " << *is << std::endl;
      std::vector<double> oldPhysData;
      std::vector<double> newPhysData;
      int nDataT, nDimT;      
      cgObj1->getSolutionDataStitched(*is, oldPhysData, nDataT, nDimT);
      int1->interpolate(nNewElm, regCntCrdsNew, oldPhysData, newPhysData);
      std::cout << "Size oldPhys = " << oldPhysData.size()
                << " Size newPhys = " << newPhysData.size()
                << " nDataT = " << nDataT
                << std::endl;
      trgVTK2->setCellDataArray((*is).c_str(), 1, newPhysData);
    }
  }
  trgVTK2->report();
  trgVTK2->write("optimizedPhys.vtu");
    
  // partition the mesh
  std::cout << " Partitioning the mesh with METIS.\n";
  //meshPartitioner* mPart = new meshPartitioner(cgObj1);
  meshPartitioner* mPart = new meshPartitioner(mesh);
  mPart->partition(nOutCgFile);
  
  // write partitin ids for the surface mesh
  std::ofstream of;
  of.open("skinMeshPart.dat", std::fstream::trunc);
  std::vector<double> elmPartIds = mPart->getPartedElm();
  for (auto it=skinElmIds.begin(); it!=skinElmIds.end(); it++)
  {
    of << elmPartIds[*it - 1] << "\n";
  } 
  of.close();

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
