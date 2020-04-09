/*
  Sample surface-to-surface mesh adaptation and solution transfer 
  utility. The mesh information (rocstar sample files) will be
  read from CGNS dumps. Partitions of the mesh first stitched
  together. All nodal, elemental and boundary condition (specific to
  rocstar in the current implementation) are read and assembled into
  a master mesh. The mesh is then registered to MAdLib/gmsh and
  a bunch of refinement actions will be performed. Finally the solutio
  data and boundary conditions will be transfered to the new mesh
  and all of the data will be written out into a vtk file. The outputs
  are:
    - stitched.vtk : only stitched mesh/grid
    - stitchedPhys.vtu : stitched mesh from input data plus all defined 
                         bcs and solution values.
    - optimized.vtk : the optimized/refined mesh.
    - optimizedPhys.vtu : optimized mesh/grid and all transfered soulution
                          and boundary condition data.
*/

// standard headers
#include <cstring>
#include <string.h>
#include <iostream>
#include <memory>

// Nemosys headers
#include "rocstarCgns.H"
#include "vtkAnalyzer.H"
#include "baseInterp.H"
#include "meshPartitioner.H"
#include "cgnsWriter.H"

// MAdLib + gmsh headers 
#include "ModelInterface.h"
#include "MAdLib.h"
#include "NodalDataManager.h"
#include "GmshEntities.h"
#include <meshBase.H>
// typedefs

/* auxiliary functions */
void getRegCenters(MAd::pMesh msh, std::vector<double>& regCntCrds);
void getFaceCenters(MAd::pMesh msh, std::vector<double>& faceCntCrds);

/*   Main Function */ 
int main(int argc, char* argv[])
{
  // check input
  // int nBCgFile = 0;
  // int nNICgFile = 0;
  
  std::string bCgFileName, niCgFileName;
  if (argc!=5) {
    std::cout << "Usage: " << argv[0] 
              << " n_b bCgFileName0 n_ni inCgFileName0" << std::endl;
    return 0;
  }
  // std::string::size_type sz;   // alias of size_t
  // nBCgFile = std::stoi(argv[1],&sz);
  // nNICgFile = std::stoi(argv[3],&sz);
  bCgFileName = argv[2];
  niCgFileName = argv[4];

  std::cout << "Reading input files #################################################\n";
  // reading cgns file1
  rocstarCgns* bCgObj = new rocstarCgns(bCgFileName);
  bCgObj->loadCgSeries(3);
  bCgObj->stitchGroup();
  // reading cgns file1
  rocstarCgns* niCgObj = new rocstarCgns(niCgFileName);
  niCgObj->loadCgSeries(4);
  niCgObj->stitchGroup();
  meshBase* stitched1 = meshBase::Create(bCgObj->getVTKMesh(), "bCgObjStitched.vtu");
  meshBase* stitched2 = meshBase::Create(niCgObj->getVTKMesh(), "niCgObjStitched.vtu");
  stitched1->write();
  stitched2->write();
  delete stitched1;
  delete stitched2;
  // stiching two files
  bCgObj->stitchMe(niCgObj); 
  
  std::cout << "Exporting mesh to MAdLib format #####################################\n";
  // exporting mesh to the MAdLib
  MAd::pGModel model = NULL;
  MAd::GM_create(&model,"");
  MAd::pMesh mesh = M_new(model);
  bCgObj->exportToMAdMesh(mesh);
  bCgObj->classifyMAdMeshOpt(mesh);
  
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
  bCgObj->getSolutionDataNames(slnNameList);  
  bCgObj->getAppendedSolutionDataName(appSlnNameList);
  slnNameList.insert(slnNameList.end(),
                     appSlnNameList.begin(), appSlnNameList.end());
  // write all data into vtk file
  for (auto is=slnNameList.begin(); is<slnNameList.end(); is++)
  {
    std::vector<double> physData;
    bCgObj->getSolutionDataStitched(*is, physData, outNData, outNDim);
    solution_type_t dt = bCgObj->getSolutionDataObj(*is)->getDataType();
    if (dt == NODAL)      
      trgVTK->setPointDataArray((*is).c_str(), 1, physData);
    else
      trgVTK->setCellDataArray((*is).c_str(), 1, physData);
  }
  trgVTK->report();
  trgVTK->write("stitchedPhys.vtu");

  
  std::cout << "Optimizing the mesh #############################################\n";
  rocstarCgns* cgObj1 = bCgObj;
  cgObj1->classifyMAdMeshBnd(mesh); // registering boundaries for the optimization step

  // perparing for element quantity interpolation
  std::vector<double> faceCntCrdsOld;
  getFaceCenters(mesh, faceCntCrdsOld);
  basicInterpolant* int1 = new basicInterpolant(3, M_numFaces(mesh), 1, faceCntCrdsOld);

  // prepare the mesh for the optimization
  std::cout << "Checking mesh sanity.\n";
  MAd::checkMesh(mesh);
  MAd::M_info(mesh);

  // defining addaptive refinement parameters
  MAd::PWLSField * sizeField = new MAd::PWLSField(mesh);
  sizeField->setCurrentSize();
  sizeField->scale(0.5);
  MAd::MeshAdapter* ma = new MAd::MeshAdapter(mesh,sizeField);
  ma->setCollapseOnBoundary(true, 1e-8);
  ma->setSwapOnBoundary(true, 1e-8);

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
  cgObj1->unclassifyMAdMeshBnd(mesh); // remove registered boundaries for proper output
  MAd::M_writeMsh (mesh, "optimizedMesh.msh", 2, NULL);
  MAd::M_info(mesh);


  // convert optimized Gmsh file to vtk
  std::cout << "Converting from optimized gmsh to vtk format.\n";
  trgGModel = new GModel("stitchedOpt"); 
  trgGModel->readMSH("optimizedMesh.msh");
  trgGModel->writeVTK("optimizedMesh.vtk", false, true);

  
  std::cout << "Transfering solution values #############################################\n";
  // write physical quantities to vtk file
  std::cout << "Writing transfered physical quantities to vtk file.\n";
  // get nodal data after refinement and write them
  vtkAnalyzer* trgVTK2;
  trgVTK2 = new vtkAnalyzer((char*)"optimizedMesh.vtk");
  trgVTK2->read();
  std::cout << "Writing cell and nodal data....\n";
  std::vector<double> faceCntCrdsNew;
  getFaceCenters(mesh, faceCntCrdsNew);
  int nNewElm = M_numFaces(mesh);
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
      std::cout << "Writing cell-based " << *is << std::endl;
      std::vector<double> oldPhysData;
      std::vector<double> newPhysData;
      int nDataT, nDimT;
      cgObj1->getSolutionDataStitched(*is, oldPhysData, nDataT, nDimT);
      int1->interpolate(nNewElm, faceCntCrdsNew, oldPhysData, newPhysData);
      trgVTK2->setCellDataArray((*is).c_str(), 1, newPhysData);      
    }
  }
  trgVTK2->report();
  trgVTK2->write("optimizedPhys.vtu");

  /*
  // partition the mesh
  //meshPartitioner* mPart = new meshPartitioner(cgObj1);
  meshPartitioner* mPart = new meshPartitioner(mesh);
  mPart->partition(nOutCgFile);
  std::vector<int> madConn = MAd::M_getConnectivities(mesh);
  std::vector<int> conn = mPart->getConns(0);
   
  // write CGNS files for the new grid
  for (int iCg=0; iCg<nOutCgFile; iCg++)
  {
     //std::ostringstream exp;
     //exp << "test_" << iCg << ".cgns";
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
   // int rCnt = 0;
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

/* get element/cell/region center coordinates */
void getFaceCenters(MAd::pMesh msh, std::vector<double>& faceCntCrds)
{
   MAd::FIter fi = M_faceIter(msh);
   // int fCnt = 0;
   while (MAd::pFace pf = FIter_next(fi)) 
   {
     double xc[3];
     MAd::F_center(pf, xc);
     faceCntCrds.push_back(xc[0]);
     faceCntCrds.push_back(xc[1]);
     faceCntCrds.push_back(xc[2]);
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
