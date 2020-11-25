/*
*/

// standard headers
#include <cstring>
#include <string.h>
#include <iostream>
#include <memory>
#include <math.h>

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
#include "MeshDataBaseIO.h"
#include <CheckOrientation.h>
#include "GmshEntities.h"

// typedefs

/* auxiliary functions */
void getRegCenters(MAd::pMesh msh, std::vector<double>& regCntCrds);
void getFaceCenters(MAd::pMesh msh, std::vector<double>& faceCntCrds);
void getVertxCoords(MAd::pMesh msh, std::vector<double>& vrtxCrds);
std::string rmPath(std::string fName);

void stlToMAdMesh(GModel* inModel, const MAd::pMesh outMesh)
{
  std::cout << "nVertex = " << inModel->getNumVertices() << std::endl;
  std::cout << "nFace = " << inModel->getNumFaces() << std::endl;
  std::cout << "nEdges = " << inModel->getNumEdges() << std::endl;
  std::cout << "nRegions = " << inModel->getNumRegions() << std::endl;
  std::cout << "nMeshElements = " << inModel->getNumMeshElements() << std::endl;
  std::cout << "nMeshVertices = " << inModel->getNumMeshVertices() << std::endl;

  // reading the vertices
  int iVert = 0;
  for (int iV=1; iV <= inModel->getNumMeshVertices(); iV++){
    iVert++;
    MVertex* vrt = inModel->getMeshVertexByTag(iV);
    //std::cout << iV << " " << vrt->x() << " " << vrt->y() << " " << vrt->z() << std::endl;
    outMesh->add_point(iV, vrt->x(), vrt->y(), vrt->z());
  }
  std::cout << "Read " << iVert << " vetices.\n";

  // reading the faces
  int nTri = 0;
  for (int iElm=1; iElm<=inModel->getNumMeshElements(); iElm++)
  {
    std::vector<MVertex*> verts;
    MElement* elm = inModel->getMeshElementByTag(iElm);
    if (elm->getType() == 3) {
      nTri++;
      MAd::pGEntity geom = (MAd::pGEntity) MAd::GM_faceByTag(outMesh->model, 0);
      elm->getVertices(verts);
      outMesh->add_triangle(verts[0]->getNum(), verts[1]->getNum(),
    		   verts[2]->getNum(), geom); 
    }
  }
  std::cout << "Read " << nTri << " triangles.\n";
  
  // classifying the mesh to make it ready for optimization
  // finding boundary faces and classifying them as 2 dimensional
  outMesh->classify_unclassified_entities();
  outMesh->destroyStandAloneEntities();
  MAd::pGEntity bnd = (MAd::pGEntity) MAd::GM_faceByTag(outMesh->model, 0);
  outMesh->classify_grid_boundaries(bnd);
}

/*   Main Function */ 
int main(int argc, char* argv[])
{
  
  std::string mshFileName;
  if (argc!=2) {
    std::cout << "Usage: " << argv[0] 
              << " mshFileName" << std::endl;
    return 0;
  }
  mshFileName = argv[1];

  // Reading gmsh file and info
  MAd::pGModel model = NULL;
  MAd::GM_create(&model,"initial");
  MAd::pMesh mesh = MAd::M_new(model);
  MAd::M_load(mesh, mshFileName.c_str());
  MAd::checkMesh(mesh);
  MAd::M_info(mesh);

  // convert to vtk
  GModel* srcGModel;
  srcGModel = new GModel("srcMsh"); 
  srcGModel->readMSH(mshFileName);
  srcGModel->writeVTK("mesh.vtk", false, true);

  // perparing for elemental quantity interpolation
  std::vector<double> regCntCrdsOld;
  getRegCenters(mesh, regCntCrdsOld);
  /*basicInterpolant* intCell = new */basicInterpolant(3, M_numRegions(mesh), 1, regCntCrdsOld);

  // perparing for nodal quantity interpolation
  std::vector<double> vrtxCrdsOld;
  getVertxCoords(mesh, vrtxCrdsOld);
  /*basicInterpolant* intNde = new */basicInterpolant(3, M_numVertices(mesh), 4, vrtxCrdsOld);

  // defining addaptive refinement parameters
  MAd::PWLSField * sizeField = new MAd::PWLSField(mesh);
  sizeField->setCurrentSize();
  sizeField->scale(1.0);
  MAd::MeshAdapter* ma = new MAd::MeshAdapter(mesh,sizeField);
  //ma->setSliverPermissionInESplit(false);
  //ma->setSliverPermissionInECollapse(false);
  //ma->setSliverQuality(0.4);
  ma->printParameters();

  // register some data on the mesh
  std::vector<double> physData;
  for (int iVrt=0; iVrt < vrtxCrdsOld.size(); iVrt+=3)
  {
    physData.push_back(std::sqrt(std::pow(vrtxCrdsOld[iVrt],2) 
                       + std::pow(vrtxCrdsOld[iVrt+1],2)
                       + std::pow(vrtxCrdsOld[iVrt+2],2)));
  }
  std::cout << "Num vertices = " << MAd::M_numVertices(mesh) << std::endl;
  std::cout << "Num physData = " << physData.size() << std::endl;
  ma->registerData("testData", physData);

  // optimize the mesh
  //std::cout << "Statistics before optimization: \n";
  //ma->printStatistics(std::cout);
  std::cout << "Optimizing the mesh ...\n";
  ma->run();
  ma->removeSlivers();
  //std::cout << "Statistics after optimization: \n";
  //ma->printStatistics(std::cout);

  // write final mesh status into a file
  //MAd::CheckFacesOrientation(mesh);

  // writing new mesh and convert to vtk
  mesh->unclassify_grid_boundaries();
  MAd::M_info(mesh);
  MAd::M_writeMsh(mesh, "transfered.msh");
  srcGModel->deleteMesh();
  srcGModel->readMSH("transfered.msh");
  srcGModel->writeVTK("transfered.vtk", false, true);

  // write physical quantities to vtk file
  std::cout << "Writing physical quantities to vtk file.\n";
  vtkAnalyzer* trgVTK;
  trgVTK = new vtkAnalyzer((char*)"transfered.vtk");
  trgVTK->read();
  trgVTK->report();

  // write all data into vtk file
  physData.clear();
  ma->getMeshData("testData", &physData);
  trgVTK->setPointDataArray("testData", 1, physData);
  trgVTK->report();
  trgVTK->write("tansferedPhys.vtu");

  /*
  // check input
  int nBCgFile = 0;
  int nNICgFile = 0;
  
  std::string bCgFileName, niCgFileName;
  if (argc!=5) {
    std::cout << "Usage: " << argv[0] 
              << " n_b bCgFileName0 n_ni inCgFileName0" << std::endl;
    return 0;
  }
  std::string::size_type sz;   // alias of size_t
  nBCgFile = std::stoi(argv[1],&sz);
  nNICgFile = std::stoi(argv[3],&sz);
  bCgFileName = argv[2];
  niCgFileName = argv[4];

  std::cout << "Reading input files #################################################\n";
  // reading cgns file1
  rocstarCgns* bCgObj = new rocstarCgns(bCgFileName);
  bCgObj->loadCgSeries(3);
  bCgObj->dummy();
  // reading cgns file1
  rocstarCgns* niCgObj = new rocstarCgns(niCgFileName);
  niCgObj->loadCgSeries(4);
  niCgObj->dummy();

  std::cout << "Stitching surface meshes ############################################\n";
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
  GModel* srcGModel;
  srcGModel = new GModel("stitched"); 
  srcGModel->readMSH("stitched.msh");
  srcGModel->writeVTK("stitched.vtk", false, true);

  
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

  // perparing for elemental quantity interpolation
  std::vector<double> faceCntCrdsOld;
  getFaceCenters(mesh, faceCntCrdsOld);
  basicInterpolant* intCell = new basicInterpolant(3, M_numFaces(mesh), 1, faceCntCrdsOld);
  // perparing for nodal quantity interpolation
  std::vector<double> vrtxCrdsOld;
  getVertxCoords(mesh, vrtxCrdsOld);
  basicInterpolant* intNde = new basicInterpolant(3, M_numVertices(mesh), 4, vrtxCrdsOld);

  // prepare the mesh for the optimization
  std::cout << "Checking mesh sanity.\n";
  MAd::checkMesh(mesh);
  MAd::M_info(mesh);

  // attach all nodal data to the mesh optimizer
  int nData, nDim;
  std::vector<std::string> cgSlnNameList;
  std::vector<std::string> cgAppSlnNameList;
  cgObj1->getSolutionDataNames(cgSlnNameList);  
  cgObj1->getAppendedSolutionDataName(cgAppSlnNameList);
  cgSlnNameList.insert(cgSlnNameList.end(),
                     cgAppSlnNameList.begin(), cgAppSlnNameList.end());

  // reading skin mesh
  std::cout << "Reading skinned mesh.\n";
  MAd::pGModel skinMdl = NULL;
  MAd::GM_create(&skinMdl,"");
  MAd::pMesh skinMsh = M_new(skinMdl);
  skinMsh->classify_unclassified_entities();
  skinMsh->destroyStandAloneEntities();
  MAd::LoadGmshMesh(skinMsh, "skinMesh.msh");
  MAd::M_info(skinMsh);

  // reading skinned mesh partition information
  std::ifstream prtFile ("skinMeshPart.dat");
  std::vector<double> elmPartIds;
  std::string line;
  if (prtFile.is_open())
  {
    while (getline(prtFile, line))
    {
      elmPartIds.push_back(std::stoi(line));
    }
  } else {
    std::cerr << "Can not file surface element partion id file skinMeshPart.dat\n";
    exit(-1);
  }
  prtFile.close();
  if (elmPartIds.size() != M_numFaces(skinMsh))
  {
    std::cout << "Incompatable partion id and skinned mesh...Aborting.\n";
    exit(-1);
  }

  // convert optimized (skinned) Gmsh file to vtk
  std::cout << "Converting from optimized (skinned) gmsh to vtk format.\n";
  srcGModel = new GModel("stitchedOpt"); 
  srcGModel->readMSH("skinMesh.msh");
  srcGModel->writeVTK("optimizedMesh.vtk", false, true);

  // quantity transfer between meshes 
  std::cout << "Transfering solution values #############################################\n";
  std::cout << "Writing transfered physical quantities to vtk file.\n";
  vtkAnalyzer* trgVTK2;
  trgVTK2 = new vtkAnalyzer((char*)"optimizedMesh.vtk");
  trgVTK2->read();
  std::cout << "Writing cell and nodal data....\n";
  std::vector<double> faceCntCrdsNew;
  getFaceCenters(skinMsh, faceCntCrdsNew);
  int nNewElm = M_numFaces(skinMsh);
  std::vector<double> vrtxCrdsNew;
  getVertxCoords(skinMsh, vrtxCrdsNew);
  int nNewVrtx = M_numVertices(skinMsh);
  for (auto is=cgSlnNameList.begin(); is<cgSlnNameList.end(); is++)
  {
    solution_type_t dt = cgObj1->getSolutionDataObj(*is)->getDataType();
    if (dt == NODAL) {
      std::cout << "Writing nodal " << *is << std::endl;
      // using MAd mesh to transfer
      //std::vector<double> physData;
      //ma->getMeshData((*is), &physData);
      //trgVTK2->setPointDataArray((*is).c_str(), 1, physData);
      // Using built-in interpolant to transfer
      std::cout << "... Using built-in transfer.\n";
      std::vector<double> oldPhysData;
      std::vector<double> newPhysData;
      int nDataT, nDimT;
      cgObj1->getSolutionDataStitched(*is, oldPhysData, nDataT, nDimT);
      intNde->interpolate(nNewVrtx, vrtxCrdsNew, oldPhysData, newPhysData);
      trgVTK2->setPointDataArray((*is).c_str(), 1, newPhysData);
      //MAd::NodalDataManagerSgl::instance().writeData((*is),((*is)+".pos").c_str());
    } else {      
      std::cout << "Writing cell-based " << *is << std::endl;
      std::vector<double> oldPhysData;
      std::vector<double> newPhysData;
      int nDataT, nDimT;
      cgObj1->getSolutionDataStitched(*is, oldPhysData, nDataT, nDimT);
      intCell->interpolate(nNewElm, faceCntCrdsNew, oldPhysData, newPhysData);
      trgVTK2->setCellDataArray((*is).c_str(), 1, newPhysData);      
    }
  }
  // DEBUG
  trgVTK2->setCellDataArray("partition", 1, elmPartIds);      
  trgVTK2->report();
  trgVTK2->write("optimizedPhys.vtu");

  
  // partition the mesh
  meshPartitioner* mPart = new meshPartitioner(skinMsh);
  mPart->setPartedElm(elmPartIds);
   
  std::cout << "Writing CGNS Files  #############################################\n";
  // write CGNS files for the new grid
  cgObj1 = niCgObj;
  int nOutCgFile = 3;
  for (int iCg=0; iCg<nOutCgFile; iCg++)
  {
     //std::ostringstream exp;
     //exp << "test_" << iCg << ".cgns";
     std::string fCgName;
     fCgName =rmPath(cgObj1->getCgFName(iCg));
     std::cout << "Writing remeshed " << fCgName << std::endl;
     // define elementary information
     cgnsWriter* cgWrtObj = new cgnsWriter(fCgName, cgObj1->getBaseName(iCg), 3, 3);
     cgWrtObj->setUnits(cgObj1->getMassUnit(), cgObj1->getLengthUnit(),
			cgObj1->getTimeUnit(), cgObj1->getTemperatureUnit(),
			cgObj1->getAngleUnit());
     cgWrtObj->setBaseItrData(cgObj1->getBaseItrName(iCg), cgObj1->getNTStep(iCg), cgObj1->getTimeStep(iCg));
     for (int iZn=1; iZn<=cgObj1->getNZone(iCg); iZn++)
     {
	cgWrtObj->setZoneItrData(cgObj1->getZoneItrName(iCg,iZn), cgObj1->getGridCrdPntr(iCg,iZn), cgObj1->getSolutionPntr(iCg,iZn));
	cgWrtObj->setZone(cgObj1->getZoneName(iCg, iZn), cgObj1->getZoneType(iCg, iZn));
	cgWrtObj->setNVrtx(mPart->getNNdePart(iCg));
	cgWrtObj->setNCell(mPart->getNElmPart(iCg));
	// define coordinates
	cgWrtObj->setGridXYZ(mPart->getCrds(iCg, MAd::M_getVrtXCrds(skinMsh)), 
			     mPart->getCrds(iCg, MAd::M_getVrtYCrds(skinMsh)), 
			     mPart->getCrds(iCg, MAd::M_getVrtZCrds(skinMsh)));
	// define connctivity
	std::cout << "sectionName = " << cgObj1->getSectionName(iCg, iZn) << std::endl;
	std::cout << "elementType = " << cgObj1->getElementType(iCg, iZn) << std::endl;
	cgWrtObj->setSection(cgObj1->getSectionName(iCg, iZn), 
			     (ElementType_t) cgObj1->getElementType(iCg, iZn), 
			     mPart->getConns(iCg));
	// define vertex and cell data 
	std::map<std::string, GridLocation_t> slnNLMap = cgObj1->getSolutionNameLocMap();
	for (auto is=slnNLMap.begin(); is!=slnNLMap.end(); is++)
	  cgWrtObj->setSolutionNode(is->first, is->second);
	// write skelleton of the file
        if (iZn == 1)
	  cgWrtObj->writeGridToFile();
        else
	  cgWrtObj->writeZoneToFile();
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
 
/* get element/cell/region center coordinates */
void getVertxCoords(MAd::pMesh msh, std::vector<double>& vrtxCrds)
{
   MAd::VIter vi = M_vertexIter(msh);
   while (MAd::pVertex pv = VIter_next(vi)) 
   {
     double xc[3];
     MAd::V_coord(pv, xc);
     vrtxCrds.push_back(xc[0]);
     vrtxCrds.push_back(xc[1]);
     vrtxCrds.push_back(xc[2]);
     /*
     std::cout << "Vertex coordinate = "
               << xc[0] << " "
               << xc[1] << " "
               << xc[2] << " "
               << std::endl;
     */
   }
}


std::string rmPath(std::string fName)
{
  std::size_t _loc = fName.find_last_of("/");
  return(fName.substr(_loc+1,fName.size()));
}
