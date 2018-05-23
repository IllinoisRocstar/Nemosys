// standard headers
#include <cstring>
#include <string.h>
#include <iostream>
#include <memory>

// Nemosys headers
#include <cgnsAnalyzer.H>
#include <meshBase.H>
#include <symmxGen.H>
#include <meshPartitioner.H>
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

  std::cout << "Reading input files #################################################\n";
  // reading cgns file1
  cgnsAnalyzer* cgObj1 = new cgnsAnalyzer(cgFileName[0]);
  cgObj1->loadGrid();
  // attaching partion id to the mesh
  std::vector<double> slnData(cgObj1->getNElement(),0);
  cgObj1->appendSolutionData("partitionOld", slnData, ELEMENTAL, cgObj1->getNElement(), 1);
  std::vector<cgnsAnalyzer*> cgObjs(nInCgFile-1);
  // adding new cgns files
  for (int iCg=1; iCg < nInCgFile; iCg++) 
  {
    cgObjs[iCg-1] = new cgnsAnalyzer(cgFileName[iCg]);
    cgObjs[iCg-1]->loadGrid();
    // appending data
    std::vector<double> slnData(cgObjs[iCg-1]->getNElement(),iCg);
    cgObjs[iCg-1]->
      appendSolutionData("partitionOld", slnData, ELEMENTAL, cgObjs[iCg-1]->getNElement(), 1);
    // stitching meshes
    cgObj1->stitchMesh(cgObjs[iCg-1], true);
  }
  if (nInCgFile > 1)
     std::cout << "Meshes stitched successfully!\n";
  
  std::cout << "Exporting mesh to VTK format ########################################\n";
  meshBase* trgVTK = meshBase::Create(cgObj1->getVTKMesh(),"stitched.vtu");
  std::cout << "Transferring physical quantities to vtk mesh ########################\n";
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
    {      
      std::cout << "Writing nodal " << *is << std::endl; 
      trgVTK->setPointDataArray((*is).c_str(), physData);
    }
    else
    {
      // gs field is 'weird' in irocstar files- we don't write it back
      if (!(*is).compare("gs"))
        continue;
      std::cout << "Writing cell-based " << *is << std::endl;
      trgVTK->setCellDataArray((*is).c_str(), physData);
    }
  }
  
  trgVTK->report();
  trgVTK->write();

  ///////////////////// REMESH WITH SIMMETRIX //////////////////////////////////////////

  std::cout << "Extracting surface mesh #############################################\n";
  meshBase* trgVTKsurf = meshBase::Create(trgVTK->extractSurface(), "skinMeshPart.vtp"); 
 
  // remeshing with simmetrix
  trgVTKsurf->write("skinMeshPart.stl");
  SymmxParams* params = new SymmxParams(); 
  params->logFName = "symmxGen.log";
  params->licFName = "/home/users/snatesh/snatesh-storage/SYMMETRIX/LICENSE/simmodsuite.lic";
  params->features = "geomsim_core,meshsim_surface,meshsim_volume,geomsim_discrete";
  meshBase* newMesh = meshBase::generateMesh("skinMeshPart.stl","simmetrix",params);
  trgVTK->setContBool(0);
  //trgVTK->transfer(newMesh,"Finite Element"); 

  ///////////////////// PARTITION THE MESH WITH METIS /////////////////////////////////

  std::cout << " Partitioning the mesh with METIS.\n";
  meshPartitioner* mPart = new meshPartitioner(newMesh);
  //meshPartitioner* mPart = new meshPartitioner(cgObj1);
  mPart->partition(nOutCgFile);
  
  //// write partitin ids for the surface mesh
  //std::ofstream of;
  //of.open("skinMeshPart.dat", std::fstream::trunc);
  //std::vector<double> elmPartIds = mPart->getPartedElm();
  //for (auto it=skinElmIds.begin(); it!=skinElmIds.end(); it++)
  //{
  //  of << elmPartIds[*it - 1] << "\n";
  //} 
  //of.close();

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
    trgVTK->transfer(tmpvtk, "Finite Element");
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
 
  if (trgVTK) delete trgVTK;
  if (trgVTKsurf) delete trgVTKsurf; 
  if (cgObj1) delete cgObj1;
  if (params) delete params;
  if (mPart) delete mPart;
  for (int i = 0; i < cgObjs.size(); ++i)
    if (cgObjs[i]) delete cgObjs[i];
  
  std::cout << "Application ended successfully!\n";

  return 0;
}
