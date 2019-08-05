// class headers for pack mesh workflow
#include "cfmeshGen.H"
#include "cfmeshParams.H"
#include "snappymeshGen.H"
#include "snappymeshParams.H"
#include "blockMeshGen.H"
#include "blockMeshParams.H"
#include "MeshManipulationFoam.H"
#include "MeshManipulationFoamParams.H"
#include "PackMeshDriver.H"
#include "meshBase.H"
#include "MeshQualityDriver.H"
#include "MeshQuality.H"
#include <boost/filesystem.hpp>

// std c++
#include <vector>
#include <tuple>
#include <memory>
#include <iostream>

// vtk
#include "ConversionDriver.H"
#include <AuxiliaryFunctions.H>
#include <vtkMesh.H>
#include <vtkCellData.h>
#include <vtkIdList.h>
#include <vtkIdList.h>
#include <vtkCellData.h>

// TODO
// 1. Some smart implementation to improve mesh if mesh is not of desired
//    quality.

// ------------------------- Pack Mesh Driver -------------------------------//

// Standard constructor controls the pack meshing process.
PackMeshDriver::PackMeshDriver(const std::string& ifname,
                              MeshManipulationFoamParams* _mparams,
                              cfmeshParams* _cfparams,
                              snappymeshParams* _snappyparams,
                              blockMeshParams* _bmparams,
                              const std::string& ofname1,
                              const std::string& ofname2)
{

  // Makes sure that constant and triSurface directories are present.
  // If directories are already present, it will not do anything.
  const char dir_path[] = "./constant";
  boost::filesystem::path dir(dir_path);
  try
  {
    boost::filesystem::create_directory(dir);
  }
  catch (boost::filesystem::filesystem_error &e)
  {
    std::cerr << "Problem in creating triSurface directory"
              << "for the snappyHexMesh" << "\n";
    std::cerr << e.what() << std::endl;
    throw;
  }

  const char dir_path1[] = "./constant/triSurface";
  boost::filesystem::path dir1(dir_path1);
  try
  {
    boost::filesystem::create_directory(dir1);
  }
  catch (boost::filesystem::filesystem_error &e)
  {
    std::cerr << "Problem in creating triSurface directory"
              << "for the snappyHexMesh" << "\n";
    std::cerr << e.what() << std::endl;
    throw;
  }


  // This object has access to all MeshManipulation utilities.
  MeshManipulationFoam* objMsh = new MeshManipulationFoam(_mparams);
  const char* nameFile = "a"; // Dummy name for input


  // RocPack output STL file lists all packs under one surface name. This step
  // renumbers STL file into n distinct surfaces and outputs number of packs
  // for mergeMeshes loop. It will overwrite the original input file
  objMsh->createControlDict();
  int nDomains = objMsh->surfSpltByTopology();

  // Second step in process is to create hexahedral mesh of packs using
  // cartesianMesh utility of CfMesh. It writes the mesh in constant/polyMesh.
  cfmeshGen* objCF = new cfmeshGen(_cfparams);
  objCF->createMeshFromSTL(nameFile);
  
  // foamToSurface utility reads mesh from constant/polyMesh, extracts the
  // surface mesh of packs and outputs them as STL file in constant/triSurface
  // directory. This output file will be used by snappyHexMesh later.
  objMsh->foamToSurface();

  // blockMesh utility takes user input for surrounding box region and 
  // generates mesh block in constant/polyMesh folder. It will overwrite the
  // previous mesh created by CfMesh. This mesh will be used as background mesh
  // by snappyHexMesh later.
  blockMeshGen* objBM = new blockMeshGen(_bmparams);
  objBM->createMeshFromSTL(nameFile);

  // snappyHexMesh reads background mesh created using blockMesh and takes
  // surface file from foamToSurface utility to snap pack surface onto back-
  // ground mesh and creates different cellZones (i.e different solids). These
  // interfaces between pack and surrounding regions are completely conformal
  // due to snappyHexMesh's unique snapping abilities.
  snappymeshGen* objSHM = new snappymeshGen(_snappyparams);
  objSHM->createMeshFromSTL(nameFile);

  // splitMeshRegions reads mesh from constant/polyMesh directory and splits
  // mesh into separate regions. Each region will represent one solid. It will
  // write bunch of domain.* directories inside constant/ and system/ folders.
  // This function does a random walking around mesh to identify different
  // regions and in that process, while naming all domains as domain.*
  // in constant folder, it skips one number for the disconnected region it
  // encounters first. This number is taken out to provide as input to merge
  // mesh. 
  int dirStat = objMsh->splitMshRegions();

  // mergeMeshes will read master domain (defined by user) from constant folder
  // and start merging other domain to it untill all slave domains are attached
  // to master domain. It loops through all domains in sequential manner and 
  // skips the missing domain.* (also skipped by splitMeshRegions) to avoid
  // runtime error.
  objMsh->mergeMeshes(dirStat, nDomains);

  // createPatch utility reads domain mesh and createPatchDict to combine all
  // different patches of multiple packs/surrounding into one patch
  // respectively
  objMsh->createPatch(dirStat);

  //Reads current mesh and write it to separate VTK/VTU files
  bool readDB = false;

  // Reads and converts pack mesh
  std::string regNme;
  if (dirStat == 1)
  	regNme = "domain2";
  else
  	regNme = "domain1";
  
  meshBase* fm = new FOAM::foamMesh(readDB);
  fm->read(regNme);
  vtkMesh* vm = new vtkMesh(fm->getDataSet(),ofname1);
  vm->report();
  vm->write();

  // Reads and converts surronding mesh
  regNme = "domain0";
  meshBase* fm2 = new FOAM::foamMesh(readDB);
  fm2->read(regNme);
  vtkMesh* vm2 = new vtkMesh(fm2->getDataSet(),ofname2);
  vm2->report();
  vm2->write();

  // Outputs useful mesh quality parameters for users
  std::string SurroundingName = "surroundingMeshQuality";
  std::string Surroundingmesh = "geom_surrounding_mesh.vtu";
  std::string PackName = "packMeshQuality";
  std::string Packmesh = "geom_pack_mesh.vtu";
  MeshQualityDriver* objSurrQ =
                new MeshQualityDriver(Surroundingmesh, SurroundingName);
  MeshQualityDriver* objPackQ =
                new MeshQualityDriver(Packmesh, PackName);

  // Cleaning up
  if (vm)
    delete vm;
  if (fm)
    delete fm;
  if (vm2)
    delete vm2;
  if (fm2)
    delete fm2;
  if (objMsh)
    delete objMsh;
  if (objSHM)
    delete objSHM;
  if (objCF)
    delete objCF;
  if (objBM)
    delete objBM;
  if (objSurrQ)
    delete objSurrQ;
  if (objPackQ)
    delete objPackQ;

  // End of workflow

  // Adds cohesive elements for hex and writes them in VTU file.
  //objMsh->addCohesiveElements(1e-13, "FinalMesh.vtu");

  //if (objMsh)
  //  delete objMsh;

  // Adds cohesive elements with thickness for hex and writes them in VTU file.
  //objMsh->addArtificialThicknessElements(1e-13, "FinalMesh.vtu");

  //if (objMsh)
  //  delete objMsh;
  
}

// Class destructor
PackMeshDriver::~PackMeshDriver()
{
  std::cout << "PackMeshDriver destroyed" << endl;
}


// Reads JSON input provided by user and json file passes down by NemDriver.
PackMeshDriver* PackMeshDriver::readJSON(const jsoncons::json inputjson)
{
  std::string ifname = inputjson["Mesh File Options"]
                          ["Input Geometry File"].as<std::string>();
  std::string ofname1 = inputjson["Mesh File Options"]
                          ["Output Pack Mesh File"].as<std::string>();
  std::string ofname2 = inputjson["Mesh File Options"]
                          ["Output Surrounding Mesh File"].as<std::string>();
  return readJSON(ifname, ofname1, ofname2, inputjson);
}

// Reads JSON file parameters and catagorizes them into separate objects
// for convenience in Pack Mesh workflow.
PackMeshDriver* PackMeshDriver::readJSON(const std::string& ifname, 
                                        const std::string& ofname1,
                                        const std::string& ofname2,
                                        const jsoncons::json inputjson)
{
  std::string meshEngine = 
      inputjson["Mesh Generation Engine"].as<std::string>();

  if (!meshEngine.compare("packmesh"))
  {
    cfmeshParams* cfparams = new cfmeshParams();
    std::string defaults1 = 
        inputjson["Meshing Parameters"]["CFMesh Parameters"].as<std::string>();

    jsoncons::json cfmparams = inputjson["Meshing Parameters"]["CFMesh Parameters"];

    // required params here
    // cad file
    if (inputjson["Mesh File Options"].has_key("Input Geometry File"))
      cfparams->geomFilePath = 
      inputjson["Mesh File Options"]["Input Geometry File"].as<std::string>();
    else
    {
      std::cerr << "A geometry file should be supplied.\n";
      throw;        
    }
        
    // mesh generator
    if (cfmparams.has_key("Generator"))
      cfparams->generator = cfmparams["Generator"].as<std::string>();
    else
    {
      std::cerr << "A mesh generation method should be selected.\n";
      std::cerr << "Options: cartesian2D tetMesh\n";
      throw;        
    }
        
    // rest of params are optional
    if (cfmparams.has_key("MaxCellSize"))
      cfparams->maxCellSize = 
        cfmparams["MaxCellSize"].as<double>();
    if (cfmparams.has_key("MinCellSize"))
      cfparams->minCellSize = 
        cfmparams["MinCellSize"].as<double>();
    if (cfmparams.has_key("BoundaryCellSize"))
      cfparams->bndryCellSize = 
        cfmparams["BoundaryCellSize"].as<double>();
    if (cfmparams.has_key("KeepCellsIntersectingBoundary"))
      cfparams->keepCellIB = 
        cfmparams["KeepCellsIntersectingBoundary"].as<double>();
    if (cfmparams.has_key("CheckForGluedMesh"))
      cfparams->chkGluMsh = 
        cfmparams["CheckForGluedMesh"].as<double>();
    if (cfmparams.has_key("AllowDisconnectedDomains"))
      cfparams->_alwDiscDomains = 
        cfmparams["AllowDisconnectedDomains"].as<bool>();

    // optional capability
    std::string cap = "BoundaryLayers";
    if (cfmparams.has_key(cap))
    {
      cfparams->_withBndLyr = true;
      cfparams->blNLyr = cfmparams[cap]["NLayers"].as<double>();
      cfparams->blThkRto = cfmparams[cap]["ThicknessRatio"].as<double>();

      if (cfmparams[cap].has_key("MaxFirstLayerThickness"))
        cfparams->maxFrstLyrThk =
          cfmparams[cap]["MaxFirstLayerThickness"].as<double>();
      if (cfmparams[cap].has_key("AllowDiscontinuity"))
        cfparams->alwDiscont =
          cfmparams[cap]["AllowDiscontinuity"].as<bool>();

      // patch boundary layers
      std::string subcap = "PatchBoundaryLayers";
      if (cfmparams[cap].has_key(subcap))
      {
        cfparams->_withBndLyrPtch = true;
        for (auto jptch : cfmparams[cap][subcap].array_range())
        {
          cfmPtchBndLyr blPatch;
          blPatch.patchName = jptch["PatchName"].as<std::string>();
          if (jptch.has_key("AllowDiscontinuity"))
            blPatch.alwDiscont = jptch["AllowDiscontinuity"].as<bool>();
          else
            blPatch.alwDiscont = false;
          if (jptch.has_key("MaxFirstLayerThickness"))
            blPatch.maxFrstLyrThk =
            jptch["MaxFirstLayerThickness"].as<int>();
          else
            blPatch.maxFrstLyrThk = -1;
          if (jptch.has_key("NLayers"))
            blPatch.blNLyr = jptch["NLayers"].as<int>();
          else
            blPatch.blNLyr = -1;
          if (jptch.has_key("ThicknessRatio"))
            blPatch.blThkRto = jptch["ThicknessRatio"].as<double>();
          else
            blPatch.blThkRto = -1.;

          (cfparams->blPatches).push_back(blPatch);
        }                          
      }
    }

    // optional capability
    cap = "SurfaceFeatureEdges";
    if (cfmparams.has_key(cap))
    {
      cfparams->_withSrfEdg = true;
      cfparams->srfEdgAng = cfmparams[cap]["Angle"].as<double>();
    }

    // optional capability
    cap = "ObjectRefinements";
    if (cfmparams.has_key(cap))
    {
      cfparams->_withObjRfn = true;
      for (auto refObj : cfmparams[cap].array_range())
      {
        cfmObjRef objRef;
        objRef.name= refObj["Name"].as<std::string>();
        for (const auto&  prm: refObj["Params"].object_range())
        {
          std::string key = std::string(prm.key());
          std::string val = prm.value().as<std::string>();
          objRef.params[key] = val;
        }
        (cfparams->objRefLst).push_back(objRef);
      }
    }

    // optional capability
    cap = "ImproveMeshQuality";
    if (cfmparams.has_key(cap))
    {
      cfparams->_withMshQlt = true;
      cfparams->qltNItr = cfmparams[cap]["NIterations"].as<int>();
      cfparams->qltNLop = cfmparams[cap]["NLoops"].as<int>();
      cfparams->qltQltThr = cfmparams[cap]["QualityThreshold"].as<double>();
      cfparams->qltNSrfItr = cfmparams[cap]["NSurfaceIterations"].as<int>();
      cfparams->qltConCelSet = 
      cfmparams[cap].get_with_default("ConstrainedCellsSet","none");
    }

    // optional capability
    cap = "LocalRefinement";
    if (cfmparams.has_key(cap))
    {
      cfparams->_withLclRef = true;
      for (auto jptch : cfmparams[cap].array_range())
      {
        cfmLclRefPatch refPatch;
        refPatch.patchName = jptch["PatchName"].as<std::string>();
        if (jptch.has_key("AdditionalRefinementLevels"))
          refPatch.aditRefLvls = jptch["AdditionalRefinementLevels"].as<int>();
        else
          refPatch.aditRefLvls = -1;
        if (jptch.has_key("CellSize"))
          refPatch.cellSize = jptch["CellSize"].as<double>();
        else
          refPatch.cellSize = -1.;

        (cfparams->refPatches).push_back(refPatch);
      }
    }

    // optional capability
    cap = "RenameBoundary";
    if (cfmparams.has_key(cap))
    {
      cfparams->_withRenBndry = true;
      cfmRenBndry renBndry;
      renBndry.defName = cfmparams[cap]["DefaultName"].as<std::string>();
      renBndry.defType = cfmparams[cap]["DefaultType"].as<std::string>();
      for (auto jnw : cfmparams[cap]["NewPatchNames"].array_range())
      {
        cfmNewPatch nwPatch = std::make_tuple( 
        jnw["Name"].as<std::string>(),
        jnw["NewName"].as<std::string>(),
        jnw["NewType"].as<std::string>() );
        renBndry.newPatches.push_back(nwPatch);
      }
      (cfparams->renBndry) = renBndry;
    }


    //Snappy
    snappymeshParams* snappyparams = new snappymeshParams();
    std::string defaults2 = 
    inputjson["Meshing Parameters"]
        ["snappyHexMesh Parameters"].as<std::string>();

    jsoncons::json shmparams =
        inputjson["Meshing Parameters"]["snappyHexMesh Parameters"];

    if (inputjson["Mesh File Options"].has_key("Input Geometry File"))
      snappyparams->geomFileName = 
      inputjson["Mesh File Options"]["Input Geometry File"].as<std::string>();
    else
    {
      std::cerr << "A geometry file should be supplied.\n";
      throw;        
    }

    if (shmparams.has_key("InputPatchName"))
      snappyparams->geomPatchName =
          shmparams["InputPatchName"].as<std::string>();
    else
    {
      std::cerr << "A patch name for input geometry must be defined.\n";
      throw;        
    }

    if (shmparams.has_key("SurfPatchName"))
      snappyparams->surfRefPatch =
          shmparams["SurfPatchName"].as<std::string>();
    else
    {
      std::cerr << "A patch name for surface refinement must be defined.\n";
      throw;        
    }

    if (shmparams.has_key("Castellated Mesh"))
        snappyparams->_withCastMesh =
              shmparams["Castellated Mesh"].as<bool>();
    else{
      std::cerr << "Please specify on/off choice using \"Castellated Mesh\""
                << "\n" << std::endl;
                throw;
    }
    if (shmparams.has_key("Snapping"))
      snappyparams->_withSnap =
        shmparams["Snapping"].as<bool>();
    else{
      std::cerr << "Please specify on/off choice using \"Snapping\""
                << "Keyword!\n" << std::endl;
                throw;
    }
    if (shmparams.has_key("Layer Addition"))
      snappyparams->_withLayers =
        shmparams["Layer Addition"].as<bool>();
    else{
      std::cerr << "Please specify on/off choice using \"Layer Addition\""
                << "Keyword!\n" << std::endl;
                throw;
    }
    if (shmparams.has_key("CellZones"))
      snappyparams->_withCellZones =
        shmparams["CellZones"].as<bool>();
    else{
      std::cerr << "Please specify on/off choice using \"CellZones\""
            << "Keyword!\n" << std::endl;
            throw;
    }
    if (shmparams.has_key("RegionRefine"))
      snappyparams->_withGeomRefReg =
        shmparams["RegionRefine"].as<bool>();
    else{
      std::cerr << "Please specify on/off choice using \"RegionRefine\""
            << "Keyword!\n" << std::endl;
            throw;
    }
    if (shmparams.has_key("SurfaceRefine"))
      snappyparams->_withSurfRefReg =
        shmparams["SurfaceRefine"].as<bool>();
    else{
      std::cerr << "Please specify on/off choice using \"SurfaceRefine\""
            << "Keyword!\n" << std::endl;
            throw;
    }
    if (shmparams.has_key("maxLocalCells"))
      snappyparams->maxLCells =
        shmparams["maxLocalCells"].as<int>();
    else{
      snappyparams->maxLCells = 2000000;
    }
    if (shmparams.has_key("maxGlobalCells"))
      snappyparams->maxGCells =
        shmparams["maxGlobalCells"].as<int>();
    else{
      snappyparams->maxGCells = 4000000;
    }
    if (shmparams.has_key("minRefCells"))
      snappyparams->minRefCells =
        shmparams["minRefCells"].as<int>();
    else{
      snappyparams->minRefCells = 0;
    }
    if (shmparams.has_key("nCellsBetweenLevels"))
      snappyparams->cellsBetnLvls =
        shmparams["nCellsBetweenLevels"].as<int>();
    else{
      snappyparams->cellsBetnLvls = 3;
    }
    if (shmparams.has_key("surfaceRefinementLvlMin"))
      snappyparams->refSurfLvlMin =
        shmparams["surfaceRefinementLvlMin"].as<int>();
    else{
      snappyparams->refSurfLvlMin = 0;
    }
    if (shmparams.has_key("surfaceRefinementLvlMax"))
      snappyparams->refSurfLvlMax =
        shmparams["surfaceRefinementLvlMax"].as<int>();
    else{
      snappyparams->refSurfLvlMax = 0;
    }
    if (shmparams.has_key("resolveFeatureAngle"))
      snappyparams->featAngle =
        shmparams["resolveFeatureAngle"].as<double>();
    else{
      snappyparams->featAngle = 60;
    }
    if (shmparams.has_key("locationInMeshX"))
      snappyparams->locMeshX =
        shmparams["locationInMeshX"].as<double>();
    else{
      std::cerr << "Location of a point in region you want to"
                << "keep cells is needed!" << std::endl;
                throw;
    }
    if (shmparams.has_key("locationInMeshY"))
        snappyparams->locMeshY =
          shmparams["locationInMeshY"].as<double>();
    else{
      std::cerr << "Location of a point in region you want to"
                << "keep cells is needed!" << std::endl;
                throw;
    }
    if (shmparams.has_key("locationInMeshZ"))
      snappyparams->locMeshZ =
        shmparams["locationInMeshZ"].as<double>();
    else{
      std::cerr << "Location of a point in region you want to"
                << "keep cells is needed!" << std::endl;
                throw;
    }
    if (shmparams.has_key("allowFreeStandingZoneFaces"))
      snappyparams->_alwFreeZone =
        shmparams["allowFreeStandingZoneFaces"].as<double>();
    else{
      snappyparams->_alwFreeZone = true;
    }
    if (shmparams.has_key("nSmoothPatch"))
      snappyparams->snapSmthPatch =
        shmparams["nSmoothPatch"].as<int>();
    else{
      snappyparams->snapSmthPatch = 4;
    }
    if (shmparams.has_key("tolerance"))
      snappyparams->snapTol =
        shmparams["tolerance"].as<double>();
    else{
      snappyparams->snapTol = 0.5;
    }
    if (shmparams.has_key("snapSolveIter"))
      snappyparams->solveSnapIter =
        shmparams["snapSolveIter"].as<int>();
    else{
      snappyparams->solveSnapIter = 200;
    }
    if (shmparams.has_key("snapRelaxIter"))
      snappyparams->relaxSnapIter =
        shmparams["snapRelaxIter"].as<int>();
    else{
      snappyparams->relaxSnapIter = 6;
    }
    if (shmparams.has_key("relativeSizes"))
      snappyparams->_relSize =
        shmparams["relativeSizes"].as<double>();
    else{
      snappyparams->_relSize = 1;
    }
    if (shmparams.has_key("expansionRatio"))
      snappyparams->expRatio =
        shmparams["expansionRatio"].as<double>();
    else{
      snappyparams->expRatio = 1.3;
    }
    if (shmparams.has_key("finalLayerThickness"))
      snappyparams->finLThick =
        shmparams["finalLayerThickness"].as<double>();
    else{
      snappyparams->finLThick = 1.0;
    }
    if (shmparams.has_key("minThickness"))
      snappyparams->minThick =
        shmparams["minThickness"].as<double>();
    else{
      snappyparams->minThick = 0.1;
    }
    if (shmparams.has_key("nGrow"))
      snappyparams->nGrow =
        shmparams["nGrow"].as<int>();
    else{
      snappyparams->nGrow = 0;
    }
    if (shmparams.has_key("featureAngle"))
      snappyparams->lyrFeatAngle =
        shmparams["featureAngle"].as<double>();
    else{
      snappyparams->lyrFeatAngle = 30;
    }
    if (shmparams.has_key("nRelaxIter"))
      snappyparams->lyrRelaxIter =
        shmparams["nRelaxIter"].as<int>();
    else{
      snappyparams->lyrRelaxIter = 3;
    }
    if (shmparams.has_key("nSmoothSurfaceNormals"))
      snappyparams->lyrSmthSurfNorm =
        shmparams["nSmoothSurfaceNormals"].as<int>();
    else{
      snappyparams->lyrSmthSurfNorm = 1;
    }
    if (shmparams.has_key("nSmoothNormals"))
      snappyparams->lyrSmthNorm =
      shmparams["nSmoothNormals"].as<int>();
    else{
      snappyparams->lyrSmthNorm = 3;
    }
    if (shmparams.has_key("nSmoothThickness"))
      snappyparams->lyrSmthThick =
        shmparams["nSmoothThickness"].as<int>();
    else{
      snappyparams->lyrSmthThick = 2;
    }
    if (shmparams.has_key("maxFaceThicknessRatio"))
      snappyparams->lyrMaxFcTR =
        shmparams["maxFaceThicknessRatio"].as<double>();
    else{
      snappyparams->lyrMaxFcTR = 0.5;
    }
    if (shmparams.has_key("maxThicknessToMedialRatio"))
      snappyparams->lyrMaxThickTMR =
        shmparams["maxThicknessToMedialRatio"].as<double>();
    else{
      snappyparams->lyrMaxThickTMR = 1.0;
    }
    if (shmparams.has_key("minMedialAxisAngle"))
      snappyparams->lyrMinMedAngl =
        shmparams["minMedialAxisAngle"].as<double>();
    else{
      snappyparams->lyrMinMedAngl = 90;
    }
    if (shmparams.has_key("nBufferCellsNoExtrude"))
      snappyparams->lyrBuffrCells =
        shmparams["nBufferCellsNoExtrude"].as<int>();
    else{
      snappyparams->lyrBuffrCells = 0;
    }
    if (shmparams.has_key("nLayerIter"))
      snappyparams->lyrIter =
        shmparams["nLayerIter"].as<int>();
    else{
      snappyparams->lyrIter = 50;
    }
    if (shmparams.has_key("maxNonOrtho"))
      snappyparams->qcMaxNOrtho =
        shmparams["maxNonOrtho"].as<int>();
    else{
      snappyparams->qcMaxNOrtho = 65;
    }
    if (shmparams.has_key("maxBoundarySkewness"))
      snappyparams->qcMaxBndrySkew =
        shmparams["maxBoundarySkewness"].as<double>();
    else{
      snappyparams->qcMaxBndrySkew = 20;
    }
    if (shmparams.has_key("maxInternalSkewness"))
      snappyparams->qcMaxIntSkew =
        shmparams["maxInternalSkewness"].as<double>();
    else{
      snappyparams->qcMaxIntSkew = 4;
    }
    if (shmparams.has_key("maxConcave"))
      snappyparams->qcMaxConc =
        shmparams["maxConcave"].as<double>();
    else{
      snappyparams->qcMaxConc = 80;
    }
    if (shmparams.has_key("minVol"))
      snappyparams->qcMinVol =
        shmparams["minVol"].as<double>();
    else{
      snappyparams->qcMinVol = 1e-13;
    }
    if (shmparams.has_key("minTetQuality"))
      snappyparams->qcMinTetQ =
        shmparams["minTetQuality"].as<double>();
    else{
      snappyparams->qcMinTetQ = 1e-15;
    }
    if (shmparams.has_key("minArea"))
      snappyparams->qcMinArea =
        shmparams["minArea"].as<double>();
    else{
      snappyparams->qcMinArea = -1;
    }
    if (shmparams.has_key("minTwist"))
      snappyparams->qcMinTwist =
        shmparams["minTwist"].as<double>();
    else{
      snappyparams->qcMinTwist = 0.02;
    }
    if (shmparams.has_key("minFaceWeight"))
      snappyparams->qcMinFaceW =
        shmparams["minFaceWeight"].as<double>();
    else{
      snappyparams->qcMinFaceW = 0.05;
    }
    if (shmparams.has_key("minVolRatio"))
      snappyparams->qcMinVolRto =
        shmparams["minVolRatio"].as<double>();
    else{
      snappyparams->qcMinVolRto = 0.01;
    }
    if (shmparams.has_key("minDeterminant"))
      snappyparams->qcMinDet =
        shmparams["minDeterminant"].as<double>();
    else{
      snappyparams->qcMinDet = 0.001;
    }
    if (shmparams.has_key("minTriangleTwist"))
      snappyparams->qcMinTrTwist =
        shmparams["minTriangleTwist"].as<double>();
    else{
      snappyparams->qcMinTrTwist = -1;
    }
    if (shmparams.has_key("qcnSmoothScale"))
      snappyparams->qcSmthScale =
        shmparams["qcnSmoothScale"].as<int>();
    else{
      snappyparams->qcSmthScale = 5;
    }
    if (shmparams.has_key("errorReduction"))
      snappyparams->qcErrRedctn =
        shmparams["errorReduction"].as<double>();
    else{
      snappyparams->qcErrRedctn = 0.75;
    }
    if (shmparams.has_key("mergeTolerance"))
      snappyparams->mergeTol =
        shmparams["mergeTolerance"].as<double>();
    else{
      snappyparams->mergeTol = 1e-06;
    }


    std::string cap2 = "GeomRefinementRegions";
    if (shmparams.has_key(cap2))
    {
      snappyparams->_withGeomRefReg = true;

      for (auto jptch2 : shmparams[cap2].array_range())
      {
        shmGeomRefine geomRef;

        if (jptch2.has_key("PatchName"))
          geomRef.patchNm = jptch2["PatchName"].as<std::string>();
        else{
          std::cerr << "Please define \"PatchName\"\n" << std::endl;
          throw;
        }
        if (jptch2.has_key("searchableShape"))
          geomRef.searchableName =
                jptch2["searchableShape"].as<std::string>();
        else{
          std::cerr << "Please define \"searchableShape\"\n" << std::endl;
          throw;
        }
        if (jptch2.has_key("shapeParams1"))
          geomRef.shapeParameters1 =
                jptch2["shapeParams1"].as<std::string>();
        else{
          std::cerr << "Insufficient Parameters\n" << std::endl;
          throw;
        }
        if (jptch2.has_key("shapeParams2"))
          geomRef.shapeParameters2 =
                jptch2["shapeParams2"].as<std::string>();
        if (jptch2.has_key("Radius"))
          geomRef.rad = jptch2["Radius"].as<double>();
        else{
          std::cerr << "Please define \"Radius\"\n" << std::endl;
          throw;
        }
        if (jptch2.has_key("Mode"))
          geomRef.mode = jptch2["Mode"].as<std::string>();
        else{
          std::cerr << "Please define \"Mode\"\n" << std::endl;
          throw;
        }
        if (jptch2.has_key("MinLevel"))
          geomRef.minLvl = jptch2["MinLevel"].as<int>();
        else{
          std::cerr << "Please define \"MinLevel\"\n" << std::endl;
          throw;
        }
        if (jptch2.has_key("MaxLevel"))
          geomRef.maxLvl = jptch2["MaxLevel"].as<int>();
        else{
              std::cerr << "Please define \"MaxLevel\"\n" << std::endl;
              throw;
        }

          (snappyparams->geomRefs).push_back(geomRef);
      }
    }


    std::string cap3 = "SurfaceRefinementRegions";
    if (shmparams.has_key(cap3))
    {
      snappyparams->_withSurfRefReg = true;

      for (auto jptch3 : shmparams[cap3].array_range())
      {

        shmRegionRef surfRef;

        if (jptch3.has_key("PatchName"))
          surfRef.refPatchNm = jptch3["PatchName"].as<std::string>();
        else{
          surfRef.refPatchNm = snappyparams->surfRefPatch;
        }

        if (jptch3.has_key("MinLevel"))
          surfRef.minLvl = jptch3["MinLevel"].as<int>();
        else{
          surfRef.minLvl = 1;
        }

        if (jptch3.has_key("MaxLevel"))
          surfRef.maxLvl = jptch3["MaxLevel"].as<int>();
        else{
          std::cerr << "Please define maximum surface refinement"
                    << "\"MaxLevel\"\n" << std::endl;
                    throw;
        }

        (snappyparams->surfRefs).push_back(surfRef);
      }
    }


    blockMeshParams* bmparams = new blockMeshParams();
    std::string defaults3 = 
     inputjson["Meshing Parameters"]["blockMesh Parameters"].as<std::string>();
    
    jsoncons::json bmshparams = inputjson["Meshing Parameters"]["blockMesh Parameters"];
    
    if (inputjson["Mesh File Options"].has_key("Input Dict File"))
      bmparams->_ownBlockMshDict = 
        inputjson["Mesh File Options"]["Input Dict File"].as<bool>();

    if (bmshparams.has_key("Block Geometry"))
      bmparams->_isBlock =
        bmshparams["Block Geometry"].as<bool>();
    else{
      std::cerr << "Define your choice of geometry using bool keys"
                << "\n" << std::endl;
                throw;
    }
    if (bmshparams.has_key("Sphere Geometry"))
      bmparams->_isSphere =
        bmshparams["Sphere Geometry"].as<bool>();
    else{
      std::cerr << "Define your choice of geometry using bool keys"
                << "\n" << std::endl;
                throw;
    }
    if (bmshparams.has_key("Cylinder/Tapered_Cone Geometry"))
      bmparams->_isCylinder_TCone =
        bmshparams["Cylinder/Tapered_Cone Geometry"].as<bool>();
    else{
      std::cerr << "Define your choice of geometry using bool keys"
                << "\n" << std::endl;
                throw;
    }
    if (bmshparams.has_key("scaleToMeters"))
      bmparams->cnvrtToMeters = 
       bmshparams["scaleToMeters"].as<double>();
    else{
      std::cerr << "Define your choice of geometry using bool keys"
                << "\n" << std::endl;
                throw;
    }

    if (bmshparams.has_key("Cell_Size")){
      bmparams->_cellSizeDefined = true;
      bmparams->cellSize = 
       bmshparams["Cell_Size"].as<double>();
    }
    else
    {
      bmparams->cellSize = -1;
      bmparams->_cellSizeDefined = false;
    }

    if (bmshparams.has_key("XdirectionCells"))
      bmparams->cellsXDir = 
       bmshparams["XdirectionCells"].as<int>();
    else{
      if (bmparams->_cellSizeDefined)
      {
        // Nothing
      }else{
      std::cerr << "Define cell numbers in X direction"
                << "\n" << std::endl;
                throw;
      }
    }
    if (bmshparams.has_key("YdirectionCells"))
      bmparams->cellsYDir = 
       bmshparams["YdirectionCells"].as<int>();
    else{
      if (bmparams->_cellSizeDefined)
      {
        // Nothing
      }else{
      std::cerr << "Define cell numbers in Y direction"
                << "\n" << std::endl;
                throw;
      }
    }
    if (bmshparams.has_key("ZdirectionCells"))
      bmparams->cellsZDir = 
       bmshparams["ZdirectionCells"].as<int>();
    else{
      if (bmparams->_cellSizeDefined)
      {
        // Nothing
      }else{
      std::cerr << "Define cell numbers in Z direction"
                << "\n" << std::endl;
                throw;
      }
    }
      
    if (bmshparams.has_key("Block Parameters"))
    {

      if (bmshparams["Block Parameters"].has_key("Auto_Generate"))
      {
        bmparams->_autoGenerateBox = true;

        if (inputjson["Mesh File Options"].has_key("Input Geometry File"))
          bmparams->packFileName = 
              inputjson["Mesh File Options"]
                      ["Input Geometry File"].as<std::string>();
        else
        {
          std::cerr << "A geometry file should be supplied.\n";
          throw;        
        }

        if (bmshparams["Block Parameters"]
              ["Auto_Generate"].has_key("Offset_XDir"))
                bmparams->offsetX = 
                  bmshparams["Block Parameters"]
                    ["Auto_Generate"]["Offset_XDir"].as<double>();
        else{
          bmparams->offsetX = 0.1;
        }
        if (bmshparams["Block Parameters"]
          ["Auto_Generate"].has_key("Offset_YDir"))
          bmparams->offsetY = 
              bmshparams["Block Parameters"]
                ["Auto_Generate"]["Offset_YDir"].as<double>();
        else{
          bmparams->offsetY = 0.1;
        }
        if (bmshparams["Block Parameters"]
          ["Auto_Generate"].has_key("Offset_ZDir"))
          bmparams->offsetZ = 
              bmshparams["Block Parameters"]
                ["Auto_Generate"]["Offset_ZDir"].as<double>();
        else{
          bmparams->offsetZ = 0.1;
        }

      }
      else{
        bmparams->_autoGenerateBox = false;
      }
      
      if (bmshparams["Block Parameters"].has_key("X1"))
        bmparams->initX = bmshparams["Block Parameters"]["X1"].as<double>();
      else{
        if (!(bmparams->_autoGenerateBox)){
          std::cerr << "Define initial point for block (X Y Z)!\n" << std::endl;
          throw;
        }
        else{
          // Nothing
        }
      }
      if (bmshparams["Block Parameters"].has_key("Y1"))
       bmparams->initY = bmshparams["Block Parameters"]["Y1"].as<double>();
      else{
        if (!(bmparams->_autoGenerateBox)){
          std::cerr << "Define initial point for block (X Y Z)!\n" << std::endl;
          throw;
        }
        else{
          // Nothing
        }
      }
      if (bmshparams["Block Parameters"].has_key("Z1"))
        bmparams->initZ = bmshparams["Block Parameters"]["Z1"].as<double>();
      else{
        if (!(bmparams->_autoGenerateBox)){
          std::cerr << "Define initial point for block (X Y Z)!\n" << std::endl;
          throw;
        }
        else{
          // Nothing
        }
      }
      if (bmshparams["Block Parameters"].has_key("LengthX"))
        bmparams->lenX = bmshparams["Block Parameters"]["LengthX"].as<double>();
      else{
        if (!(bmparams->_autoGenerateBox)){
          std::cerr << "Define initial point for block (X Y Z)!\n" << std::endl;
          throw;
        }
        else{
          // Nothing
        }
      }
      if (bmshparams["Block Parameters"].has_key("LengthY"))
        bmparams->lenY = bmshparams["Block Parameters"]["LengthY"].as<double>();
      else{
        if (!(bmparams->_autoGenerateBox)){
          std::cerr << "Define initial point for block (X Y Z)!\n" << std::endl;
          throw;
        }
        else{
          // Nothing
        }
      }
      if (bmshparams["Block Parameters"].has_key("LengthZ"))
        bmparams->lenZ = bmshparams["Block Parameters"]["LengthZ"].as<double>();
      else{
        if (!(bmparams->_autoGenerateBox)){
          std::cerr << "Define initial point for block (X Y Z)!\n" << std::endl;
          throw;
        }
        else{
          // Nothing
        }
      }
      if (bmshparams["Block Parameters"].has_key("GradingXdir"))
        bmparams->smplGradingX = 
        bmshparams["Block Parameters"]["GradingXdir"].as<double>();
        else{
          bmparams->smplGradingX = 1;
        }
      if (bmshparams["Block Parameters"].has_key("GradingYdir"))
        bmparams->smplGradingY = 
        bmshparams["Block Parameters"]["GradingYdir"].as<double>();
        else{
          bmparams->smplGradingY = 1;
        }
      if (bmshparams["Block Parameters"].has_key("GradingZdir"))
        bmparams->smplGradingZ = 
        bmshparams["Block Parameters"]["GradingZdir"].as<double>();
        else{
          bmparams->smplGradingZ = 1;
        }
    }
    
    if (bmshparams.has_key("Sphere Parameters"))
    {
      
      if (bmshparams["Sphere Parameters"].has_key("Center X"))
        bmparams->centerX =
            bmshparams["Sphere Parameters"]["Center X"].as<double>();
      else{
        std::cerr << "Define sphere center!\n" << std::endl;
        throw;
      }
      if (bmshparams["Sphere Parameters"].has_key("Center Y"))
        bmparams->centerY =
            bmshparams["Sphere Parameters"]["Center Y"].as<double>();
      else{
        std::cerr << "Define sphere center!\n" << std::endl;
        throw;
      }
      if (bmshparams["Sphere Parameters"].has_key("Center Z"))
        bmparams->centerZ =
            bmshparams["Sphere Parameters"]["Center Z"].as<double>();
      else{
        std::cerr << "Define sphere center!\n" << std::endl;
        throw;
      }
      if (bmshparams["Sphere Parameters"].has_key("Radius"))
        bmparams->radius =
            bmshparams["Sphere Parameters"]["Radius"].as<double>();
      else{
        std::cerr << "Define sphere radius!\n" << endl;
        throw;
      }
      if (bmshparams["Sphere Parameters"].has_key("GradingXdir"))
        bmparams->sphrGradingX =
            bmshparams["Sphere Parameters"]["GradingXdir"].as<double>();
      else{
        bmparams->sphrGradingX = 1;
      }
      if (bmshparams["Sphere Parameters"].has_key("GradingYdir"))
        bmparams->sphrGradingY =
            bmshparams["Sphere Parameters"]["GradingYdir"].as<double>();
      else{
        bmparams->sphrGradingY = 1;
      }
      if (bmshparams["Sphere Parameters"].has_key("GradingXdir"))
        bmparams->sphrGradingZ =
            bmshparams["Sphere Parameters"]["GradingZdir"].as<double>();
      else{
        bmparams->sphrGradingZ = 1;
      }
    }
    
    if (bmshparams.has_key("Cylinder/Tapered_Cone Parameters"))
    {
      if (bmshparams["Cylinder/Tapered_Cone Parameters"].has_key("Center X"))
        bmparams->centerCyl[0] =
              bmshparams["Cylinder/Tapered_Cone Parameters"]
                                                ["Center X"].as<double>();
      else{
        std::cerr << "Define center point for cylinder axis\n" << std::endl;
        throw;
      }
      if (bmshparams["Cylinder/Tapered_Cone Parameters"].has_key("Center Y"))
        bmparams->centerCyl[1] =
              bmshparams["Cylinder/Tapered_Cone Parameters"]
                                                ["Center Y"].as<double>();
      else{
        std::cerr << "Define center point for cylinder axis\n" << std::endl;
        throw;
      }
      if (bmshparams["Cylinder/Tapered_Cone Parameters"].has_key("Center Z"))
        bmparams->centerCyl[2] =
              bmshparams["Cylinder/Tapered_Cone Parameters"]
                                                ["Center Z"].as<double>();
      else{
        std::cerr << "Define center point for cylinder axis\n" << std::endl;
        throw;
      }
      if (bmshparams["Cylinder/Tapered_Cone Parameters"].has_key("Radius1"))
        bmparams->radius1 =
              bmshparams["Cylinder/Tapered_Cone Parameters"]
                                                ["Radius1"].as<double>();
      
      if (bmshparams["Cylinder/Tapered_Cone Parameters"].has_key("Radius2")){
        bmparams->radius2 =
                bmshparams["Cylinder/Tapered_Cone Parameters"]
                                                  ["Radius2"].as<double>();
      }
      else
      {
        bmparams->radius2 =
                bmshparams["Cylinder/Tapered_Cone Parameters"]
                  ["Radius1"].as<double>();
      }

      if (bmshparams["Cylinder/Tapered_Cone Parameters"].has_key("GradingXdir"))
        bmparams->cylGrading[0] =
            bmshparams["Cylinder/Tapered_Cone Parameters"]
                  ["GradingXdir"].as<double>();
      else{
        bmparams->cylGrading[0] = 1;
      }
      if (bmshparams["Cylinder/Tapered_Cone Parameters"].has_key("GradingYdir"))
        bmparams->cylGrading[1] =
            bmshparams["Cylinder/Tapered_Cone Parameters"]
                  ["GradingYdir"].as<double>();
      else{
        bmparams->cylGrading[1] = 1;
      }
      if (bmshparams["Cylinder/Tapered_Cone Parameters"].has_key("GradingXdir"))
        bmparams->cylGrading[2] =
            bmshparams["Cylinder/Tapered_Cone Parameters"]
                  ["GradingZdir"].as<double>();
      else{
        bmparams->cylGrading[2] = 1;
      }
      
      if (bmshparams["Cylinder/Tapered_Cone Parameters"].has_key("Height"))
        bmparams->height =
            bmshparams["Cylinder/Tapered_Cone Parameters"]
                ["Height"].as<double>();
      else{
        std::cerr << "Define cylinder height\n" << endl;
      }
    }


    MeshManipulationFoamParams* mparams = new MeshManipulationFoamParams();
    std::string defaults4 = 
      inputjson["Meshing Parameters"]
          ["MeshManipulation Parameters"].as<std::string>();
    
    jsoncons::json pmshparams =
      inputjson["Meshing Parameters"]["MeshManipulation Parameters"];

    
    // Parameter parsing starts here
    if (pmshparams.has_key("Enable SurfLambdaMuSmooth"))
          mparams->_doSurfaceLMSmth =
            pmshparams["Enable SurfLambdaMuSmooth"].as<bool>();
    if (pmshparams.has_key("Enable splitMeshRegions"))
          mparams->_doSplitMshRegs =
            pmshparams["Enable splitMeshRegions"].as<bool>();
    if (pmshparams.has_key("Enable MergeMeshes"))
          mparams->_doMergeMsh = 
            pmshparams["Enable MergeMeshes"].as<bool>();
    if (pmshparams.has_key("Enable CreatePatch"))
          mparams->_doCreatePtchs = 
            pmshparams["Enable CreatePatch"].as<bool>();
    if (pmshparams.has_key("Enable foamToSurface"))
          mparams->_doFoam2Surf = 
            pmshparams["Enable foamToSurface"].as<bool>();
    if (pmshparams.has_key("Enable surfaceSplitByTopology"))
          mparams->_doSurfSplit = 
            pmshparams["Enable surfaceSplitByTopology"].as<bool>();

    if (pmshparams.has_key("SurfLambdaMuSmooth Parameters"))
    {
      if (pmshparams["SurfLambdaMuSmooth Parameters"].has_key(
              "AddFeatureFile?"))
          mparams->_addFeatureFile =
            pmshparams["SurfLambdaMuSmooth Parameters"]
              ["AddFeatureFile?"].as<bool>();

      if (pmshparams["SurfLambdaMuSmooth Parameters"].has_key(
              "Input STL File"))
          mparams->slmssurfaceFile =
            pmshparams["SurfLambdaMuSmooth Parameters"]
              ["Input STL File"].as<std::string>();

      if (pmshparams["SurfLambdaMuSmooth Parameters"].has_key(
              "Output STL File"))
          mparams->slmsoutputFile =
            pmshparams["SurfLambdaMuSmooth Parameters"]
              ["Output STL File"].as<std::string>();

      if (pmshparams["SurfLambdaMuSmooth Parameters"].has_key(
              "Lambda"))
          mparams->lambda =
            pmshparams["SurfLambdaMuSmooth Parameters"]
              ["Lambda"].as<double>();

      if (pmshparams["SurfLambdaMuSmooth Parameters"].has_key(
              "Mu"))
          mparams->mu =
            pmshparams["SurfLambdaMuSmooth Parameters"]
              ["Mu"].as<double>();

      if (pmshparams["SurfLambdaMuSmooth Parameters"].has_key(
              "Smoothing Interations"))
          mparams->slmsIterations =
            pmshparams["SurfLambdaMuSmooth Parameters"]
              ["Smoothing Interations"].as<int>();           
    }


    if (pmshparams.has_key("splitMeshRegions Parameters"))
    {
      if (pmshparams["splitMeshRegions Parameters"].has_key("overwrite?"))
        mparams->_overwriteMsh = 
        pmshparams["splitMeshRegions Parameters"]["overwrite?"].as<bool>();

      if (pmshparams["splitMeshRegions Parameters"].has_key("usecellZones?"))
        mparams->_cellZones = 
        pmshparams["splitMeshRegions Parameters"]["usecellZones?"].as<bool>();
    }


    if (pmshparams.has_key("mergeMeshes Parameters"))
    {
      if (pmshparams["mergeMeshes Parameters"].has_key("Master Region"))
        mparams->masterCase =
          pmshparams["mergeMeshes Parameters"]
            ["Master Region"].as<std::string>();

      if (pmshparams["mergeMeshes Parameters"].has_key("Add Region"))
        mparams->addCase = 
          pmshparams["mergeMeshes Parameters"]
            ["Add Region"].as<std::string>();

      if (pmshparams["mergeMeshes Parameters"].has_key("overwrite?"))
        mparams->_overwriteMergeMsh =
          pmshparams["mergeMeshes Parameters"]
            ["overwrite?"].as<bool>();

      if (pmshparams["mergeMeshes Parameters"].has_key("Master Region Path"))
        mparams->masterCasePath =
          pmshparams["mergeMeshes Parameters"]
            ["Master Region Path"].as<std::string>();

      if (pmshparams["mergeMeshes Parameters"].has_key("Add Region Path"))
        mparams->addCasePath = 
          pmshparams["mergeMeshes Parameters"]
            ["Add Region Path"].as<std::string>();

      if (pmshparams["mergeMeshes Parameters"].has_key("Number of Domains"))
        mparams->numDomains = 
          pmshparams["mergeMeshes Parameters"]
            ["Number of Domains"].as<int>();
      else{
        mparams->numDomains = -1;
      }
    }


    if (pmshparams.has_key("createPatch Parameters"))
    {
     if (pmshparams["createPatch Parameters"].has_key("Surrounding PatchName"))
        mparams->surroundingName = 
          pmshparams["createPatch Parameters"]
            ["Surrounding PatchName"].as<std::string>();

      if (pmshparams["createPatch Parameters"].has_key("Packs PatchName"))
        mparams->packsName = 
          pmshparams["createPatch Parameters"]
            ["Packs PatchName"].as<std::string>();

     if (pmshparams["createPatch Parameters"].has_key("Surrounding PatchType"))
        mparams->srrndngPatchType = 
          pmshparams["createPatch Parameters"]
            ["Surrounding PatchType"].as<std::string>();

      if (pmshparams["createPatch Parameters"].has_key("Packs PatchType"))
        mparams->packsPatchType = 
          pmshparams["createPatch Parameters"]
            ["Packs PatchType"].as<std::string>();

      if (pmshparams["createPatch Parameters"].has_key("overwrite?"))
        mparams->_overwritecpMsh = 
          pmshparams["createPatch Parameters"]
            ["overwrite?"].as<bool>();

      if (pmshparams["createPatch Parameters"].has_key("Packs Path"))
        mparams->pathPacks = 
          pmshparams["createPatch Parameters"]
            ["Packs Path"].as<std::string>();

      if (pmshparams["createPatch Parameters"].has_key("Surrounding Path"))
        mparams->pathSurrounding = 
          pmshparams["createPatch Parameters"]
            ["Surrounding Path"].as<std::string>();

    }
  
    if (pmshparams.has_key("foamToSurface Parameters"))
    {
      if (pmshparams["foamToSurface Parameters"].has_key("Output File Path"))
        mparams->outSurfName = 
          pmshparams["foamToSurface Parameters"]
            ["Output File Path"].as<std::string>();
    }

   if (pmshparams.has_key("surfaceSplitByTopology Parameters"))
   {
    if (pmshparams["surfaceSplitByTopology Parameters"].has_key("Input File"))
        mparams->surfFile = 
          pmshparams["surfaceSplitByTopology Parameters"]
            ["Input File"].as<std::string>();
    else{
      mparams->surfFile = cfparams->geomFilePath;
    }
    if (pmshparams["surfaceSplitByTopology Parameters"].has_key("Output File"))
        mparams->outSurfFile = 
          pmshparams["surfaceSplitByTopology Parameters"]
            ["Output File"].as<std::string>();
    else{
      mparams->outSurfFile = cfparams->geomFilePath;
    }
   }
   else{
      mparams->surfFile = cfparams->geomFilePath;
      mparams->outSurfFile = cfparams->geomFilePath;
   }

    
    // Give all data to PackMesh Driver
    PackMeshDriver* pckmshdrvobj = new PackMeshDriver(ifname, mparams, 
                                                      cfparams, snappyparams,
                                                      bmparams, ofname1,
                                                      ofname2);
    return pckmshdrvobj;
  } 
  else
  {
    std::cout << "Mesh generation engine " << meshEngine 
      << " is not supported" << std::endl;
    exit(1);
  }
}