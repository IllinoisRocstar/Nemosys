#ifdef HAVE_CFMSH
#  include <boost/filesystem.hpp>

#  include "MeshManipulationFoam.H"
#  include "MeshManipulationFoamParams.H"
#  include "MeshQuality.H"
#  include "MeshQualityDriver.H"
#  include "blockMeshGen.H"
#  include "blockMeshParams.H"
#  include "snappymeshGen.H"
#  include "snappymeshParams.H"
#endif

#include "PackMeshDriver.H"
#include "meshBase.H"
#include "rocPack.H"

// std c++
#include <iostream>
#include <memory>
#include <tuple>
#include <vector>

// vtk
#include <AuxiliaryFunctions.H>
#include <vtkCellData.h>
#include <vtkIdList.h>
#include <vtkMesh.H>
#include <vtkGeoMesh.H>

#include "ConversionDriver.H"

// TODO
// 1. Some implementation to improve mesh if mesh is not of desired
//    quality.

// ------------------------- Pack Mesh Driver -------------------------------//

#ifdef HAVE_CFMSH
// Standard constructor controls the pack meshing process.
PackMeshDriver::PackMeshDriver(
    const std::string &ifname, MeshManipulationFoamParams *_mparams,
    snappymeshParams *_snappyparams,
    blockMeshParams *_bmparams, const std::string &ofname_pack,
    const std::string &ofname_surrndng, const std::string &ofname_merged,
    const bool &useRocpack, const double &locAdjust) {
  std::cout << "PackMeshDriver constructed!" << std::endl;

  // Makes sure that constant and triSurface directories are present.
  // If directories are already present, it will not do anything.
  const char dir_path[] = "./constant";
  boost::filesystem::path dir(dir_path);
  try {
    boost::filesystem::create_directory(dir);
  } catch (boost::filesystem::filesystem_error &e) {
    std::cerr << "Problem in creating triSurface directory"
              << "for the snappyHexMesh"
              << "\n";
    std::cerr << e.what() << std::endl;
    throw;
  }

  const char dir_path1[] = "./constant/triSurface";
  boost::filesystem::path dir1(dir_path1);
  try {
    boost::filesystem::create_directory(dir1);
  } catch (boost::filesystem::filesystem_error &e) {
    std::cerr << "Problem in creating triSurface directory"
              << "for the snappyHexMesh"
              << "\n";
    std::cerr << e.what() << std::endl;
    throw;
  }

  // This object has access to all MeshManipulation utilities.
  MeshManipulationFoam *objMsh = new MeshManipulationFoam(_mparams);
  const char *nameFile = "a";  // Dummy name for input

  // A rocpack method that creates stl and then moves it to triSurface using
  // boost.
  if (useRocpack) {
    std::string hexOutSTL = ifname + ".stl";
    auto *objrocPck = new NEM::GEO::rocPack(ifname, hexOutSTL);

    objrocPck->removeBoundaryVolumes();
    objrocPck->rocPack2Surf();

    if (objrocPck) delete objrocPck;

    const std::string dir_path1 = hexOutSTL;
    boost::filesystem::path dir1(dir_path1);

    const std::string dir_path2 = "./constant/triSurface/" + ifname + ".stl";
    boost::filesystem::path dir2(dir_path2);

    boost::filesystem::copy_file(
        dir1, dir2, boost::filesystem::copy_option::overwrite_if_exists);
  } else {
    const std::string dir_path1 = ifname;
    boost::filesystem::path dir1(dir_path1);

    const std::string dir_path2 = "./constant/triSurface/" + ifname;
    boost::filesystem::path dir2(dir_path2);

    boost::filesystem::copy_file(
        dir1, dir2, boost::filesystem::copy_option::overwrite_if_exists);
  }

  // blockMesh utility takes user input for surrounding box region and
  // generates mesh block in constant/polyMesh folder. It will overwrite the
  // previous mesh created by CfMesh. This mesh will be used as background mesh
  // by snappyHexMesh later.
  blockMeshGen *objBM = new blockMeshGen(_bmparams);
  objBM->createMeshFromSTL(nameFile);

  // Updating location for next process
  if (locAdjust > 0) {
    _snappyparams->locMeshX = _bmparams->coordsBox.first[0] + 0.001 + locAdjust;
    _snappyparams->locMeshY = _bmparams->coordsBox.first[1] + 0.001 + locAdjust;
    _snappyparams->locMeshZ = _bmparams->coordsBox.first[2] + 0.001 + locAdjust;
  } else {
    _snappyparams->locMeshX = _bmparams->coordsBox.first[0] + 0.001;
    _snappyparams->locMeshY = _bmparams->coordsBox.first[1] + 0.001;
    _snappyparams->locMeshZ = _bmparams->coordsBox.first[2] + 0.001;
  }

  // snappyHexMesh reads background mesh created using blockMesh and takes
  // surface file from foamToSurface utility to snap pack surface onto back-
  // ground mesh and creates different cellZones (i.e different solids). These
  // interfaces between pack and surrounding regions are completely conformal
  // due to snappyHexMesh's unique snapping abilities.
  snappymeshGen *objSHM = new snappymeshGen(_snappyparams);
  objSHM->createMeshFromSTL(nameFile);

  // splitMeshRegions reads mesh from constant/polyMesh directory and splits
  // mesh into separate regions. Each region will represent one solid. It will
  // write bunch of domain.* directories inside constant/ and system/ folders.
  // This function does a random walking around mesh to identify different
  // regions and in that process, while naming all domains as domain.*
  // in constant folder, it skips one number for the disconnected region it
  // encounters first. This number is taken out to provide as input to merge
  // mesh.
  std::pair<std::vector<int>, std::string> dirStat = objMsh->splitMshRegions();

  int skippedDir = dirStat.first[0];
  int totalRegs = dirStat.first[1] - 1;
  std::string surroundingRegion = dirStat.second;

  // mergeMeshes will read master domain (defined by user) from constant folder
  // and start merging other domain to it untill all slave domains are attached
  // to master domain. It loops through all domains in sequential manner and
  // skips the missing domain.* (also skipped by splitMeshRegions) to avoid
  // runtime error.
  std::cout << "Total # of domains are = " << totalRegs << std::endl;
  if (totalRegs == 1) {
    // Nothing
  } else {
    objMsh->mergeMeshes(skippedDir, totalRegs);
  }

  // createPatch utility reads domain mesh and createPatchDict to combine all
  // different patches of multiple packs/surrounding into one patch
  // respectively
  if (totalRegs == 1) {
    // Nothing
  } else {
    objMsh->createPatch(skippedDir);
  }

  // Reads current mesh and write it to separate VTK/VTU files
  bool readDB = false;

  // Reads and converts pack mesh
  std::string regNme;
  if (skippedDir == 1)
    regNme = "domain2";
  else
    regNme = "domain1";

  if (totalRegs == 1) regNme = _snappyparams->singleSolidPatch;

  meshBase *fm = new FOAM::foamMesh(readDB);
  fm->read(regNme);
  vtkMesh *vm2 = new vtkMesh(fm->getDataSet(), ofname_pack);
  std::vector<double> physIdPack 
    = std::vector<double>(fm->getNumberOfCells(),1);
  vm2->setCellDataArray("PhysGrpId",physIdPack);
  vm2->write();

  // Reads and converts surronding mesh
  regNme = surroundingRegion;

  if (totalRegs == 1) regNme = surroundingRegion;
  meshBase *fm2 = new FOAM::foamMesh(readDB);
  fm2->read(regNme);
  std::vector<double> physIdSurrounding 
    = std::vector<double>(fm2->getNumberOfCells(),0);
  vtkMesh *vm3 = new vtkMesh(fm2->getDataSet(), ofname_surrndng);
  vm3->setCellDataArray("PhysGrpId",physIdSurrounding);
  vm3->write();

  // Merge meshes and write
  vtkMesh *vm = new vtkMesh(vm2->getDataSet(),ofname_merged);
  vm->merge(vm3->getDataSet());
  vm->report();
  vm->write();

  // // Outputs useful mesh quality parameters for users
  // std::string SurroundingName = "surroundingMeshQuality";
  // std::string Surroundingmesh = "geom_surrounding_mesh.vtu";
  // std::string PackName = "packMeshQuality";
  // std::string Packmesh = "geom_pack_mesh.vtu";
  // MeshQualityDriver *objSurrQ =
  //     new MeshQualityDriver(Surroundingmesh, SurroundingName);
  // MeshQualityDriver *objPackQ = new MeshQualityDriver(Packmesh, PackName);

  // Cleaning up
  if (vm) delete vm;
  if (fm) delete fm;
  if (vm2) delete vm2;
  if (vm3) delete vm3;
  if (fm2) delete fm2;
  if (objMsh) delete objMsh;
  if (objSHM) delete objSHM;
  if (objBM) delete objBM;
  // if (objSurrQ) delete objSurrQ;
  // if (objPackQ) delete objPackQ;
  // End of workflow
}
#endif

PackMeshDriver::PackMeshDriver(
    const std::string &ifname, const std::string &ofname,
    const bool &scaleVolumes, const double &scaleValue, const double &meshSize,
    const bool &rmvBndryPacks, const bool &setPeriodicGeom,
    const bool &setPeriodicMesh, const bool &enable2PhysGrps,
    const bool &enableMultiPhysGrps, const bool &wantGeometryOnly,
    const bool &createCohesive, const bool &enablePatches,
    const std::vector<double> &transferMesh, const bool &customDomain,
    const std::vector<double> &domainBounds, const int &mshAlgorithm,
    const bool &enableDefaultOutput, const bool &enablePhysGrpPerShape,
    const int &refineLevel, const double &upperThreshold,
    const double &lowerThreshold, const bool &preserveSize,
    const int &elemOrder) {
  auto *objrocPck = new NEM::GEO::rocPack(ifname, ofname);

  if ((enable2PhysGrps && enableMultiPhysGrps) ||
      (enablePhysGrpPerShape && enable2PhysGrps) ||
      (enablePhysGrpPerShape && enableMultiPhysGrps)) {
    std::cerr << "Please select only one of the options for physical groups"
              << std::endl;
    throw;
  }

  if ((customDomain && setPeriodicGeom) || (customDomain && setPeriodicMesh)) {
    std::cerr
        << "WARNING!! -> Cannot make geometry periodic using custom bounds."
        << " using remove geometries on boundary option instead for meshing!"
        << std::endl;
  }

  if (upperThreshold > 0.0 || lowerThreshold > 0.0)
    objrocPck->applyFilter(upperThreshold, lowerThreshold);

  objrocPck->translateAll(transferMesh[0], transferMesh[1], transferMesh[2]);
  objrocPck->setMeshingAlgorithm(mshAlgorithm);

  if (elemOrder > 2 || elemOrder <= 0) {
    std::cerr << "Only element orders 1 and 2 are supported!" << std::endl;
    throw;
  }

  objrocPck->setElementOrder(elemOrder);

  if (createCohesive) {
    objrocPck->sanityCheckOn();
    objrocPck->enableCohesiveElements();
  }

  if (preserveSize) objrocPck->setSizePreservation();

  if (refineLevel > 0) objrocPck->assignRefinement(refineLevel);

  if (customDomain) objrocPck->setCustomDomain(domainBounds);

  if (rmvBndryPacks) objrocPck->removeBoundaryVolumes();

  if (scaleVolumes) objrocPck->shrinkVolumes(scaleValue);

  if (meshSize != 0) objrocPck->setMeshSize(meshSize);

  if (setPeriodicGeom) objrocPck->setPeriodicGeometry();

  if (enableDefaultOutput) objrocPck->enableDefOuts();

  if (wantGeometryOnly) objrocPck->rocPack2Surf();

  if (enableMultiPhysGrps) objrocPck->enablePhysicalGrps();

  if (enable2PhysGrps) objrocPck->enableTwoPhysGrps();

  if (enablePhysGrpPerShape) objrocPck->enablePhysicalGroupsPerShape();

  if (enablePatches) objrocPck->enableSurfacePatches();

  if (setPeriodicMesh) {
    if (customDomain) {
      objrocPck->removeBoundaryVolumes();
      objrocPck->setPeriodicMesh();
      objrocPck->rocPack2Periodic3D();
    } else {
      objrocPck->setPeriodicGeometry();
      objrocPck->setPeriodicMesh();
      objrocPck->rocPack2Periodic3D();
    }
  }

  if (objrocPck) delete objrocPck;
}

// Class destructor
PackMeshDriver::~PackMeshDriver() {
  std::cout << "PackMeshDriver destroyed" << endl;
}

// Reads JSON input provided by user and json file passes down by NemDriver.
PackMeshDriver *PackMeshDriver::readJSON(const jsoncons::json inputjson) {
  std::string meshType = inputjson["Mesh Type"].as<std::string>();

#ifndef HAVE_CFMSH
  if (meshType == "hexahedral") {
    std::cerr << "Build NEMoSys with -DENABLE_CFMSH=ON for hexahedral meshing!"
              << std::endl;
    throw;
  }
#endif

  if (meshType == "tetrahedral") {
    std::string ifname =
        inputjson["Mesh File Options"].contains("Input Rocpack File")
            ? inputjson["Mesh File Options"]["Input Rocpack File"]
                  .as<std::string>()
            : "empty";

    std::string outFile =
        inputjson["Mesh File Options"].contains("Output Mesh File")
            ? inputjson["Mesh File Options"]["Output Mesh File"]
                  .as<std::string>()
            : "empty";

    if (ifname == "empty" || outFile == "empty") {
      std::cout << "Please define Input Rocpack File and Output Mesh File!"
                << std::endl;
      throw;
    }

    return readJSON(ifname, outFile, inputjson);
  } else if (meshType == "geometry") {
    std::string ifname =
        inputjson["Mesh File Options"].contains("Input Rocpack File")
            ? inputjson["Mesh File Options"]["Input Rocpack File"]
                  .as<std::string>()
            : "empty";

    std::string outFile =
        inputjson["Mesh File Options"].contains("Output Mesh File")
            ? inputjson["Mesh File Options"]["Output Mesh File"]
                  .as<std::string>()
            : "empty";

    // wantGeometryOnly = true;

    if (ifname == "empty" || outFile == "empty") {
      std::cout << "Please define Input Rocpack File and Output Mesh File!"
                << std::endl;
      throw;
    }
    return readJSON(ifname, outFile, inputjson);
  }
#ifdef HAVE_CFMSH
  else if (meshType == "hexahedral") {
    bool useRocPack = false;

    if ((inputjson["Mesh File Options"].contains("Input Rocpack File")) && 
       (inputjson["Mesh File Options"].contains("Input Geometry File"))) {
      std::cerr << "Please define only on of the following options, not both!" 
                << std::endl << " 1. \"Input Rocpack File\"" << std::endl
                << " 2. \"Input Geometry File\"" << std::endl;
      throw;
    }

    std::string ifname =
        inputjson["Mesh File Options"].contains("Input Rocpack File")
            ? inputjson["Mesh File Options"]["Input Rocpack File"]
                  .as<std::string>()
            : "empty";

    if (ifname == "empty") {
      ifname = inputjson["Mesh File Options"].contains("Input Geometry File")
                   ? inputjson["Mesh File Options"]["Input Geometry File"]
                         .as<std::string>()
                   : "empty";
    } else {
      useRocPack = true;
    }

    if (ifname == "empty") {
      std::cout << "Please define Input Rocpack File or Input Geometry File!"
                << std::endl;
      throw;
    }

    std::string ofname_pack =
        inputjson["Mesh File Options"]["Output Pack Mesh File"]
            .as<std::string>();
    std::string ofname_surrndng =
        inputjson["Mesh File Options"]["Output Surrounding Mesh File"]
            .as<std::string>();

    std::string ofname_merged = 
            inputjson["Mesh File Options"].contains("Output Combined Mesh File")
                   ? inputjson["Mesh File Options"]["Output Combined Mesh File"]
                         .as<std::string>()
                   : "PackMesh.vtu";

    return readJSON(ifname, ofname_pack, ofname_surrndng, ofname_merged, 
                    useRocPack, inputjson);
  }
#endif
  else {
    std::cout << "Define Keywords" << std::endl;
    throw;
  }
}

#ifdef HAVE_CFMSH
// Reads JSON file parameters and catagorizes them into separate objects
// for convenience in Pack Mesh workflow.
PackMeshDriver *PackMeshDriver::readJSON(const std::string &ifname,
                                         const std::string &ofname_pack,
                                         const std::string &ofname_surrndng,
                                         const std::string &ofname_merged,
                                         const bool useRocpack,
                                         const jsoncons::json inputjson) {
  std::string meshEngine =
      inputjson["Mesh Generation Engine"].as<std::string>();

  if (!meshEngine.compare("packmesh")) {
    // Snappy hex mesh params
    snappymeshParams *params = new snappymeshParams();
    std::string defaults =
        inputjson["Meshing Parameters"]["snappyHexMesh Parameters"]
            .as<std::string>();
    double locAdjust = -1.0;
    if (!defaults.compare("default")) {
      // MeshGenDriver* mshgndrvobj =
      //     new MeshGenDriver(ifname, meshEngine, params, ofname);
      // return mshgndrvobj;
    } else {
      jsoncons::json shmparams =
          inputjson["Meshing Parameters"]["snappyHexMesh Parameters"];
          
      jsoncons::json geomParams = shmparams["Geometry Definition"];

      if (inputjson["Mesh File Options"].contains("Input Geometry File"))
        params->geomFileName =
            inputjson["Mesh File Options"]["Input Geometry File"]
                .as<std::string>();
      else {
        params->geomFileName = ifname + ".stl";
      }

      // General booleans
      if (shmparams.contains("Castellated Mesh"))
        params->_withCastMesh = shmparams["Castellated Mesh"].as<bool>();
      else {
        std::cerr << "Please specify on/off choice using \"Castellated Mesh\""
                  << "\n"
                  << std::endl;
        throw;
      }
      if (shmparams.contains("Snapping"))
        params->_withSnap = shmparams["Snapping"].as<bool>();
      else {
        std::cerr << "Please specify on/off choice using \"Snapping\""
                  << "Keyword!\n"
                  << std::endl;
        throw;
      }
      if (shmparams.contains("Layer Addition"))
        params->_withLayers = shmparams["Layer Addition"].as<bool>();
      else {
        std::cerr << "Please specify on/off choice using \"Layer Addition\""
                  << "Keyword!\n"
                  << std::endl;
        throw;
      }

      // Geometry definition with patch defining ability
      bool multiPatches = false;
      if (geomParams.contains("Enable Multi Patches"))
        multiPatches = geomParams["Enable Multi Patches"].as<bool>();
      else {
        std::cerr << "Please provide boolean \"Enable Multi Patches\".\n";
        throw;
      }

      if (!multiPatches) {
        if (geomParams.contains("InputPatchName"))
          params->singleSolidPatch =
              geomParams["InputPatchName"].as<std::string>();
        else {
          std::cerr << "A patch name for input geometry must be defined.\n";
          throw;
        }
      } else {
        // Geometry STL Patch Definition
        std::string cap4 = "Geometry Patches";
        if (geomParams.contains(cap4)) {
          params->_withMultiPatches = true;

          for (auto jptch3 : geomParams[cap4].array_range()) {
            shmSTLDefinition stlPatches;

            if (jptch3.contains("Geometry Patch Name"))
              stlPatches.STLPatchName =
                  jptch3["Geometry Patch Name"].as<std::string>();
            else {
              std::cerr
                  << "Please provide Output Patch Name for STL file provided"
                  << std::endl;
              throw;
            }

            if (jptch3.contains("Output Patch Name"))
              stlPatches.snappyPatchName =
                  jptch3["Output Patch Name"].as<std::string>();
            else {
              std::cerr
                  << "Please provide Output Patch Name for STL file provided"
                  << std::endl;
              throw;
            }

            (params->stlPatchDefs).push_back(stlPatches);
          }
        } else {
          std::cerr
              << "User has selected to define multiple geometry patches.\n"
              << "Please define them using \"Geometry Patches\" keyword"
              << std::endl;
          throw;
        }
      }

      // Define searchable shapes with whatever patch name user wants
      std::string cap5 = "Custom Patches";
      if (geomParams.contains(cap5)) {
        for (auto jptch3 : geomParams[cap5].array_range()) {
          shmSearchableShapes SSS;
          if (jptch3.contains("Custom Patch Name"))
            SSS.patchNm = jptch3["Custom Patch Name"].as<std::string>();
          else {
            std::cerr << "Please provide Custom Patch Name" << std::endl;
            throw;
          }

          if (jptch3.contains("Searchable Shape"))
            SSS.searchableName = jptch3["Searchable Shape"].as<std::string>();
          else {
            std::cerr << "Please provide Searchable Shape Name" << std::endl;
            throw;
          }

          if (SSS.searchableName == "searchableBox") {
            if (jptch3.contains("minimum bound"))
              SSS.shapeParameters1 = jptch3["minimum bound"].as<std::string>();
            else {
              std::cerr << "Please provide minimum bound (i.e (-1 -1 -1))"
                        << std::endl;
              throw;
            }
            if (jptch3.contains("maximum bound"))
              SSS.shapeParameters2 = jptch3["maximum bound"].as<std::string>();
            else {
              std::cerr << "Please provide maximum bound (i.e (1 1 1))"
                        << std::endl;
              throw;
            }
          } else if (SSS.searchableName == "searchableCylinder") {
            if (jptch3.contains("Axis Point 1"))
              SSS.shapeParameters1 = jptch3["Axis Point 1"].as<std::string>();
            else {
              std::cerr << "Please provide Axis Point 1 (i.e (-1 -1 -1))"
                        << std::endl;
              throw;
            }
            if (jptch3.contains("Axis Point 2"))
              SSS.shapeParameters2 = jptch3["Axis Point 2"].as<std::string>();
            else {
              std::cerr << "Please provide Axis Point 2 (i.e (1 1 1))"
                        << std::endl;
              throw;
            }
            if (jptch3.contains("Radius"))
              SSS.rad = jptch3["Radius"].as<double>();
            else {
              std::cerr << "Please provide Radius" << std::endl;
              throw;
            }
          } else if (SSS.searchableName == "searchableSphere") {
            if (jptch3.contains("Center"))
              SSS.shapeParameters1 = jptch3["Center"].as<std::string>();
            else {
              std::cerr << "Please provide Center (i.e (1 1 1))" << std::endl;
              throw;
            }
            if (jptch3.contains("Radius"))
              SSS.rad = jptch3["Radius"].as<double>();
            else {
              std::cerr << "Please provide Radius" << std::endl;
              throw;
            }
          } else {
            std::cerr << SSS.searchableName << " is not supported yet!"
                      << std::endl;
            throw;
          }
          (params->srchShape).push_back(SSS);
        }
      }

      // Castellated Mesh Controls Inputs
      jsoncons::json castMeshParams = shmparams["Castellated Mesh Controls"];

      if (castMeshParams.contains("CellZones"))
        params->_withCellZones = castMeshParams["CellZones"].as<bool>();
      else {
        std::cerr << "Please specify on/off choice using \"CellZones\""
                  << "Keyword!\n"
                  << std::endl;
        throw;
      }
      if (castMeshParams.contains("RegionRefine"))
        params->_withGeomRefReg = castMeshParams["RegionRefine"].as<bool>();
      else {
        std::cerr << "Please specify on/off choice using \"RegionRefine\""
                  << "Keyword!\n"
                  << std::endl;
        throw;
      }
      if (castMeshParams.contains("SurfaceRefine"))
        params->_withSurfRefReg = castMeshParams["SurfaceRefine"].as<bool>();
      else {
        std::cerr << "Please specify on/off choice using \"SurfaceRefine\""
                  << "Keyword!\n"
                  << std::endl;
        throw;
      }
      if (castMeshParams.contains("GeneralGapLevelIncrement"))
        params->castMeshGpLvl =
            castMeshParams["GeneralGapLevelIncrement"].as<int>();
      else {
        params->castMeshGpLvl = 1;
      }
      if (castMeshParams.contains("locInMesh_adjust"))
        locAdjust = castMeshParams["locInMesh_adjust"].as<double>();
      else
        locAdjust = -1.0;
      if (castMeshParams.contains("maxLocalCells"))
        params->maxLCells = castMeshParams["maxLocalCells"].as<int>();
      else {
        params->maxLCells = 2000000;
      }
      if (castMeshParams.contains("maxGlobalCells"))
        params->maxGCells = castMeshParams["maxGlobalCells"].as<int>();
      else {
        params->maxGCells = 4000000;
      }
      if (castMeshParams.contains("minRefCells"))
        params->minRefCells = castMeshParams["minRefCells"].as<int>();
      else {
        params->minRefCells = 0;
      }
      if (castMeshParams.contains("nCellsBetweenLevels"))
        params->cellsBetnLvls = castMeshParams["nCellsBetweenLevels"].as<int>();
      else {
        params->cellsBetnLvls = 3;
      }
      if (castMeshParams.contains("surfaceRefinementLvlMin"))
        params->refSurfLvlMin =
            castMeshParams["surfaceRefinementLvlMin"].as<int>();
      else {
        params->refSurfLvlMin = 0;
      }
      if (castMeshParams.contains("surfaceRefinementLvlMax"))
        params->refSurfLvlMax =
            castMeshParams["surfaceRefinementLvlMax"].as<int>();
      else {
        params->refSurfLvlMax = 0;
      }
      if (castMeshParams.contains("resolveFeatureAngle"))
        params->featAngle = castMeshParams["resolveFeatureAngle"].as<double>();
      else {
        params->featAngle = 60;
      }
      if (castMeshParams.contains("gapLevelIncrement"))
        params->gPLvlInc = castMeshParams["gapLevelIncrement"].as<int>();
      else
        params->gPLvlInc = 1;
      if (castMeshParams.contains("planarAngle"))
        params->planarAngle = castMeshParams["planarAngle"].as<int>();
      else
        params->planarAngle = 1;
      if (castMeshParams.contains("locationInMeshX"))
        params->locMeshX = castMeshParams["locationInMeshX"].as<double>();
      else {
        params->locMeshX = -1;
      }
      if (castMeshParams.contains("locationInMeshY"))
        params->locMeshY = castMeshParams["locationInMeshY"].as<double>();
      else {
        params->locMeshY = -1;
      }
      if (castMeshParams.contains("locationInMeshZ"))
        params->locMeshZ = castMeshParams["locationInMeshZ"].as<double>();
      else {
        params->locMeshZ = -1;
      }
      if (castMeshParams.contains("allowFreeStandingZoneFaces"))
        params->_alwFreeZone =
            castMeshParams["allowFreeStandingZoneFaces"].as<bool>();
      else {
        params->_alwFreeZone = true;
      }

      // Refinement Surfaces
      std::string cap3 = "SurfaceRefinementRegions";
      if (castMeshParams.contains(cap3)) {
        for (auto jptch3 : castMeshParams[cap3].array_range()) {
          shmSurfRefine surfRef;

          if (jptch3.contains("Patch Name"))
            surfRef.refPatchNm = jptch3["Patch Name"].as<std::string>();
          else {
            surfRef.refPatchNm = params->singleSolidPatch;
          }

          if (jptch3.contains("Patch Type"))
            surfRef.patchType = jptch3["Patch Type"].as<int>();
          else {
            surfRef.patchType = "NO";
          }

          if (jptch3.contains("MinLevel"))
            surfRef.minLvl = jptch3["MinLevel"].as<double>();
          else {
            surfRef.minLvl = 1;
          }

          if (jptch3.contains("MaxLevel"))
            surfRef.maxLvl = jptch3["MaxLevel"].as<int>();
          else {
            std::cerr << "Please define maximum surface refinement"
                      << "\"MaxLevel\"\n"
                      << std::endl;
            throw;
          }
          (params->surfRefs).push_back(surfRef);
        }
      }

      // Region Refinement
      std::string cap2 = "GeomRefinementRegions";
      if (shmparams.contains(cap2)) {
        for (auto jptch2 : shmparams[cap2].array_range()) {
          shmRegionRefine geomRef;

          if (jptch2.contains("Patch Name"))
            geomRef.patchNm = jptch2["Patch Name"].as<std::string>();
          else {
            std::cerr << "Please define \"PatchName\"\n" << std::endl;
            throw;
          }
          if (jptch2.contains("Mode"))
            geomRef.mode = jptch2["Mode"].as<std::string>();
          else {
            std::cerr << "Please define \"Mode\"\n" << std::endl;
            throw;
          }
          if (jptch2.contains("MinLevel"))
            geomRef.minLvl = jptch2["MinLevel"].as<double>();
          else {
            std::cerr << "Please define \"MinLevel\"\n" << std::endl;
            throw;
          }
          if (jptch2.contains("MaxLevel"))
            geomRef.maxLvl = jptch2["MaxLevel"].as<int>();
          else {
            std::cerr << "Please define \"MaxLevel\"\n" << std::endl;
            throw;
          }

          (params->geomRefs).push_back(geomRef);
        }
      }

      // Features
      std::string cap6 = "Feature File";
      if (castMeshParams.contains(cap6)) {
        params->_withFeatureEdgeFile = true;
        for (auto jptch3 : castMeshParams[cap6].array_range()) {
          shmFeatureEdgeRef ftrOne;
          if (jptch3.contains("File Name"))
            ftrOne.fileName = jptch3["File Name"].as<int>();
          else {
            std::cerr << "Please provide feature file name" << std::endl;
            throw;
          }

          if (jptch3.contains("MinLevel"))
            ftrOne.minLvl = jptch3["MinLevel"].as<double>();
          else {
            ftrOne.minLvl = 1;
          }

          if (jptch3.contains("MaxLevel"))
            ftrOne.maxLvl = jptch3["MaxLevel"].as<int>();
          else {
            std::cerr << "Please define maximum surface refinement"
                      << "\"MaxLevel\"\n"
                      << std::endl;
            throw;
          }
          (params->ftrEdge).push_back(ftrOne);
        }
      }

      // Snapping Controls
      if (inputjson["Meshing Parameters"].contains("Snapping Controls")) {
        jsoncons::json snapParams = shmparams["Snapping Controls"];
      
        if (snapParams.contains("nSmoothPatch"))
          params->snapSmthPatch = snapParams["nSmoothPatch"].as<int>();
        else {
          params->snapSmthPatch = 4;
        }
        if (snapParams.contains("tolerance"))
          params->snapTol = snapParams["tolerance"].as<double>();
        else {
          params->snapTol = 0.5;
        }
        if (snapParams.contains("snapSolveIter"))
          params->solveSnapIter = snapParams["snapSolveIter"].as<int>();
        else {
          params->solveSnapIter = 200;
        }
        if (snapParams.contains("snapRelaxIter"))
          params->relaxSnapIter = snapParams["snapRelaxIter"].as<int>();
        else {
          params->relaxSnapIter = 6;
        }
        if (snapParams.contains("nFeatureSnapIter"))
          params->nFeatureSnapIter = snapParams["nFeatureSnapIter"].as<int>();
        else
          params->nFeatureSnapIter = 10;
        if (snapParams.contains("implicitFeatureSnap"))
          params->implicitFeatureSnap =
              snapParams["implicitFeatureSnap"].as<bool>();
        else
          params->implicitFeatureSnap = false;
        if (snapParams.contains("explicitFeatureSnap"))
          params->explicitFeatureSnap =
              snapParams["explicitFeatureSnap"].as<bool>();
        else
          params->explicitFeatureSnap = true;
        if (snapParams.contains("multiRegionFeatureSnap"))
          params->multiRegionFeatureSnap =
              snapParams["multiRegionFeatureSnap"].as<bool>();
        else
          params->multiRegionFeatureSnap = false;
      }

      // Layer controls
      if (inputjson["Meshing Parameters"].contains("Mesh Layers Controls")) {
        jsoncons::json layerParams = shmparams["Mesh Layers Controls"];
      
        if (layerParams.contains("relativeSizes"))
          params->_relSize = layerParams["relativeSizes"].as<bool>();
        else {
          params->_relSize = 1;
        }
        if (layerParams.contains("expansionRatio"))
          params->expRatio = layerParams["expansionRatio"].as<double>();
        else {
          params->expRatio = 1.3;
        }
        if (layerParams.contains("finalLayerThickness"))
          params->finLThick = layerParams["finalLayerThickness"].as<double>();
        else {
          params->finLThick = 1.0;
        }
        if (layerParams.contains("minThickness"))
          params->minThick = layerParams["minThickness"].as<double>();
        else {
          params->minThick = 0.1;
        }
        if (layerParams.contains("firstLayerThickness"))
          params->firstLyrThickness =
              layerParams["firstLayerThickness"].as<double>();
        else
          params->firstLyrThickness = -1.0;
        if (layerParams.contains("thickness"))
          params->thickness = layerParams["thickness"].as<double>();
        else
          params->thickness = -1.0;
        if (layerParams.contains("nGrow"))
          params->nGrow = layerParams["nGrow"].as<int>();
        else {
          params->nGrow = 0;
        }
        if (layerParams.contains("featureAngle"))
          params->lyrFeatAngle = layerParams["featureAngle"].as<double>();
        else {
          params->lyrFeatAngle = 30;
        }
        if (layerParams.contains("nRelaxIter"))
          params->lyrRelaxIter = layerParams["nRelaxIter"].as<int>();
        else {
          params->lyrRelaxIter = 3;
        }
        if (layerParams.contains("nSmoothSurfaceNormals"))
          params->lyrSmthSurfNorm =
              layerParams["nSmoothSurfaceNormals"].as<int>();
        else {
          params->lyrSmthSurfNorm = 1;
        }
        if (layerParams.contains("nSmoothNormals"))
          params->lyrSmthNorm = layerParams["nSmoothNormals"].as<int>();
        else {
          params->lyrSmthNorm = 3;
        }
        if (layerParams.contains("nSmoothThickness"))
          params->lyrSmthThick = layerParams["nSmoothThickness"].as<int>();
        else {
          params->lyrSmthThick = 2;
        }
        if (layerParams.contains("maxFaceThicknessRatio"))
          params->lyrMaxFcTR = layerParams["maxFaceThicknessRatio"].as<double>();
        else {
          params->lyrMaxFcTR = 0.5;
        }
        if (layerParams.contains("maxThicknessToMedialRatio"))
          params->lyrMaxThickTMR =
              layerParams["maxThicknessToMedialRatio"].as<double>();
        else {
          params->lyrMaxThickTMR = 1.0;
        }
        if (layerParams.contains("minMedialAxisAngle"))
          params->lyrMinMedAngl = layerParams["minMedialAxisAngle"].as<double>();
        else {
          params->lyrMinMedAngl = 90;
        }
        if (layerParams.contains("nBufferCellsNoExtrude"))
          params->lyrBuffrCells = layerParams["nBufferCellsNoExtrude"].as<int>();
        else {
          params->lyrBuffrCells = 0;
        }
        if (layerParams.contains("nLayerIter"))
          params->lyrIter = layerParams["nLayerIter"].as<int>();
        else {
          params->lyrIter = 50;
        }
        if (layerParams.contains("nRelaxedIter"))
          params->nRelaxedIter = layerParams["nRelaxedIter"].as<int>();
        else {
          params->nRelaxedIter = 20;
        }
        if (layerParams.contains("slipFeatureAngle"))
          params->slipFeatureAngle = layerParams["slipFeatureAngle"].as<int>();
        else {
          params->slipFeatureAngle = 20;
        }
        if (layerParams.contains("nMedialAxisIter"))
          params->nMedialAxisIter = layerParams["nMedialAxisIter"].as<int>();
        else
          params->nMedialAxisIter = -1;
        if (layerParams.contains("nSmoothDisplacement"))
          params->nSmoothDisplacement =
              layerParams["nSmoothDisplacement"].as<int>();
        else
          params->nSmoothDisplacement = -1;

        // Layers
        std::string cap7 = "Layers";
        if (layerParams.contains(cap7)) {
          for (auto jptch3 : layerParams[cap7].array_range()) {
            shmLayers lyrOne;
            if (jptch3.contains("Patch Name"))
              lyrOne.patchName = jptch3["File Name"].as<int>();
            else {
              std::cerr << "Please provide patch name for layers" << std::endl;
              throw;
            }

            if (jptch3.contains("nSurfaceLayers"))
              lyrOne.nSurfaceLayers = jptch3["nSurfaceLayers"].as<int>();
            else {
              lyrOne.nSurfaceLayers = 1;
            }

            if (jptch3.contains("expansionRatio"))
              lyrOne.expansionRatio = jptch3["expansionRatio"].as<int>();
            else {
              lyrOne.expansionRatio = 1;
            }

            if (jptch3.contains("finalLayerThickness"))
              lyrOne.finalLayerThickness =
                  jptch3["finalLayerThickness"].as<int>();
            else {
              lyrOne.finalLayerThickness = 1;
            }

            if (jptch3.contains("firstLayerThickness"))
              lyrOne.firstLyrThickness =
                  jptch3["firstLayerThickness"].as<double>();
            else {
              lyrOne.firstLyrThickness = -1.0;
            }

            if (jptch3.contains("thickness"))
              lyrOne.finalLayerThickness = jptch3["thickness"].as<double>();
            else {
              lyrOne.thickness = -1.0;
            }

            if (jptch3.contains("minThickness"))
              lyrOne.minThickness = jptch3["minThickness"].as<int>();
            else {
              lyrOne.minThickness = 1;
            }
            params->layerVec.push_back(lyrOne);
          }
        }
      }

      // Mesh Quality Controls
      if (inputjson["Meshing Parameters"].contains("Mesh Quality Controls")) {
        jsoncons::json qcMeshParams = shmparams["Mesh Quality Controls"];

        if (qcMeshParams.contains("maxNonOrtho"))
          params->qcMaxNOrtho = qcMeshParams["maxNonOrtho"].as<int>();
        else {
          params->qcMaxNOrtho = 65;
        }
        if (qcMeshParams.contains("maxBoundarySkewness"))
          params->qcMaxBndrySkew =
              qcMeshParams["maxBoundarySkewness"].as<double>();
        else {
          params->qcMaxBndrySkew = 20;
        }
        if (qcMeshParams.contains("maxInternalSkewness"))
          params->qcMaxIntSkew = qcMeshParams["maxInternalSkewness"].as<double>();
        else {
          params->qcMaxIntSkew = 4;
        }
        if (qcMeshParams.contains("maxConcave"))
          params->qcMaxConc = qcMeshParams["maxConcave"].as<double>();
        else {
          params->qcMaxConc = 80;
        }
        if (qcMeshParams.contains("minVol"))
          params->qcMinVol = qcMeshParams["minVol"].as<double>();
        else {
          params->qcMinVol = 1e-13;
        }
        if (qcMeshParams.contains("minTetQuality"))
          params->qcMinTetQ = qcMeshParams["minTetQuality"].as<double>();
        else {
          params->qcMinTetQ = 1e-15;
        }
        if (qcMeshParams.contains("minArea"))
          params->qcMinArea = qcMeshParams["minArea"].as<double>();
        else {
          params->qcMinArea = -1;
        }
        if (qcMeshParams.contains("minTwist"))
          params->qcMinTwist = qcMeshParams["minTwist"].as<double>();
        else {
          params->qcMinTwist = 0.02;
        }
        if (qcMeshParams.contains("minFaceWeight"))
          params->qcMinFaceW = qcMeshParams["minFaceWeight"].as<double>();
        else {
          params->qcMinFaceW = 0.05;
        }
        if (qcMeshParams.contains("minVolRatio"))
          params->qcMinVolRto = qcMeshParams["minVolRatio"].as<double>();
        else {
          params->qcMinVolRto = 0.01;
        }
        if (qcMeshParams.contains("minDeterminant"))
          params->qcMinDet = qcMeshParams["minDeterminant"].as<double>();
        else {
          params->qcMinDet = 0.001;
        }
        if (qcMeshParams.contains("minTriangleTwist"))
          params->qcMinTrTwist = qcMeshParams["minTriangleTwist"].as<double>();
        else {
          params->qcMinTrTwist = -1;
        }
        if (qcMeshParams.contains("qcnSmoothScale"))
          params->qcSmthScale = qcMeshParams["qcnSmoothScale"].as<int>();
        else {
          params->qcSmthScale = 5;
        }
        if (qcMeshParams.contains("errorReduction"))
          params->qcErrRedctn = qcMeshParams["errorReduction"].as<double>();
        else {
          params->qcErrRedctn = 0.75;
        }
        if (shmparams.contains("mergeTolerance"))
          params->mergeTol = shmparams["mergeTolerance"].as<double>();
        else {
          params->mergeTol = 1e-06;
        }
      }
    }

    blockMeshParams *bmparams = new blockMeshParams();
    std::string defaults3 =
        inputjson["Meshing Parameters"]["blockMesh Parameters"]
            .as<std::string>();

    jsoncons::json bmshparams =
        inputjson["Meshing Parameters"]["blockMesh Parameters"];

    if (inputjson["Mesh File Options"].contains("Input Dict File"))
      bmparams->_ownBlockMshDict =
          inputjson["Mesh File Options"]["Input Dict File"].as<bool>();

    if (bmshparams.contains("Block Geometry"))
      bmparams->_isBlock = bmshparams["Block Geometry"].as<bool>();
    else {
      std::cerr << "Define your choice of geometry using bool keys"
                << "\n"
                << std::endl;
      throw;
    }
    if (bmshparams.contains("Sphere Geometry"))
      bmparams->_isSphere = bmshparams["Sphere Geometry"].as<bool>();
    else {
      std::cerr << "Define your choice of geometry using bool keys"
                << "\n"
                << std::endl;
      throw;
    }
    if (bmshparams.contains("Cylinder/Tapered_Cone Geometry"))
      bmparams->_isCylinder_TCone =
          bmshparams["Cylinder/Tapered_Cone Geometry"].as<bool>();
    else {
      std::cerr << "Define your choice of geometry using bool keys"
                << "\n"
                << std::endl;
      throw;
    }
    if (bmshparams.contains("scaleToMeters"))
      bmparams->cnvrtToMeters = bmshparams["scaleToMeters"].as<double>();
    else {
      std::cerr << "Define your choice of geometry using bool keys"
                << "\n"
                << std::endl;
      throw;
    }

    if (bmshparams.contains("Cell_Size")) {
      bmparams->_cellSizeDefined = true;
      bmparams->cellSize = bmshparams["Cell_Size"].as<double>();
    } else {
      bmparams->cellSize = -1;
      bmparams->_cellSizeDefined = false;
    }

    if (bmshparams.contains("XdirectionCells"))
      bmparams->cellsXDir = bmshparams["XdirectionCells"].as<int>();
    else {
      if (bmparams->_cellSizeDefined) {
        // Nothing
      } else {
        std::cerr << "Define cell numbers in X direction"
                  << "\n"
                  << std::endl;
        throw;
      }
    }
    if (bmshparams.contains("YdirectionCells"))
      bmparams->cellsYDir = bmshparams["YdirectionCells"].as<int>();
    else {
      if (bmparams->_cellSizeDefined) {
        // Nothing
      } else {
        std::cerr << "Define cell numbers in Y direction"
                  << "\n"
                  << std::endl;
        throw;
      }
    }
    if (bmshparams.contains("ZdirectionCells"))
      bmparams->cellsZDir = bmshparams["ZdirectionCells"].as<int>();
    else {
      if (bmparams->_cellSizeDefined) {
        // Nothing
      } else {
        std::cerr << "Define cell numbers in Z direction"
                  << "\n"
                  << std::endl;
        throw;
      }
    }

    if (bmshparams.contains("Block Parameters")) {
      if (bmshparams["Block Parameters"].contains("Auto_Generate")) {
        bmparams->_autoGenerateBox = true;

        if (inputjson["Mesh File Options"].contains("Input Geometry File"))
          bmparams->packFileName =
              inputjson["Mesh File Options"]["Input Geometry File"]
                  .as<std::string>();
        else {
          bmparams->packFileName = ifname + ".stl";
        }

        if (bmshparams["Block Parameters"]["Auto_Generate"].contains(
                "Offset_XDir"))
          bmparams->offsetX =
              bmshparams["Block Parameters"]["Auto_Generate"]["Offset_XDir"]
                  .as<double>();
        else {
          bmparams->offsetX = 0.1;
        }
        if (bmshparams["Block Parameters"]["Auto_Generate"].contains(
                "Offset_YDir"))
          bmparams->offsetY =
              bmshparams["Block Parameters"]["Auto_Generate"]["Offset_YDir"]
                  .as<double>();
        else {
          bmparams->offsetY = 0.1;
        }
        if (bmshparams["Block Parameters"]["Auto_Generate"].contains(
                "Offset_ZDir"))
          bmparams->offsetZ =
              bmshparams["Block Parameters"]["Auto_Generate"]["Offset_ZDir"]
                  .as<double>();
        else {
          bmparams->offsetZ = 0.1;
        }

      } else {
        bmparams->_autoGenerateBox = false;
      }

      if (bmshparams["Block Parameters"].contains("X1"))
        bmparams->initX = bmshparams["Block Parameters"]["X1"].as<double>();
      else {
        if (!(bmparams->_autoGenerateBox)) {
          std::cerr << "Define initial point for block (X Y Z)!\n" << std::endl;
          throw;
        } else {
          // Nothing
        }
      }
      if (bmshparams["Block Parameters"].contains("Y1"))
        bmparams->initY = bmshparams["Block Parameters"]["Y1"].as<double>();
      else {
        if (!(bmparams->_autoGenerateBox)) {
          std::cerr << "Define initial point for block (X Y Z)!\n" << std::endl;
          throw;
        } else {
          // Nothing
        }
      }
      if (bmshparams["Block Parameters"].contains("Z1"))
        bmparams->initZ = bmshparams["Block Parameters"]["Z1"].as<double>();
      else {
        if (!(bmparams->_autoGenerateBox)) {
          std::cerr << "Define initial point for block (X Y Z)!\n" << std::endl;
          throw;
        } else {
          // Nothing
        }
      }
      if (bmshparams["Block Parameters"].contains("LengthX"))
        bmparams->lenX = bmshparams["Block Parameters"]["LengthX"].as<double>();
      else {
        if (!(bmparams->_autoGenerateBox)) {
          std::cerr << "Define initial point for block (X Y Z)!\n" << std::endl;
          throw;
        } else {
          // Nothing
        }
      }
      if (bmshparams["Block Parameters"].contains("LengthY"))
        bmparams->lenY = bmshparams["Block Parameters"]["LengthY"].as<double>();
      else {
        if (!(bmparams->_autoGenerateBox)) {
          std::cerr << "Define initial point for block (X Y Z)!\n" << std::endl;
          throw;
        } else {
          // Nothing
        }
      }
      if (bmshparams["Block Parameters"].contains("LengthZ"))
        bmparams->lenZ = bmshparams["Block Parameters"]["LengthZ"].as<double>();
      else {
        if (!(bmparams->_autoGenerateBox)) {
          std::cerr << "Define initial point for block (X Y Z)!\n" << std::endl;
          throw;
        } else {
          // Nothing
        }
      }
      if (bmshparams["Block Parameters"].contains("GradingXdir"))
        bmparams->smplGradingX =
            bmshparams["Block Parameters"]["GradingXdir"].as<double>();
      else {
        bmparams->smplGradingX = 1;
      }
      if (bmshparams["Block Parameters"].contains("GradingYdir"))
        bmparams->smplGradingY =
            bmshparams["Block Parameters"]["GradingYdir"].as<double>();
      else {
        bmparams->smplGradingY = 1;
      }
      if (bmshparams["Block Parameters"].contains("GradingZdir"))
        bmparams->smplGradingZ =
            bmshparams["Block Parameters"]["GradingZdir"].as<double>();
      else {
        bmparams->smplGradingZ = 1;
      }
    }

    if (bmshparams.contains("Sphere Parameters")) {
      if (bmshparams["Sphere Parameters"].contains("Center X"))
        bmparams->centerX =
            bmshparams["Sphere Parameters"]["Center X"].as<double>();
      else {
        std::cerr << "Define sphere center!\n" << std::endl;
        throw;
      }
      if (bmshparams["Sphere Parameters"].contains("Center Y"))
        bmparams->centerY =
            bmshparams["Sphere Parameters"]["Center Y"].as<double>();
      else {
        std::cerr << "Define sphere center!\n" << std::endl;
        throw;
      }
      if (bmshparams["Sphere Parameters"].contains("Center Z"))
        bmparams->centerZ =
            bmshparams["Sphere Parameters"]["Center Z"].as<double>();
      else {
        std::cerr << "Define sphere center!\n" << std::endl;
        throw;
      }
      if (bmshparams["Sphere Parameters"].contains("Radius"))
        bmparams->radius =
            bmshparams["Sphere Parameters"]["Radius"].as<double>();
      else {
        std::cerr << "Define sphere radius!\n" << endl;
        throw;
      }
      if (bmshparams["Sphere Parameters"].contains("GradingXdir"))
        bmparams->sphrGradingX =
            bmshparams["Sphere Parameters"]["GradingXdir"].as<double>();
      else {
        bmparams->sphrGradingX = 1;
      }
      if (bmshparams["Sphere Parameters"].contains("GradingYdir"))
        bmparams->sphrGradingY =
            bmshparams["Sphere Parameters"]["GradingYdir"].as<double>();
      else {
        bmparams->sphrGradingY = 1;
      }
      if (bmshparams["Sphere Parameters"].contains("GradingXdir"))
        bmparams->sphrGradingZ =
            bmshparams["Sphere Parameters"]["GradingZdir"].as<double>();
      else {
        bmparams->sphrGradingZ = 1;
      }
    }

    if (bmshparams.contains("Cylinder/Tapered_Cone Parameters")) {
      if (bmshparams["Cylinder/Tapered_Cone Parameters"].contains("Center X"))
        bmparams->centerCyl[0] =
            bmshparams["Cylinder/Tapered_Cone Parameters"]["Center X"]
                .as<double>();
      else {
        std::cerr << "Define center point for cylinder axis\n" << std::endl;
        throw;
      }
      if (bmshparams["Cylinder/Tapered_Cone Parameters"].contains("Center Y"))
        bmparams->centerCyl[1] =
            bmshparams["Cylinder/Tapered_Cone Parameters"]["Center Y"]
                .as<double>();
      else {
        std::cerr << "Define center point for cylinder axis\n" << std::endl;
        throw;
      }
      if (bmshparams["Cylinder/Tapered_Cone Parameters"].contains("Center Z"))
        bmparams->centerCyl[2] =
            bmshparams["Cylinder/Tapered_Cone Parameters"]["Center Z"]
                .as<double>();
      else {
        std::cerr << "Define center point for cylinder axis\n" << std::endl;
        throw;
      }
      if (bmshparams["Cylinder/Tapered_Cone Parameters"].contains("Radius1"))
        bmparams->radius1 =
            bmshparams["Cylinder/Tapered_Cone Parameters"]["Radius1"]
                .as<double>();

      if (bmshparams["Cylinder/Tapered_Cone Parameters"].contains("Radius2")) {
        bmparams->radius2 =
            bmshparams["Cylinder/Tapered_Cone Parameters"]["Radius2"]
                .as<double>();
      } else {
        bmparams->radius2 =
            bmshparams["Cylinder/Tapered_Cone Parameters"]["Radius1"]
                .as<double>();
      }

      if (bmshparams["Cylinder/Tapered_Cone Parameters"].contains(
              "GradingXdir"))
        bmparams->cylGrading[0] =
            bmshparams["Cylinder/Tapered_Cone Parameters"]["GradingXdir"]
                .as<double>();
      else {
        bmparams->cylGrading[0] = 1;
      }
      if (bmshparams["Cylinder/Tapered_Cone Parameters"].contains(
              "GradingYdir"))
        bmparams->cylGrading[1] =
            bmshparams["Cylinder/Tapered_Cone Parameters"]["GradingYdir"]
                .as<double>();
      else {
        bmparams->cylGrading[1] = 1;
      }
      if (bmshparams["Cylinder/Tapered_Cone Parameters"].contains(
              "GradingXdir"))
        bmparams->cylGrading[2] =
            bmshparams["Cylinder/Tapered_Cone Parameters"]["GradingZdir"]
                .as<double>();
      else {
        bmparams->cylGrading[2] = 1;
      }

      if (bmshparams["Cylinder/Tapered_Cone Parameters"].contains("Height"))
        bmparams->height =
            bmshparams["Cylinder/Tapered_Cone Parameters"]["Height"]
                .as<double>();
      else {
        std::cerr << "Define cylinder height\n" << endl;
      }
    }

    MeshManipulationFoamParams *mparams = new MeshManipulationFoamParams();
    // std::string defaults4 =
    //     inputjson["Meshing Parameters"]["MeshManipulation Parameters"]
    //         .as<std::string>();

    jsoncons::json pmshparams;
    if (inputjson["Meshing Parameters"].contains("MeshManipulation Parameters")) {
      pmshparams = inputjson["Meshing Parameters"]["MeshManipulation Parameters"];

      // Parameter parsing starts here
      if (pmshparams.contains("Enable SurfLambdaMuSmooth"))
        mparams->_doSurfaceLMSmth =
            pmshparams["Enable SurfLambdaMuSmooth"].as<bool>();
      if (pmshparams.contains("Enable splitMeshRegions"))
        mparams->_doSplitMshRegs =
            pmshparams["Enable splitMeshRegions"].as<bool>();
      if (pmshparams.contains("Enable MergeMeshes"))
        mparams->_doMergeMsh = pmshparams["Enable MergeMeshes"].as<bool>();
      if (pmshparams.contains("Enable CreatePatch"))
        mparams->_doCreatePtchs = pmshparams["Enable CreatePatch"].as<bool>();
      if (pmshparams.contains("Enable foamToSurface"))
        mparams->_doFoam2Surf = pmshparams["Enable foamToSurface"].as<bool>();
      if (pmshparams.contains("Enable surfaceSplitByTopology"))
        mparams->_doSurfSplit =
            pmshparams["Enable surfaceSplitByTopology"].as<bool>();

      if (pmshparams.contains("SurfLambdaMuSmooth Parameters")) {
        if (pmshparams["SurfLambdaMuSmooth Parameters"].contains(
                "AddFeatureFile?"))
          mparams->_addFeatureFile =
              pmshparams["SurfLambdaMuSmooth Parameters"]["AddFeatureFile?"]
                  .as<bool>();

        if (pmshparams["SurfLambdaMuSmooth Parameters"].contains(
                "Input STL File"))
          mparams->slmssurfaceFile =
              pmshparams["SurfLambdaMuSmooth Parameters"]["Input STL File"]
                  .as<std::string>();

        if (pmshparams["SurfLambdaMuSmooth Parameters"].contains(
                "Output STL File"))
          mparams->slmsoutputFile =
              pmshparams["SurfLambdaMuSmooth Parameters"]["Output STL File"]
                  .as<std::string>();

        if (pmshparams["SurfLambdaMuSmooth Parameters"].contains("Lambda"))
          mparams->lambda =
              pmshparams["SurfLambdaMuSmooth Parameters"]["Lambda"].as<double>();

        if (pmshparams["SurfLambdaMuSmooth Parameters"].contains("Mu"))
          mparams->mu =
              pmshparams["SurfLambdaMuSmooth Parameters"]["Mu"].as<double>();

        if (pmshparams["SurfLambdaMuSmooth Parameters"].contains(
                "Smoothing Interations"))
          mparams->slmsIterations =
              pmshparams["SurfLambdaMuSmooth Parameters"]["Smoothing Interations"]
                  .as<int>();
      }

      if (pmshparams.contains("splitMeshRegions Parameters")) {
        if (pmshparams["splitMeshRegions Parameters"].contains("overwrite?"))
          mparams->_overwriteMsh =
              pmshparams["splitMeshRegions Parameters"]["overwrite?"].as<bool>();
        else
          mparams->_overwriteMsh = true;

        if (pmshparams["splitMeshRegions Parameters"].contains("usecellZones?"))
          mparams->_cellZones =
              pmshparams["splitMeshRegions Parameters"]["usecellZones?"]
                  .as<bool>();
        else
          mparams->_cellZones = true;
      } else {
        mparams->_overwriteMsh = true;
        mparams->_cellZones = true;
      }

      mparams->addCase = params->singleSolidPatch;
      mparams->masterCase = "domain1";  // Not Used
      if (pmshparams.contains("mergeMeshes Parameters")) {
        if (pmshparams["mergeMeshes Parameters"].contains("overwrite?"))
          mparams->_overwriteMergeMsh =
              pmshparams["mergeMeshes Parameters"]["overwrite?"].as<bool>();
        else
          mparams->_overwriteMergeMsh = true;

        if (pmshparams["mergeMeshes Parameters"].contains("Master Region Path"))
          mparams->masterCasePath =
              pmshparams["mergeMeshes Parameters"]["Master Region Path"]
                  .as<std::string>();
        else
          mparams->masterCasePath = ".";

        if (pmshparams["mergeMeshes Parameters"].contains("Add Region Path"))
          mparams->addCasePath =
              pmshparams["mergeMeshes Parameters"]["Add Region Path"]
                  .as<std::string>();
        else
          mparams->addCasePath = ".";
      } else {
        mparams->masterCase = "domain1";  // Not Used
        mparams->addCase = params->singleSolidPatch;
        mparams->_overwriteMergeMsh = true;
        mparams->masterCasePath = ".";
        mparams->addCasePath = ".";
      }

      if (pmshparams.contains("createPatch Parameters")) {
        if (pmshparams["createPatch Parameters"].contains(
                "Surrounding PatchName"))
          mparams->surroundingName =
              pmshparams["createPatch Parameters"]["Surrounding PatchName"]
                  .as<std::string>();
        else
          mparams->surroundingName = "Soil";

        if (pmshparams["createPatch Parameters"].contains("Packs PatchName"))
          mparams->packsName =
              pmshparams["createPatch Parameters"]["Packs PatchName"]
                  .as<std::string>();
        else
          mparams->packsName = "Rocks";

        if (pmshparams["createPatch Parameters"].contains(
                "Surrounding PatchType"))
          mparams->srrndngPatchType =
              pmshparams["createPatch Parameters"]["Surrounding PatchType"]
                  .as<std::string>();
        else
          mparams->srrndngPatchType = "wall";

        if (pmshparams["createPatch Parameters"].contains("Packs PatchType"))
          mparams->packsPatchType =
              pmshparams["createPatch Parameters"]["Packs PatchType"]
                  .as<std::string>();
        else
          mparams->packsPatchType = "wall";

        if (pmshparams["createPatch Parameters"].contains("overwrite?"))
          mparams->_overwritecpMsh =
              pmshparams["createPatch Parameters"]["overwrite?"].as<bool>();
        else
          mparams->_overwritecpMsh = true;
      } else {
        mparams->surroundingName = "Soil";
        mparams->packsName = "Rocks";
        mparams->srrndngPatchType = "wall";
        mparams->packsPatchType = "wall";
        mparams->_overwritecpMsh = true;
      }

      if (pmshparams.contains("foamToSurface Parameters")) {
        if (pmshparams["foamToSurface Parameters"].contains("Output File Path"))
          mparams->outSurfName =
              pmshparams["foamToSurface Parameters"]["Output File Path"]
                  .as<std::string>();
      }

      if (pmshparams.contains("surfaceSplitByTopology Parameters")) {
        if (pmshparams["surfaceSplitByTopology Parameters"].contains(
                "Input File"))
          mparams->surfFile =
              pmshparams["surfaceSplitByTopology Parameters"]["Input File"]
                  .as<std::string>();
        else {
          mparams->surfFile = params->geomFileName;
        }
        if (pmshparams["surfaceSplitByTopology Parameters"].contains(
                "Output File"))
          mparams->outSurfFile =
              pmshparams["surfaceSplitByTopology Parameters"]["Output File"]
                  .as<std::string>();
        else {
          mparams->outSurfFile = "SurfaceSplitOut.stl";
        }
      } else {
        mparams->surfFile = params->geomFileName;
      }
    } else {
      // All default options here?
      mparams->_overwriteMsh = true;
      mparams->_cellZones = true;
      mparams->masterCase = "domain1";
      mparams->addCase = params->singleSolidPatch;
      mparams->_overwriteMergeMsh = true;
      mparams->masterCasePath = ".";
      mparams->addCasePath = ".";
      mparams->surroundingName = "Soil";
      mparams->packsName = "Rocks";
      mparams->srrndngPatchType = "wall";
      mparams->packsPatchType = "wall";
      mparams->_overwritecpMsh = true;
      mparams->surfFile = params->geomFileName;
      mparams->outSurfFile = "SurfaceSplitOut.stl";
    }

    // Give all data to PackMesh Driver
    PackMeshDriver *pckmshdrvobj =
        new PackMeshDriver(ifname, mparams, params, bmparams, 
                           ofname_pack, ofname_surrndng, ofname_merged, 
                           useRocpack,locAdjust);
    return pckmshdrvobj;
  } else {
    std::cout << "Mesh generation engine " << meshEngine << " is not supported"
              << std::endl;
    exit(1);
  }
}
#endif

PackMeshDriver *PackMeshDriver::readJSON(const std::string &ifname,
                                         const std::string &ofname,
                                         const jsoncons::json inputjson) {
  std::string meshEngine =
      inputjson["Mesh Generation Engine"].as<std::string>();

  if (!meshEngine.compare("packmesh")) {
    bool wantGeometryOnly = false;
    std::string meshType = inputjson["Mesh Type"].as<std::string>();
    if (meshType == "geometry") wantGeometryOnly = true;

    if (wantGeometryOnly) {
      bool setPeriodicMesh = false;

      std::vector<double> transferMesh;
      transferMesh.resize(3);

      std::vector<double> domainBounds;
      domainBounds.resize(6);

      bool customDomain = false;

      double meshSize =
          inputjson["Meshing Parameters"].get_with_default("Mesh Size", 0.);
      int mshAlgorithm =
          inputjson["Meshing Parameters"].get_with_default("Mesh Algorithm", 1);
      bool scaleVolumes = inputjson["Meshing Parameters"].get_with_default(
          "Scale Geometries", false);
      double scaleValue =
          inputjson["Meshing Parameters"].get_with_default("Scale Value", 1.);
      bool rmvBndryPacks = inputjson["Meshing Parameters"].get_with_default(
          "Remove geometries on boundary", false);
      bool enable2PhysGrps = inputjson["Meshing Parameters"].get_with_default(
          "Enable two physical groups", false);
      bool enableMultiPhysGrps =
          inputjson["Meshing Parameters"].get_with_default(
              "Enable multi physical groups", false);
      bool createCohesive = inputjson["Meshing Parameters"].get_with_default(
          "Create cohesive elements", false);
      bool enablePatches = inputjson["Meshing Parameters"].get_with_default(
          "Enable Patches", false);
      bool setPeriodicGeom = inputjson["Meshing Parameters"].get_with_default(
          "Set Periodic Geometry", false);
      bool enableOutBool = inputjson["Meshing Parameters"].get_with_default(
          "Enable Default Outputs", false);
      bool enablePhysGrpPerShape =
          inputjson["Meshing Parameters"].get_with_default(
              "Enable physical group per shape", false);
      bool preserveSize = inputjson["Meshing Parameters"].get_with_default(
          "Enable Size Preservation", false);
      int refineLevel = inputjson["Meshing Parameters"].get_with_default(
          "Refinement Levels", 0);
      double upperThreshold = inputjson["Meshing Parameters"].get_with_default(
          "Upper Threshold", 0.);
      double lowerThreshold = inputjson["Meshing Parameters"].get_with_default(
          "Lower Threshold", 0.);
      int elemOrder =
          inputjson["Meshing Parameters"].get_with_default("Element Order", 1);

      if (inputjson["Meshing Parameters"].contains("Custom Domain")) {
        customDomain = true;
        if (inputjson["Meshing Parameters"]["Custom Domain"].contains(
                "Initial_X"))
          domainBounds[0] =
              inputjson["Meshing Parameters"]["Custom Domain"]["Initial_X"]
                  .as<double>();

        if (inputjson["Meshing Parameters"]["Custom Domain"].contains(
                "Initial_Y"))
          domainBounds[1] =
              inputjson["Meshing Parameters"]["Custom Domain"]["Initial_Y"]
                  .as<double>();

        if (inputjson["Meshing Parameters"]["Custom Domain"].contains(
                "Initial_Z"))
          domainBounds[2] =
              inputjson["Meshing Parameters"]["Custom Domain"]["Initial_Z"]
                  .as<double>();

        if (inputjson["Meshing Parameters"]["Custom Domain"].contains(
                "Length_X"))
          domainBounds[3] =
              inputjson["Meshing Parameters"]["Custom Domain"]["Length_X"]
                  .as<double>();

        if (inputjson["Meshing Parameters"]["Custom Domain"].contains(
                "Length_Y"))
          domainBounds[4] =
              inputjson["Meshing Parameters"]["Custom Domain"]["Length_Y"]
                  .as<double>();

        if (inputjson["Meshing Parameters"]["Custom Domain"].contains(
                "Length_Z"))
          domainBounds[5] =
              inputjson["Meshing Parameters"]["Custom Domain"]["Length_Z"]
                  .as<double>();
      } else
        customDomain = false;

      if (inputjson["Meshing Parameters"].contains("TransferMesh_X"))
        transferMesh[0] =
            inputjson["Meshing Parameters"]["TransferMesh_X"].as<double>();

      if (inputjson["Meshing Parameters"].contains("TransferMesh_Y"))
        transferMesh[1] =
            inputjson["Meshing Parameters"]["TransferMesh_Y"].as<double>();

      if (inputjson["Meshing Parameters"].contains("TransferMesh_Z"))
        transferMesh[2] =
            inputjson["Meshing Parameters"]["TransferMesh_Z"].as<double>();

      if (mshAlgorithm == 1 || mshAlgorithm == 2 || mshAlgorithm == 5 ||
          mshAlgorithm == 6 || mshAlgorithm == 7 || mshAlgorithm == 8 ||
          mshAlgorithm == 9) {
      } else {
        std::cerr
            << "Valid choices for 2D meshing algorithm are "
            << " 1 (MeshAdapt), 2 (Automatic),  5 (Delunay), "
            << " 6 (Frontal Delunay), 7 (BAMG), 8 (Frontal Delunay for Quads)"
            << " 9 (Packing of Parallelograms)." << std::endl;
        throw;
      }

      PackMeshDriver *pckmshdrvobj = new PackMeshDriver(
          ifname, ofname, scaleVolumes, scaleValue, meshSize, rmvBndryPacks,
          setPeriodicGeom, setPeriodicMesh, enable2PhysGrps,
          enableMultiPhysGrps, wantGeometryOnly, createCohesive, enablePatches,
          transferMesh, customDomain, domainBounds, mshAlgorithm, enableOutBool,
          enablePhysGrpPerShape, refineLevel, upperThreshold, lowerThreshold,
          preserveSize, elemOrder);
      return pckmshdrvobj;
    } else {
      bool setPeriodicMesh = true;

      std::vector<double> transferMesh;
      transferMesh.resize(3);

      std::vector<double> domainBounds;
      domainBounds.resize(6);

      bool customDomain = false;

      double meshSize =
          inputjson["Meshing Parameters"].get_with_default("Mesh Size", 0.);
      int mshAlgorithm =
          inputjson["Meshing Parameters"].get_with_default("Mesh Algorithm", 1);
      bool scaleVolumes = inputjson["Meshing Parameters"].get_with_default(
          "Scale Geometries", false);
      double scaleValue =
          inputjson["Meshing Parameters"].get_with_default("Scale Value", 1.);
      bool rmvBndryPacks = inputjson["Meshing Parameters"].get_with_default(
          "Remove geometries on boundary", false);
      bool enable2PhysGrps = inputjson["Meshing Parameters"].get_with_default(
          "Enable two physical groups", false);
      bool enableMultiPhysGrps =
          inputjson["Meshing Parameters"].get_with_default(
              "Enable multi physical groups", false);
      bool createCohesive = inputjson["Meshing Parameters"].get_with_default(
          "Create cohesive elements", false);
      bool enablePatches = inputjson["Meshing Parameters"].get_with_default(
          "Enable Patches", false);
      bool setPeriodicGeom = inputjson["Meshing Parameters"].get_with_default(
          "Set Periodic Geometry", false);
      bool enableOutBool = inputjson["Meshing Parameters"].get_with_default(
          "Enable Default Outputs", false);
      bool enablePhysGrpPerShape =
          inputjson["Meshing Parameters"].get_with_default(
              "Enable physical group per shape", false);
      bool preserveSize = inputjson["Meshing Parameters"].get_with_default(
          "Enable Size Preservation", false);
      int refineLevel = inputjson["Meshing Parameters"].get_with_default(
          "Refinement Levels", 0);
      double upperThreshold = inputjson["Meshing Parameters"].get_with_default(
          "Upper Threshold", 0.);
      double lowerThreshold = inputjson["Meshing Parameters"].get_with_default(
          "Lower Threshold", 0.);
      int elemOrder =
          inputjson["Meshing Parameters"].get_with_default("Element Order", 1);

      if (inputjson["Meshing Parameters"].contains("Custom Domain")) {
        customDomain = true;
        if (inputjson["Meshing Parameters"]["Custom Domain"].contains(
                "Initial_X"))
          domainBounds[0] =
              inputjson["Meshing Parameters"]["Custom Domain"]["Initial_X"]
                  .as<double>();

        if (inputjson["Meshing Parameters"]["Custom Domain"].contains(
                "Initial_Y"))
          domainBounds[1] =
              inputjson["Meshing Parameters"]["Custom Domain"]["Initial_Y"]
                  .as<double>();

        if (inputjson["Meshing Parameters"]["Custom Domain"].contains(
                "Initial_Z"))
          domainBounds[2] =
              inputjson["Meshing Parameters"]["Custom Domain"]["Initial_Z"]
                  .as<double>();

        if (inputjson["Meshing Parameters"]["Custom Domain"].contains(
                "Length_X"))
          domainBounds[3] =
              inputjson["Meshing Parameters"]["Custom Domain"]["Length_X"]
                  .as<double>();

        if (inputjson["Meshing Parameters"]["Custom Domain"].contains(
                "Length_Y"))
          domainBounds[4] =
              inputjson["Meshing Parameters"]["Custom Domain"]["Length_Y"]
                  .as<double>();

        if (inputjson["Meshing Parameters"]["Custom Domain"].contains(
                "Length_Z"))
          domainBounds[5] =
              inputjson["Meshing Parameters"]["Custom Domain"]["Length_Z"]
                  .as<double>();
      } else
        customDomain = false;

      if (inputjson["Meshing Parameters"].contains("TransferMesh_X"))
        transferMesh[0] =
            inputjson["Meshing Parameters"]["TransferMesh_X"].as<double>();

      if (inputjson["Meshing Parameters"].contains("TransferMesh_Y"))
        transferMesh[1] =
            inputjson["Meshing Parameters"]["TransferMesh_Y"].as<double>();

      if (inputjson["Meshing Parameters"].contains("TransferMesh_Z"))
        transferMesh[2] =
            inputjson["Meshing Parameters"]["TransferMesh_Z"].as<double>();

      if (mshAlgorithm == 1 || mshAlgorithm == 4 || mshAlgorithm == 7 ||
          mshAlgorithm == 9 || mshAlgorithm == 10) {
      } else {
        std::cerr
            << "Valid choices for 3D meshing algorithm are "
            << " 1 (Delunay), 4 (Frontal), 7 (MMG3D), 9 (R-Tree), 10 (HXT)"
            << std::endl;
        throw;
      }

      PackMeshDriver *pckmshdrvobj = new PackMeshDriver(
          ifname, ofname, scaleVolumes, scaleValue, meshSize, rmvBndryPacks,
          setPeriodicGeom, setPeriodicMesh, enable2PhysGrps,
          enableMultiPhysGrps, wantGeometryOnly, createCohesive, enablePatches,
          transferMesh, customDomain, domainBounds, mshAlgorithm, enableOutBool,
          enablePhysGrpPerShape, refineLevel, upperThreshold, lowerThreshold,
          preserveSize, elemOrder);
      return pckmshdrvobj;
    }
  } else {
    std::cout << "Mesh generation engine " << meshEngine << " is not supported"
              << std::endl;
    exit(1);
  }
}
