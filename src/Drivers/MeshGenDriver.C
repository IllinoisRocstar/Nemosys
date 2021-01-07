#include "MeshGenDriver.H"

#include "AuxiliaryFunctions.H"

#include "gmshGen.H"
#include "geoMeshFactory.H"
#include "gmshParams.H"

#ifdef HAVE_NGEN
#  include "netgenGen.H"
#  include "netgenParams.H"
#endif
#ifdef HAVE_SIMMETRIX
#  include "simmetrixGen.H"
#  include "simmetrixParams.H"
#endif
#ifdef HAVE_CFMSH
#  include "blockMeshGen.H"
#  include "blockMeshParams.H"
#  include "cfmeshGen.H"
#  include "cfmeshParams.H"
#  include "snappymeshGen.H"
#  include "snappymeshParams.H"
#endif

// std c++
#include <iostream>
#include <memory>
#include <tuple>
#include <vector>

namespace NEM {
namespace DRV {

// ----------------------------- MeshGen Driver
// -----------------------------------//

MeshGenDriver::MeshGenDriver(const std::string &ifname,
                             const std::string &meshEngine,
                             meshingParams *_params,
                             const std::string &ofname) {
  std::cout << "MeshGenDriver created" << std::endl;
  params = _params;

  meshGen *generator = nullptr;

  if (meshEngine == "netgen") {
#ifdef HAVE_NGEN
    generator = new netgenGen(dynamic_cast<netgenParams *>(params));
    std::string newname = nemAux::trim_fname(ifname, ".vol");
#else
    std::cerr << "NETGEN is not enabled during build."
              << " Build NEMoSys with ENABLE_NETGEN to use this method."
              << std::endl;
    exit(1);
#endif
  } else if (meshEngine == "gmsh") {
    generator =
        new NEM::GEN::gmshGen(dynamic_cast<NEM::GEN::gmshParams *>(params));
    std::string newname = nemAux::trim_fname(ifname, ".msh");
#ifdef HAVE_SIMMETRIX
  } else if (meshEngine == "simmetrix") {
    generator = new simmetrixGen(dynamic_cast<simmetrixParams *>(params));
    std::string newname = nemAux::trim_fname(ifname, ".vtu");
#endif
#ifdef HAVE_CFMSH
  } else if (meshEngine == "cfmesh") {
    generator = new cfmeshGen(dynamic_cast<cfmeshParams *>(params));
    std::string newname = nemAux::trim_fname(ifname, ".vtu");
  } else if (meshEngine == "snappyHexMesh") {
    generator = new snappymeshGen(dynamic_cast<snappymeshParams *>(params));
    std::string newname = nemAux::trim_fname(ifname, ".vtu");
  } else if (meshEngine == "blockMesh") {
    generator = new blockMeshGen(dynamic_cast<blockMeshParams *>(params));
    std::string newname = nemAux::trim_fname(ifname, ".vtu");
#endif
  } else {
    std::cerr << meshEngine << " is not a supported meshing engine"
              << std::endl;
    exit(1);
  }

  if (!generator) {
    std::cerr << "Meshing with engine " << meshEngine << " failed. Aborting."
              << std::endl;
    exit(1);
  }

  int status = generator->createMeshFromSTL(ifname.c_str());

  if (status) {
    std::cerr << "Mesh Engine " << meshEngine << " not recognized" << std::endl;
    exit(1);
  }

  // output file type
  std::string outputType = nemAux::find_ext(ofname);
  std::shared_ptr<meshBase> mesh;
  if (outputType == ".msh" && meshEngine == "gmsh") {
    // createMeshFromSTL (in this case) outputs a .msh file by default
    return;
  }

  // otherwise, resort to meshBase

  if (meshEngine == "netgen") {
    std::string newname = nemAux::trim_fname(ifname, ".vol");
    mesh = meshBase::CreateShared(meshBase::exportVolToVtk(newname));
  } else if (meshEngine == "gmsh") {
    std::string newname = nemAux::trim_fname(ifname, ".msh");
    // conversion from STL using gmsh engine writes a ".msh" file by default
    if (outputType == ".msh") return;
    mesh = meshBase::CreateShared(meshBase::exportGmshToVtk(newname));
  } else {
    // default
    std::string newname = nemAux::trim_fname(ifname, ".vtu");
    mesh = meshBase::CreateShared(
        meshBase::Create(generator->getDataSet(), newname));
  }

  mesh->setFileName(ofname);
  mesh->report();
  mesh->write();
}

std::shared_ptr<meshBase> MeshGenDriver::getNewMesh() const {
  return mesh ? mesh : nullptr;
}

MeshGenDriver::~MeshGenDriver() {
  std::cout << "MeshGenDriver destroyed" << std::endl;
}

MeshGenDriver *MeshGenDriver::readJSON(const jsoncons::json &inputjson) {
  std::string ifname =
      inputjson["Mesh File Options"]["Input Geometry File"].as<std::string>();
  std::string ofname =
      inputjson["Mesh File Options"]["Output Mesh File"].as<std::string>();
  return readJSON(ifname, ofname, inputjson);
}

MeshGenDriver *MeshGenDriver::readJSON(const std::string &ifname,
                                       const std::string &ofname,
                                       const jsoncons::json &inputjson) {
  std::string meshEngine =
      inputjson["Mesh Generation Engine"].as<std::string>();
  if (meshEngine == "netgen") {
    std::string defaults =
        inputjson["Meshing Parameters"]["Netgen Parameters"].as<std::string>();
    if (defaults == "default") {
      auto *params = new netgenParams();
      auto *mshgndrvobj = new MeshGenDriver(ifname, meshEngine, params, ofname);
      return mshgndrvobj;
    } else {
      jsoncons::json ngparams =
          inputjson["Meshing Parameters"]["Netgen Parameters"];

      auto *params = new netgenParams();

      if (ngparams.contains("uselocalh"))
        params->uselocalh = ngparams["uselocalh"].as<bool>();
      if (ngparams.contains("maxh"))
        params->maxh = ngparams["maxh"].as<double>();
      if (ngparams.contains("fineness"))
        params->fineness = ngparams["fineness"].as<double>();
      if (ngparams.contains("grading"))
        params->grading = ngparams["grading"].as<double>();
      if (ngparams.contains("elementsperedge"))
        params->elementsperedge = ngparams["elementsperedge"].as<double>();
      if (ngparams.contains("elementspercurve"))
        params->elementspercurve = ngparams["elementspercurve"].as<double>();
      if (ngparams.contains("closeedgeenable"))
        params->closeedgeenable = ngparams["closeedgeenable"].as<bool>();
      if (ngparams.contains("closeedgefact"))
        params->closeedgefact = ngparams["closeedgefact"].as<double>();
      if (ngparams.contains("second_order"))
        params->second_order = ngparams["second_order"].as<bool>();
      if (ngparams.contains("meshsize_filename"))
        params->meshsize_filename =
            ngparams["meshsize_filename"].as<std::string>();
      if (ngparams.contains("quad_dominated"))
        params->quad_dominated = ngparams["quad_dominated"].as<bool>();
      if (ngparams.contains("optvolmeshenable"))
        params->optvolmeshenable = ngparams["optvolmeshenable"].as<bool>();
      if (ngparams.contains("optsteps_2d"))
        params->optsteps_2d = ngparams["optsteps_2d"].as<int>();
      if (ngparams.contains("optsteps_3d"))
        params->optsteps_3d = ngparams["optsteps_3d"].as<int>();
      if (ngparams.contains("invert_tets"))
        params->invert_tets = ngparams["invert_tets"].as<bool>();
      if (ngparams.contains("invert_trigs"))
        params->invert_trigs = ngparams["invert_trigs"].as<bool>();
      if (ngparams.contains("check_overlap"))
        params->check_overlap = ngparams["check_overlap"].as<bool>();
      if (ngparams.contains("check_overlapping_boundary"))
        params->check_overlapping_boundary =
            ngparams["check_overlapping_boundary"].as<bool>();
      if (ngparams.contains("refine_with_geometry_adaptation"))
        params->refine_with_geom =
            ngparams["refine_with_geometry_adaptation"].as<bool>();
      if (ngparams.contains("refine_without_geometry_adaptation"))
        params->refine_without_geom =
            ngparams["refine_without_geometry_adaptation"].as<bool>();

      auto *mshgndrvobj = new MeshGenDriver(ifname, meshEngine, params, ofname);
      return mshgndrvobj;
    }
  } else if (meshEngine == "gmsh") {
    std::cout << "Gmsh mesh engine selected" << std::endl;

    std::string defaults;
    if (inputjson.contains("Meshing Parameters")) {
      if (inputjson["Meshing Parameters"].contains("Gmsh Parameters")) {
        defaults =
            inputjson["Meshing Parameters"]["Gmsh Parameters"].as_string();
      } else {
        std::cerr << "Error: 'Gmsh Parameters' not found in JSON" << std::endl;
        exit(-1);
      }
    } else {
      std::cerr << "Error: 'Meshing Parameters' not found in JSON" << std::endl;
      exit(-1);
    }
    if (defaults == "default") {
      auto *params = new NEM::GEN::gmshParams();
      auto *mshgndrvobj = new MeshGenDriver(ifname, meshEngine, params, ofname);
      return mshgndrvobj;
    } else {
      jsoncons::json gparams =
          inputjson["Meshing Parameters"]["Gmsh Parameters"];

      auto params = new NEM::GEN::gmshParams();

      // Get the output file name for later
      std::string ofname =
          inputjson["Mesh File Options"]["Output Mesh File"].as_string();
      params->ofname = ofname;

      if (gparams.contains("minSize"))
        params->minSize = gparams["minSize"].as_double();
      if (gparams.contains("maxSize"))
        params->maxSize = gparams["maxSize"].as_double();
      if (gparams.contains("surfaceAlgorithm"))
        params->algo2D = gparams["surfaceAlgorithm"].as_string();
      if (gparams.contains("volumeAlgorithm"))
        params->algo3D = gparams["volumeAlgorithm"].as_string();
      if (gparams.contains("extendSizeFromBoundary"))
        params->extSizeFromBoundary =
            gparams["extendSizeFromBoundary"].as<bool>();
      if (gparams.contains("sizeFromCurvature"))
        params->sizeFromCurvature = gparams["sizeFromCurvature"].as<bool>();
      if (gparams.contains("minElementsPer2Pi"))
        params->minElePer2Pi = gparams["minElementsPer2Pi"].as<int>();
      if (gparams.contains("optimize"))
        params->optimize = gparams["optimize"].as<bool>();
      if (gparams.contains("optimizeThreshold"))
        params->optimizeThreshold = gparams["optimizeThreshold"].as_double();
      if (gparams.contains("elementOrder"))
        params->elementOrder = gparams["elementOrder"].as<int>();
      if (gparams.contains("subdivisionAlgorithm"))
        params->subdivisionAlg = gparams["subdivisionAlgorithm"].as<int>();
      if (gparams.contains("saveAll"))
        params->saveAll = gparams["saveAll"].as<bool>();
      if (gparams.contains("fragment"))
        params->fragmentAll = gparams["fragment"].as<bool>();

      if (gparams.contains("ColorMap")) {
        std::map<std::string, std::string> color2groupMap;
        for (const auto &cm : gparams["ColorMap"].array_range()) {
          std::string color;
          std::string name;
          if (cm.contains("Color")) {
            int r, g, b;
            r = cm["Color"].at(0).as<int>();
            g = cm["Color"].at(1).as<int>();
            b = cm["Color"].at(2).as<int>();
            color = std::to_string(r) + "," + std::to_string(g) + "," +
                    std::to_string(b);
          } else {
            std::cerr << "Error: Keyword 'Color' for Color Map not found."
                      << std::endl;
            throw;
          }
          if (cm.contains("Group")) {
            name = cm["Group"].as<std::string>();
          } else {
            std::cerr << "Error: Keyword 'Group' for Color Map not found."
                      << std::endl;
            throw;
          }
          color2groupMap[color] = name;
        }
        params->mColorMap = true;
        params->color2groupMap = color2groupMap;
      }

      if (gparams.contains("TransfiniteBlocks")) {
        for (const auto &tb : gparams["TransfiniteBlocks"].array_range()) {
          NEM::GEN::TransfiniteBlock block;
          if (tb.contains("Volume")) {
            block.id = tb["Volume"].as<int>();
          } else {
            std::cerr << "Error: Keyword 'id' for transfinite block not found."
                      << std::endl;
            throw;
          }
          if (tb.contains("Axis")) {
            // normalize and store each direction vector
            for (int i = 0; i < 3; ++i) {
              double x = tb["Axis"].at(i).at(0).as<double>();
              double y = tb["Axis"].at(i).at(1).as<double>();
              double z = tb["Axis"].at(i).at(2).as<double>();
              double axis_len = std::sqrt(x * x + y * y + z * z);
              if (axis_len < 1e-3) {
                std::cerr << "Axis " << i << " for block " << block.id
                          << " is too small. Please prescribe an axis with length > 1e-3"
                          << std::endl;
                exit(1);
              }
              block.axis[i][0] = x / axis_len;
              block.axis[i][1] = y / axis_len;
              block.axis[i][2] = z / axis_len;
            }
          } else {
            std::cerr
                << "Error: Keyword 'axis' for transfinite block not found."
                << std::endl;
            throw;
          }

          auto setTransfiniteCurve = [&tb,
                                      &block](const std::string &axis) -> void {
            int idx = -1;
            if (axis == "x") idx = 0;
            if (axis == "y") idx = 1;
            if (axis == "z") idx = 2;
            if (idx == -1) {
              std::cerr << "Invalid axis name '" << axis
                        << "'. Valid axis names are : "
                        << "'x', 'y', and 'z'. Aborting." << std::endl;
              exit(1);
            }
            if (tb.contains(axis)) {
              block.vert[idx] = tb[axis]["Vertices"].as<int>();
              if (tb[axis].contains("Bump")) {
                block.type[idx] = "Bump";
                block.coef[idx] = tb[axis]["Bump"].as<double>();
              } else if (tb[axis].contains("Progression")) {
                block.type[idx] = "Progression";
                block.coef[idx] = tb[axis]["Progression"].as<double>();
              }
            } else {
              std::cerr << "Error : Keyword '" << axis
                        << "' for transfinite block not found." << std::endl;
              throw;
            }
          };

          setTransfiniteCurve("x");
          setTransfiniteCurve("y");
          setTransfiniteCurve("z");

          params->transfiniteBlocks[block.id] = block;
        }
        params->mTransfiniteVolumes = true;
      }

      // Size Fields parsing
      std::string cap = "SizeFields";
      if (gparams.contains(cap)) {
        std::cout << "Parsing SizeFields ...";
        params->mSizeField = true;
        for (const auto &sf : gparams[cap].array_range()) {
          // Instantiate volSizeField struct object
          NEM::GEN::volSizeField sizeField;

          // Get the size field Type
          if (sf.contains("Type"))
            sizeField.type = sf["Type"].as_string();
          else {
            std::cerr << "Error: Keyword 'Type' for Size Field not found."
                      << std::endl;
            throw;
          }

          // Get the size field ID
          if (sf.contains("ID"))
            sizeField.id = sf["ID"].as<int>();
          else {
            std::cerr << "Error: Keyword 'ID' for Size Field not found."
                      << std::endl;
            throw;
          }

          // Get the size field Params
          if (sf.contains("Params")) {
            std::string key;
            double val;
            std::pair<std::string, double> p;

            for (const auto &prm : sf["Params"].object_range()) {
              key = std::string(prm.key());

              // Check what type of data each object has
              if (prm.value().is_array()) {
                if (prm.value()[0].is_string()) {
                  std::pair<std::string, std::vector<std::string>> p_strg;
                  p_strg.first = key;
                  for (int i = 0; i < prm.value().size(); i++) {
                    p_strg.second.push_back(prm.value()[i].as_string());
                    // std::cout << prm.value()[i].as_string() << std::endl;
                  }
                  sizeField.strg_list_params.push_back(p_strg);
                } else {
                  std::pair<std::string, std::vector<double>> p_num;
                  p_num.first = key;
                  for (int i = 0; i < prm.value().size(); i++) {
                    p_num.second.push_back(prm.value()[i].as_double());
                    // std::cout << prm.value()[i].as<int>() << std::endl;
                  }
                  sizeField.num_list_params.push_back(p_num);
                }
              } else {
                val = prm.value().as_double();
                p = {key, val};
                sizeField.params.push_back(p);
              }
            }
            (params->sizeFields).push_back(sizeField);
          } else {
            std::cerr << "Error: Size Field of Type " << sizeField.type
                      << " and ID " << sizeField.id
                      << " has no 'Params' keyword" << std::endl;
            throw;
          }
        }
      }
      // Background Field specification
      if (gparams.contains("BackgroundField")) {
        params->bgField = gparams["BackgroundField"].as<int>();
      } else if (params->mSizeField && !gparams.contains("BackgroundField") &&
                 gparams[cap].size() > 1) {
        std::cout << "Warning: Mesh Background Field not specified."
                  << " Using size field with highest ID." << std::endl;
      }
      std::cout << " Done." << std::endl;

      auto* gmshGenerator =
          new MeshGenDriver(ifname, meshEngine, params, ofname);
      return gmshGenerator;
    }
  } else if (meshEngine == "simmetrix") {
#ifndef HAVE_SIMMETRIX
    std::cerr << "NEMoSys must be recompiled with Simmetrix support."
              << std::endl;
    exit(1);
#else
    if (!inputjson.contains("License File")) {
      std::cerr << "Simmetrix license file must be specified in JSON"
                << std::endl;
      exit(1);
    }

    auto *params = new simmetrixParams();
    params->licFName = inputjson["License File"].as<std::string>();
    params->features = inputjson["Features"].as<std::string>();
    params->logFName = inputjson["Log File"].as<std::string>();

    std::string defaults =
        inputjson["Meshing Parameters"]["Simmetrix Parameters"]
            .as<std::string>();
    if (defaults == "default") {
      auto *mshgndrvobj = new MeshGenDriver(ifname, meshEngine, params, ofname);
      return mshgndrvobj;
    } else {
      jsoncons::json simmxParams =
          inputjson["Meshing Parameters"]["Simmetrix Parameters"];
      if (simmxParams.contains("Mesh Size"))
        params->meshSize = simmxParams["Mesh Size"].as<double>();
      if (simmxParams.contains("Anisotropic Curvature Refinement"))
        params->anisoMeshCurv =
            simmxParams["Anisotropic Curvature Refinement"].as<double>();
      if (simmxParams.contains("Global Gradation Rate"))
        params->glbSizeGradRate =
            simmxParams["Global Gradation Rate"].as<double>();
      if (simmxParams.contains("Surface Mesh Improver Gradation Rate"))
        params->surfMshImprovGradRate =
            simmxParams["Surface Mesh Improver Gradation Rate"].as<double>();
      if (simmxParams.contains("Surface Mesh Improver Min Size"))
        params->surfMshImprovMinSize =
            simmxParams["Surface Mesh Improver Min Size"].as<double>();

      auto *mshgndrvobj = new MeshGenDriver(ifname, meshEngine, params, ofname);
      return mshgndrvobj;
    }
#endif
  } else if (meshEngine == "cfmesh") {
#ifndef HAVE_CFMSH
    std::cerr << "NEMoSys must be recompiled with cfMesh support." << std::endl;
    exit(1);
#else
    auto *params = new cfmeshParams();
    std::string defaults =
        inputjson["Meshing Parameters"]["CFMesh Parameters"].as<std::string>();
    if (defaults == "default") {
      auto *mshgndrvobj = new MeshGenDriver(ifname, meshEngine, params, ofname);
      return mshgndrvobj;
    } else {
      jsoncons::json cfmparams =
          inputjson["Meshing Parameters"]["CFMesh Parameters"];

      // required params here
      // cad file
      if (inputjson["Mesh File Options"].contains("Input Geometry File"))
        params->geomFilePath =
            inputjson["Mesh File Options"]["Input Geometry File"]
                .as<std::string>();
      else {
        std::cerr << "A geometry file should be supplied.\n";
        throw;
      }

      // mesh generator
      if (cfmparams.contains("Generator"))
        params->generator = cfmparams["Generator"].as<std::string>();
      else {
        std::cerr << "A mesh generation method should be selected.\n";
        std::cerr << "Options: cartesian2D tetMesh\n";
        throw;
      }

      // rest of params are optional
      if (cfmparams.contains("MaxCellSize"))
        params->maxCellSize = cfmparams["MaxCellSize"].as<double>();
      if (cfmparams.contains("MinCellSize"))
        params->minCellSize = cfmparams["MinCellSize"].as<double>();
      if (cfmparams.contains("BoundaryCellSize"))
        params->bndryCellSize = cfmparams["BoundaryCellSize"].as<double>();
      if (cfmparams.contains("BoundaryCellSizeRefinementThickness"))
        params->bndryCellSizeRefThk =
            cfmparams["BoundaryCellSizeRefinementThickness"].as<double>();
      if (cfmparams.contains("KeepCellsIntersectingBoundary"))
        params->keepCellIB =
            cfmparams["KeepCellsIntersectingBoundary"].as<bool>();
      if (cfmparams.contains("CheckForGluedMesh"))
        params->chkGluMsh = cfmparams["CheckForGluedMesh"].as<bool>();
      if (cfmparams.contains("AllowDisconnectedDomains"))
        params->_alwDiscDomains =
            cfmparams["AllowDisconnectedDomains"].as<bool>();
      else
        params->_alwDiscDomains = false;

      // optional capability
      std::string cap = "BoundaryLayers";
      if (cfmparams.contains(cap)) {
        params->_withBndLyr = true;
        params->blNLyr = cfmparams[cap]["NLayers"].as<int>();
        params->blThkRto = cfmparams[cap]["ThicknessRatio"].as<double>();
        if (cfmparams[cap].contains("MaxFirstLayerThickness"))
          params->maxFrstLyrThk =
              cfmparams[cap]["MaxFirstLayerThickness"].as<double>();
        if (cfmparams[cap].contains("AllowDiscontinuity"))
          params->alwDiscont = cfmparams[cap]["AllowDiscontinuity"].as<bool>();

        // patch boundary layers
        std::string subcap = "PatchBoundaryLayers";
        if (cfmparams[cap].contains(subcap)) {
          params->_withBndLyrPtch = true;
          for (const auto &jptch : cfmparams[cap][subcap].array_range()) {
            cfmPtchBndLyr blPatch;
            blPatch.patchName = jptch["PatchName"].as<std::string>();
            if (jptch.contains("AllowDiscontinuity"))
              blPatch.alwDiscont = jptch["AllowDiscontinuity"].as<bool>();
            else
              blPatch.alwDiscont = false;
            if (jptch.contains("MaxFirstLayerThickness"))
              blPatch.maxFrstLyrThk = jptch["MaxFirstLayerThickness"].as<int>();
            else
              blPatch.maxFrstLyrThk = -1;
            if (jptch.contains("NLayers"))
              blPatch.blNLyr = jptch["NLayers"].as<int>();
            else
              blPatch.blNLyr = -1;
            if (jptch.contains("ThicknessRatio"))
              blPatch.blThkRto = jptch["ThicknessRatio"].as<double>();
            else
              blPatch.blThkRto = -1.;
            (params->blPatches).push_back(blPatch);
          }
        }
      }

      // optional capability
      cap = "SurfaceFeatureEdges";
      if (cfmparams.contains(cap)) {
        params->_withSrfEdg = true;
        params->srfEdgAng = cfmparams[cap]["Angle"].as<double>();
      }

      // optional capability
      cap = "ObjectRefinements";
      if (cfmparams.contains(cap)) {
        params->_withObjRfn = true;
        for (const auto &refObj : cfmparams[cap].array_range()) {
          cfmObjRef objRef;
          objRef.name = refObj["Name"].as<std::string>();
          for (const auto &prm : refObj["Params"].object_range()) {
            std::string key = std::string(prm.key());
            std::string val = prm.value().as<std::string>();
            objRef.params[key] = val;
          }
          (params->objRefLst).push_back(objRef);
        }
      }

      // optional capability
      cap = "ImproveMeshQuality";
      if (cfmparams.contains(cap)) {
        params->_withMshQlt = true;
        params->qltNItr = cfmparams[cap]["NIterations"].as<int>();
        params->qltNLop = cfmparams[cap]["NLoops"].as<int>();
        params->qltQltThr = cfmparams[cap]["QualityThreshold"].as<double>();
        params->qltNSrfItr = cfmparams[cap]["NSurfaceIterations"].as<int>();
        params->qltConCelSet =
            cfmparams[cap].get_with_default("ConstrainedCellsSet", "none");
      }

      // optional capability
      cap = "LocalRefinement";
      if (cfmparams.contains(cap)) {
        params->_withLclRef = true;
        for (const auto &jptch : cfmparams[cap].array_range()) {
          cfmLclRefPatch refPatch;
          refPatch.patchName = jptch["PatchName"].as<std::string>();
          if (jptch.contains("AdditionalRefinementLevels"))
            refPatch.aditRefLvls =
                jptch["AdditionalRefinementLevels"].as<int>();
          else
            refPatch.aditRefLvls = -1;
          if (jptch.contains("RefinementThickness"))
            refPatch.refThickness = jptch["RefinementThickness"].as<double>();
          else
            refPatch.refThickness = -1;
          if (jptch.contains("CellSize"))
            refPatch.cellSize = jptch["CellSize"].as<double>();
          else
            refPatch.cellSize = -1.;
          (params->refPatches).push_back(refPatch);
        }
      }

      // optional capability
      cap = "RenameBoundary";
      if (cfmparams.contains(cap)) {
        params->_withRenBndry = true;
        cfmRenBndry renBndry;
        renBndry.defName = cfmparams[cap]["DefaultName"].as<std::string>();
        renBndry.defType = cfmparams[cap]["DefaultType"].as<std::string>();
        for (const auto &jnw : cfmparams[cap]["NewPatchNames"].array_range()) {
          cfmNewPatch nwPatch = std::make_tuple(
              jnw["Name"].as<std::string>(), jnw["NewName"].as<std::string>(),
              jnw["NewType"].as<std::string>());
          renBndry.newPatches.push_back(nwPatch);
        }
        (params->renBndry) = renBndry;
      }

      auto *mshgndrvobj = new MeshGenDriver(ifname, meshEngine, params, ofname);
      return mshgndrvobj;
    }
#endif
  } else if (!meshEngine.compare("snappyHexMesh")) {
#ifndef HAVE_CFMSH
    std::cerr << "Nemosys must be recompiled with cfMesh support" << std::endl;
    exit(1);
#else
    snappymeshParams *params = new snappymeshParams();
    std::string defaults =
        inputjson["Meshing Parameters"]["snappyHexMesh Parameters"]
            .as<std::string>();
    if (!defaults.compare("default")) {
      MeshGenDriver *mshgndrvobj =
          new MeshGenDriver(ifname, meshEngine, params, ofname);
      return mshgndrvobj;
    } else {
      jsoncons::json shmparams =
          inputjson["Meshing Parameters"]["snappyHexMesh Parameters"];

      jsoncons::json geomParams = shmparams["Geometry Definition"];
      jsoncons::json castMeshParams = shmparams["Castellated Mesh Controls"];
      jsoncons::json snapParams = shmparams["Snapping Controls"];
      jsoncons::json layerParams = shmparams["Mesh Layers Controls"];
      jsoncons::json qcMeshParams = shmparams["Mesh Quality Controls"];

      if (inputjson["Mesh File Options"].contains("Input Geometry File"))
        params->geomFileName =
            inputjson["Mesh File Options"]["Input Geometry File"]
                .as<std::string>();
      else {
        std::cerr << "A geometry file should be supplied.\n";
        throw;
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
        std::cerr << "Location of a point in region you want to"
                  << "keep cells is needed!" << std::endl;
        throw;
      }
      if (castMeshParams.contains("locationInMeshY"))
        params->locMeshY = castMeshParams["locationInMeshY"].as<double>();
      else {
        std::cerr << "Location of a point in region you want to"
                  << "keep cells is needed!" << std::endl;
        throw;
      }
      if (castMeshParams.contains("locationInMeshZ"))
        params->locMeshZ = castMeshParams["locationInMeshZ"].as<double>();
      else {
        std::cerr << "Location of a point in region you want to"
                  << "keep cells is needed!" << std::endl;
        throw;
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
        // params->_withSurfRefReg = true;

        for (auto jptch3 : castMeshParams[cap3].array_range()) {
          shmSurfRefine surfRef;

          if (jptch3.contains("Patch Name"))
            surfRef.refPatchNm = jptch3["Patch Name"].as<std::string>();
          else {
            surfRef.refPatchNm = params->singleSolidPatch;
          }

          if (jptch3.contains("Patch Type"))
            surfRef.patchType = jptch3["Patch Type"].as<std::string>();
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
      if (castMeshParams.contains(cap2)) {
        // params->_withGeomRefReg = true;

        for (auto jptch2 : castMeshParams[cap2].array_range()) {
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
            ftrOne.fileName = jptch3["File Name"].as<std::string>();
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

      // Layer controls
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
            lyrOne.patchName = jptch3["Patch Name"].as<std::string>();
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
            lyrOne.expansionRatio = jptch3["expansionRatio"].as<double>();
          else {
            lyrOne.expansionRatio = 1;
          }

          if (jptch3.contains("finalLayerThickness"))
            lyrOne.finalLayerThickness =
                jptch3["finalLayerThickness"].as<double>();
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
            lyrOne.minThickness = jptch3["minThickness"].as<double>();
          else {
            lyrOne.minThickness = 1;
          }
          params->layerVec.push_back(lyrOne);
        }
      }

      // Mesh Quality Controls
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

      auto *mshgndrvobj = new MeshGenDriver(ifname, meshEngine, params, ofname);
      return mshgndrvobj;
    }
#endif
  } else if (!meshEngine.compare("blockMesh")) {
#ifndef HAVE_CFMSH
    std::cerr << "Nemosys must be recompiled with cfMesh support" << std::endl;
    exit(1);
#else
    blockMeshParams *params = new blockMeshParams();
    std::string defaults =
        inputjson["Meshing Parameters"]["blockMesh Parameters"]
            .as<std::string>();
    if (!defaults.compare("default")) {
      MeshGenDriver *mshgndrvobj =
          new MeshGenDriver(ifname, meshEngine, params, ofname);
      return mshgndrvobj;
    } else {
      jsoncons::json bmshparams =
          inputjson["Meshing Parameters"]["blockMesh Parameters"];

      if (inputjson["Mesh File Options"].contains("Input Dict File"))
        params->_ownBlockMshDict =
            inputjson["Mesh File Options"]["Input Dict File"].as<bool>();

      // Parameter parsing starts here
      if (bmshparams.contains("Block Geometry"))
        params->_isBlock = bmshparams["Block Geometry"].as<bool>();
      else {
        std::cerr << "Define your choice of geometry using bool keys"
                  << "\n"
                  << std::endl;
        throw;
      }
      if (bmshparams.contains("Sphere Geometry"))
        params->_isSphere = bmshparams["Sphere Geometry"].as<bool>();
      else {
        std::cerr << "Define your choice of geometry using bool keys"
                  << "\n"
                  << std::endl;
        throw;
      }
      if (bmshparams.contains("Cylinder/Tapered_Cone Geometry"))
        params->_isCylinder_TCone =
            bmshparams["Cylinder/Tapered_Cone Geometry"].as<bool>();
      else {
        std::cerr << "Define your choice of geometry using bool keys"
                  << "\n"
                  << std::endl;
        throw;
      }
      if (bmshparams.contains("scaleToMeters"))
        params->cnvrtToMeters = bmshparams["scaleToMeters"].as<double>();
      else {
        std::cerr << "Define your choice of geometry using bool keys"
                  << "\n"
                  << std::endl;
        throw;
      }

      if (bmshparams.contains("Cell_Size")) {
        params->_cellSizeDefined = true;
        params->cellSize = bmshparams["Cell_Size"].as<double>();
      } else {
        params->cellSize = -1;
        params->_cellSizeDefined = false;
      }

      if (bmshparams.contains("XdirectionCells"))
        params->cellsXDir = bmshparams["XdirectionCells"].as<int>();
      else {
        if (params->_cellSizeDefined) {
          // Nothing
        } else {
          std::cerr << "Define cell numbers in X direction"
                    << "\n"
                    << std::endl;
          throw;
        }
      }
      if (bmshparams.contains("YdirectionCells"))
        params->cellsYDir = bmshparams["YdirectionCells"].as<int>();
      else {
        if (params->_cellSizeDefined) {
          // Nothing
        } else {
          std::cerr << "Define cell numbers in Y direction"
                    << "\n"
                    << std::endl;
          throw;
        }
      }
      if (bmshparams.contains("ZdirectionCells"))
        params->cellsZDir = bmshparams["ZdirectionCells"].as<int>();
      else {
        if (params->_cellSizeDefined) {
          // Nothing
        } else {
          std::cerr << "Define cell size for mesh"
                    << "\n"
                    << std::endl;
          throw;
        }
      }

      if (bmshparams.contains("Block Parameters")) {
        if (bmshparams["Block Parameters"].contains("Auto_Generate")) {
          params->_autoGenerateBox = true;

          if (inputjson["Mesh File Options"].contains("Input Geometry File"))
            params->packFileName =
                inputjson["Mesh File Options"]["Input Geometry File"]
                    .as<std::string>();
          else {
            std::cerr << "A geometry file should be supplied.\n";
            throw;
          }

          if (bmshparams["Block Parameters"]["Auto_Generate"].contains(
                  "Offset_XDir"))
            params->offsetX =
                bmshparams["Block Parameters"]["Auto_Generate"]["Offset_XDir"]
                    .as<double>();
          else {
            params->offsetX = 0.1;
          }
          if (bmshparams["Block Parameters"]["Auto_Generate"].contains(
                  "Offset_YDir"))
            params->offsetY =
                bmshparams["Block Parameters"]["Auto_Generate"]["Offset_YDir"]
                    .as<double>();
          else {
            params->offsetY = 0.1;
          }
          if (bmshparams["Block Parameters"]["Auto_Generate"].contains(
                  "Offset_ZDir"))
            params->offsetZ =
                bmshparams["Block Parameters"]["Auto_Generate"]["Offset_ZDir"]
                    .as<double>();
          else {
            params->offsetZ = 0.1;
          }

        } else {
          params->_autoGenerateBox = false;
        }

        if (bmshparams["Block Parameters"].contains("X1"))
          params->initX = bmshparams["Block Parameters"]["X1"].as<double>();
        else {
          if (!(params->_autoGenerateBox)) {
            std::cerr << "Define initial point for block (X Y Z)!\n"
                      << std::endl;
            throw;
          } else {
            std::cout << "Box will be generated automatically" << std::endl;
          }
        }
        if (bmshparams["Block Parameters"].contains("Y1"))
          params->initY = bmshparams["Block Parameters"]["Y1"].as<double>();
        else {
          if (!(params->_autoGenerateBox)) {
            std::cerr << "Define initial point for block (X Y Z)!\n"
                      << std::endl;
            throw;
          } else {
            std::cout << "Box will be generated automatically" << std::endl;
          }
        }
        if (bmshparams["Block Parameters"].contains("Z1"))
          params->initZ = bmshparams["Block Parameters"]["Z1"].as<double>();
        else {
          if (!(params->_autoGenerateBox)) {
            std::cerr << "Define initial point for block (X Y Z)!\n"
                      << std::endl;
            throw;
          } else {
            std::cout << "Box will be generated automatically" << std::endl;
          }
        }
        if (bmshparams["Block Parameters"].contains("LengthX"))
          params->lenX = bmshparams["Block Parameters"]["LengthX"].as<double>();
        else {
          if (!(params->_autoGenerateBox)) {
            std::cerr << "Define desired box length in X,Y,Z direction\n"
                      << std::endl;
            throw;
          } else {
            std::cout << "Box will be generated automatically" << std::endl;
          }
        }
        if (bmshparams["Block Parameters"].contains("LengthY"))
          params->lenY = bmshparams["Block Parameters"]["LengthY"].as<double>();
        else {
          if (!(params->_autoGenerateBox)) {
            std::cerr << "Define desired box length in X,Y,Z direction\n"
                      << std::endl;
            throw;
          } else {
            std::cout << "Box will be generated automatically" << std::endl;
          }
        }
        if (bmshparams["Block Parameters"].contains("LengthZ"))
          params->lenZ = bmshparams["Block Parameters"]["LengthZ"].as<double>();
        else {
          if (!(params->_autoGenerateBox)) {
            std::cerr << "Define desired box length in X,Y,Z direction\n"
                      << std::endl;
            throw;
          } else {
            std::cout << "Box will be generated automatically" << std::endl;
          }
        }
        if (bmshparams["Block Parameters"].contains("GradingXdir"))
          params->smplGradingX =
              bmshparams["Block Parameters"]["GradingXdir"].as<double>();
        else {
          params->smplGradingX = 1;
        }
        if (bmshparams["Block Parameters"].contains("GradingYdir"))
          params->smplGradingY =
              bmshparams["Block Parameters"]["GradingYdir"].as<double>();
        else {
          params->smplGradingY = 1;
        }
        if (bmshparams["Block Parameters"].contains("GradingZdir"))
          params->smplGradingZ =
              bmshparams["Block Parameters"]["GradingZdir"].as<double>();
        else {
          params->smplGradingZ = 1;
        }
      }

      if (bmshparams.contains("Sphere Parameters")) {
        if (bmshparams["Sphere Parameters"].contains("Center X"))
          params->centerX =
              bmshparams["Sphere Parameters"]["Center X"].as<double>();
        else {
          std::cerr << "Define sphere center!\n" << std::endl;
          throw;
        }
        if (bmshparams["Sphere Parameters"].contains("Center Y"))
          params->centerY =
              bmshparams["Sphere Parameters"]["Center Y"].as<double>();
        else {
          std::cerr << "Define sphere center!\n" << std::endl;
          throw;
        }
        if (bmshparams["Sphere Parameters"].contains("Center Z"))
          params->centerZ =
              bmshparams["Sphere Parameters"]["Center Z"].as<double>();
        else {
          std::cerr << "Define sphere center!\n" << std::endl;
          throw;
        }
        if (bmshparams["Sphere Parameters"].contains("Radius"))
          params->radius =
              bmshparams["Sphere Parameters"]["Radius"].as<double>();
        else {
          std::cerr << "Define sphere radius!\n" << endl;
          throw;
        }
        if (bmshparams["Sphere Parameters"].contains("GradingXdir"))
          params->sphrGradingX =
              bmshparams["Sphere Parameters"]["GradingXdir"].as<double>();
        else {
          params->sphrGradingX = 1;
        }
        if (bmshparams["Sphere Parameters"].contains("GradingYdir"))
          params->sphrGradingY =
              bmshparams["Sphere Parameters"]["GradingYdir"].as<double>();
        else {
          params->sphrGradingY = 1;
        }
        if (bmshparams["Sphere Parameters"].contains("GradingXdir"))
          params->sphrGradingZ =
              bmshparams["Sphere Parameters"]["GradingZdir"].as<double>();
        else {
          params->sphrGradingZ = 1;
        }
      }

      if (bmshparams.contains("Cylinder/Tapered_Cone Parameters")) {
        if (bmshparams["Cylinder/Tapered_Cone Parameters"].contains("Center X"))
          params->centerCyl[0] =
              bmshparams["Cylinder/Tapered_Cone Parameters"]["Center X"]
                  .as<double>();
        else {
          std::cerr << "Define center point for cylinder axis\n" << std::endl;
          throw;
        }
        if (bmshparams["Cylinder/Tapered_Cone Parameters"].contains("Center Y"))
          params->centerCyl[1] =
              bmshparams["Cylinder/Tapered_Cone Parameters"]["Center Y"]
                  .as<double>();
        else {
          std::cerr << "Define center point for cylinder axis\n" << std::endl;
          throw;
        }
        if (bmshparams["Cylinder/Tapered_Cone Parameters"].contains("Center Z"))
          params->centerCyl[2] =
              bmshparams["Cylinder/Tapered_Cone Parameters"]["Center Z"]
                  .as<double>();
        else {
          std::cerr << "Define center point for cylinder axis\n" << std::endl;
          throw;
        }
        if (bmshparams["Cylinder/Tapered_Cone Parameters"].contains("Radius1"))
          params->radius1 =
              bmshparams["Cylinder/Tapered_Cone Parameters"]["Radius1"]
                  .as<double>();

        if (bmshparams["Cylinder/Tapered_Cone Parameters"].contains(
                "Radius2")) {
          params->radius2 =
              bmshparams["Cylinder/Tapered_Cone Parameters"]["Radius2"]
                  .as<double>();
        } else {
          params->radius2 =
              bmshparams["Cylinder/Tapered_Cone Parameters"]["Radius1"]
                  .as<double>();
        }

        if (bmshparams["Cylinder/Tapered_Cone Parameters"].contains(
                "GradingXdir"))
          params->cylGrading[0] =
              bmshparams["Cylinder/Tapered_Cone Parameters"]["GradingXdir"]
                  .as<double>();
        else {
          params->cylGrading[0] = 1;
        }
        if (bmshparams["Cylinder/Tapered_Cone Parameters"].contains(
                "GradingYdir"))
          params->cylGrading[1] =
              bmshparams["Cylinder/Tapered_Cone Parameters"]["GradingYdir"]
                  .as<double>();
        else {
          params->cylGrading[1] = 1;
        }
        if (bmshparams["Cylinder/Tapered_Cone Parameters"].contains(
                "GradingXdir"))
          params->cylGrading[2] =
              bmshparams["Cylinder/Tapered_Cone Parameters"]["GradingZdir"]
                  .as<double>();
        else {
          params->cylGrading[2] = 1;
        }

        if (bmshparams["Cylinder/Tapered_Cone Parameters"].contains("Height"))
          params->height =
              bmshparams["Cylinder/Tapered_Cone Parameters"]["Height"]
                  .as<double>();
        else {
          std::cerr << "Define cylinder height\n" << endl;
        }
      }

      auto *mshgndrvobj = new MeshGenDriver(ifname, meshEngine, params, ofname);
      return mshgndrvobj;
    }

#endif
  }

  else {
    std::cerr << "Mesh generation engine " << meshEngine << " is not supported"
              << std::endl;
    exit(1);
  }
}

} // namespace DRV
} // namespace NEM