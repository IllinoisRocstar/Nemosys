#include <MeshGenDriver.H>
#include <netgenGen.H>
#include <netgenParams.H>
#ifdef HAVE_SYMMX
  #include <symmxGen.H>
  #include <symmxParams.H>
#endif
#ifdef HAVE_CFMSH
  #include <cfmeshGen.H>
  #include <cfmeshParams.H>
#endif

// std c++
#include <vector>
#include <tuple>

// ----------------------------- MeshGen Driver -----------------------------------//

MeshGenDriver::MeshGenDriver(const std::string& ifname, const std::string& meshEngine, 
                             meshingParams* _params, const std::string& ofname)
{
  params = _params;
  mesh = meshBase::CreateShared(meshBase::generateMesh(ifname, meshEngine, params));
  mesh->setFileName(ofname);
  mesh->report();
  mesh->write();
  std::cout << "MeshGenDriver created" << std::endl;
}

std::shared_ptr<meshBase> MeshGenDriver::getNewMesh()
{
  if (mesh)
    return mesh;
}

MeshGenDriver::~MeshGenDriver()
{
  std::cout << "MeshGenDriver destroyed" << std::endl;
}

MeshGenDriver* MeshGenDriver::readJSON(const json& inputjson)
{
  std::string ifname = inputjson["Mesh File Options"]
                                ["Input Geometry File"].as<std::string>();
  std::string ofname = inputjson["Mesh File Options"]
                                ["Output Mesh File"].as<std::string>();
  return readJSON(ifname, ofname, inputjson);
}

MeshGenDriver* MeshGenDriver::readJSON(const std::string& ifname, 
                                       const std::string& ofname,
                                       const json& inputjson)
{
  std::string meshEngine = inputjson["Mesh Generation Engine"].as<std::string>();
  if (!meshEngine.compare("netgen"))
  {
    std::string defaults = inputjson["Meshing Parameters"]
                                    ["Netgen Parameters"].as<std::string>();
    if (!defaults.compare("default"))
    {
      netgenParams* params = new netgenParams();
      MeshGenDriver* mshgndrvobj = new MeshGenDriver(ifname, meshEngine, params, ofname);
      return mshgndrvobj;
    }
    else
    {
      json ngparams = inputjson["Meshing Parameters"]["Netgen Parameters"];
      
      netgenParams* params = new netgenParams();  
     
      if (ngparams.has_key("uselocalh")) 
        params->uselocalh = ngparams["uselocalh"].as<bool>();
      if (ngparams.has_key("maxh"))
        params->maxh = ngparams["maxh"].as<double>();
      if (ngparams.has_key("fineness"))
        params->fineness = ngparams["fineness"].as<double>();
      if (ngparams.has_key("grading"))
        params->grading = ngparams["grading"].as<double>();
      if (ngparams.has_key("elementsperedge"))
        params->elementsperedge = ngparams["elementsperedge"].as<double>();
      if (ngparams.has_key("elementspercurve"))
        params->elementspercurve = ngparams["elementspercurve"].as<double>();
      if (ngparams.has_key("closeedgeenable"))
        params->closeedgeenable = ngparams["closeedgeenable"].as<bool>();
      if (ngparams.has_key("closeedgefact"))
        params->closeedgefact  = ngparams["closeedgefact"].as<double>();
      if (ngparams.has_key("second_order"))
        params->second_order  = ngparams["second_order"].as<bool>();
      if (ngparams.has_key("meshsize_filename"))
        params->meshsize_filename = ngparams["meshsize_filename"].as<std::string>();
      if (ngparams.has_key("quad_dominated"))
        params->quad_dominated  = ngparams["quad_dominated"].as<bool>();
      if (ngparams.has_key("optvolmeshenable"))
        params->optvolmeshenable = ngparams["optvolmeshenable"].as<bool>();
      if (ngparams.has_key("optsteps_2d"))
        params->optsteps_2d = ngparams["optsteps_2d"].as<int>();
      if (ngparams.has_key("optsteps_3d"))
        params->optsteps_3d = ngparams["optsteps_3d"].as<int>();
      if (ngparams.has_key("invert_tets"))
        params->invert_tets = ngparams["invert_tets"].as<bool>();
      if (ngparams.has_key("invert_trigs"))
        params->invert_trigs = ngparams["invert_trigs"].as<bool>();
      if (ngparams.has_key("check_overlap"))
        params->check_overlap = ngparams["check_overlap"].as<bool>();
      if (ngparams.has_key("check_overlapping_boundary"))
        params->check_overlapping_boundary = ngparams["check_overlapping_boundary"].as<bool>();
      if (ngparams.has_key("refine_with_geometry_adaptation"))
        params->refine_with_geom = ngparams["refine_with_geometry_adaptation"].as<bool>();
      if (ngparams.has_key("refine_without_geometry_adaptation"))
        params->refine_without_geom = ngparams["refine_without_geometry_adaptation"].as<bool>(); 

      MeshGenDriver* mshgndrvobj = new MeshGenDriver(ifname, meshEngine, params, ofname);
      return mshgndrvobj;
    }
  }
  else if (!meshEngine.compare("simmetrix"))
  {
    #ifndef HAVE_SYMMX
      std::cerr << "Nemosys must be recompiled with simmetrix support" << std::endl;
      exit(1);
    #else
      if (!inputjson.has_key("License File"))
      {
        std::cerr << "Simmetrix License file must be specified in json" << std::endl;
        exit(1);
      }
      
      symmxParams* params = new symmxParams();
      params->licFName = inputjson["License File"].as<std::string>();
      params->features = inputjson["Features"].as<std::string>();
      params->logFName = inputjson["Log File"].as<std::string>();
      
      std::string defaults = inputjson["Meshing Parameters"]["Simmetrix Parameters"].as<std::string>();
      if (!defaults.compare("default"))
      {
        MeshGenDriver* mshgndrvobj = new MeshGenDriver(ifname, meshEngine, params, ofname);
        return mshgndrvobj; 
      }
      else
      {
        json symmxparams = inputjson["Meshing Parameters"]["Simmetrix Parameters"];
        if (symmxparams.has_key("Mesh Size"))
          params->meshSize = symmxparams["Mesh Size"].as<double>();
        if (symmxparams.has_key("Anisotropic Curvature Refinement"))
          params->anisoMeshCurv = symmxparams["Anisotropic Curvature Refinement"].as<double>();
        if (symmxparams.has_key("Global Gradation Rate"))
          params->glbSizeGradRate = symmxparams["Global Gradation Rate"].as<double>();
        if (symmxparams.has_key("Surface Mesh Improver Gradation Rate"))
          params->surfMshImprovGradRate = symmxparams["Surface Mesh Improver Gradation Rate"].as<double>();
        if (symmxparams.has_key("Surface Mesh Improver Min Size"))
          params->surfMshImprovMinSize = symmxparams["Surface Mesh Improver Min Size"].as<double>();
      
        MeshGenDriver* mshgndrvobj = new MeshGenDriver(ifname, meshEngine, params, ofname);
        return mshgndrvobj;
      } 
    #endif
  }

  else if (!meshEngine.compare("cfmesh"))
  {
    #ifndef HAVE_CFMSH
      std::cerr << "Nemosys must be recompiled with cfMesh support" << std::endl;
      exit(1);
    #else
      cfmeshParams* params = new cfmeshParams();
      std::string defaults = 
          inputjson["Meshing Parameters"]["CFMesh Parameters"].as<std::string>();
      if (!defaults.compare("default"))
      {
        MeshGenDriver* mshgndrvobj = 
            new MeshGenDriver(ifname, meshEngine, params, ofname);
        return mshgndrvobj;
      }
      else
      {
        json cfmparams = inputjson["Meshing Parameters"]["CFMesh Parameters"];

        // required params here
        // cad file
        if (inputjson["Mesh File Options"].has_key("Input Geometry File"))
            params->geomFilePath = 
                inputjson["Mesh File Options"]["Input Geometry File"].as<std::string>();
        else
        {
            std::cerr << "A geometry file should be supplied.\n";
            throw;        
        }
        
        // mesh generator
        if (cfmparams.has_key("Generator"))
            params->generator = cfmparams["Generator"].as<std::string>();
        else
        {
            std::cerr << "A mesh generation method should be selected.\n";
            std::cerr << "Options: cartesian2D tetMesh\n";
            throw;        
        }
        
        // rest of params are optional
        if (cfmparams.has_key("MaxCellSize"))
          params->maxCellSize = 
              cfmparams["MaxCellSize"].as<double>();
        if (cfmparams.has_key("MinCellSize"))
          params->minCellSize = 
              cfmparams["MinCellSize"].as<double>();
        if (cfmparams.has_key("BoundaryCellSize"))
          params->bndryCellSize = 
              cfmparams["BoundaryCellSize"].as<double>();
        if (cfmparams.has_key("KeepCellsIntersectingBoundary"))
          params->keepCellIB = 
              cfmparams["KeepCellsIntersectingBoundary"].as<double>();
        if (cfmparams.has_key("CheckForGluedMesh"))
          params->chkGluMsh = 
              cfmparams["CheckForGluedMesh"].as<double>();

        // optional capability
        std::string cap = "BoundaryLayers";
        if (cfmparams.has_key(cap))
        {
            params->_withBndLyr = true;
            params->blNLyr = cfmparams[cap]["NLayers"].as<double>();
            params->blThkRto = cfmparams[cap]["ThicknessRatio"].as<double>();
            if (cfmparams[cap].has_key("MaxFirstLayerThickness"))
                params->maxFrstLyrThk = cfmparams[cap]["MaxFirstLayerThickness"].as<double>();
            if (cfmparams[cap].has_key("AllowDiscontinuity"))
                params->alwDiscont = cfmparams[cap]["AllowDiscontinuity"].as<bool>();

            // patch boundary layers
            std::string subcap = "PatchBoundaryLayers";
            if (cfmparams[cap].has_key(subcap))
            {
                params->_withBndLyrPtch = true;
                for (auto jptch : cfmparams[cap][subcap].array_range())
                {
                    cfmPtchBndLyr blPatch;
                    blPatch.patchName = jptch["PatchName"].as<std::string>();
                    if (jptch.has_key("AllowDiscontinuity"))
                        blPatch.alwDiscont = jptch["AllowDiscontinuity"].as<bool>();
                    else
                        blPatch.alwDiscont = false;
                    if (jptch.has_key("MaxFirstLayerThickness"))
                        blPatch.maxFrstLyrThk = jptch["MaxFirstLayerThickness"].as<int>();
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
                    (params->blPatches).push_back(blPatch);
                }
                            
            }
        }

        // optional capability
        cap = "SurfaceFeatureEdges";
        if (cfmparams.has_key(cap))
        {
            params->_withSrfEdg = true;
            params->srfEdgAng = cfmparams[cap]["Angle"].as<double>();
        }

        // optional capability
        cap = "ObjectRefinements";
        if (cfmparams.has_key(cap))
        {
            params->_withObjRfn = true;
            for (auto refObj : cfmparams[cap].array_range())
            {
                cfmObjRef objRef;
                objRef.name= refObj["Name"].as<std::string>();
                for (const auto&  prm: refObj["Params"].object_range())
                {
                    std::string key = prm.key();
                    std::string val = prm.value().as<std::string>();
                    objRef.params[key] = val;
                }
                (params->objRefLst).push_back(objRef);
            }
        }

        // optional capability
        cap = "ImproveMeshQuality";
        if (cfmparams.has_key(cap))
        {
            params->_withMshQlt = true;
            params->qltNItr = cfmparams[cap]["NIterations"].as<int>();
            params->qltNLop = cfmparams[cap]["NLoops"].as<int>();
            params->qltQltThr = cfmparams[cap]["QualityThreshold"].as<double>();
            params->qltNSrfItr = cfmparams[cap]["NSurfaceIterations"].as<int>();
            params->qltConCelSet = 
                cfmparams[cap].get_with_default("ConstrainedCellsSet","none");
        }

        // optional capability
        cap = "LocalRefinement";
        if (cfmparams.has_key(cap))
        {
            params->_withLclRef = true;
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
                (params->refPatches).push_back(refPatch);
            }
        }

        // optional capability
        cap = "RenameBoundary";
        if (cfmparams.has_key(cap))
        {
            params->_withRenBndry = true;
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
            (params->renBndry) = renBndry;
        }

        MeshGenDriver* mshgndrvobj = new MeshGenDriver(ifname, meshEngine, 
                params, ofname);
        return mshgndrvobj;
      } 

    #endif
  }
  else
  {
    std::cout << "Mesh generation engine " << meshEngine << " is not supported" << std::endl;
    exit(1);
  }
}
