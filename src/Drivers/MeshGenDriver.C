#include <MeshGenDriver.H>
#include <netgenGen.H>
#include <netgenParams.H>
#ifdef HAVE_SYMMX
  #include <symmxGen.H>
  #include <symmxParams.H>
#endif
// ----------------------------- MeshGen Driver -----------------------------------//

MeshGenDriver::MeshGenDriver(std::string ifname, std::string meshEngine, 
                             meshingParams* _params, std::string ofname)
{
  params = _params;
  mesh = meshBase::generateMesh(ifname, meshEngine, params);
  mesh->setFileName(ofname);
  mesh->report();
  mesh->write();
  std::cout << "MeshGenDriver created" << std::endl;
}

MeshGenDriver::~MeshGenDriver()
{
  if (mesh)
  {
    delete mesh;
    mesh = 0;
  }
  if (params)
  {
    delete params;
    params = 0;
  } 
  std::cout << "MeshGenDriver destroyed" << std::endl;
}

MeshGenDriver* MeshGenDriver::readJSON(json inputjson)
{
  std::string meshEngine = inputjson["Mesh Generation Engine"].as<std::string>();
  std::string ifname = inputjson["Mesh File Options"]
                                ["Input Geometry File"].as<std::string>();
  std::string ofname = inputjson["Mesh File Options"]
                                ["Output Mesh File"].as<std::string>();
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

  else
  {
    std::cout << "Mesh generation engine " << meshEngine << " is not supported" << std::endl;
    exit(1);
  }
}
