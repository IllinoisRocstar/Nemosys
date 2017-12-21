#include <MeshGenDriver.H>
#include <meshGen.H>
// ----------------------------- MeshGen Driver -----------------------------------//

MeshGenDriver::MeshGenDriver(std::string ifname, std::string meshEngine, 
                             meshingParams* _params, std::string ofname)
{
  params = _params;
  mesh = new meshBase();
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
      NetgenParams* params = new NetgenParams();
      MeshGenDriver* mshgndrvobj = new MeshGenDriver(ifname, meshEngine, params, ofname);
      return mshgndrvobj;
    }
    else
    {
      json ngparams = inputjson["Meshing Parameters"]["Netgen Parameters"];
      
      NetgenParams* params = new NetgenParams();  
     
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

      MeshGenDriver* mshgndrvobj = new MeshGenDriver(ifname, meshEngine, params, ofname);
      return mshgndrvobj;
    }


  }
  else
  {
    std::cout << "Mesh generation engine " << meshEngine << " is not supported" << std::endl;
    exit(1);
  }
}
