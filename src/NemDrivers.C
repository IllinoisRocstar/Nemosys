#include <NemDrivers.H>
#include <meshGen.H>

//------------------------------ Factory of Drivers ----------------------------------------//
NemDriver* NemDriver::readJSON(std::string ifname)
{
  std::ifstream inputStream(ifname);
  if (!inputStream.good() || find_ext(ifname) != ".json")
  {
    std::cout << "Error opening file " << ifname << std::endl;
    exit(1);
  }
  if (find_ext(ifname) != ".json")
  {
    std::cout << "Input File must be in .json format" << std::endl;
    exit(1);
  }

  json inputjson;
  inputStream >> inputjson;
  std::string program_type = inputjson["Program Type"].as<std::string>();
  if (!program_type.compare("Transfer"))
  {
    return TransferDriver::readJSON(inputjson); 
  }
  else if (!program_type.compare("Refinement"))
  {
    return RefineDriver::readJSON(inputjson);
  }
  else if (!program_type.compare("Mesh Generation"))
  {
    return MeshGenDriver::readJSON(inputjson);
  }
  else
  {
    std::cout << "Program Type " << program_type 
              << " is not supported by  Nemosys" << std::endl;
    exit(1);
  }
  
}

//----------------------- Transfer Driver -----------------------------------------//
TransferDriver::TransferDriver(std::string srcmsh, std::string trgmsh,
                               std::string method, std::string ofname)
{
  source = new meshUser(srcmsh);
  target = new meshUser(trgmsh);
  Timer T;
  T.start();
  source->transfer(target, method);
  T.stop();
  std::cout << "Time spent transferring data (ms) " << T.elapsed() << std::endl;
  target->write(ofname); 
}


TransferDriver::TransferDriver(std::string srcmsh, std::string trgmsh, std::string method,
                               std::vector<std::string> arrayNames, std::string ofname)
{
  source = new meshUser(srcmsh);
  target = new meshUser(trgmsh);
  Timer T;
  T.start();
  source->transfer(target, method, arrayNames);
  T.stop();
  std::cout << "Time spent transferring data (ms) " << T.elapsed() << std::endl;
  target->write(ofname); 
}

TransferDriver::~TransferDriver()
{
  if (source)
  {
    delete source;
    source = 0;
  }
  if (target)
  {
    delete target;
    target = 0;
  }
}

TransferDriver* TransferDriver::readJSON(json inputjson)
{
  std::string srcmsh; 
  std::string trgmsh; 
  std::string outmsh; 
  std::string method;
  std::string transferAll;
  bool transferall = 1;
  std::vector<std::string> arrayNames;

  srcmsh = inputjson["Mesh File Options"]
                    ["Input Mesh Files"]
                    ["Source Mesh"].as<std::string>();
  trgmsh = inputjson["Mesh File Options"]
                    ["Input Mesh Files"]
                    ["Target Mesh"].as<std::string>();
  outmsh = inputjson["Mesh File Options"]
                    ["Output Mesh File"].as<std::string>();
  method = inputjson["Transfer Options"]
                    ["Method"].as<std::string>(); 

  transferAll = inputjson["Transfer Options"]
                         ["Transfer All Arrays"].as<std::string>();
  
  if (!transferAll.compare("False") || !transferAll.compare("false"))
  {
    arrayNames = inputjson["Transfer Options"]
                          ["Array Names"].as<std::vector<std::string>>();
    transferall = 0;
  }
 
  TransferDriver* trnsdrvobj;
  if (transferall)
  {
    trnsdrvobj = new TransferDriver(srcmsh, trgmsh, method, outmsh);
  } 
  else
  {
    trnsdrvobj = new TransferDriver(srcmsh, trgmsh, method, arrayNames, outmsh); 
  }
  
  return trnsdrvobj;

}

TransferDriver* TransferDriver::readJSON(std::string ifname)
{
  std::ifstream inputStream(ifname);
  if (!inputStream.good() || find_ext(ifname) != ".json")
  {
    std::cout << "Error opening file " << ifname << std::endl;
    exit(1);
  }
  if (find_ext(ifname) != ".json")
  {
    std::cout << "Input File must be in .json format" << std::endl;
    exit(1);
  }

  json inputjson;
  inputStream >> inputjson;
  return TransferDriver::readJSON(inputjson);
    
}

// -------------------------------- Refine Driver -------------------------------------//
RefineDriver::RefineDriver(std::string _mesh, std::string method, std::string arrayName,
                           double dev_mult, bool maxIsmin, double edgescale, std::string ofname)
{
  mesh = new meshUser(_mesh);
  std::cout << std::endl;
  mesh->report();
  std::cout << std::endl;
  mesh->refineMesh(method, arrayName, dev_mult, maxIsmin, edgescale, ofname);  
  
  meshUser* refinedmesh = new meshUser(ofname);
  std::cout << std::endl;
  refinedmesh->report();
  std::cout << std::endl;
  if (refinedmesh)
  {
    delete refinedmesh;
    refinedmesh = 0;
  }
  
}

RefineDriver::RefineDriver(std::string _mesh, std::string method, double edgescale, std::string ofname)
{
  mesh = new meshUser(_mesh);
  std::cout << std::endl;
  mesh->report();
  std::cout << std::endl;
  mesh->refineMesh(method, edgescale, ofname);
  
  meshUser* refinedmesh = new meshUser(ofname);
  std::cout << std::endl;
  refinedmesh->report();
  std::cout << std::endl;
  if (refinedmesh)
  {
    delete refinedmesh;
    refinedmesh = 0;
  }
}

RefineDriver::~RefineDriver()
{
  if (mesh)
  {
    delete mesh;
    mesh = 0;
  }
}

RefineDriver* RefineDriver::readJSON(json inputjson)
{
  RefineDriver* refdrvobj;
  std::string _mesh;
  std::string ofname;
  std::string method;
  std::string arrayName;
  double dev_mult;
  bool maxIsmin;
  double edgescale = 0;
  _mesh = inputjson["Mesh File Options"]["Input Mesh File"].as<std::string>();
  ofname = inputjson["Mesh File Options"]["Output Mesh File"].as<std::string>();
  method = inputjson["Refinement Options"]["Refinement Method"].as<std::string>();
  if (!method.compare("uniform"))
  {
    edgescale = inputjson["Refinement Options"]["Edge Scaling"].as<double>();
    refdrvobj = new RefineDriver(_mesh, method, edgescale, ofname);   
  }
  else
  {
    arrayName = inputjson["Refinement Options"]["Array Name"].as<std::string>();
    dev_mult = inputjson["Refinement Options"]["StdDev Multiplier"].as<double>();
    maxIsmin = inputjson["Refinement Options"]["Max Is Min for Scaling"].as<bool>();
    refdrvobj = new RefineDriver(_mesh, method, arrayName, dev_mult, maxIsmin, 
                                 edgescale, ofname);
  }
  return refdrvobj;  
}

RefineDriver* RefineDriver::readJSON(std::string ifname)
{
  std::ifstream inputStream(ifname);
  if (!inputStream.good() || find_ext(ifname) != ".json")
  {
    std::cout << "Error opening file " << ifname << std::endl;
    exit(1);
  }
  if (find_ext(ifname) != ".json")
  {
    std::cout << "Input File must be in .json format" << std::endl;
    exit(1);
  }


  json inputjson;
  inputStream >> inputjson;

  return RefineDriver::readJSON(inputjson);
}

// ----------------------------- MeshGen Driver -----------------------------------//
MeshGenDriver::MeshGenDriver(std::string ifname, std::string meshEngine, 
                             std::string ofname)
{
  mesh = new meshUser();
  mesh->generateMesh(ifname, meshEngine);
  mesh->setFileName(ofname);
  mesh->report();
  mesh->write();
}

MeshGenDriver::MeshGenDriver(std::string ifname, std::string meshEngine, 
                             meshingParams* _params, std::string ofname)
{
  params = _params;
  mesh = new meshUser();
  mesh->generateMesh(ifname, meshEngine, params);
  mesh->setFileName(ofname);
  mesh->report();
  mesh->write();
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
      MeshGenDriver* mshgndrvobj = new MeshGenDriver(ifname, meshEngine, ofname);
      return mshgndrvobj;
    }
    else
    {
      
      NetgenParams* params = new NetgenParams();  
      
      params->uselocalh = inputjson["Meshing Parameters"]["Netgen Parameters"]
                                   ["uselocalh"].as<bool>();
      params->maxh = inputjson["Meshing Parameters"]["Netgen Parameters"]
                              ["maxh"].as<double>();
      params->fineness = inputjson["Meshing Parameters"]["Netgen Parameters"]
                                  ["fineness"].as<double>();
      params->grading = inputjson["Meshing Parameters"]["Netgen Parameters"]  
                                  ["grading"].as<double>();
      params->elementsperedge = inputjson["Meshing Parameters"]["Netgen Parameters"]
                                         ["elementsperedge"].as<double>();
      params->elementspercurve = inputjson["Meshing Parameters"]["Netgen Parameters"]
                                          ["elementspercurve"].as<double>();
      params->closeedgeenable = inputjson["Meshing Parameters"]["Netgen Parameters"]
                                         ["closeedgeenable"].as<bool>();
      params->closeedgefact  = inputjson["Meshing Parameters"]["Netgen Parameters"]
                                        ["closeedgefact"].as<double>();
      params->second_order  = inputjson["Meshing Parameters"]["Netgen Parameters"]
                                       ["second_order"].as<bool>();
      params->meshsize_filename = inputjson["Meshing Parameters"]["Netgen Parameters"]
                                           ["meshsize_filename"].as<std::string>();
      params->quad_dominated  = inputjson["Meshing Parameters"]["Netgen Parameters"]
                                         ["quad_dominated"].as<bool>();
      params->optvolmeshenable = inputjson["Meshing Parameters"]["Netgen Parameters"]
                                          ["optvolmeshenable"].as<bool>();
      params->optsteps_2d = inputjson["Meshing Parameters"]["Netgen Parameters"]
                                     ["optsteps_2d"].as<int>();
      params->optsteps_3d = inputjson["Meshing Parameters"]["Netgen Parameters"]
                                     ["optsteps_3d"].as<int>();
      params->invert_tets = inputjson["Meshing Parameters"]["Netgen Parameters"]
                                     ["invert_tets"].as<bool>();
      params->invert_trigs = inputjson["Meshing Parameters"]["Netgen Parameters"]
                                      ["invert_trigs"].as<bool>();
      params->check_overlap = inputjson["Meshing Parameters"]["Netgen Parameters"]
                                       ["check_overlap"].as<bool>();
      params->check_overlapping_boundary = inputjson["Meshing Parameters"]["Netgen Parameters"]
                                                    ["check_overlapping_boundary"].as<bool>();
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

