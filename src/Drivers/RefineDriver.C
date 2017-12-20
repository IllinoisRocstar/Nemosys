#include <RefineDriver.H>

// -------------------------------- Refine Driver -------------------------------------//
RefineDriver::RefineDriver(std::string _mesh, std::string method, std::string arrayName,
                           double dev_mult, bool maxIsmin, double edgescale, std::string ofname)
{
  mesh = meshBase::Create(_mesh);
  std::cout << std::endl;
  mesh->report();
  std::cout << std::endl;
  mesh->refineMesh(method, arrayName, dev_mult, maxIsmin, edgescale, ofname);  
  
  meshBase* refinedmesh = meshBase::Create(ofname);
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
  mesh = meshBase::Create(_mesh);
  std::cout << std::endl;
  mesh->report();
  std::cout << std::endl;
  mesh->refineMesh(method, edgescale, ofname);
  
  meshBase* refinedmesh = meshBase::Create(ofname);
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

