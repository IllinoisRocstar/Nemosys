#include <NemDrivers.H>


TransferDriver::TransferDriver(std::string srcmsh, std::string trgmsh,
                               std::string method, std::string ofname)
{
  source = new meshUser(srcmsh);
  target = new meshUser(trgmsh);
  Timer T;
  T.start();
  source->transfer(target, method);
  T.stop();
  std::cout << "Time spent transferring data (ms)" << T.elapsed() << std::endl;
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
  std::cout << "Time spent transferring data (ms)" << T.elapsed() << std::endl;
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

  RefineDriver* refdrvobj;
  json inputjson;
  inputStream >> inputjson;
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



