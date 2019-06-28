#include "RefineDriver.H"
#include "AuxiliaryFunctions.H"

#include <iostream>
#include <fstream>

// -------------------------------- Refine Driver -------------------------------------//
RefineDriver::RefineDriver(const std::string &_mesh,
                           const std::string &method,
                           const std::string &arrayName,
                           double dev_mult,
                           bool maxIsmin,
                           double edgescale,
                           const std::string &ofname,
                           bool transferData,
                           double sizeFactor)
{
  std::cout << "RefineDriver created" << std::endl;
  std::cout << "Size Factor = " << sizeFactor << std::endl;
  mesh = meshBase::Create(_mesh);
  std::cout << "\n";
  mesh->report();
  std::cout << "\n";
  mesh->refineMesh(method, arrayName, dev_mult, maxIsmin, edgescale, ofname,
                   transferData, sizeFactor);
}

RefineDriver::RefineDriver(const std::string &_mesh,
                           const std::string &method,
                           double edgescale,
                           const std::string &ofname,
                           bool transferData)
{
  std::cout << "RefineDriver created" << std::endl;
  mesh = meshBase::Create(_mesh);
  std::cout << "\n";
  mesh->report();
  std::cout << "\n";
  mesh->refineMesh(method, edgescale, ofname, transferData);
}

RefineDriver::RefineDriver(const std::string &_mesh,
                           const std::string &method,
                           const std::string &arrayName,
                           int order,
                           const std::string &ofname,
                           bool transferData)
{
  mesh = meshBase::Create(_mesh);
  std::cout << "\n";
  mesh->report();
  std::cout << "\n";
  mesh->refineMesh(method, arrayName, order, ofname, transferData);
  std::cout << "RefineDriver created" << std::endl;
}

RefineDriver::~RefineDriver()
{
  delete mesh;
  std::cout << "RefineDriver destroyed" << std::endl;
}

RefineDriver *RefineDriver::readJSON(const jsoncons::json &inputjson)
{
  RefineDriver *refdrvobj;
  std::string _mesh;
  std::string ofname;
  std::string method;
  std::string arrayName;
  double dev_mult;
  bool maxIsmin, transferData;
  double edgescale = 0;
  _mesh = inputjson["Mesh File Options"]["Input Mesh File"].as<std::string>();
  ofname = inputjson["Mesh File Options"]["Output Mesh File"].as<std::string>();
  method = inputjson["Refinement Options"]["Refinement Method"].as<std::string>();
  transferData = inputjson["Refinement Options"]["Transfer Data"].as<bool>();
  if (method == "uniform")
  {
    edgescale = inputjson["Refinement Options"]["Edge Scaling"].as<double>();
    refdrvobj = new RefineDriver(_mesh, method, edgescale, ofname,
                                 transferData);
  }
  else if (method == "Z2 Error Estimator")
  {
    arrayName = inputjson["Refinement Options"]["Array Name"].as<std::string>();
    int order = inputjson["Refinement Options"]["Shape Function Order"].as<int>();
    refdrvobj = new RefineDriver(_mesh, method, arrayName, order, ofname,
                                 transferData);
  }
  else
  {
    arrayName = inputjson["Refinement Options"]["Array Name"].as<std::string>();
    dev_mult = inputjson["Refinement Options"]["StdDev Multiplier"].as<double>();
    maxIsmin = inputjson["Refinement Options"]["Max Is Min for Scaling"].as<bool>();
    double sizeFactor;
    sizeFactor = inputjson["Refinement Options"].has_key("Size Factor")
                 ? inputjson["Refinement Options"]["Size Factor"].as<double>()
                 : 1.0;
    refdrvobj = new RefineDriver(_mesh, method, arrayName, dev_mult, maxIsmin,
                                 edgescale, ofname, transferData, sizeFactor);
  }
  return refdrvobj;
}

RefineDriver *RefineDriver::readJSON(const std::string &ifname)
{
  std::ifstream inputStream(ifname);
  if (!inputStream.good() || nemAux::find_ext(ifname) != ".json")
  {
    std::cerr << "Error opening file " << ifname << std::endl;
    exit(1);
  }
  if (nemAux::find_ext(ifname) != ".json")
  {
    std::cerr << "Input File must be in .json format" << std::endl;
    exit(1);
  }

  jsoncons::json inputjson;
  inputStream >> inputjson;

  // checking if array
  if (inputjson.is_array())
  {
    std::cerr
        << "Warning: Input is an array. Only first element will be processed\n";
    return RefineDriver::readJSON(inputjson[0]);
  }
  else
  {
    return RefineDriver::readJSON(inputjson);
  }
}
