#include <MeshQualityDriver.H>

MeshQualityDriver::MeshQualityDriver(std::string _mesh, std::string ofname)
{
  mesh = meshBase::Create(_mesh);
  mesh->checkMesh(ofname);
  std::cout << "MeshQualityDriver created" << std::endl;
}

MeshQualityDriver::~MeshQualityDriver()
{
  if (mesh)
  {
    delete mesh;
    mesh = 0;
  }
  std::cout << "MeshQualityDriver destroyed" << std::endl;
}

MeshQualityDriver* MeshQualityDriver::readJSON(json inputjson)
{

  MeshQualityDriver* qualdrvobj;
  std::string _mesh;
  std::string ofname;
  _mesh = inputjson["Input Mesh File"].as<std::string>();
  ofname = inputjson["Output File"].as<std::string>();
  
  qualdrvobj = new MeshQualityDriver(_mesh, ofname); 
  return qualdrvobj;  
}
