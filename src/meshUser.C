#include <meshUser.H>


/*int meshUser::generateMesh(std::string filename, std::string meshEngine)
{
  if (filename.find(".stl") == -1)
  {
    std::cout << "Only CAD files in STL format are supported" << std::endl;
    exit(1);
  }
  mesh = meshBase::generateMesh(filename,meshEngine);
  write_ext.assign(".vtu");
  fname.assign(trim_fname(filename,".vol"));
  std::cout << "user constructed" << std::endl;
  
  return 0;

}*/
int meshUser::generateMesh(std::string filename, std::string meshEngine, meshingParams* params)
{
  if (filename.find(".stl") == -1)
  {
    std::cout << "Only CAD files in STL format are supported" << std::endl;
    exit(1);
  }
  mesh = meshBase::generateMesh(filename,meshEngine, params);
  write_ext.assign(".vtu");
  fname.assign(trim_fname(filename,".vol"));
  std::cout << "user constructed" << std::endl;
  
  return 0;
}

// get number of points in mesh
int meshUser::getNumberOfPoints() 
{ 
  return mesh->getNumberOfPoints(); 
}

// get number of cells in mesh
int meshUser::getNumberOfCells() 
{ 
  return mesh->getNumberOfCells(); 
}

void meshUser::report()
{
  mesh->report(&fname[0u]);
}

// get point with id
std::vector<double> meshUser::getPoint(int id) 
{ 
  return mesh->getPoint(id);
}

// get cell with id : returns point indices and respective coordinates
std::map<int,std::vector<double>> meshUser::getCell(int id) 
{ 
  return mesh->getCell(id);
}

// print coordinates of point with id
void meshUser::printPoint(int id)
{ 
  //std::vector<double> coord = mesh->getPoint(id);
  printVec(mesh->getPoint(id));
}

// print point ids and coordinates of cell with id
void meshUser::printCell(int id)
{
  std::map<int, std::vector<double>> cell = mesh->getCell(id);        
  std::map<int, std::vector<double>>::iterator it = cell.begin();
  while(it != cell.end())
  {
    std::cout << it->first << " ";
    printVec(it->second);
    it++;
  }
}

// write mesh to file
void meshUser::write()
{
  mesh->write(&(trim_fname(fname,write_ext))[0u], write_ext);
}

void meshUser::write(std::string fname)
{
  mesh->write(fname,write_ext);
}

// transfer point data with given id from this user to target
int meshUser::transfer(meshUser* target, std::string method, int arrayID)
{
  return mesh->transfer(target->getMesh(),method,arrayID);
}

// transfer point data with given ids from this user to target
int meshUser::transfer(meshUser* target, std::string method, 
                       const std::vector<int>& arrayIDs)
{
  return mesh->transfer(target->getMesh(),method,arrayIDs);
}

// transfer all point data from this user to target
int meshUser::transfer(meshUser* target, std::string method)
{
  return mesh->transfer(target->getMesh(), method);
}

// TODO: Probably can delete this function
// transfer point data with given name from this user to target
int meshUser::transfer(meshUser* target, std::string method, std::string arrayName)
{
  int arrayID = mesh->IsArrayName(arrayName);
  if (arrayID == -1)
  {
    std::cout << "Array " << arrayName 
              << " not found in set of point data arrays" << std::endl;
    exit(1);
  }
  
  return mesh->transfer(target->getMesh(),method);
}

// transfer point data with given names fro mthis user to target
int meshUser::transfer(meshUser* target, std::string method, 
             const std::vector<string>& arrayNames)
{
  std::vector<int> arrayIDs(arrayNames.size());
  for (int i = 0; i < arrayNames.size(); ++i)
  {
    int id = mesh->IsArrayName(arrayNames[i]);
    if (id == -1)
    {
      std::cout << "Array " << arrayNames[i] 
                << " not found in set of point data arrays" << std::endl;
      exit(1);
    }
    arrayIDs[i] = id;
  }
  return mesh->transfer(target->getMesh(),method,arrayIDs);
}

// size field generation
void meshUser::generateSizeField(std::string method, int arrayID, double dev_mult, bool maxIsmin)
{
  mesh->generateSizeField(method, arrayID, dev_mult, maxIsmin);
}

// convert to gmsh format without data (only volume elements)
void meshUser::writeMSH(std::string fname)
{
  mesh->writeMSH(fname);
}

// convert to gmsh format with specified point or cell data
void meshUser::writeMSH(std::string fname, std::string pointOrCell, int arrayID)
{
  mesh->writeMSH(fname, pointOrCell, arrayID);
}

void meshUser::refineMesh(std::string method, int arrayID, 
                          double dev_mult, bool maxIsmin, 
                          double edge_scale, std::string ofname)
{
  mesh->refineMesh(method, arrayID, dev_mult, maxIsmin, edge_scale, ofname);
}

void meshUser::refineMesh(std::string method, std::string arrayName, 
                          double dev_mult, bool maxIsmin, 
                          double edge_scale, std::string ofname)
{
  int arrayID = mesh->IsArrayName(arrayName);
  if (arrayID == -1)
  {
    std::cout << "Array " << arrayName
              << " not fuond in set of point data arrays" << std::endl;
    exit(1);
  }
  mesh->refineMesh(method, arrayID, dev_mult, maxIsmin, edge_scale, ofname);
}
void meshUser::refineMesh(std::string method, double edge_scale, std::string ofname)
{
  mesh->refineMesh(method, 0, 0, 0, edge_scale, ofname);
}
