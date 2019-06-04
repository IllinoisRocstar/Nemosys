#include <StlToVtk.H>

#include "AuxiliaryFunctions.H"

void exportStlToVtk(std::string ifname)
{
  std::ifstream stl(ifname.c_str());
  if (!stl.good())
  {
    std::cout << "error reading " << ifname << std::endl;
    exit(1);
  }

  std::map<std::vector<double>,int> point_map;
  typedef std::pair<std::map<std::vector<double>,int>::iterator, bool> insert_return;
  insert_return ret;

  std::vector<triangle> triangles;

  std::string line;
  double tmp;
  int i = 1;
  std::vector<double> p1,p2,p3;
  int j,k,l;
  j = -1;
  while(getline(stl,line))
  {
    std::size_t found = line.find("vertex");
    if (found != std::string::npos)
    { 
      std::string data = line.substr(found+6,std::string::npos);
      std::stringstream ss(data);
      while(ss >> tmp)
      {
        if (i == 1)
          p1.push_back(tmp);
        if (i == 2)
          p2.push_back(tmp);
        if (i == 3)
          p3.push_back(tmp);
      }  
      i+=1;
    }
    if (i == 4)
    {
      triangles.push_back(triangle(p1,p2,p3));

      ret = point_map.insert(std::pair<std::vector<double>,int> (p1,j));
      if (ret.second == true)
      {
        j+=1;
        point_map[p1] = j;
      }
      ret = point_map.insert(std::pair<std::vector<double>,int> (p2,j));
      if (ret.second == true)
      {
        j+=1;
        point_map[p2] = j;

      }
      ret = point_map.insert(std::pair<std::vector<double>,int> (p3,j));
      if (ret.second == true)
      {
        j+=1;
        point_map[p3] = j;
      }
      
      p1.clear();
      p2.clear();
      p3.clear();
      i = 1;
    }
  }

  std::multimap<int, std::vector<double>> flipped = nemAux::flip_map(point_map);

  std::ofstream vtk(nemAux::trim_fname(ifname, ".vtk"));
  if (!vtk.good()) {
    std::cout << "Error opening file " << nemAux::trim_fname(ifname, ".vtk")
              << std::endl;
    exit(1);
  }

  vtk << "# vtk DataFile Version 2.0" << std::endl 
      << "Converted From Netgen" << std::endl
      << "ASCII" << std::endl
      << "DATASET UNSTRUCTURED_GRID" << std::endl 
      << "POINTS " << point_map.size() << " double" << std::endl;
  
  std::multimap<int,std::vector<double>>::iterator flip_it = flipped.begin();
  while (flip_it != flipped.end())
  {
    for (int i = 0; i < 3; ++i)
      vtk << flip_it->second[i] << " ";
    vtk << std::endl; 
    ++flip_it;
  }
 
  vtk << "\nCELLS " << triangles.size() << " " << triangles.size()*4 << std::endl;
  for (int i = 0; i < triangles.size(); ++i)
  {
    vtk << 3 << " "
        << point_map[triangles[i].p1] << " "
        << point_map[triangles[i].p2] << " "
        << point_map[triangles[i].p3] << std::endl;
  }

  vtk << "\nCELL_TYPES " << triangles.size() << std::endl;
  for (int i = 0; i < triangles.size(); ++i)
    vtk << 5 << std::endl; 
}





