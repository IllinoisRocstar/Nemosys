#include <meshGen.H>
#include <meshingParams.H>
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

meshGen* meshGen::Create(std::string fname, std::string meshEngine)
{
  if (!meshEngine.compare("netgen"))
  {
    netgenGen* generator = new netgenGen();
    return generator;
  }
  #ifdef HAVE_SYMMX 
  else if (!meshEngine.compare("simmetrix"))
  {
    symmxGen* generator = new symmxGen();
    return generator;  
  }
  #endif
  #ifdef HAVE_CFMSH 
  else if (!meshEngine.compare("cfmesh"))
  {
    cfmeshGen* generator = new cfmeshGen();
    return generator;
  }
  #endif
  else
  {
    std::cout << meshEngine << " is not a supported meshing engine" << std::endl;
    exit(1);
  }
}

meshGen* meshGen::Create(std::string fname, std::string meshEngine, meshingParams* params)
{
  if (!meshEngine.compare("netgen"))
  {
    netgenGen* generator = new netgenGen(dynamic_cast<netgenParams*>(params));
    return generator;
  }
  #ifdef HAVE_SYMMX
  else if (!meshEngine.compare("simmetrix"))
  {
    symmxGen* generator = new symmxGen(dynamic_cast<symmxParams*>(params));
    return generator;
  }
  #endif
  #ifdef HAVE_CFMSH
  else if (!meshEngine.compare("cfmesh"))
  {
    cfmeshGen* generator = new cfmeshGen(dynamic_cast<cfmeshParams*>(params));
    return generator;
  }
  #endif
  else
  {
    std::cout << meshEngine << " is not a supported meshing engine" << std::endl;
    exit(1);
  }
}

vtkSmartPointer<vtkDataSet> meshGen::getDataSet()
{
  return dataSet;
} 
