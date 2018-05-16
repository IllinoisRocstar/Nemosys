#include <meshGen.H>
#include <netgenGen.H>
#ifdef HAVE_SYMMX
  #include <symmxGen.H>
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
    netgenGen* generator = new netgenGen(dynamic_cast<NetgenParams*>(params));
    return generator;
  }
  #ifdef HAVE_SYMMX
  else if (!meshEngine.compare("simmetrix"))
  {
    symmxGen* generator = new symmxGen(dynamic_cast<SymmxParams*>(params));
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
