#include <meshGen.H>
#include <meshingParams.H>
#include <netgenGen.H>
#include <netgenParams.H>
#ifdef HAVE_SIMMETRIX
  #include <simmetrixGen.H>
  #include <simmetrixParams.H>
#endif

meshGen* meshGen::Create(std::string fname, std::string meshEngine)
{
  if (!meshEngine.compare("netgen"))
  {
    netgenGen* generator = new netgenGen();
    return generator;
  }
  #ifdef HAVE_SIMMETRIX 
  else if (!meshEngine.compare("simmetrix"))
  {
    simmetrixGen* generator = new simmetrixGen();
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
  #ifdef HAVE_SIMMETRIX
  else if (!meshEngine.compare("simmetrix"))
  {
    simmetrixGen* generator = new simmetrixGen(dynamic_cast<simmetrixParams*>(params));
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
