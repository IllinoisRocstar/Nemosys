#include "meshGen.H"

#include "meshingParams.H"
#ifdef HAVE_NGEN
  #include "netgenGen.H"
  #include "netgenParams.H"
#endif
#ifdef HAVE_SIMMETRIX
  #include "simmetrixGen.H"
  #include "simmetrixParams.H"
#endif
#ifdef HAVE_CFMSH
  #include "cfmeshGen.H"
  #include "cfmeshParams.H"
#endif


meshGen *meshGen::Create(const std::string &fname,
                         const std::string &meshEngine)
{
  if (meshEngine == "netgen")
  {
#ifdef HAVE_NGEN
    auto *generator = new netgenGen();
    return generator;
#else
    std::cerr << "NETGEN is not enabled during build."
              << " Build NEMoSys with ENABLE_NETGEN to use this method."
              << std::endl;
    exit(1);
#endif
  }
#ifdef HAVE_SIMMETRIX
  else if (meshEngine == "simmetrix")
  {
    auto *generator = new simmetrixGen();
    return generator;
  }
#endif
#ifdef HAVE_CFMSH
  else if (meshEngine == "cfmesh")
  {
    auto *generator = new cfmeshGen();
    return generator;
  }
#endif
  else
  {
    std::cerr << meshEngine << " is not a supported meshing engine"
              << std::endl;
    exit(1);
  }
}

meshGen *meshGen::Create(const std::string &fname,
                         const std::string &meshEngine,
                         meshingParams *params)
{
  if (meshEngine == "netgen")
  {
#ifdef HAVE_NGEN
    auto *generator = new netgenGen(dynamic_cast<netgenParams *>(params));
    return generator;
#else
    std::cerr << "NETGEN is not enabled during build."
              << " Build NEMoSys with ENABLE_NETGEN to use this method."
              << std::endl;
    exit(1);
#endif
  }
#ifdef HAVE_SIMMETRIX
  else if (meshEngine == "simmetrix")
  {
    auto *generator = new simmetrixGen(
        dynamic_cast<simmetrixParams *>(params));
    return generator;
  }
#endif
#ifdef HAVE_CFMSH
  else if (meshEngine == "cfmesh")
  {
    auto *generator = new cfmeshGen(dynamic_cast<cfmeshParams *>(params));
    return generator;
  }
#endif
  else
  {
    std::cerr << meshEngine << " is not a supported meshing engine"
              << std::endl;
    exit(1);
  }
}
