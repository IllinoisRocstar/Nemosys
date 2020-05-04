#include "meshGen.H"

#include "gmshGen.H"
#include "gmshParams.H"
#include "meshingParams.H"
#ifdef HAVE_NGEN
#  include "netgenGen.H"
#  include "netgenParams.H"
#endif
#ifdef HAVE_SIMMETRIX
#  include "simmetrixGen.H"
#  include "simmetrixParams.H"
#endif
#ifdef HAVE_CFMSH
#  include "blockMeshGen.H"
#  include "blockMeshParams.H"
#  include "cfmeshGen.H"
#  include "cfmeshParams.H"
#  include "snappymeshGen.H"
#  include "snappymeshParams.H"
#endif

meshGen *meshGen::Create(const std::string &fname,
                         const std::string &meshEngine) {
  if (meshEngine == "netgen") {
#ifdef HAVE_NGEN
    auto *generator = new netgenGen();
    return generator;
#else
    std::cerr << "NETGEN is not enabled during build."
              << " Build NEMoSys with ENABLE_NETGEN to use this method."
              << std::endl;
    exit(1);
#endif
  } else if (meshEngine == "gmsh") {
    auto *generator = new NEM::GEN::gmshGen();
    return generator;
  }
#ifdef HAVE_SIMMETRIX
  else if (meshEngine == "simmetrix") {
    auto *generator = new simmetrixGen();
    return generator;
  }
#endif
#ifdef HAVE_CFMSH
  else if (meshEngine == "cfmesh") {
    auto *generator = new cfmeshGen();
    return generator;
  } else if (meshEngine == "snappyHexMesh") {
    auto *generator = new snappymeshGen();
    return generator;
  } else if (meshEngine == "blockMesh") {
    auto *generator = new blockMeshGen();
    return generator;
  }
#endif
  else {
    std::cerr << meshEngine << " is not a supported meshing engine"
              << std::endl;
    exit(1);
  }
}

meshGen *meshGen::Create(const std::string &fname,
                         const std::string &meshEngine, meshingParams *params) {
  if (meshEngine == "netgen") {
#ifdef HAVE_NGEN
    auto *generator = new netgenGen(dynamic_cast<netgenParams *>(params));
    return generator;
#else
    std::cerr << "NETGEN is not enabled during build."
              << " Build NEMoSys with ENABLE_NETGEN to use this method."
              << std::endl;
    exit(1);
#endif
  } else if (meshEngine == "gmsh") {
    auto *generator =
        new NEM::GEN::gmshGen(dynamic_cast<NEM::GEN::gmshParams *>(params));
    return generator;
  }
#ifdef HAVE_SIMMETRIX
  else if (meshEngine == "simmetrix") {
    auto *generator = new simmetrixGen(dynamic_cast<simmetrixParams *>(params));
    return generator;
  }
#endif
#ifdef HAVE_CFMSH
  else if (meshEngine == "cfmesh") {
    auto *generator = new cfmeshGen(dynamic_cast<cfmeshParams *>(params));
    return generator;
  } else if (meshEngine == "snappyHexMesh") {
    auto *generator =
        new snappymeshGen(dynamic_cast<snappymeshParams *>(params));
    return generator;
  } else if (meshEngine == "blockMesh") {
    auto *generator = new blockMeshGen(dynamic_cast<blockMeshParams *>(params));
    return generator;
  }
#endif
  else {
    std::cerr << meshEngine << " is not a supported meshing engine"
              << std::endl;
    exit(1);
  }
}
