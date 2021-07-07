#ifdef HAVE_NGEN

#include "MeshGeneration/netgenGen.H"
#include "MeshGeneration/netgenParams.H"

#include "AuxiliaryFunctions.H"


netgenGen::netgenGen()
    : mp(), refine_with_geom(), refine_without_geom()
{
  std::cout << "initializing netgen mesh generator" << std::endl;
  nglib::Ng_Init();
  mesh = nglib::Ng_NewMesh();
}

netgenGen::netgenGen(const netgenParams *params)
    : refine_with_geom(), refine_without_geom()
{
  std::cout << "initializing netgen mesh generator" << std::endl;
  nglib::Ng_Init();
  set_mp(params);
  mesh = nglib::Ng_NewMesh();
}

netgenGen::~netgenGen()
{
  std::cout << "finalizing netgen mesh generator" << std::endl;
  if (mesh)
  {
    nglib::Ng_DeleteMesh(mesh);
    mesh = nullptr;
  }
  nglib::Ng_Exit();
}


void netgenGen::set_mp(const netgenParams* params)
{
  mp.uselocalh                   = params->uselocalh;
  mp.maxh                        = params->maxh;
  mp.fineness                    = params->fineness;
  mp.grading                     = params->grading;
  mp.elementsperedge             = params->elementsperedge;
  mp.elementspercurve            = params->elementspercurve;
  mp.closeedgeenable             = params->closeedgeenable;
  mp.closeedgefact               = params->closeedgefact;
  mp.second_order                = params->second_order;
  mp.meshsize_filename = const_cast<char *>(params->meshsize_filename.c_str());
  mp.quad_dominated              = params->quad_dominated;
  mp.optvolmeshenable            = params->optvolmeshenable;
  mp.optsteps_2d                 = params->optsteps_2d;
  mp.optsteps_3d                 = params->optsteps_3d;
  mp.invert_tets                 = params->invert_tets;
  mp.invert_trigs                = params->invert_trigs;
  mp.check_overlap               = params->check_overlap;
  mp.check_overlapping_boundary  = params->check_overlapping_boundary;
  refine_without_geom            = params->refine_without_geom;
  refine_with_geom               = params->refine_with_geom;
}

int netgenGen::createMeshFromSTL(const char *fname)
{
  // Define pointer to STL Geometry
  nglib::Ng_STL_Geometry *stl_geom;

  // Result of Netgen Operations
  nglib::Ng_Result ng_res;

  int np, ne;

  // Check that the STL file exists
  std::ifstream inSTL(fname);
  if (!inSTL.good()) {
    std::cerr << "Error reading in STL file: " << fname << std::endl;
    return 1;
  }

  // Read in the STL File
  stl_geom = nglib::Ng_STL_LoadGeometry(fname);

  // This object seems to always be defined. If the file does
  // not exits, the process continues with no triangle facets
  // loaded.
  if (!stl_geom)
  {
    std::cerr << "Error reading in STL File: " << fname << std::endl;
    return 1;
  }
  std::cout << "Successfully loaded STL File: " << fname << std::endl;

  std::cout << "Initialise the STL Geometry structure...." << std::endl;
  ng_res = nglib::Ng_STL_InitSTLGeometry(stl_geom);
  if (ng_res != nglib::NG_OK)
  {
    std::cerr << "Error Initialising the STL Geometry....Aborting!!" << std::endl;
    return 1;
  }

  std::cout << "Start Edge Meshing...." << std::endl;
  ng_res = nglib::Ng_STL_MakeEdges(stl_geom, mesh, &mp);
  if (ng_res != nglib::NG_OK)
  {
    std::cerr << "Error in Edge Meshing....Aborting!!" << std::endl;
    return 1;
  }

  std::cout << "Start Surface Meshing...." << std::endl;
  ng_res = nglib::Ng_STL_GenerateSurfaceMesh(stl_geom, mesh, &mp);
  if (ng_res != nglib::NG_OK)
  {
    std::cerr << "Error in Surface Meshing....Aborting!!" << std::endl;
    return 1;
  }

  std::cout << "Start Volume Meshing...." << std::endl;
  ng_res = nglib::Ng_GenerateVolumeMesh(mesh, &mp);
  if (ng_res != nglib::NG_OK)
  {
    std::cerr << "Error in Volume Meshing....Aborting!!" << std::endl;
    return 1;
  }

  std::cout << "Meshing successfully completed....!!" << std::endl;

  // volume mesh output
  np = nglib::Ng_GetNP(mesh);
  std::cout << "Points: " << np << std::endl;

  ne = nglib::Ng_GetNE(mesh);
  std::cout << "Elements: " << ne << std::endl;

  // refinement without or with geometry adaption:
  if (refine_without_geom)
    nglib::Ng_Uniform_Refinement(mesh);
  else if (refine_with_geom)
    nglib::Ng_STL_Uniform_Refinement(stl_geom, mesh);

  std::cout << "elements after refinement: " << nglib::Ng_GetNE(mesh)
            << std::endl;
  std::cout << "points   after refinement: " << nglib::Ng_GetNP(mesh)
            << std::endl;

  std::cout << "Saving Mesh in VOL Format...." << std::endl;
  std::string newfname(fname);
  nglib::Ng_SaveMesh(mesh, nemAux::trim_fname(newfname, ".vol").c_str());
  return 0;
}

#endif
