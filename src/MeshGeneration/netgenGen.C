#include <netgenGen.H>
#include <netgenParams.H>
#include <AuxiliaryFunctions.H>

using namespace nglib;   


netgenGen::netgenGen():mp()
{
  std::cout << "initializing netgen mesh generator" << std::endl;
  nglib::Ng_Init();
  mesh = nglib::Ng_NewMesh();
}

netgenGen::netgenGen(netgenParams* params)
{
  std::cout << "initializing netgen mesh generator" << std::endl;
  nglib::Ng_Init();
  set_mp(params);
  mesh = nglib::Ng_NewMesh();
}   

netgenGen::~netgenGen()
{
  std::cout << "finalizing netgen mesh generator" << std::endl;
  if(mesh)
  {
    nglib::Ng_DeleteMesh(mesh);
    mesh = 0; 
  }
  nglib::Ng_Exit();
}

void netgenGen::set_mp(netgenParams* params)
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
  mp.meshsize_filename           = &(params->meshsize_filename)[0u];          
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

using std::cout; using std::endl;

int netgenGen::createMeshFromSTL(const char* fname)
{
  // Define pointer to STL Geometry
  Ng_STL_Geometry *stl_geom;

  // Result of Netgen Operations
  Ng_Result ng_res;

  int np, ne; 

  // Read in the STL File
  stl_geom = Ng_STL_LoadGeometry(fname);
  if(!stl_geom)
  {   
     cout << "Error reading in STL File: " << fname << endl;
   return 1;
  }   
  cout << "Successfully loaded STL File: " << fname << endl;

  cout << "Initialise the STL Geometry structure...." << endl;
  ng_res = Ng_STL_InitSTLGeometry(stl_geom);
  if(ng_res != NG_OK)
  {   
     cout << "Error Initialising the STL Geometry....Aborting!!" << endl;
    return 1;
  }   

  cout << "Start Edge Meshing...." << endl;
  ng_res = Ng_STL_MakeEdges(stl_geom, mesh, &mp);
  if(ng_res != NG_OK)
  {   
     cout << "Error in Edge Meshing....Aborting!!" << endl;
    return 1;
  }   

  cout << "Start Surface Meshing...." << endl;
  ng_res = Ng_STL_GenerateSurfaceMesh(stl_geom, mesh, &mp);
  if(ng_res != NG_OK)
  {   
     cout << "Error in Surface Meshing....Aborting!!" << endl;
    return 1;
  }   
  
  cout << "Start Volume Meshing...." << endl;
  ng_res = Ng_GenerateVolumeMesh (mesh, &mp);
  if(ng_res != NG_OK)
  {   
     cout << "Error in Volume Meshing....Aborting!!" << endl;
   return 1;
  }   
  
  cout << "Meshing successfully completed....!!" << endl;

  // volume mesh output
  np = Ng_GetNP(mesh);
  cout << "Points: " << np << endl;

  ne = Ng_GetNE(mesh);
  cout << "Elements: " << ne << endl;

  // refinement without or with geomety adaption:
  if (refine_without_geom)
    Ng_Uniform_Refinement (mesh);
  else if (refine_with_geom)
    Ng_STL_Uniform_Refinement (stl_geom, mesh);

  cout << "elements after refinement: " << Ng_GetNE(mesh) << endl;
  cout << "points   after refinement: " << Ng_GetNP(mesh) << endl;

  cout << "Saving Mesh in VOL Format...." << endl;
  std::string newfname(fname);
  Ng_SaveMesh(mesh,&(trim_fname(newfname,".vol"))[0u]);
  return 0;
}

