#include <meshGen.H>

using namespace nglib;   
int meshNetgen::createMeshFromSTL(char* fname)
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

  // Set the Meshing Parameters to be used
  mp.maxh = 1.0e+6;
  mp.fineness = 0.4;
  mp.second_order = 0;

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




  // refinement without geomety adaption:
  // Ng_Uniform_Refinement (mesh);

  // refinement with geomety adaption:   
  Ng_STL_Uniform_Refinement (stl_geom, mesh);

  cout << "elements after refinement: " << Ng_GetNE(mesh) << endl;
  cout << "points   after refinement: " << Ng_GetNP(mesh) << endl;

  cout << "Saving Mesh in VOL Format...." << endl;
  std::string newfname(fname);
  Ng_SaveMesh(mesh,&(trim_fname(newfname,".vol"))[0u]);
  return 0;
}


  
