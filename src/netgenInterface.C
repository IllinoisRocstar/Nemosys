#include <netgenInterface.H>

using namespace nglib;   
int netgenInterface::createMeshFromSTL(char* fname)
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
  Ng_Meshing_Parameters mp; 
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

  //cout << "Saving Mesh in VOL Format...." << endl;
  //Ng_SaveMesh(mesh,"test.vol");


  // refinement without geomety adaption:
  // Ng_Uniform_Refinement (mesh);

  // refinement with geomety adaption:   
  Ng_STL_Uniform_Refinement (stl_geom, mesh);

  cout << "elements after refinement: " << Ng_GetNE(mesh) << endl;
  cout << "points   after refinement: " << Ng_GetNP(mesh) << endl;

  //Ng_SaveMesh(mesh,"test_ref.vol");
  numPoints = Ng_GetNP(mesh);
  numCells = Ng_GetNE(mesh);
  return 0;
}

int netgenInterface::exportToVTK(char* fname)
{
  if (!(numPoints && numCells))
  {
    std::cout << "Mesh has no points or no cells!" << std::endl;
    exit(1);
  }
    
  std::ofstream vtk(fname);

  if (!vtk.good())
  {
    std::cout << "Error opening: " << fname << std::endl;
    exit(1);
  }

  vtk << "# vtk DataFile Version 2.0" << std::endl 
      << "Converted From Netgen" << std::endl
      << "ASCII" << std::endl
      << "DATASET UNSTRUCTURED_GRID" << std::endl 
      << "POINTS " << numPoints << " double" << std::endl;

  // netgen is 1-based index
  for (int i = 1; i <= numPoints; ++i)
  {
    double point[3];
    Ng_GetPoint(mesh, i, point);
    for (int j = 0; j < 3; ++j)
      vtk << point[j] << " ";
    vtk << std::endl;         
  }

  vtk << "\nCELLS " << numCells << " " << numCells*5 << std::endl;
  for (int i = 1; i <= numCells; ++i)
  { 
    int tet[4];
    Ng_GetVolumeElement(mesh, i, tet);
    
    vtk << 4 << " ";
    for (int j = 0; j < 4; ++j)
      vtk << tet[j]-1 << " ";
    vtk << std::endl;
  }

  vtk << "\nCELL_TYPES " << numCells << std::endl;
  for (int i = 0; i < numCells; ++i)
    vtk << VTK_TETRA << std::endl; 
 
}



