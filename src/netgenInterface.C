#include <netgenInterface.H>
#include <meshing.hpp>
#include <string.h>

//TODO:      Add writing surface elements to exportToVTK

using namespace nglib;   
int netgenInterface::createMeshFromSTL(char* fname)
{

  if (numPoints || numCells)
  {
    std::cout << "mesh is already populated!\ndelete and proceed" << std::endl;
    exit(1);      
  }
  
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

int netgenInterface::importFromVol(char* fname)
{
  if (numPoints || numCells)
  {
    std::cout << "mesh is already populated!\ndelete and proceed" << std::endl;
    exit(1);      
  }
  
  int status = Ng_MergeMesh(mesh, fname);
  numPoints = Ng_GetNP(mesh);
  numCells = Ng_GetNE(mesh);

  if (!(numPoints && numCells))
  {
    std::cout << "Error loading file " << fname 
              << "\nMesh has no points or no cells!\n";
    exit(1);
  }
    
  return status;
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

  int numSurfCells = Ng_GetNSE(mesh);

  std::cout << "\tNumber Of Points: " << numPoints << std::endl;
  std::cout << "\tNumber Of Surface Elements: " << numSurfCells << std::endl;
  std::cout << "\tNumber Of Volume Elements: " << numCells << std::endl;


  vtk << "\nCELLS " << numCells+numSurfCells << " " << numCells*5 + numSurfCells*4 << std::endl;
  for (int i = 1; i <= numSurfCells; ++i)
  {
    int tri[3];
    Ng_GetSurfaceElement(mesh, i, tri);
    vtk << 3 << " ";
    for (int j = 0; j < 3; ++j)
      vtk << tri[j]-1 << " ";
    vtk << std::endl;
  }

  for (int i = 1; i <= numCells; ++i)
  { 
    int tet[4];
    Ng_GetVolumeElement(mesh, i, tet);
    vtk << 4 << " ";
    for (int j = 0; j < 4; ++j)
      vtk << tet[j]-1 << " ";
    vtk << std::endl;
  }

  vtk << "\nCELL_TYPES " << numCells+numSurfCells << std::endl;
  for (int i = 0; i < numSurfCells; ++i)
    vtk << 5 << std::endl; 
  for (int i = 0; i < numCells; ++i)
    vtk << 10 << std::endl;
  return 0;
}

/*int netgenInterface::importFromVTK(char* fname)
{
  if (numPoints || numCells)
  {
    std::cout << "mesh is already populated!\ndelete and proceed" << std::endl;
    exit(1);      
  }

  vtkAnalyzer* vtkMesh = new vtkAnalyzer(fname);
  vtkMesh->report();
  vtkDataSet* dataSet = vtkMesh->getDataSet();
  
  // add points
  for (int i = 0; i < vtkMesh->numberOfPoints; ++i)
  {
    double* pntCrd =  vtkMesh->getPointCoords(i);
    Ng_AddPoint(mesh, pntCrd);
  }

  // add surface elements (if they exist) and volume elements
  bool has_surf = 0;
  for (int i = 0; i < vtkMesh->numberOfCells; ++i)
  {
    vtkIdList* vtkpntIds = dataSet->GetCell(i)->GetPointIds();
    int numIds = vtkpntIds->GetNumberOfIds();
    int pntIds[numIds];
    for (int j = 0; j < numIds; ++j)
      pntIds[j] = vtkpntIds->GetId((numIds-1)-j)+1; // netgen is 1-based index
    if (numIds == 3)
    {
      Ng_AddSurfaceElement(mesh,NG_TRIG,pntIds);
      has_surf = 1;
    }
    if (numIds == 4)
      Ng_AddVolumeElement(mesh,NG_TET,pntIds);  
  }

  if (!has_surf)
  {
    std::cout << "\tNo surface elements detected in " << fname << std::endl
              << "\tDetecting and adding boundary elements to mesh ...." << std::endl;
    std::multimap<int, std::vector<int> > boundaries = vtkMesh->findBoundaryFaces();  
    std::multimap<int,std::vector<int> >::iterator it;
    for (it = boundaries.begin(); it!=boundaries.end(); ++it)
    {
      int pntIds[3];
      for (int i = 0; i < 3; ++i)
        pntIds[i] = it->second[2-i]+1; // netgen is 1-based index
      Ng_AddSurfaceElement(mesh,NG_TRIG,pntIds); 
    }
    std::cout << "\tAdded " << boundaries.size() << " boundary elements" << std::endl;
  }
  


 
  // Set the Meshing Parameters to be used
  Ng_Meshing_Parameters mp; 
  mp.maxh = 1;
  mp.fineness = 1;
  mp.uselocalh = 1; 
  mp.optsteps_2d = 1;
  mp.optsteps_3d = 1;
  mp.grading = .1;
  Ng_OptimizeVolume(mesh, &mp);
   
  //Ng_Uniform_Refinement (mesh);

  // testing mesh import
  std::string old_fname(fname);
  std::string new_fname = trim_fname(old_fname,".vol");
  Ng_SaveMesh(mesh, &new_fname[0u]);
  numPoints = Ng_GetNP(mesh);
  numCells = Ng_GetNE(mesh);

  delete vtkMesh;  
  return 0;
}*/

  /* Default constructor for the Mesh Parameters class

   Note: This constructor initialises the variables in the 
   class with the following default values
   - #uselocalh: 1
   - #maxh: 1000.0
   - #fineness: 0.5
   - #grading: 0.3
   - #elementsperedge: 2.0
   - #elementspercurve: 2.0
   - #closeedgeenable: 0
   - #closeedgefact: 2.0
   - #secondorder: 0
   - #meshsize_filename: null
   - #quad_dominated: 0
   - #optsurfmeshenable: 1
   - #optvolmeshenable: 1
   - #optsteps_2d: 3
   - #optsteps_3d: 3
   - #invert_tets: 0
   - #invert_trigs:0 
   - #check_overlap: 1
   - #check_overlapping_boundary: 1 */
