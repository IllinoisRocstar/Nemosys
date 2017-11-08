#ifndef netgenInterface_H
#define netgenInterface_H

// Nemosys
#include <meshPhys.H>

// third party
namespace nglib
{
  #include <nglib.h>
}

/*
  TODO: This class may do the following:
    1) Read in stl geometry
    2) Generate netgen mesh
    3) Write to vtk format
    4) Use refinement params from simulation with meshPhys 
    5) Use netgen to perform h refinement
    6) Write back to vtk
    7) Interpolate data to new mesh and write solution vtk

    8) Also take point cloud of vtk and add points and do volume meshing 
*/

class netgenInterface
{
  // constructor and destructor
  public:
    netgenInterface():numPoints(0),numCells(0)
    {
      nglib::Ng_Init();
      mesh = nglib::Ng_NewMesh();
    }
    
    ~netgenInterface()
    {
      nglib::Ng_DeleteMesh(mesh);
      nglib::Ng_Exit();
    }      

  // netgen mesh creation
  public:
    int createMeshFromSTL(char* fname);
    int exportToVTK(char* fname);
 
  private:
//    nglib::Ng_Meshing_Parameters meshParams; // params for netgen meshing  
    nglib::Ng_Mesh* mesh; // netgen mesh object
    int numPoints, numCells; // number of points and cells in mesh
};

/*
   Default constructor for the Mesh Parameters class

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
   - #check_overlapping_boundary: 1
*/


#endif
