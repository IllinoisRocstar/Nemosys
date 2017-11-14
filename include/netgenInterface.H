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
      if(mesh) nglib::Ng_DeleteMesh(mesh);
      nglib::Ng_Exit();
    }      

  // netgen mesh creation
  public:
    int createMeshFromSTL(char* fname);
    int importFromVol(char* fname);
    int exportToVTK(char* fname);
    int importFromVTK(char* fname); 
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

//   int uselocalh;                      //!< Switch to enable / disable usage of local mesh size modifiers
//
//   double maxh;                        //!< Maximum global mesh size allowed
//   double minh;                        //!< Minimum global mesh size allowed
//
//   double fineness;                    //!< Mesh density: 0...1 (0 => coarse; 1 => fine)
//   double grading;                     //!< Mesh grading: 0...1 (0 => uniform mesh; 1 => aggressive local grading)
//
//   double elementsperedge;             //!< Number of elements to generate per edge of the geometry
//   double elementspercurve;            //!< Elements to generate per curvature radius
//
//   int closeedgeenable;                //!< Enable / Disable mesh refinement at close edges
//   double closeedgefact;               //!< Factor to use for refinement at close edges (larger => finer)
//
//   int minedgelenenable;         //!< Enable / Disable user defined minimum edge length for edge subdivision
//   double minedgelen;                  //!< Minimum edge length to use while subdividing the edges (default = 1e-4)
//
//   int second_order;                   //!< Generate second-order surface and volume elements
//   int quad_dominated;                 //!< Creates a Quad-dominated mesh 
//
//   char * meshsize_filename;           //!< Optional external mesh size file 
//
//   int optsurfmeshenable;              //!< Enable / Disable automatic surface mesh optimization
//   int optvolmeshenable;               //!< Enable / Disable automatic volume mesh optimization
//
//   int optsteps_3d;                     //!< Number of optimize steps to use for 3-D mesh optimization
//   int optsteps_2d;                     //!< Number of optimize steps to use for 2-D mesh optimization
//
//   // Philippose - 13/09/2010
//   // Added a couple more parameters into the meshing parameters list 
//   // from Netgen into Nglib
//   int invert_tets;                    //!< Invert all the volume elements
//   int invert_trigs;                   //!< Invert all the surface triangle elements
//
//   int check_overlap;                  //!< Check for overlapping surfaces during Surface meshing
//   int check_overlapping_boundary;     //!< Check for overlapping surface elements before volume meshing
#endif
