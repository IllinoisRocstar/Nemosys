#ifndef BLOCKMESHPARAMS_H
#define BLOCKMESHPARAMS_H

#ifdef HAVE_CFMSH

#include "meshingParams.H"
#include <string>
#include <vector>
#include <map>
#include <tuple>

//***********************************************************************************************//

// TODO
// 1. Let user define patch names and types
// 2. Let user select how they want to define Box parameters
//    (i.e. 2 points, 1 point and lengths, center and lengths)

// In Development
/*struct bmshPatchDefine
{
  bool _isPatch;      // if given patch is a patch
  bool _isWall;     // if given patch is a wall
  bool _isSymmPlane;    // if given patch is a symmetric plane
  bool _isEmptyPtch;    // if given patch is an empty patch
  std::string ptchName;  // Patch Name
};*/


/** 
    @brief blockMeshParams contains the parameters important for automatic
          meshing using blockMeshGen class. These parameters are assigned a 
          value during user input pasing through JSON in MeshGenDriver.
**/
class blockMeshParams : public meshingParams
{
  // --- access
  public:
  /** @brief blockMeshParams standard constructor
  **/
  blockMeshParams(); // constructor

  /** @brief blockMeshParams standard destructor
  **/
  ~blockMeshParams(){}; //destructor
  
  // --- Booleans for shape preferences
  /** @brief Enables option to use user provided blockMeshDict
  **/
  bool _ownBlockMshDict;  // User provided blockMeshDict
  
  /** @brief Enables meshing of 3D block
  **/
  bool _isBlock;    // Structured Block mesh
  
  /** @brief Enables meshing of 3D sphere
  **/
  bool _isSphere;   // Structured Sphere mesh
  
  /** @brief Enables meshing of 3D cylinder or tapered cone
  **/
  bool _isCylinder_TCone; // Structured Cylinder or Tapered Cone Mesh
  
  /** @brief when enabled, allows user to define mesh cell size instead of
            number of cells in perticular direction.
  **/
  bool _cellSizeDefined;  // User defined mesh density
  
  // --- Some general options
  /** @brief Defines mesh scale to meters
  **/
  double cnvrtToMeters; // Scale to meters
  
  /** @brief Defines number of cells in X direction
  **/
  int cellsXDir;  // # of cells in X Direction
  
  /** @brief Defines number of cells in Y direction
  **/
  int cellsYDir;  // # of cells in Y Direction
  
  /** @brief Defines number of cells in Z direction
  **/
  int cellsZDir;  // # of cells Z Direction

  // --- All blockMesh parameters are accepted in meters
  // --- Automatic Box Generation for Pack Mesh
  /** @brief Input surface file for generating block automatically
            surrounding whole geometry domain with 10% offset in all direction 
  **/
  std::string packFileName;

  /** @brief Enables generation of block around geometry automatically
  **/
  bool _autoGenerateBox;  // just for pack meshing
  
  /** @brief offset in X direction for auto box generation
  **/
  double offsetX;     // offset for auto box generation
  
  /** @brief offset in Y direction for auto box generation
  **/
  double offsetY;     // offset for auto box generation
  
  /** @brief offset in Z direction for auto box generation
  **/
  double offsetZ;     // offset for auto box generation
  
  /** @brief If _cellSizeDefined option is enabled, user can define mesh cell
            size here.
  **/
  double cellSize;    // mesh density defined

  // --- Options for Box
  /** @brief location of initial point
  **/
  double initX;   // location of initial point
  
  /** @brief location of initial point
  **/
  double initY;   // location of initial point
  
  /** @brief location of initial point
  **/
  double initZ;   // location of initial point
  
  /** @brief length in respective direction
  **/
  double lenX;    // length in respective direction
  
  /** @brief length in respective direction
  **/
  double lenY;    // length in respective direction
  
  /** @brief length in respective direction
  **/
  double lenZ;    // length in respective direction
  
  /** @brief defines grading in hex block (Default 1)
  **/
  double smplGradingX;  // defines grading in hex block (Default 1)
  
  /** @brief defines grading in hex block (Default 1)
  **/
  double smplGradingY;  // defines grading in hex block (Default 1)
  
  /** @brief defines grading in hex block (Default 1)
  **/
  double smplGradingZ;  // defines grading in hex block (Default 1)

  // -- Options for sphere mesh
  /** @brief
  **/
  double centerX;   // Sphere center
  
  /** @brief Defines X coordinate of sphere center
  **/
  double centerY;   // Sphere center
  
  /** @brief Defines Y coordinate of sphere center
  **/
  double centerZ;   // Sphere center
  
  /** @brief Defines Z coordinate of sphere center
  **/
  double radius;    // Radius of sphere
  
  /** @brief User define grading in X direction (Default 1)
  **/
  double sphrGradingX;  // User define grading (Default 1)
  
  /** @brief User define grading in Y direction (Default 1)
  **/
  double sphrGradingY;  // User define grading (Default 1)
  
  /** @brief User define grading in Z direction (Default 1)
  **/
  double sphrGradingZ;  // User define grading (Default 1)

  // -- Options for cylinder and tapered cone mesh
  /** @brief Vector for storing cylinder center coordinates
  **/
  std::vector<double> centerCyl = std::vector<double>(3); //Center location
  
  /** @brief Radius on one end of cylinder
  **/
  double radius1;   // radius 1
  
  /** @brief Radius on the other end of cylinder
  **/
  double radius2;   // radius 2
  
  /** @brief Height of cylinder
  **/
  double height;    // height of cylinder
  
  /** @brief Vector of mesh grading in X, Y, and Z direction respectively
  **/
  std::vector<double> cylGrading = std::vector<double>(3);  //defines grading

};

#endif

#endif