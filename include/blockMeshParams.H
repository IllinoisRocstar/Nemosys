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


class blockMeshParams : public meshingParams
{
  public:
  // Default initialization
  blockMeshParams(); // constructor
  ~blockMeshParams(){}; //destructor
  
  // Booleans for shape preferences
  bool _ownBlockMshDict;  // User provided blockMeshDict
  bool _isBlock;    // Structured Block mesh
  bool _isSphere;   // Structured Sphere mesh
  bool _isCylinder_TCone; // Structured Cylinder or Tapered Cone Mesh
  bool _cellSizeDefined;  // User defined mesh density


  // Some general options
  double cnvrtToMeters; // Scale to meters
  int cellsXDir;  // # of cells in X Direction
  int cellsYDir;  // # of cells in Y Direction
  int cellsZDir;  // # of cells Z Direction


  // All blockMesh parameters are accepted in meters

  // Automatic Box Generation for Pack Mesh
  std::string packFileName;
  bool _autoGenerateBox;  // just for pack meshing
  double offsetX;     // offset for auto box generation
  double offsetY;     // offset for auto box generation
  double offsetZ;     // offset for auto box generation
  double cellSize;    // mesh density defined

  // Options for Box
  double initX;   // location of initial point
  double initY;   // location of initial point
  double initZ;   // location of initial point
  double lenX;    // length in respective direction
  double lenY;    // length in respective direction
  double lenZ;    // length in respective direction
  double smplGradingX;  // defines grading in hex block (Default 1)
  double smplGradingY;  // defines grading in hex block (Default 1)
  double smplGradingZ;  // defines grading in hex block (Default 1)


  // Options for sphere mesh
  double centerX;   // Sphere center
  double centerY;   // Sphere center
  double centerZ;   // Sphere center
  double radius;    // Radius of sphere
  double sphrGradingX;  // User define grading (Default 1)
  double sphrGradingY;  // User define grading (Default 1)
  double sphrGradingZ;  // User define grading (Default 1)


  // Options for cylinder and tapered cone mesh
  std::vector<double> centerCyl = std::vector<double>(3); //Center location
  double radius1;   // radius 1
  double radius2;   // radius 2
  double height;    // height of cylinder
  std::vector<double> cylGrading = std::vector<double>(3);  //defines grading


};

#endif

#endif