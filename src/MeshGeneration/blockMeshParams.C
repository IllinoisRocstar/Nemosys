#include "blockMeshParams.H"

// All the default parameters

blockMeshParams::blockMeshParams() // contructor body
{
  // Default geometry options
  _isBlock = true;
  _ownBlockMshDict = false;
  _isSphere = false;
  _isCylinder_TCone = false;

  // General options
  cnvrtToMeters = 1;  // Scale to meters
  cellsXDir = 80; // Cells in X direction
  cellsYDir = 80; // Cells in Y direction
  cellsZDir = 80; // Cells in Z direction

  // Block default options
  initX = 0;  // Initial point
  initY = 0;  // Initial point
  initZ = 0;  // Initial point
  lenX = 1; // Length from initial point
  lenY = 1; // Length from initial point
  lenZ = 1; // Length from initial point
  smplGradingX = 1; // Defines grading in each direction
  smplGradingY = 1; // Defines grading in each direction
  smplGradingZ = 1; // Defines grading in each direction

  // Sphere default options
  centerX = 0;  // Sphere center
  centerY = 0;  // Sphere center
  centerZ = 0;  // Sphere center
  sphrGradingX = 1; // Grading in each direction
  sphrGradingY = 1; // Grading in each direction
  sphrGradingZ = 1; // Grading in each direction

  std::vector<double> minPt = std::vector<double>(3,0.0);
  std::vector<double> maxPt = std::vector<double>(3,0.0);
  coordsBox = std::make_pair(minPt,maxPt);

}
