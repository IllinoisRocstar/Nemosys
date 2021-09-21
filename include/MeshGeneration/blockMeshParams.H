#ifndef NEMOSYS_BLOCKMESHPARAMS_H_
#define NEMOSYS_BLOCKMESHPARAMS_H_

#ifdef HAVE_CFMSH

#include "nemosys_export.h"
#include "MeshGeneration/meshingParams.H"

#include <string>
#include <vector>
#include <map>
#include <tuple>
#include <jsoncons/config/jsoncons_config.hpp>  // for jsoncons::optional

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

// --- All blockMesh parameters are accepted in meters

// --- Base class for blockmesh shape options
struct NEMOSYS_EXPORT bmShape {
  virtual ~bmShape() = default;

 protected:
  bmShape() = default;
  bmShape(const bmShape &) = default;
  bmShape(bmShape &&) = default;
  bmShape &operator=(const bmShape &) = default;
  bmShape &operator=(bmShape &&) = default;
};

// --- Automatic Box Generation for Pack Mesh
struct NEMOSYS_EXPORT bmAutoGenBox {
  /** @brief Input surface file for generating block automatically
            surrounding whole geometry domain with 10% offset in all direction
  **/
  std::string packFileName;
  /** @brief offset in X/Y/Z direction for auto box generation
   **/
  std::array<double, 3> offset{0.1, 0.1, 0.1};
};

// --- Options for Box
struct NEMOSYS_EXPORT bmBox : public bmShape {
  /** @brief location of initial point
   **/
  std::array<double, 3> init{0., 0., 0.};  // location of initial point

  /** @brief length in X/Y/Z direction
   **/
  std::array<double, 3> len{1., 1., 1.};

  /** @brief defines grading in hex block (Default 1)
   **/
  std::array<double, 3> smplGrading{1., 1., 1.};

  /** @brief Pair of min and max box locations for autogenerate
   **/
  std::pair<std::array<double, 3>, std::array<double, 3>> coordsBox{
      {0., 0., 0.}, {0., 0., 0.}};

  jsoncons::optional<bmAutoGenBox> autoGenerate;
};

// -- Options for sphere mesh
struct NEMOSYS_EXPORT bmSphere : public bmShape {
  /** @brief Defines X/Y/Z coordinate of sphere center
   **/
  std::array<double, 3> center{0., 0., 0.};  // Sphere center

  /** @brief Defines radius of sphere
   **/
  double radius;  // Radius of sphere

  /** @brief User defined grading in X/Y/Z direction (Default 1)
   **/
  std::array<double, 3> sphrGrading{1., 1., 1.};
};

// -- Options for cylinder and tapered cone mesh
struct NEMOSYS_EXPORT bmCylTaperedCone : public bmShape {
  /** @brief Vector for storing cylinder center coordinates
   **/
  std::array<double, 3> centerCyl{0., 0., 0.};  // Center location

  /** @brief Radius on one end of cylinder
   **/
  double radius1;  // radius 1

  /** @brief Radius on the other end of cylinder; if not set, radius1 used
   **/
  jsoncons::optional<double> radius2;  // radius 2

  /** @brief Height of cylinder
   **/
  double height;  // height of cylinder

  /** @brief Vector of mesh grading in X, Y, and Z direction respectively
   **/
  std::array<double, 3> cylGrading{0., 0., 0.};  // defines grading
};

/**
    @brief blockMeshParams contains the parameters important for automatic
          meshing using blockMeshGen class. These parameters are assigned a
          value during user input pasing through JSON in MeshGenDriver.
**/
class NEMOSYS_EXPORT blockMeshParams : public meshingParams {
  // --- access
 public:
  /** @brief blockMeshParams standard constructor
   **/
  blockMeshParams() = default;  // constructor

  // --- Booleans for shape preferences
  /** @brief Enables option to use user provided blockMeshDict
   **/
  // TODO: Unused
  bool ownBlockMshDict{false};  // User provided blockMeshDict

  /** @brief 3D block, sphere, cylinder, or tapered cone
   **/
  std::shared_ptr<bmShape> shape{};

  // --- Some general options
  /** @brief Defines mesh scale to meters
   **/
  double cnvrtToMeters{1};  // Scale to meters

  /** @brief Defines number of cells in X/Y/Z direction; unused if cellSize has
      a value
   **/
  std::array<int, 3> nCells{40, 40, 40};

  /** @brief user can define mesh cell size here (if not set, nCells used
      instead).
   **/
  jsoncons::optional<double> cellSize{};  // mesh density defined

  /** @brief Boolean for if the operation is packmesh
   **/
  bool isPackMesh{true};
};

#endif

#endif  // NEMOSYS_BLOCKMESHPARAMS_H_
