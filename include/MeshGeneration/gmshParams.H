#ifndef NEMOSYS_GMSHPARAMS_H_
#define NEMOSYS_GMSHPARAMS_H_

#include "nemosys_export.h"
#include "MeshGeneration/meshingParams.H"

#include <array>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>

namespace NEM {

namespace GEN {

/**
  * @brief A structure for defining volumetric mesh size fields.
  Size fields supported in this structure are Ball, Cylinder, Box, and Frustum.
**/
struct NEMOSYS_EXPORT volSizeField {
  std::string type; /**< Type of size field, eg. Ball,Cylinder,Box */
  int id;           /**< Size Field ID */
  /**< Vector of pairs to store size field parameters */
  std::vector<std::pair<std::string, double>> params;

  /**< Vector of pairs to store list of numbers */
  std::vector<std::pair<std::string, std::vector<double>>> num_list_params;

  /**< Vector of pairs to store list of strings */
  std::vector<std::pair<std::string, std::vector<std::string>>>
      strg_list_params;
};

/**
 @brief A struct for defining hexahedral transfinite volumes.
 */
struct NEMOSYS_EXPORT TransfiniteBlock {
  // specified gmsh id (tag) of volume to be treated as transfinite
  int id;

  // 3 direction vectors specifying axis of block (hexahedron)
  std::array<std::array<double, 3>, 3> axis;

  // number of vertices in each axial direction
  std::array<int, 3> vert;

  // each axis has either "Progression" or "Bump" type
  // default is "Progression" with coef 1
  std::array<std::string, 3> type{"Progression", "Progression", "Progression"};

  // value used in "Progression" or "Bump" along given axis
  std::array<double, 3> coef{1.0, 1.0, 1.0};

  bool operator<(const TransfiniteBlock &rhs) const { return id < rhs.id; }
};

/**
  @brief gmshParams contains all parameters essential for mesh generation
      using gmshGen class methods. These parameters are assigned a value
      at time of user input parsing through JSON in meshGenDriver.
**/
class NEMOSYS_EXPORT gmshParams : public meshingParams {
 public:
  /** @brief gmshParams standard constructor
   **/
  gmshParams() = default;

  /** @brief Output mesh file name
   **/
  std::string ofname{};

  /** @brief Minimum global mesh size
   **/
  double minSize{0.01};

  /** @brief Maximum global mesh size
   **/
  double maxSize{50.0};

  /** @brief Surface meshing algorithm
   **/
  std::string algo2D{"Frontal"};

  /** @brief Volume meshing algorithm
   **/
  std::string algo3D{"HXT"};

  /** @brief Extend mesh size from boundary option
   **/
  bool extSizeFromBoundary{true};

  /** @brief Mesh size based on curvature option
   **/
  bool sizeFromCurvature{false};

  /** @brief Minimum number of mesh elements per two Pi
   **/
  int minElePer2Pi{6};

  /** @brief Whether to optimize mesh or not
   **/
  bool optimize{false};

  /** @brief Mesh optimization threshold, between 0 and 1
   **/
  double optimizeThreshold{0.3};

  /** @brief Size field ID to use as background field
   **/
  int bgField{-1};

  /**
   * @brief Element order
   */
  int elementOrder{1};

  /**
   * @brief Subdivision algorithm 0: none, 1: all quads, 2: all hexas
   */
  int subdivisionAlg{1};

  /**
   * @brief Save all elements
   */
  bool saveAll{true};

  /**
   * @brief Whether to call boolean fragments on all volumes
   */
  bool fragmentAll{false};

  /** @brief Vector for volSizeField struct
   **/
  std::vector<volSizeField> sizeFields;

  /** @brief Map from RGB to physical name
   **/
  std::map<std::array<int, 3>, std::string> color2groupMap;

  /** @brief Map from volume id to transfinite hexahedron information.
   **/
  std::set<TransfiniteBlock> transfiniteBlocks;

  /**
   * @brief Get list of file extensions supported by gmsh
   * @return supported gmsh mesh file extensions
   */
  static const std::array<const char *, 6> &getMeshExtensions() {
    static constexpr std::array<const char *, 6> meshExtensions{
        ".inp", ".unv", ".p3d", ".stl", ".vtk", ".su2"};
    return meshExtensions;
  }
};

}  // namespace GEN

}  // namespace NEM

#endif  // NEMOSYS_GMSHPARAMS_H_
