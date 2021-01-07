#ifndef GMSHPARAMS_H
#define GMSHPARAMS_H

#include "meshingParams.H"

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
struct volSizeField {
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
struct TransfiniteBlock {
  // specified gmsh id (tag) of volume to be treated as transfinite
  size_t id;

  // 3 direction vectors specifying axis of block (hexahedron)
  double axis[3][3];

  // number of vertices in each axial direction
  int vert[3];

  // each axis has either "Progression" or "Bump" type
  // default is "Progression" with coef 1
  std::string type[3] = { "Progression", "Progression", "Progression" };

  // value used in "Progression" or "Bump" along given axis
  double coef[3] = { 1.0, 1.0, 1.0 };
};

/**
  @brief gmshParams contains all parameters essential for mesh generation
      using gmshGen class methods. These parameters are assigned a value
      at time of user input parsing through JSON in meshGenDriver.
**/
class gmshParams : public meshingParams {
 public:
  /** @brief gmshParams standard constructor
   **/
  gmshParams();

  /** @brief gmshParams standard destructor
   **/
  ~gmshParams(){};

  /** @brief Output mesh file name
   **/
  std::string ofname;

  /** @brief Minimum global mesh size
   **/
  double minSize;

  /** @brief Maximum global mesh size
   **/
  double maxSize;

  /** @brief Surface meshing algorithm
   **/
  std::string algo2D;

  /** @brief Volume meshing algorithm
   **/
  std::string algo3D;

  /** @brief Extend mesh size from boundary option
   **/
  bool extSizeFromBoundary;

  /** @brief Mesh size based on curvature option
   **/
  bool sizeFromCurvature;

  /** @brief Minimum number of mesh elements per two Pi
   **/
  int minElePer2Pi;

  /** @brief Whether to optimize mesh or not
   **/
  bool optimize;

  /** @brief Mesh optimization threshold, between 0 and 1
   **/
  double optimizeThreshold;

  /** @brief Mesh size field bool
   **/
  bool mSizeField;

  /** @brief Mesh color map bool
   **/
  bool mColorMap;

  /** @brief Mesh transfinite volume bool
   **/
  bool mTransfiniteVolumes;

  /** @brief Size field ID to use as background field
   **/
  int bgField;

  /**
   * @brief Element order
   */
  int elementOrder;

  /**
   * @brief Subdivision algorithm 0: none, 1: all quads, 2: all hexas
   */
  int subdivisionAlg;

  /**
   * @brief Save all elements
   */
  bool saveAll;

  /**
   * @brief Whether to call boolean fragments on all volumes
   */
  bool fragmentAll;

  /** @brief Vector for volSizeField struct
   **/
  std::vector<volSizeField> sizeFields;

  /** @brief Vector of supported gmsh mesh file extensions
   **/
  std::set<std::string> meshExtensions;

  /** @brief Map from RGB float tuple (cast to string) to physical name
   **/
  std::map<std::string, std::string> color2groupMap;

  /** @brief Map from volume id to transfinite hexahedron information.
   **/
  std::map<size_t, TransfiniteBlock> transfiniteBlocks;
};

}  // namespace GEN

}  // namespace NEM

#endif
