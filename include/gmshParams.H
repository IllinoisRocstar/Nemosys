#ifndef GMSHPARAMS_H
#define GMSHPARAMS_H

#include "meshingParams.H"

#include <set>
#include <string>
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

  /** @brief Size field ID to use as background field
   **/
  int bgField;

  /** @brief Vector for volSizeField struct
   **/
  std::vector<volSizeField> sizeFields;

  /** @brief Vector of supported gmsh mesh file extensions
   **/
  std::set<std::string> meshExtensions;
};

}  // namespace GEN

}  // namespace NEM

#endif
