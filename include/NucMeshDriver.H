#ifndef NEMOSYS_NUCMESHDRIVER_H_
#define NEMOSYS_NUCMESHDRIVER_H_

#include "nemosys_export.h"
#include "AuxiliaryFunctions.H"
#include "NemDriver.H"
#include "shape.H"

/**
 * @brief Driver class to parse input JSON file and create mixed-element meshes.
 */
class NEMOSYS_EXPORT NucMeshDriver : public NemDriver {
 public:
  /**
   * @brief Constructs the geometry and generates mesh from input JSON
   * @param shapes Geometry and Mesh JSON object from input JSON file
   */
  explicit NucMeshDriver(jsoncons::json shapes);

  ~NucMeshDriver() override;

  /**
   * @brief high level driver method
   * @param inputjson input JSON file
   * @return NucMesh constructor
   */
  static NucMeshDriver *readJSON(const jsoncons::json &inputjson);

 public:
  /**
   * @brief Parses the Geometry and Mesh data of JSON
   * @param shapes Geometry and Mesh JSON object
   */
  void parseGeomAndMesh(jsoncons::json shapes);

  /**
   * @brief Parses the global options data of JSON
   * @param opts Global Options JSON object
   */
  void parseOptions(jsoncons::json opts);

  /**
   * @brief Parses the saved objects data of JSON
   * @param savedObj Saved Objects JSON object
   */
  void parseSavedObjects(jsoncons::json savedObj);

  /**
   * @brief Augments the shapes radii so that the meshed area is conserved
   * @param obj JSON object
   * @param shapetype The shape type (circle,polygon) being evaluated
   * @param tolerance The iterative solver tolerance
   */
  jsoncons::json correctMeshArea(jsoncons::json obj, std::string shapetype,
                                 double tolerance = 1e-8);

  /**
   * @brief Parses circle data of JSON and constructs Circles
   * @param circ Circle JSON object
   * @param conserving Boolean for whether mesh conservation is being done for
   * that circle
   */
  void makeCircles(jsoncons::json circ, bool conserving = false);

  /**
   * @brief Parses circleInPolys data of JSON and constructs CirclesInPolys
   * @param circInPoly CircleInPoly JSON object
   * @param conserving Boolean for whether mesh conservation is being done for
   * that circle in the polygon
   */
  void makeCirclesInPolys(jsoncons::json circInPoly, bool conserving = false);

  /**
   * @brief Parses polygon data of JSON and constructs Polygon
   * @param poly Polygon JSON object
   * @param ns number of sides of polygon
   */
  void makePolygons(jsoncons::json poly, int ns = 0);

  /**
   * @brief Parses array data of JSON and constructs Array
   * @param arr Array JSON object
   */
  void makeArray(jsoncons::json arr);

  /**
   * @brief Parses rectangular array data of JSON and constructs Array
   * @param arr Rectangular array JSON object
   */
  void makeRectangularArray(jsoncons::json arr);

  /**
   * @brief Parses cartesian array data of JSON and constructs Array
   * @param arr Cartesian array JSON object
   */
  void makeCartesianArray(jsoncons::json arr);

  /**
   * @brief Parses polar array data of JSON and constructs Array
   * @param arr Polar array JSON object
   */
  void makePolarArray(jsoncons::json arr);

  /**
   * @brief Parses hexagonal array data of JSON and constructs Array
   * @param arr Hexagonal array JSON object
   */
  void makeHexagonalArray(jsoncons::json arr);

  /**
   * @brief Checks that polygon radii are larger than last circle radius
   * @param circle_radii Circle radii
   * @param poly_radii Polygon radii
   * @param nSides Number of polygon sides
   */
  double checkInscribedCircle(std::vector<double> circle_radii,
                              std::vector<double> poly_radii, int nSides);

 private:
  std::map<int, NEM::GEO::shape *>
      shape_map; /**< map for shape index and shape
                  * @note used for identifying surface ids
                  * that change due to Boolean operations */
  std::map<std::string, int>
      phystag_map;       /**< map for physical entity tags and physical names */
  int id;                /**< shape object index
                          * @note used in shape_map */
  int physTag;           /**< physical entity tag id index
                          * @note used in phystag_map */
  std::string extension; /**< output file name extension */
  std::string ofname;    /**< output file name */
  std::map<std::string, std::pair<std::string, jsoncons::json>>
      savedobj_map; /**< map for Saved Objects declared in input JSON file
                     * @note contents of map are Alias, Shape Type, and Shape
                     * JSON object */
  bool extrude;     /**< bool to construct 3D mesh */
  bool gui;         /**< bool to open Gmsh gui */
  bool skipAll;     /**< bool to skip all break statements */
  std::vector<int> layers;     /** vector of 3D mesh layers */
  std::vector<double> heights; /** vector of 3D mesh layer heights */

  double paused_time; /** accumulates the paused time while GUI is open */

  int parsingCount; /** counter to make sure parseGeomAndMesh doesn't go to
                       infinity */

  /**
   * @brief Opens the Gmsh GUI
   */
  void openGUI();
};

#endif  // NEMOSYS_NUCMESHDRIVER_H_
