#ifndef NEMOSYS_NUCMESHDRIVER_H_
#define NEMOSYS_NUCMESHDRIVER_H_

#include "nemosys_export.h"
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
   * @brief Parses the global options data of JSON
   * @param opts Global Options JSON object
   */
  static void parseOptions(jsoncons::json opts);

  /**
   * @brief Parses the saved objects data of JSON
   * @param savedObj Saved Objects JSON object
   */
  void parseSavedObjects(jsoncons::json savedObj);

  /**
   * @brief Parses circle data of JSON and constructs Circles
   * @param circ Circle JSON object
   */
  void makeCircles(jsoncons::json circ);

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
};

#endif  // NEMOSYS_NUCMESHDRIVER_H_
