#ifndef NEMOSYS_TEMPLATEMESHDRIVER_H_
#define NEMOSYS_TEMPLATEMESHDRIVER_H_

#include "AuxiliaryFunctions.H"
#include "NemDriver.H"
#include "nemosys_export.h"

/**
 * @brief Driver class to generate meshes according to a template.
 */
class NEMOSYS_EXPORT TemplateMeshDriver : public NemDriver {
public:
  /**
   * @brief Generates mesh for template and parameters defined in JSON file
   * @param inputjson Input JSON file
   */
  explicit TemplateMeshDriver(jsoncons::json inputjson);

  ~TemplateMeshDriver() override;

  /**
   * @brief high level driver method
   * @param inputjson input JSON file
   * @return NucMesh constructor
   */
  static TemplateMeshDriver *readJSON(const jsoncons::json &inputjson);

private:
  /**
   * @brief Executes the template file
   * @param tplName Template file to be executed
   * @param outName Output mesh file name
   * @param params Parameters file name
   */
  void executeTemplate(std::string tplName, std::string outName,
                       std::string params);

  /**
   * @brief Reads the parameters for Spiral Tape geometry/mesh from the input
   * JSON file and generates the appropriate parameters file
   * @param inputjson Input JSON object
   * @return The date and time stamped parameters file name
   */
  std::string spiralTapePipe(jsoncons::json inputjson);

  /**
   * @brief Inserts the parameters file name into the template file at the
   * appropriate line
   * @param file Template file to be executed
   */
  void insertParams(std::string file);

  /**
   * @brief Encrypts the template file
   * @param inFile Input file name
   * @param outFile Output file name
   */
  void encrypt(std::string inFile, std::string outFile);

  /**
   * @brief Decrypts the template file
   * @param filename Input file name
   * @return Returns the decrypted file name
   */
  std::string decrypt(std::string filename);

  std::string params_name; /**< time stamped parameters file name */
};

#endif // NEMOSYS_TEMPLATEMESHDRIVER_H_
