#ifndef NEMOSYS_TEMPLATEMESHDRIVER_H_
#define NEMOSYS_TEMPLATEMESHDRIVER_H_

#include "nemosys_export.h"
#include "Drivers/NemDriver.H"

namespace NEM {
namespace DRV {

/**
 * @brief Driver class to generate meshes according to a template.
 */
class NEMOSYS_EXPORT TemplateMeshDriver : public NemDriver {
 public:
  using Files = DriverOutFile;

  struct NEMOSYS_EXPORT Opts {
    explicit Opts(std::string templateName);
    std::string templateName{};
    jsoncons::json templateParams{};
    JSONCONS_TYPE_TRAITS_FRIEND
   private:
    Opts() = default;
  };

  explicit TemplateMeshDriver(Files file, Opts opts);

  const Files &getFiles() const;
  void setFiles(Files file);
  const Opts &getOpts() const;
  void setOpts(Opts opts);
  void execute() const override;

  ~TemplateMeshDriver() override;

  JSONCONS_TYPE_TRAITS_FRIEND

 private:
  TemplateMeshDriver();

  Files file_;
  Opts opts_;

  static constexpr const char *programType = "Template Mesh Generation";
  jsoncons::string_view getProgramType() const override;

  /**
   * @brief Executes the template file
   * @param tplName Template file to be executed
   * @param outName Output mesh file name
   * @param params Parameters file name
   */
  static void executeTemplate(std::string tplName, std::string outName,
                              std::string params);

  /**
   * @brief Reads the parameters for Spiral Tape geometry/mesh from the input
   * JSON file and generates the appropriate parameters file
   * @param inputjson Input JSON object
   * @return The date and time stamped parameters file name
   */
  static std::string spiralTapePipe(jsoncons::json inputjson);

  /**
   * @brief Inserts the parameters file name into the template file at the
   * appropriate line
   * @param file Template file to be executed
   * @param params_name Parameters file name
   */
  static void insertParams(std::string file, std::string params_name);

  /**
   * @brief Encrypts the template file
   * @param inFile Input file name
   * @param outFile Output file name
   */
  static void encrypt(std::string inFile, std::string outFile);

  /**
   * @brief Decrypts the template file
   * @param filename Input file name
   * @return Returns the decrypted file name
   */
  static std::string decrypt(std::string filename);
};

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_TEMPLATEMESHDRIVER_H_
