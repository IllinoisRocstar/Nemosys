#ifndef NEMOSYS_NUCMESHDRIVER_H_
#define NEMOSYS_NUCMESHDRIVER_H_

#include "nemosys_export.h"
#include "Drivers/NemDriver.H"

#include "shape.H"

namespace NEM {
namespace DRV {

/**
 * @brief Driver class to parse input JSON file and create mixed-element meshes.
 */
class NEMOSYS_EXPORT NucMeshDriver : public NemDriver {
 public:
  using Files = DriverOutFile;
  NucMeshDriver(Files file, jsoncons::json geometryAndMesh);

  const Files &getFiles() const;
  void setFiles(Files file);
  const jsoncons::json &getGeometryAndMesh() const;
  void setGeometryAndMesh(jsoncons::json geometryAndMesh);
  void execute() const override;

  JSONCONS_TYPE_TRAITS_FRIEND

 private:
  NucMeshDriver();

  Files file_;
  jsoncons::json geometryAndMesh_;

  static constexpr const char *programType = "NucMesh Generation";
  jsoncons::string_view getProgramType() const override;
};

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_NUCMESHDRIVER_H_
