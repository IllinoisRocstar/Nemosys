#ifndef NEMOSYS_FUELELEMENT_H_
#define NEMOSYS_FUELELEMENT_H_

#include "nemosys_export.h"
#include "NucMesh/ShapeBase.H"

#include <array>
#include <map>
#include <utility>
#include <vector>

namespace NEM {
namespace NUCMESH {

struct Plate;

/**
 * @brief The ATR Fuel Element
 */
class NEMOSYS_EXPORT FuelElement : public ShapeBase {
 public:
  /**
   * Container to describe the set of fuel plates
   * @param startRadius The starting radius of the first plate
   * @param angle Angle arc of fuel element
   * @param plates Radius and meshing description of each plate
   * @param center Center of shapes
   */
  FuelElement(double startRadius, double angle, std::vector<Plate> plates,
              const std::array<double, 3> &center = {0, 0, 0});

  /**
   * Gets the starting radius of the fuel element
   * @return The starting radius
   */
  double getStartRadius() const;

  /**
   * Sets the starting radius of the fuel element
   * @param radius The input starting radius
   */
  void setStartRadius(double radius);

  /**
   * Gets the angle of the fuel element
   * @return The angle
   */
  double getAngle() const;

  /**
   * Sets the angle of the fuel element
   * @param angle The input angle
   */
  void setAngle(double angle);

  /**
   * Gets the vector of plates
   * @return Vector of plates
   */
  const std::vector<Plate> &getPlates() const;

  /**
   * Sets the vector of plates
   * @param plates Vector of plates
   */
  void setPlates(std::vector<Plate> plates);

  /**
   * Adds a plate to the container of plates
   * @param plate Plate to be added
   */
  void addPlate(const Plate &plate);

 protected:
  NEM::GEO::GeoManager createGeo() const override;

 private:
  double startRadius_;         // starting radius of fuel element
  double angle_;               // angle of fuel element
  std::vector<Plate> plates_;  // vector of plates
};

struct NEMOSYS_EXPORT PlateMeshOption {
 public:
  enum class MeshingType {
    TRI,
    QUAD,
    STRUCT,
  };

  static PlateMeshOption ApplyTriMesh() {
    return PlateMeshOption{MeshingType::TRI};
  };
  static PlateMeshOption ApplyQuadMesh() {
    return PlateMeshOption{MeshingType::QUAD};
  };
  static PlateMeshOption ApplyStructuredMesh(
      const std::array<int, 2> &numElems) {
    return PlateMeshOption{numElems};
  }

  MeshingType meshingType;  // Mesh type struct
  // Structured mesh number of segments (radial,circum)
  std::array<int, 2> numSegments;

 private:
  explicit PlateMeshOption(MeshingType type);
  explicit PlateMeshOption(const std::array<int, 2> &numElems);
};

struct NEMOSYS_EXPORT Plate {
  /**
   * Describe a set of faces in a @p Fuel Element
   * @param thickness Plate thickness
   * @param fuelThickness Fuel meat thickness
   * @param sideWidth Width of side section of plate
   * @param endWidth Width of end section of plate
   * @param fuelOffset Offset to fuel
   * @param meshType_map Map of options for meshing
   * @param materials_map Map of materials
   * @param sideset_map Name of sideset (optional)
   */
  Plate(double thickness, double fuelThickness, double sideWidth,
        double endWidth, double fuelOffset,
        const std::map<std::string, PlateMeshOption> &meshType_map,
        const std::map<std::string, std::string> &materials_map,
        const std::map<std::string, std::string> &sideset_map = {})
      : thickness(thickness),
        fuelThickness(fuelThickness),
        sideWidth(sideWidth),
        endWidth(endWidth),
        fuelOffset(fuelOffset),
        meshType_map(meshType_map),
        materials_map(materials_map),
        sideset_map(sideset_map) {}
  double thickness, fuelThickness, sideWidth, endWidth, fuelOffset;
  std::map<std::string, PlateMeshOption> meshType_map;
  std::map<std::string, std::string> materials_map;
  std::map<std::string, std::string> sideset_map;
};

}  // namespace NUCMESH
}  // namespace NEM

#endif  // NEMOSYS_FUELELEMENT_H_
