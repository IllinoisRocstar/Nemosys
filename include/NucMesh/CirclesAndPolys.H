/*******************************************************************************
* Promesh                                                                      *
* Copyright (C) 2022, IllinoisRocstar LLC. All rights reserved.                *
*                                                                              *
* Promesh is the property of IllinoisRocstar LLC.                              *
*                                                                              *
* IllinoisRocstar LLC                                                          *
* Champaign, IL                                                                *
* www.illinoisrocstar.com                                                      *
* promesh@illinoisrocstar.com                                                  *
*******************************************************************************/
/*******************************************************************************
* This file is part of Promesh                                                 *
*                                                                              *
* This version of Promesh is free software: you can redistribute it and/or     *
* modify it under the terms of the GNU Lesser General Public License as        *
* published by the Free Software Foundation, either version 3 of the License,  *
* or (at your option) any later version.                                       *
*                                                                              *
* Promesh is distributed in the hope that it will be useful, but WITHOUT ANY   *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS    *
* FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more *
* details.                                                                     *
*                                                                              *
* You should have received a copy of the GNU Lesser General Public License     *
* along with this program. If not, see <https://www.gnu.org/licenses/>.        *
*                                                                              *
*******************************************************************************/
#ifndef NEMOSYS_CIRCLESANDPOLYS_H_
#define NEMOSYS_CIRCLESANDPOLYS_H_

#include "nemosys_export.h"
#include "NucMesh/ShapeBase.H"

#include <array>
#include <utility>
#include <vector>

namespace NEM {
namespace NUCMESH {

struct PolyRing;
struct Ring;

/**
 * @brief A set of concentric circles and polygons and the faces they enclose
 */
class NEMOSYS_EXPORT CirclesAndPolys : public ShapeBase {
 public:
  /**
   * Container to describe a set of concentric circles/polygons
   * @param numSides Number of sides for polygons (note circles are treated as
   * @p numSides arcs)
   * @param rings Radius and meshing description of each ring
   * @param center Center of shapes
   */
  CirclesAndPolys(int numSides, std::vector<PolyRing> rings,
                  const std::array<double, 3> &center = {0, 0, 0});

  int getNumSides() const;
  void setNumSides(int numSides);
  const std::vector<PolyRing> &getRings() const;
  void setRings(std::vector<PolyRing> rings);
  void addRing(const PolyRing &ring);

 protected:
  NEM::GEO::GeoManager createGeo() const override;

 private:
  int numSides_;
  std::vector<PolyRing> rings_;
};

/**
 * @brief A set of concentric circles and the faces they enclose
 */
class NEMOSYS_EXPORT Circles : public ShapeBase {
 public:
  /**
   * Container to describe a set of concentric circles (note circles treated as
   * two 180 degree arcs)
   * @param rings Radius and meshing description of each ring
   * @param center Center of shapes
   */
  explicit Circles(std::vector<Ring> rings,
                   const std::array<double, 3> &center = {0, 0, 0});

  const std::vector<Ring> &getRings() const;
  void setRings(std::vector<Ring> rings);
  void addRing(const Ring &ring);

 protected:
  NEM::GEO::GeoManager createGeo() const override;

 private:
  std::vector<Ring> rings_;
};

struct NEMOSYS_EXPORT RingMeshOption {
 public:
  enum class MeshingType {
    TRI,
    QUAD,
    STRUCT,
  };

  static RingMeshOption ApplyTriMesh() {
    return RingMeshOption{MeshingType::TRI};
  };
  static RingMeshOption ApplyQuadMesh() {
    return RingMeshOption{MeshingType::QUAD};
  };
  static RingMeshOption ApplyStructuredMesh(
      const std::array<int, 2> &numElems) {
    return RingMeshOption{numElems};
  }

  MeshingType meshingType; // Mesh type struct
  // Structured mesh number of segments (radial,circum)
  std::array<int, 2> numSegments;

 private:
  explicit RingMeshOption(MeshingType type);
  explicit RingMeshOption(const std::array<int, 2> &numElems);
};

struct NEMOSYS_EXPORT Ring {
  /**
   * Describe a set of faces in a @p Circles
   * @param radius Outer radius
   * @param meshType Options for meshing
   * @param material Name of material
   * @param sideset Name of sideset
   */
  Ring(double radius, const RingMeshOption &meshType, std::string material = {},
       std::string sideset = {})
      : radius(radius),
        meshType(meshType),
        material(std::move(material)),
        sideset(std::move(sideset)) {}
  double radius;
  RingMeshOption meshType;
  std::string material;
  std::string sideset;
};

struct NEMOSYS_EXPORT PolyRing : public Ring {
  enum class ShapeType { CIRCLE, POLYGON };
  /**
   * Describe a set of faces in a @p CirclesAndPolys
   * @param shapeType Circular or Polygon ring
   * @param radius Outer radius
   * @param meshType Options for meshing
   * @param rotation Rotation of outer edges
   * @param material Material name string
   * @param sideset Sideset name string
   */
  PolyRing(ShapeType shapeType, double radius, const RingMeshOption &meshType,
           double rotation = 0, std::string material = {},
           std::string sideset = {})
      : Ring(radius, meshType, std::move(material), std::move(sideset)),
        shapeType(shapeType),
        rotation(rotation) {}
  ShapeType shapeType;
  double rotation;
};

}  // namespace NUCMESH
}  // namespace NEM

#endif  // NEMOSYS_CIRCLESANDPOLYS_H_
