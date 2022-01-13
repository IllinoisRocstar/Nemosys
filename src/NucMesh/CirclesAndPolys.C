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
#include "NucMesh/CirclesAndPolys.H"

#include <cassert>

#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <gce_MakeCirc.hxx>
#include <gce_MakePln.hxx>
#include <gp_Circ.hxx>
#include <gp_Pnt.hxx>

#include "Geometry/GeoManager.H"
#include "NucMesh/NucMeshShapeData.H"

namespace NEM {
namespace NUCMESH {

namespace {

constexpr double increment(int numSides) { return 360. / numSides; }

constexpr double baseRotation(int numSides) {
  return -1. * (180 - increment(numSides)) / 2.;
}

std::array<double, 3> getRingPoint(const std::array<double, 3> &center,
                                   int numSides, double radius, double rotation,
                                   int pointIdx) {
  return ShapeBase::getRotatedPoint(
      center, {radius, (pointIdx * increment(numSides) + rotation +
                        baseRotation(numSides))});
}

double modRotation(int numSides, double rotation) {
  return std::fmod(rotation, increment(numSides));
}

gp_Pnt coords2Point(const std::array<double, 3> &points) {
  return {points[0], points[1], points[2]};
}

std::pair<std::vector<TopoDS_Vertex>, std::vector<TopoDS_Edge>>
getRingVerticesAndCircEdges(
    const std::array<double, 3> &center, const int numSides,
    const double radius, double rotation,
    PolyRing::ShapeType shapeType = PolyRing::ShapeType::CIRCLE) {
  std::vector<TopoDS_Vertex> verts;
  std::vector<TopoDS_Edge> edges;
  assert(verts.empty() && edges.empty());
  verts.reserve(numSides);
  rotation = modRotation(numSides, rotation);
  for (int i = 0; i < numSides; ++i) {
    auto coords = getRingPoint(center, numSides, radius, rotation, i);
    verts.emplace_back(BRepBuilderAPI_MakeVertex{coords2Point(coords)});
  }
  if (shapeType == PolyRing::ShapeType::CIRCLE) {
    gp_Pnt centerPt{center[0], center[1], center[2]};
    auto circlePoint1 = getRingPoint(center, 2, radius, 0, 0);
    auto circlePoint2 = getRingPoint(center, 2, radius, 90., 0);
    gp_Circ circle{
        gce_MakeCirc{centerPt,
                     gce_MakePln{coords2Point(circlePoint1),
                                 coords2Point(circlePoint2), centerPt}
                         .Value(),
                     radius}};
    for (int i = 0; i < numSides; ++i) {
      edges.emplace_back(BRepBuilderAPI_MakeEdge{
          circle, verts[i], verts[i == numSides - 1 ? 0 : i + 1]});
    }
  } else {
    assert(shapeType == PolyRing::ShapeType::POLYGON);
    for (int i = 0; i < numSides; ++i) {
      edges.emplace_back(BRepBuilderAPI_MakeEdge{
          verts[i], verts[i == numSides - 1 ? 0 : i + 1]});
    }
  }
  return {std::move(verts), std::move(edges)};
}

std::vector<TopoDS_Edge> getRadialEdges(
    const std::vector<TopoDS_Vertex> &innerVerts,
    const std::vector<TopoDS_Vertex> &outerVerts) {
  assert(innerVerts.size() == outerVerts.size());
  std::vector<TopoDS_Edge> radialEdges;
  radialEdges.reserve(innerVerts.size());
  for (std::size_t i = 0; i < innerVerts.size(); ++i) {
    radialEdges.emplace_back(
        BRepBuilderAPI_MakeEdge{innerVerts[i], outerVerts[i]});
  }
  return radialEdges;
}

std::vector<TopoDS_Face> getFaces(const std::vector<TopoDS_Edge> &innerEdges,
                                  const std::vector<TopoDS_Edge> &outerEdges,
                                  const std::vector<TopoDS_Edge> &radialEdges) {
  assert(innerEdges.size() == outerEdges.size() &&
         innerEdges.size() == radialEdges.size());
  std::vector<TopoDS_Face> faces;
  faces.reserve(innerEdges.size());
  for (std::size_t i = 0; i < innerEdges.size(); ++i) {
    // Does BRepBuilderAPI_MakeWire care about orientation?
    faces.emplace_back(BRepBuilderAPI_MakeFace{BRepBuilderAPI_MakeWire{
        outerEdges[i], radialEdges[i == innerEdges.size() - 1 ? 0 : i + 1],
        innerEdges[i], radialEdges[i]}});
  }
  return faces;
}

}  // namespace

RingMeshOption::RingMeshOption(const std::array<int, 2> &numElems)
    : meshingType(MeshingType::STRUCT), numSegments(numElems) {}

RingMeshOption::RingMeshOption(MeshingType type)
    : meshingType(type), numSegments() {}

CirclesAndPolys::CirclesAndPolys(int numSides, std::vector<PolyRing> rings,
                                 const std::array<double, 3> &center)
    : ShapeBase(center), numSides_(numSides), rings_(std::move(rings)) {}

int CirclesAndPolys::getNumSides() const { return numSides_; }

void CirclesAndPolys::setNumSides(int numSides) { this->numSides_ = numSides; }

const std::vector<PolyRing> &CirclesAndPolys::getRings() const {
  return rings_;
}

void CirclesAndPolys::setRings(std::vector<PolyRing> rings) {
  this->rings_ = std::move(rings);
}

void CirclesAndPolys::addRing(const PolyRing &ring) {
  this->rings_.emplace_back(ring);
}

template <typename RingT, typename FVertCirc>
NEM::GEO::GeoManager drawNested(int numSides, const std::vector<RingT> &rings,
                                const std::array<double, 3> &center,
                                FVertCirc funcVertsCircEdges) {
  NEM::GEO::GeoManager output{2};
  if (rings.empty()) { return output; }

  TopoDS_Vertex centerVert{BRepBuilderAPI_MakeVertex{coords2Point(center)}};

  std::vector<std::vector<TopoDS_Vertex>> vertices{};  // along the rings
  vertices.reserve(rings.size());
  std::vector<std::vector<TopoDS_Edge>> circumferentialEdges{};
  circumferentialEdges.reserve(rings.size());
  for (const auto &ring : rings) {
    auto vertsAndEdges = funcVertsCircEdges(ring);
    vertices.emplace_back(std::move(vertsAndEdges.first));
    circumferentialEdges.emplace_back(std::move(vertsAndEdges.second));
  }

  std::vector<std::vector<TopoDS_Edge>> radialEdges{};
  radialEdges.reserve(rings.size());
  std::vector<std::vector<TopoDS_Face>> faces{};
  faces.reserve(rings.size());
  {
    auto prevVerts = vertices.begin();
    auto iterVerts = rings.size() >= 2 ? prevVerts + 1 : vertices.end();
    auto prevCircEdges = circumferentialEdges.begin();
    auto iterCircEdges =
        rings.size() >= 2 ? prevCircEdges + 1 : circumferentialEdges.end();
    {
      radialEdges.emplace_back();
      faces.emplace_back();
      auto &innerEdges = circumferentialEdges[0];
      BRepBuilderAPI_MakeWire wireBuilder{};
      for (const auto &edge : innerEdges) { wireBuilder.Add(edge); }
      faces.back().emplace_back(BRepBuilderAPI_MakeFace{wireBuilder});
    }
    for (; iterVerts != vertices.end();
         ++prevVerts, ++iterVerts, ++prevCircEdges, ++iterCircEdges) {
      radialEdges.emplace_back(getRadialEdges(*prevVerts, *iterVerts));
      faces.emplace_back(
          getFaces(*prevCircEdges, *iterCircEdges, radialEdges.back()));
    }
  }

  {
    {  // innermost ring
      switch (rings[0].meshType.meshingType) {
        case RingMeshOption::MeshingType::STRUCT:
        case RingMeshOption::MeshingType::QUAD:
          output.insertConstruct<QuadMeshSurface>(faces[0][0],
                                                  rings[0].material);
          break;
        case RingMeshOption::MeshingType::TRI:
        default:
          output.insertConstruct<TriMeshSurface>(faces[0][0],
                                                 rings[0].material);
          break;
      }
      if (!rings[0].sideset.empty()) {
        auto &ssetName = rings[0].sideset;
        for (int i = 0; i < numSides; ++i) {
          output.insertConstruct<SideSetEdge>(circumferentialEdges[0][i],
                                              ssetName);
        }
      }
    }
    auto prevCircEdges = circumferentialEdges.begin();
    auto iterCircEdges =
        rings.size() >= 2 ? prevCircEdges + 1 : circumferentialEdges.end();
    auto iterRadialEdges = radialEdges.begin() + 1;
    auto iterFaces = faces.begin() + 1;
    auto prevRing = rings.begin();
    auto iterRing = rings.size() >= 2 ? prevRing + 1 : rings.end();
    for (; iterRing != rings.end(); prevCircEdges = iterCircEdges,
                                    ++iterCircEdges, ++iterRadialEdges,
                                    ++iterFaces, ++iterRing, ++prevRing) {
      auto &meshingParams = iterRing->meshType;
      auto meshingType = meshingParams.meshingType;
      auto &prevMeshingParams = prevRing->meshType;
      for (int i = 0; i < numSides; ++i) {
        switch (meshingType) {
          case RingMeshOption::MeshingType::STRUCT:
          case RingMeshOption::MeshingType::QUAD:
            output.insertConstruct<QuadMeshSurface>(iterFaces->at(i),
                                                    iterRing->material);
            break;
          case RingMeshOption::MeshingType::TRI:
          default:
            output.insertConstruct<TriMeshSurface>(iterFaces->at(i),
                                                   iterRing->material);
        }
      }
      if (meshingType == RingMeshOption::MeshingType::STRUCT) {
        if (prevMeshingParams.meshingType ==
                RingMeshOption::MeshingType::STRUCT &&
            prevMeshingParams.numSegments[1] != meshingParams.numSegments[1]) {
          std::cerr << "Different number of nodes in circumferential direction "
                       "in adjacent rings may cause mesh errors.\n";
        }
        for (int i = 0; i < numSides; ++i) {
          // If previous ring is NOT struct mesh, ensure
          if (prevMeshingParams.meshingType !=
              RingMeshOption::MeshingType::STRUCT) {
            // In case added to map previously, override it
            output.getMap()[prevCircEdges->at(i)].reset(new EdgeSegments{
                prevRing->sideset, meshingParams.numSegments[1]});
          }
          output.insertConstruct<EdgeSegments>(iterCircEdges->at(i),
                                               iterRing->sideset,
                                               meshingParams.numSegments[1]);
          output.insertConstruct<EdgeSegments>(iterRadialEdges->at(i), "",
                                               meshingParams.numSegments[0]);
        }
      } else {
        if (!iterRing->sideset.empty()) {
          auto &ssetName = iterRing->sideset;
          for (int i = 0; i < numSides; ++i) {
            output.insertConstruct<SideSetEdge>(iterCircEdges->at(i), ssetName);
          }
        }
      }
    }
  }
  return output;
}

NEM::GEO::GeoManager CirclesAndPolys::createGeo() const {
  return drawNested(
      numSides_, rings_, this->getCenter(), [this](const PolyRing &ring) {
        return getRingVerticesAndCircEdges(this->getCenter(), this->numSides_,
                                           ring.radius, ring.rotation,
                                           ring.shapeType);
      });
}

Circles::Circles(std::vector<Ring> rings, const std::array<double, 3> &center)
    : ShapeBase(center), rings_(std::move(rings)) {}

const std::vector<Ring> &Circles::getRings() const { return rings_; }

void Circles::setRings(std::vector<Ring> rings) {
  this->rings_ = std::move(rings);
}

void Circles::addRing(const Ring &ring) { this->rings_.emplace_back(ring); }

NEM::GEO::GeoManager Circles::createGeo() const {
  return drawNested(2, rings_, this->getCenter(), [this](const Ring &ring) {
    return getRingVerticesAndCircEdges(this->getCenter(), 2, ring.radius, 0);
  });
}

}  // namespace NUCMESH
}  // namespace NEM