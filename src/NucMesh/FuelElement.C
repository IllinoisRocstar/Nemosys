#include "NucMesh/FuelElement.H"

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

#include <STEPControl_Writer.hxx>

namespace NEM {
namespace NUCMESH {

namespace {

constexpr double increment() { return 180.; }

std::array<double, 3> getRingPoint(const std::array<double, 3> &center,
                                   double radius, double rotation) {
  return ShapeBase::getRotatedPoint(center, {radius, rotation});
}

double getAngleFromOffset(double offset, double radius) {
  double theta = 2. * asin(offset / (2. * radius));  // radians
  return theta;
}

gp_Pnt coords2Point(const std::array<double, 3> &points) {
  return {points[0], points[1], points[2]};
}

std::array<double, 3> getPlatePoint(const std::array<double, 3> &center,
                                    double radius, double offset,
                                    double angle) {
  double theta;
  angle *= M_PI / 180.;
  if (offset >= 0.)
    theta = getAngleFromOffset(offset, radius);
  else if (offset < 0. && offset > -1.e10)
    theta = angle + getAngleFromOffset(offset, radius);
  else
    theta = angle;
  return {center[0] + radius * std::cos(theta),
          center[1] + radius * std::sin(theta), center[2]};
}

std::pair<std::vector<TopoDS_Vertex>, std::vector<TopoDS_Edge>>
getPlateVerticesAndCircEdges(const std::array<double, 3> &center, double radius,
                             double angle, int plateIndex, double fuelThickness,
                             double thickness, double endWidth,
                             double sideWidth, double fuelOffset) {
  // NOTE: radius here is the starting radius of this plate
  std::vector<TopoDS_Vertex> verts;
  std::vector<TopoDS_Edge> edges;
  std::vector<double> offsets = {0.,
                                 endWidth,
                                 endWidth + sideWidth,
                                 endWidth + sideWidth + fuelOffset,
                                 -1. * (endWidth + sideWidth + fuelOffset),
                                 -1. * (endWidth + sideWidth),
                                 -1. * endWidth,
                                 -1.e30};
  std::array<double, 3> circlePoint1{};
  std::array<double, 3> circlePoint2{};
  gp_Pnt centerPt;
  gp_Circ circle;
  // Make vertices & circumferential edges
  // If fuel meat exists
  if (fuelThickness > 0.) {
    centerPt = {center[0], center[1], center[2]};
    int vn = plateIndex == 0 ? 20 : 14;
    int en = plateIndex == 0 ? 16 : 10;
    verts.reserve(vn);
    edges.reserve(en);
    if (plateIndex == 0) {
      // Verts 0-2
      for (int i = 0; i < 3; ++i) {
        auto coords = getPlatePoint(center, radius, offsets[i], angle);
        auto vert = BRepBuilderAPI_MakeVertex{coords2Point(coords)};
        verts.emplace_back(vert);
      }
      // Verts 3-5
      for (int i = 5; i < 8; ++i) {
        auto coords = getPlatePoint(center, radius, offsets[i], angle);
        verts.emplace_back(BRepBuilderAPI_MakeVertex{coords2Point(coords)});
      }
      circlePoint1 = getRingPoint(center, radius, 0.);
      circlePoint2 = getRingPoint(center, radius, 90.);
      circle = {gce_MakeCirc{centerPt,
                             gce_MakePln{coords2Point(circlePoint1),
                                         coords2Point(circlePoint2), centerPt}
                                 .Value(),
                             radius}};
      for (int i = 0; i < 5; ++i) {
        edges.emplace_back(
            BRepBuilderAPI_MakeEdge{circle, verts[i], verts[i + 1]});
      }
    }
    // Verts 6-7
    radius += (thickness - fuelThickness) / 2.;
    for (int i = 2; i < 4; ++i) {
      auto coords = getPlatePoint(center, radius, offsets[i], angle);
      verts.emplace_back(BRepBuilderAPI_MakeVertex{coords2Point(coords)});
    }
    // Verts 8-9
    for (int i = 4; i < 6; ++i) {
      auto coords = getPlatePoint(center, radius, offsets[i], angle);
      verts.emplace_back(BRepBuilderAPI_MakeVertex{coords2Point(coords)});
    }
    circlePoint1 = getRingPoint(center, radius, 0.);
    circlePoint2 = getRingPoint(center, radius, 90.);
    circle = gce_MakeCirc{centerPt,
                          gce_MakePln{coords2Point(circlePoint1),
                                      coords2Point(circlePoint2), centerPt}
                              .Value(),
                          radius};
    int start = plateIndex > 0 ? 0 : 6;
    int end = plateIndex > 0 ? 3 : start + 3;
    for (int i = start; i < end; ++i) {
      edges.emplace_back(
          BRepBuilderAPI_MakeEdge{circle, verts[i], verts[i + 1]});
    }

    // Verts 10-11
    radius += fuelThickness;
    for (int i = 2; i < 4; ++i) {
      auto coords = getPlatePoint(center, radius, offsets[i], angle);
      verts.emplace_back(BRepBuilderAPI_MakeVertex{coords2Point(coords)});
    }
    // Verts 12-13
    for (int i = 4; i < 6; ++i) {
      auto coords = getPlatePoint(center, radius, offsets[i], angle);
      verts.emplace_back(BRepBuilderAPI_MakeVertex{coords2Point(coords)});
    }
    circlePoint1 = getRingPoint(center, radius, 0.);
    circlePoint2 = getRingPoint(center, radius, 90.);
    circle = gce_MakeCirc{centerPt,
                          gce_MakePln{coords2Point(circlePoint1),
                                      coords2Point(circlePoint2), centerPt}
                              .Value(),
                          radius};
    for (int i = start + 4; i < end + 4; ++i) {
      edges.emplace_back(
          BRepBuilderAPI_MakeEdge{circle, verts[i], verts[i + 1]});
    }
    // Verts 14-16
    radius += (thickness - fuelThickness) / 2.;
    for (int i = 0; i < 3; ++i) {
      auto coords = getPlatePoint(center, radius, offsets[i], angle);
      verts.emplace_back(BRepBuilderAPI_MakeVertex{coords2Point(coords)});
    }
    // Verts 17-19
    for (int i = 5; i < 8; ++i) {
      auto coords = getPlatePoint(center, radius, offsets[i], angle);
      verts.emplace_back(BRepBuilderAPI_MakeVertex{coords2Point(coords)});
    }
    circlePoint1 = getRingPoint(center, radius, 0.);
    circlePoint2 = getRingPoint(center, radius, 90.);
    circle = gce_MakeCirc{centerPt,
                          gce_MakePln{coords2Point(circlePoint1),
                                      coords2Point(circlePoint2), centerPt}
                              .Value(),
                          radius};
    // 14 19
    for (int i = start + 8; i < end + 10; ++i) {
      edges.emplace_back(
          BRepBuilderAPI_MakeEdge{circle, verts[i], verts[i + 1]});
    }
    // Else its a space plate
  } else {
    //  NOTE: If plateIndex > 0, don't make first row vertices or edges.
    //        Use the previous vertices and edges to make the next plate.
    int vn = plateIndex == 0 ? 12 : 6;
    int en = plateIndex == 0 ? 10 : 5;
    verts.reserve(vn);
    edges.reserve(en);
    if (plateIndex == 0) {
      // Verts 0-2
      for (int i = 0; i < 3; ++i) {
        auto coords = getPlatePoint(center, radius, offsets[i], angle);
        auto vert = BRepBuilderAPI_MakeVertex{coords2Point(coords)};
        verts.emplace_back(vert);
      }
      // Verts 3-5
      for (int i = 5; i < 8; ++i) {
        auto coords = getPlatePoint(center, radius, offsets[i], angle);
        verts.emplace_back(BRepBuilderAPI_MakeVertex{coords2Point(coords)});
      }
      centerPt = {center[0], center[1], center[2]};
      circlePoint1 = getRingPoint(center, radius, 0.);
      circlePoint2 = getRingPoint(center, radius, 90.);
      circle = {gce_MakeCirc{centerPt,
                             gce_MakePln{coords2Point(circlePoint1),
                                         coords2Point(circlePoint2), centerPt}
                                 .Value(),
                             radius}};
      for (int i = 0; i < 5; ++i) {
        edges.emplace_back(
            BRepBuilderAPI_MakeEdge{circle, verts[i], verts[i + 1]});
      }
    }
    // Verts 6-8
    radius += thickness;
    for (int i = 0; i < 3; ++i) {
      auto coords = getPlatePoint(center, radius, offsets[i], angle);
      verts.emplace_back(BRepBuilderAPI_MakeVertex{coords2Point(coords)});
    }
    // Verts 9-11
    for (int i = 5; i < 8; ++i) {
      auto coords = getPlatePoint(center, radius, offsets[i], angle);
      verts.emplace_back(BRepBuilderAPI_MakeVertex{coords2Point(coords)});
    }
    circlePoint1 = getRingPoint(center, radius, 0.);
    circlePoint2 = getRingPoint(center, radius, 90.);
    circle = gce_MakeCirc{centerPt,
                          gce_MakePln{coords2Point(circlePoint1),
                                      coords2Point(circlePoint2), centerPt}
                              .Value(),
                          radius};
    int start = plateIndex > 0 ? 0 : 6;
    int end = plateIndex > 0 ? 5 : start + 5;
    for (int i = start; i < end; ++i) {
      edges.emplace_back(
          BRepBuilderAPI_MakeEdge{circle, verts[i], verts[i + 1]});
    }
  }
  return {std::move(verts), std::move(edges)};
}

// Method to create first plate radial edges
std::vector<TopoDS_Edge> getFirstRadialEdges(
    const std::vector<TopoDS_Vertex> &innerVerts) {
  std::vector<TopoDS_Edge> radialEdges;
  // If a fuel plate
  if (innerVerts.size() > 6) {
    radialEdges.reserve(12);
    radialEdges.emplace_back(
        BRepBuilderAPI_MakeEdge{innerVerts[0], innerVerts[14]});
    radialEdges.emplace_back(
        BRepBuilderAPI_MakeEdge{innerVerts[1], innerVerts[15]});
    radialEdges.emplace_back(
        BRepBuilderAPI_MakeEdge{innerVerts[4], innerVerts[18]});
    radialEdges.emplace_back(
        BRepBuilderAPI_MakeEdge{innerVerts[5], innerVerts[19]});

    radialEdges.emplace_back(
        BRepBuilderAPI_MakeEdge{innerVerts[2], innerVerts[6]});
    radialEdges.emplace_back(
        BRepBuilderAPI_MakeEdge{innerVerts[3], innerVerts[9]});
    radialEdges.emplace_back(
        BRepBuilderAPI_MakeEdge{innerVerts[10], innerVerts[16]});
    radialEdges.emplace_back(
        BRepBuilderAPI_MakeEdge{innerVerts[13], innerVerts[17]});

    for (int i = 6; i < 10; i++) {
      radialEdges.emplace_back(
          BRepBuilderAPI_MakeEdge{innerVerts[i], innerVerts[i + 4]});
    }
    // Space Plate
  } else {
    radialEdges.reserve(5);
    for (int i = 0; i < 6; ++i) {
      radialEdges.emplace_back(
          BRepBuilderAPI_MakeEdge{innerVerts[i], innerVerts[i + 6]});
    }
  }
  return radialEdges;
}

// Get radial edges for plates 2 and up
std::vector<TopoDS_Edge> getRadialEdges(
    const std::vector<TopoDS_Vertex> &innerVerts,
    const std::vector<TopoDS_Vertex> &outerVerts, int plateIndex) {
  std::vector<TopoDS_Edge> radialEdges;
  // Have to make custom connections
  if (plateIndex == 1) {
    // Fuel after fuel
    if (innerVerts.size() > 6 && outerVerts.size() > 6) {
      radialEdges.reserve(12);
      radialEdges.emplace_back(
          BRepBuilderAPI_MakeEdge{innerVerts[14], outerVerts[8]});
      radialEdges.emplace_back(
          BRepBuilderAPI_MakeEdge{innerVerts[15], outerVerts[9]});
      radialEdges.emplace_back(
          BRepBuilderAPI_MakeEdge{innerVerts[18], outerVerts[12]});
      radialEdges.emplace_back(
          BRepBuilderAPI_MakeEdge{innerVerts[19], outerVerts[13]});

      radialEdges.emplace_back(
          BRepBuilderAPI_MakeEdge{innerVerts[16], outerVerts[0]});
      radialEdges.emplace_back(
          BRepBuilderAPI_MakeEdge{innerVerts[17], outerVerts[3]});
      radialEdges.emplace_back(
          BRepBuilderAPI_MakeEdge{outerVerts[4], outerVerts[10]});
      radialEdges.emplace_back(
          BRepBuilderAPI_MakeEdge{outerVerts[7], outerVerts[11]});

      for (int i = 0; i < 4; ++i) {
        radialEdges.emplace_back(
            BRepBuilderAPI_MakeEdge{outerVerts[i], outerVerts[i + 4]});
      }
      // Space after fuel
    } else if (innerVerts.size() > 6 && outerVerts.size() == 6) {
      radialEdges.reserve(6);
      for (int i = 0; i < 6; ++i) {
        radialEdges.emplace_back(
            BRepBuilderAPI_MakeEdge{innerVerts[i + 14], outerVerts[i]});
      }
      // Fuel after Space
    } else if (innerVerts.size() == 6 && outerVerts.size() > 6) {
      radialEdges.reserve(12);
      radialEdges.emplace_back(
          BRepBuilderAPI_MakeEdge{innerVerts[0], outerVerts[8]});
      radialEdges.emplace_back(
          BRepBuilderAPI_MakeEdge{innerVerts[1], outerVerts[9]});
      radialEdges.emplace_back(
          BRepBuilderAPI_MakeEdge{innerVerts[4], outerVerts[12]});
      radialEdges.emplace_back(
          BRepBuilderAPI_MakeEdge{innerVerts[5], outerVerts[13]});

      radialEdges.emplace_back(
          BRepBuilderAPI_MakeEdge{innerVerts[2], outerVerts[0]});
      radialEdges.emplace_back(
          BRepBuilderAPI_MakeEdge{innerVerts[9], outerVerts[3]});
      radialEdges.emplace_back(
          BRepBuilderAPI_MakeEdge{outerVerts[4], outerVerts[10]});
      radialEdges.emplace_back(
          BRepBuilderAPI_MakeEdge{outerVerts[7], outerVerts[11]});

      for (int i = 0; i < 4; ++i) {
        radialEdges.emplace_back(
            BRepBuilderAPI_MakeEdge{outerVerts[i], outerVerts[i + 4]});
      }
    }
  } else {
    // Fuel after fuel
    if (innerVerts.size() > 6 && outerVerts.size() > 6) {
      radialEdges.reserve(12);
      radialEdges.emplace_back(
          BRepBuilderAPI_MakeEdge{innerVerts[8], outerVerts[8]});
      radialEdges.emplace_back(
          BRepBuilderAPI_MakeEdge{innerVerts[9], outerVerts[9]});
      radialEdges.emplace_back(
          BRepBuilderAPI_MakeEdge{innerVerts[12], outerVerts[12]});
      radialEdges.emplace_back(
          BRepBuilderAPI_MakeEdge{innerVerts[13], outerVerts[13]});

      radialEdges.emplace_back(
          BRepBuilderAPI_MakeEdge{innerVerts[10], outerVerts[0]});
      radialEdges.emplace_back(
          BRepBuilderAPI_MakeEdge{innerVerts[11], outerVerts[3]});
      radialEdges.emplace_back(
          BRepBuilderAPI_MakeEdge{outerVerts[4], outerVerts[10]});
      radialEdges.emplace_back(
          BRepBuilderAPI_MakeEdge{outerVerts[7], outerVerts[11]});

      for (int i = 0; i < 4; ++i) {
        radialEdges.emplace_back(
            BRepBuilderAPI_MakeEdge{outerVerts[i], outerVerts[i + 4]});
      }
    }
    // Space after fuel
    else if (innerVerts.size() > 6 && outerVerts.size() == 6) {
      radialEdges.reserve(6);
      for (int i = 0; i < 6; ++i) {
        radialEdges.emplace_back(
            BRepBuilderAPI_MakeEdge{innerVerts[i + 8], outerVerts[i]});
      }
      // Fuel after Space
    } else if (innerVerts.size() == 6 && outerVerts.size() > 6) {
      radialEdges.reserve(12);
      radialEdges.emplace_back(
          BRepBuilderAPI_MakeEdge{innerVerts[0], outerVerts[8]});
      radialEdges.emplace_back(
          BRepBuilderAPI_MakeEdge{innerVerts[1], outerVerts[9]});
      radialEdges.emplace_back(
          BRepBuilderAPI_MakeEdge{innerVerts[4], outerVerts[12]});
      radialEdges.emplace_back(
          BRepBuilderAPI_MakeEdge{innerVerts[5], outerVerts[13]});

      radialEdges.emplace_back(
          BRepBuilderAPI_MakeEdge{innerVerts[2], outerVerts[0]});
      radialEdges.emplace_back(
          BRepBuilderAPI_MakeEdge{innerVerts[3], outerVerts[3]});
      radialEdges.emplace_back(
          BRepBuilderAPI_MakeEdge{outerVerts[4], outerVerts[10]});
      radialEdges.emplace_back(
          BRepBuilderAPI_MakeEdge{outerVerts[7], outerVerts[11]});

      for (int i = 0; i < 4; ++i) {
        radialEdges.emplace_back(
            BRepBuilderAPI_MakeEdge{outerVerts[i], outerVerts[i + 4]});
      }
    }
    // Space after Space
    else if (innerVerts.size() == 6 && outerVerts.size() == 6) {
      radialEdges.reserve(6);
      for (int i = 0; i < 6; ++i) {
        radialEdges.emplace_back(
            BRepBuilderAPI_MakeEdge{innerVerts[i], outerVerts[i]});
      }
    }
  }
  return radialEdges;
}

// Create the faces of the first plate
std::vector<TopoDS_Face> getFirstFaces(
    const std::vector<TopoDS_Edge> &innerEdges,
    const std::vector<TopoDS_Edge> &radialEdges) {
  std::vector<TopoDS_Face> faces;
  // Make faces in specific order
  // If Fuel Plate
  if (innerEdges.size() > 5) {
    faces.reserve(9);
    {
      faces.emplace_back(BRepBuilderAPI_MakeFace{BRepBuilderAPI_MakeWire{
          innerEdges[0], radialEdges[0], innerEdges[11], radialEdges[1]}});

      // Make wire with 6 edges, have to use Add method
      auto wire = BRepBuilderAPI_MakeWire{innerEdges[1], radialEdges[1],
                                          innerEdges[12], radialEdges[6]};
      wire.Add(radialEdges[8]);
      wire.Add(radialEdges[4]);
      faces.emplace_back(BRepBuilderAPI_MakeFace{wire});

      wire = BRepBuilderAPI_MakeWire{innerEdges[3], radialEdges[2],
                                     innerEdges[14], radialEdges[7]};
      wire.Add(radialEdges[11]);
      wire.Add(radialEdges[5]);
      faces.emplace_back(BRepBuilderAPI_MakeFace{wire});

      faces.emplace_back(BRepBuilderAPI_MakeFace{BRepBuilderAPI_MakeWire{
          innerEdges[4], radialEdges[3], innerEdges[15], radialEdges[2]}});
    }
    {
      auto wire = BRepBuilderAPI_MakeWire{innerEdges[2], radialEdges[4],
                                          innerEdges[5], innerEdges[6]};
      wire.Add(innerEdges[7]);
      wire.Add(radialEdges[5]);
      faces.emplace_back(BRepBuilderAPI_MakeFace{wire});

      wire = BRepBuilderAPI_MakeWire{innerEdges[8], radialEdges[6],
                                     innerEdges[13], radialEdges[7]};
      wire.Add(innerEdges[10]);
      wire.Add(innerEdges[9]);
      faces.emplace_back(BRepBuilderAPI_MakeFace{wire});
    }
    {
      faces.emplace_back(BRepBuilderAPI_MakeFace{BRepBuilderAPI_MakeWire{
          innerEdges[5], radialEdges[8], innerEdges[8], radialEdges[9]}});
      faces.emplace_back(BRepBuilderAPI_MakeFace{BRepBuilderAPI_MakeWire{
          innerEdges[6], radialEdges[9], innerEdges[9], radialEdges[10]}});
      faces.emplace_back(BRepBuilderAPI_MakeFace{BRepBuilderAPI_MakeWire{
          innerEdges[7], radialEdges[10], innerEdges[10], radialEdges[11]}});
    }
  }
  // Space Plate
  else {
    faces.reserve(5);
    faces.emplace_back(BRepBuilderAPI_MakeFace{BRepBuilderAPI_MakeWire{
        innerEdges[0], radialEdges[0], innerEdges[5], radialEdges[1]}});
    faces.emplace_back(BRepBuilderAPI_MakeFace{BRepBuilderAPI_MakeWire{
        innerEdges[1], radialEdges[1], innerEdges[6], radialEdges[2]}});
    faces.emplace_back(BRepBuilderAPI_MakeFace{BRepBuilderAPI_MakeWire{
        innerEdges[3], radialEdges[3], innerEdges[8], radialEdges[4]}});
    faces.emplace_back(BRepBuilderAPI_MakeFace{BRepBuilderAPI_MakeWire{
        innerEdges[4], radialEdges[4], innerEdges[9], radialEdges[5]}});
    faces.emplace_back(BRepBuilderAPI_MakeFace{BRepBuilderAPI_MakeWire{
        innerEdges[2], radialEdges[2], innerEdges[7], radialEdges[3]}});
  }
  return faces;
}

// Get the faces for plates 2 and up
std::vector<TopoDS_Face> getFaces(const std::vector<TopoDS_Edge> &innerEdges,
                                  const std::vector<TopoDS_Edge> &outerEdges,
                                  const std::vector<TopoDS_Edge> &radialEdges,
                                  int plateIndex) {
  std::vector<TopoDS_Face> faces;
  if (plateIndex == 1) {
    // Fuel after fuel
    if (innerEdges.size() > 5 && outerEdges.size() > 5) {
      faces.reserve(9);
      // ends and sides
      {
        faces.emplace_back(BRepBuilderAPI_MakeFace{BRepBuilderAPI_MakeWire{
            innerEdges[11], radialEdges[0], outerEdges[6], radialEdges[1]}});

        auto wire = BRepBuilderAPI_MakeWire{innerEdges[12], radialEdges[1],
                                            outerEdges[7], radialEdges[6]};
        wire.Add(radialEdges[8]);
        wire.Add(radialEdges[4]);
        faces.emplace_back(BRepBuilderAPI_MakeFace{wire});

        wire = BRepBuilderAPI_MakeWire{innerEdges[14], radialEdges[2],
                                       outerEdges[9], radialEdges[7]};
        wire.Add(radialEdges[11]);
        wire.Add(radialEdges[5]);
        faces.emplace_back(BRepBuilderAPI_MakeFace{wire});

        faces.emplace_back(BRepBuilderAPI_MakeFace{BRepBuilderAPI_MakeWire{
            innerEdges[15], radialEdges[2], outerEdges[10], radialEdges[3]}});
      }
      // plate
      {
        auto wire = BRepBuilderAPI_MakeWire{innerEdges[13], radialEdges[4],
                                            outerEdges[0], outerEdges[1]};
        wire.Add(outerEdges[2]);
        wire.Add(radialEdges[5]);
        faces.emplace_back(BRepBuilderAPI_MakeFace{wire});

        wire = BRepBuilderAPI_MakeWire{outerEdges[3], radialEdges[6],
                                       outerEdges[8], radialEdges[7]};
        wire.Add(outerEdges[5]);
        wire.Add(outerEdges[4]);
        faces.emplace_back(BRepBuilderAPI_MakeFace{wire});
      }
      // fuel meat
      {
        for (int i = 0; i < 3; ++i) {
          faces.emplace_back(BRepBuilderAPI_MakeFace{
              BRepBuilderAPI_MakeWire{outerEdges[i], radialEdges[i + 8],
                                      outerEdges[i + 3], radialEdges[i + 9]}});
        }
      }
    }  // Space after fuel
    else if (innerEdges.size() > 5 && outerEdges.size() == 5) {
      faces.reserve(5);
      faces.emplace_back(BRepBuilderAPI_MakeFace{BRepBuilderAPI_MakeWire{
          innerEdges[11], radialEdges[0], outerEdges[0], radialEdges[1]}});
      faces.emplace_back(BRepBuilderAPI_MakeFace{BRepBuilderAPI_MakeWire{
          innerEdges[12], radialEdges[1], outerEdges[1], radialEdges[2]}});
      faces.emplace_back(BRepBuilderAPI_MakeFace{BRepBuilderAPI_MakeWire{
          innerEdges[14], radialEdges[3], outerEdges[3], radialEdges[4]}});
      faces.emplace_back(BRepBuilderAPI_MakeFace{BRepBuilderAPI_MakeWire{
          innerEdges[15], radialEdges[4], outerEdges[4], radialEdges[5]}});
      // middle
      faces.emplace_back(BRepBuilderAPI_MakeFace{BRepBuilderAPI_MakeWire{
          innerEdges[13], radialEdges[2], outerEdges[2], radialEdges[3]}});

      // Fuel after Space
    } else if (innerEdges.size() == 5 && outerEdges.size() > 5) {
      faces.reserve(9);
      // ends and sides
      {
        faces.emplace_back(BRepBuilderAPI_MakeFace{BRepBuilderAPI_MakeWire{
            innerEdges[5], radialEdges[0], outerEdges[6], radialEdges[1]}});

        auto wire = BRepBuilderAPI_MakeWire{innerEdges[6], radialEdges[1],
                                            outerEdges[7], radialEdges[6]};
        wire.Add(radialEdges[8]);
        wire.Add(radialEdges[4]);
        faces.emplace_back(BRepBuilderAPI_MakeFace{wire});

        wire = BRepBuilderAPI_MakeWire{innerEdges[8], radialEdges[2],
                                       outerEdges[9], radialEdges[7]};
        wire.Add(radialEdges[11]);
        wire.Add(radialEdges[5]);
        faces.emplace_back(BRepBuilderAPI_MakeFace{wire});

        faces.emplace_back(BRepBuilderAPI_MakeFace{BRepBuilderAPI_MakeWire{
            innerEdges[9], radialEdges[2], outerEdges[10], radialEdges[3]}});
      }
      // plate
      {
        auto wire = BRepBuilderAPI_MakeWire{innerEdges[7], radialEdges[4],
                                            outerEdges[0], outerEdges[1]};
        wire.Add(outerEdges[2]);
        wire.Add(radialEdges[5]);
        faces.emplace_back(BRepBuilderAPI_MakeFace{wire});

        wire = BRepBuilderAPI_MakeWire{outerEdges[3], radialEdges[6],
                                       outerEdges[8], radialEdges[7]};
        wire.Add(outerEdges[5]);
        wire.Add(outerEdges[4]);
        faces.emplace_back(BRepBuilderAPI_MakeFace{wire});
      }
      // fuel meat
      {
        for (int i = 0; i < 3; ++i) {
          faces.emplace_back(BRepBuilderAPI_MakeFace{
              BRepBuilderAPI_MakeWire{outerEdges[i], radialEdges[i + 8],
                                      outerEdges[i + 3], radialEdges[i + 9]}});
        }
      }
      // Space after Space
    } else if (innerEdges.size() == 5 && outerEdges.size() == 5) {
      faces.reserve(5);
      faces.emplace_back(BRepBuilderAPI_MakeFace{BRepBuilderAPI_MakeWire{
          innerEdges[5], radialEdges[0], outerEdges[0], radialEdges[1]}});
      faces.emplace_back(BRepBuilderAPI_MakeFace{BRepBuilderAPI_MakeWire{
          innerEdges[6], radialEdges[1], outerEdges[1], radialEdges[2]}});
      faces.emplace_back(BRepBuilderAPI_MakeFace{BRepBuilderAPI_MakeWire{
          innerEdges[8], radialEdges[3], outerEdges[3], radialEdges[4]}});
      faces.emplace_back(BRepBuilderAPI_MakeFace{BRepBuilderAPI_MakeWire{
          innerEdges[9], radialEdges[4], outerEdges[4], radialEdges[5]}});
      // middle
      faces.emplace_back(BRepBuilderAPI_MakeFace{BRepBuilderAPI_MakeWire{
          innerEdges[7], radialEdges[2], outerEdges[2], radialEdges[3]}});
    }
  } else {
    // Fuel after fuel
    if (innerEdges.size() > 5 && outerEdges.size() > 5) {
      faces.reserve(9);
      // ends and sides
      {
        faces.emplace_back(BRepBuilderAPI_MakeFace{BRepBuilderAPI_MakeWire{
            innerEdges[6], radialEdges[1], outerEdges[6], radialEdges[0]}});

        auto wire = BRepBuilderAPI_MakeWire{innerEdges[7], radialEdges[1],
                                            outerEdges[7], radialEdges[6]};
        wire.Add(radialEdges[8]);
        wire.Add(radialEdges[4]);
        faces.emplace_back(BRepBuilderAPI_MakeFace{wire});

        wire = BRepBuilderAPI_MakeWire{innerEdges[9], radialEdges[2],
                                       outerEdges[9], radialEdges[7]};
        wire.Add(radialEdges[11]);
        wire.Add(radialEdges[5]);
        faces.emplace_back(BRepBuilderAPI_MakeFace{wire});

        faces.emplace_back(BRepBuilderAPI_MakeFace{BRepBuilderAPI_MakeWire{
            innerEdges[10], radialEdges[2], outerEdges[10], radialEdges[3]}});
      }
      // plate
      {
        auto wire = BRepBuilderAPI_MakeWire{innerEdges[8], radialEdges[4],
                                            outerEdges[0], outerEdges[1]};
        wire.Add(outerEdges[2]);
        wire.Add(radialEdges[5]);
        faces.emplace_back(BRepBuilderAPI_MakeFace{wire});

        wire = BRepBuilderAPI_MakeWire{outerEdges[3], radialEdges[6],
                                       outerEdges[8], radialEdges[7]};
        wire.Add(outerEdges[5]);
        wire.Add(outerEdges[4]);
        faces.emplace_back(BRepBuilderAPI_MakeFace{wire});
      }
      // fuel meat
      {
        for (int i = 0; i < 3; ++i) {
          faces.emplace_back(BRepBuilderAPI_MakeFace{
              BRepBuilderAPI_MakeWire{outerEdges[i], radialEdges[i + 8],
                                      outerEdges[i + 3], radialEdges[i + 9]}});
        }
      }
    }  // Space after fuel
    else if (innerEdges.size() > 5 && outerEdges.size() == 5) {
      faces.reserve(5);
      faces.emplace_back(BRepBuilderAPI_MakeFace{BRepBuilderAPI_MakeWire{
          innerEdges[6], radialEdges[0], outerEdges[0], radialEdges[1]}});
      faces.emplace_back(BRepBuilderAPI_MakeFace{BRepBuilderAPI_MakeWire{
          innerEdges[7], radialEdges[1], outerEdges[1], radialEdges[2]}});
      faces.emplace_back(BRepBuilderAPI_MakeFace{BRepBuilderAPI_MakeWire{
          innerEdges[9], radialEdges[3], outerEdges[3], radialEdges[4]}});
      faces.emplace_back(BRepBuilderAPI_MakeFace{BRepBuilderAPI_MakeWire{
          innerEdges[10], radialEdges[4], outerEdges[4], radialEdges[5]}});
      // middle
      faces.emplace_back(BRepBuilderAPI_MakeFace{BRepBuilderAPI_MakeWire{
          innerEdges[8], radialEdges[2], outerEdges[2], radialEdges[3]}});

      // Fuel after Space
    } else if (innerEdges.size() == 5 && outerEdges.size() > 5) {
      faces.reserve(9);
      // ends and sides
      {
        faces.emplace_back(BRepBuilderAPI_MakeFace{BRepBuilderAPI_MakeWire{
            innerEdges[0], radialEdges[0], outerEdges[6], radialEdges[1]}});

        auto wire = BRepBuilderAPI_MakeWire{innerEdges[1], radialEdges[1],
                                            outerEdges[7], radialEdges[6]};
        wire.Add(radialEdges[8]);
        wire.Add(radialEdges[4]);
        faces.emplace_back(BRepBuilderAPI_MakeFace{wire});

        wire = BRepBuilderAPI_MakeWire{innerEdges[3], radialEdges[2],
                                       outerEdges[9], radialEdges[7]};
        wire.Add(radialEdges[11]);
        wire.Add(radialEdges[5]);
        faces.emplace_back(BRepBuilderAPI_MakeFace{wire});

        faces.emplace_back(BRepBuilderAPI_MakeFace{BRepBuilderAPI_MakeWire{
            innerEdges[4], radialEdges[2], outerEdges[10], radialEdges[3]}});
      }
      // plate
      {
        auto wire = BRepBuilderAPI_MakeWire{innerEdges[2], radialEdges[4],
                                            outerEdges[0], outerEdges[1]};
        wire.Add(outerEdges[2]);
        wire.Add(radialEdges[5]);
        faces.emplace_back(BRepBuilderAPI_MakeFace{wire});

        wire = BRepBuilderAPI_MakeWire{outerEdges[3], radialEdges[6],
                                       outerEdges[8], radialEdges[7]};
        wire.Add(outerEdges[5]);
        wire.Add(outerEdges[4]);
        faces.emplace_back(BRepBuilderAPI_MakeFace{wire});
      }
      // fuel meat
      {
        for (int i = 0; i < 3; ++i) {
          faces.emplace_back(BRepBuilderAPI_MakeFace{
              BRepBuilderAPI_MakeWire{outerEdges[i], radialEdges[i + 8],
                                      outerEdges[i + 3], radialEdges[i + 9]}});
        }
      }
      // Space after Space
    } else if (innerEdges.size() == 5 && outerEdges.size() == 5) {
      faces.reserve(5);
      faces.emplace_back(BRepBuilderAPI_MakeFace{BRepBuilderAPI_MakeWire{
          innerEdges[0], radialEdges[0], outerEdges[0], radialEdges[1]}});
      faces.emplace_back(BRepBuilderAPI_MakeFace{BRepBuilderAPI_MakeWire{
          innerEdges[1], radialEdges[1], outerEdges[1], radialEdges[2]}});
      faces.emplace_back(BRepBuilderAPI_MakeFace{BRepBuilderAPI_MakeWire{
          innerEdges[3], radialEdges[3], outerEdges[3], radialEdges[4]}});
      faces.emplace_back(BRepBuilderAPI_MakeFace{BRepBuilderAPI_MakeWire{
          innerEdges[4], radialEdges[4], outerEdges[4], radialEdges[5]}});
      // middle
      faces.emplace_back(BRepBuilderAPI_MakeFace{BRepBuilderAPI_MakeWire{
          innerEdges[2], radialEdges[2], outerEdges[2], radialEdges[3]}});
    }
  }
  return faces;
}

void getMaterialAndSideset(const std::string &zone,
                           const std::map<std::string, std::string> &materials,
                           const std::map<std::string, std::string> &sideset,
                           std::string &mat, std::string &ssName) {
  // Find the material name for the zone
  auto find = materials.find(zone);
  if (find == materials.end())
    std::cerr << "Material for 'Side' not found" << std::endl;
  else
    mat = find->second;

  // Find the sideset name if it exits
  if (!sideset.empty()) {
    auto srch = sideset.find(zone);
    if (find != sideset.end()) {
      ssName = srch->second;
    }
  }
}

}  // namespace

PlateMeshOption::PlateMeshOption(const std::array<int, 2> &numElems)
    : meshingType(MeshingType::STRUCT), numSegments(numElems) {}

PlateMeshOption::PlateMeshOption(MeshingType type)
    : meshingType(type), numSegments() {}

FuelElement::FuelElement(double startRadius, double angle,
                         std::vector<Plate> plates,
                         const std::array<double, 3> &center)
    : ShapeBase(center),
      startRadius_(startRadius),
      angle_(angle),
      plates_(std::move(plates)) {}

double FuelElement::getStartRadius() const { return startRadius_; }

void FuelElement::setStartRadius(double radius) { this->startRadius_ = radius; }

double FuelElement::getAngle() const { return angle_; }

void FuelElement::setAngle(double angle) { this->angle_ = angle; }

const std::vector<Plate> &FuelElement::getPlates() const { return plates_; }

void FuelElement::setPlates(std::vector<Plate> plates) {
  this->plates_ = std::move(plates);
}

void FuelElement::addPlate(const Plate &plate) {
  this->plates_.emplace_back(plate);
}

template <typename PlateT, typename FVertCirc>
NEM::GEO::GeoManager drawNested(double startRadius, double angle,
                                const std::vector<PlateT> &plates,
                                const std::array<double, 3> &center,
                                FVertCirc funcVertsCircEdges) {
  // Initialize a GeoManager of dim 2
  NEM::GEO::GeoManager output{2};
  if (plates.empty()) {
    std::cerr << "plates is empty" << std::endl;
    return output;
  }

  // Make the center vertex
  TopoDS_Vertex centerVert{BRepBuilderAPI_MakeVertex{coords2Point(center)}};
  // Container of vertices in each plate
  std::vector<std::vector<TopoDS_Vertex>> vertices{};  // along the plates
  vertices.reserve(plates.size());
  //  Conatiner of circumferential edges in each plate
  std::vector<std::vector<TopoDS_Edge>> circumferentialEdges{};
  int plateCount = 0;
  auto previousPlate = plates.begin();
  double currentRadius = startRadius;
  for (const auto &plate : plates) {
    if (plateCount != 0) {
      currentRadius += previousPlate->thickness;
      previousPlate++;
    }
    auto vertsAndEdges = funcVertsCircEdges(plate, currentRadius, plateCount);
    // Increment plate count
    plateCount++;
    vertices.emplace_back(std::move(vertsAndEdges.first));
    circumferentialEdges.emplace_back(std::move(vertsAndEdges.second));
  }

  std::vector<std::vector<TopoDS_Edge>> radialEdges{};
  radialEdges.reserve(plates.size());
  std::vector<std::vector<TopoDS_Face>> faces{};
  faces.reserve(plates.size());
  {
    // First plate radial edges
    auto prevVerts = vertices.begin();
    radialEdges.emplace_back(getFirstRadialEdges(*prevVerts));
    auto prevCircEdges = circumferentialEdges.begin();
    faces.emplace_back(getFirstFaces(*prevCircEdges, radialEdges.back()));

    // Get subsequent plate radial edges
    auto iterVerts = plates.size() >= 2 ? prevVerts + 1 : vertices.end();
    auto iterCircEdges =
        plates.size() >= 2 ? prevCircEdges + 1 : circumferentialEdges.end();
    int count = 1;
    for (; iterVerts != vertices.end();
         ++prevVerts, ++iterVerts, ++prevCircEdges, ++iterCircEdges, ++count) {
      radialEdges.emplace_back(getRadialEdges(*prevVerts, *iterVerts, count));
      faces.emplace_back(
          getFaces(*prevCircEdges, *iterCircEdges, radialEdges.back(), count));
    }
  }

  // Start adding mesh/material properties
  // For the first plate
  {
    // Current plate meshing params map
    auto &meshingParams = plates[0].meshType_map;
    auto &materials = plates[0].materials_map;
    auto &sideset = plates[0].sideset_map;
    for (const auto &param : meshingParams) {
      // Get the meshing type
      auto meshingType = param.second.meshingType;
      if (param.first == "End") {
        std::string mat{};
        std::string ssName{};
        getMaterialAndSideset(param.first, materials, sideset, mat, ssName);
        switch (meshingType) {
          case PlateMeshOption::MeshingType::STRUCT:
          case PlateMeshOption::MeshingType::QUAD:
            output.insertConstruct<QuadMeshSurface>(faces[0].at(0), mat);
            output.insertConstruct<QuadMeshSurface>(faces[0].at(3), mat);
            break;
          case PlateMeshOption::MeshingType::TRI:
          default:
            output.insertConstruct<TriMeshSurface>(faces[0].at(0), mat);
            output.insertConstruct<TriMeshSurface>(faces[0].at(3), mat);
        }
        if (meshingType == PlateMeshOption::MeshingType::STRUCT) {
          output.insertConstruct<EdgeSegments>(circumferentialEdges[0].at(0),
                                               ssName,
                                               param.second.numSegments[1]);
          output.insertConstruct<EdgeSegments>(
              circumferentialEdges[0].at(plates[0].fuelThickness > 0 ? 11 : 5),
              ssName, param.second.numSegments[1]);
          output.insertConstruct<EdgeSegments>(circumferentialEdges[0].at(4),
                                               ssName,
                                               param.second.numSegments[1]);
          output.insertConstruct<EdgeSegments>(
              circumferentialEdges[0].at(plates[0].fuelThickness > 0 ? 15 : 9),
              ssName, param.second.numSegments[1]);
          if (plates[0].fuelThickness > 0) {
            for (int i = 0; i < 4; ++i) {
              output.insertConstruct<EdgeSegments>(radialEdges[0].at(i), "",
                                                   param.second.numSegments[0]);
            }
          } else {
            for (int i = 0; i < 6; ++i) {
              if (i == 2 || i == 3)
                continue;
              else
                output.insertConstruct<EdgeSegments>(
                    radialEdges[0].at(i), "", param.second.numSegments[0]);
            }
          }
        } else {
          if (!ssName.empty()) {
            output.insertConstruct<SideSetEdge>(circumferentialEdges[0].at(0),
                                                ssName);
            output.insertConstruct<SideSetEdge>(
                circumferentialEdges[0].at(plates[0].fuelThickness > 0 ? 11
                                                                       : 5),
                ssName);
            output.insertConstruct<SideSetEdge>(circumferentialEdges[0].at(4),
                                                ssName);
            output.insertConstruct<SideSetEdge>(
                circumferentialEdges[0].at(plates[0].fuelThickness > 0 ? 15
                                                                       : 9),
                ssName);
          }
        }
      } else if (param.first == "Side") {
        std::string mat{};
        std::string ssName{};
        getMaterialAndSideset(param.first, materials, sideset, mat, ssName);
        switch (meshingType) {
          case PlateMeshOption::MeshingType::STRUCT:
          case PlateMeshOption::MeshingType::QUAD:
            output.insertConstruct<QuadMeshSurfaceComposite>(faces[0].at(1),
                                                             mat);
            output.insertConstruct<QuadMeshSurfaceComposite>(faces[0].at(2),
                                                             mat);
            break;
          case PlateMeshOption::MeshingType::TRI:
          default:
            output.insertConstruct<TriMeshSurface>(faces[0].at(1), mat);
            output.insertConstruct<TriMeshSurface>(faces[0].at(2), mat);
        }
        if (meshingType == PlateMeshOption::MeshingType::STRUCT) {
          output.insertConstruct<EdgeSegments>(circumferentialEdges[0].at(1),
                                               ssName,
                                               param.second.numSegments[1]);
          output.insertConstruct<EdgeSegments>(
              circumferentialEdges[0].at(plates[0].fuelThickness > 0 ? 12 : 6),
              ssName, param.second.numSegments[1]);
          output.insertConstruct<EdgeSegments>(circumferentialEdges[0].at(3),
                                               ssName,
                                               param.second.numSegments[1]);
          output.insertConstruct<EdgeSegments>(
              circumferentialEdges[0].at(plates[0].fuelThickness > 0 ? 14 : 8),
              ssName, param.second.numSegments[1]);

          // Specifying only 2 radial edges even though these faces have 4
          // radial edges total. Leaving the others to be defined by the
          // 'Plate' and 'Fuel' mesh options.
          output.insertConstruct<EdgeSegments>(radialEdges[0].at(1), "",
                                               param.second.numSegments[0]);
          output.insertConstruct<EdgeSegments>(
              radialEdges[0].at(plates[0].fuelThickness > 0 ? 2 : 4), "",
              param.second.numSegments[0]);
        } else {
          if (!ssName.empty()) {
            output.insertConstruct<SideSetEdge>(circumferentialEdges[0].at(1),
                                                ssName);
            output.insertConstruct<SideSetEdge>(
                circumferentialEdges[0].at(plates[0].fuelThickness > 0 ? 12
                                                                       : 6),
                ssName);
            output.insertConstruct<SideSetEdge>(circumferentialEdges[0].at(2),
                                                ssName);
            output.insertConstruct<SideSetEdge>(
                circumferentialEdges[0].at(plates[0].fuelThickness > 0 ? 14
                                                                       : 8),
                ssName);
          }
        }
      } else if (param.first == "Plate") {
        std::string mat{};
        std::string ssName{};
        getMaterialAndSideset(param.first, materials, sideset, mat, ssName);
        // Also get mesh params for Fuel
        int fuelCircSegments = 0;
        int plateSideSegments = 1;
        auto srch2 = meshingParams.find("Fuel");
        if (srch2 == meshingParams.end() && plates[0].fuelThickness > 0)
          std::cerr << "Mesh params for 'Fuel' not found" << std::endl;
        else {
          fuelCircSegments = srch2->second.numSegments[1];
          if (param.second.numSegments[1] > fuelCircSegments) {
            plateSideSegments =
                (param.second.numSegments[1] - fuelCircSegments) / 2;
          }
        }
        switch (meshingType) {
          case PlateMeshOption::MeshingType::STRUCT:
          case PlateMeshOption::MeshingType::QUAD:
            if (plates[0].fuelThickness > 0) {
              output.insertConstruct<QuadMeshSurfaceComposite>(faces[0].at(4),
                                                               mat);
              output.insertConstruct<QuadMeshSurfaceComposite>(faces[0].at(5),
                                                               mat);
              output.insertConstruct<QuadMeshSurfaceComposite>(faces[0].at(6),
                                                               mat);
              output.insertConstruct<QuadMeshSurfaceComposite>(faces[0].at(8),
                                                               mat);
            } else
              output.insertConstruct<QuadMeshSurfaceComposite>(faces[0].at(4),
                                                               mat);
            break;
          case PlateMeshOption::MeshingType::TRI:
          default:
            if (plates[0].fuelThickness > 0) {
              output.insertConstruct<TriMeshSurface>(faces[0].at(4), mat);
              output.insertConstruct<TriMeshSurface>(faces[0].at(5), mat);
              output.insertConstruct<TriMeshSurface>(faces[0].at(6), mat);
              output.insertConstruct<TriMeshSurface>(faces[0].at(8), mat);
            } else
              output.insertConstruct<TriMeshSurface>(faces[0].at(4), mat);
        }
        if (meshingType == PlateMeshOption::MeshingType::STRUCT) {
          // Specifying only 1 circum edge even though these faces have 4
          // circum edges total. Leaving the others to be defined by the
          // 'Fuel' mesh options.
          output.insertConstruct<EdgeSegmentsComposite>(
              circumferentialEdges[0].at(2), ssName,
              param.second.numSegments[1]);
          output.insertConstruct<EdgeSegmentsComposite>(
              circumferentialEdges[0].at(plates[0].fuelThickness > 0 ? 13 : 7),
              ssName, param.second.numSegments[1]);

          if (plates[0].fuelThickness > 0) {
            output.insertConstruct<EdgeSegmentsComposite>(
                circumferentialEdges[0].at(5), ssName, plateSideSegments);
            output.insertConstruct<EdgeSegmentsComposite>(
                circumferentialEdges[0].at(7), ssName, plateSideSegments);
            output.insertConstruct<EdgeSegmentsComposite>(
                circumferentialEdges[0].at(8), ssName, plateSideSegments);
            output.insertConstruct<EdgeSegmentsComposite>(
                circumferentialEdges[0].at(10), ssName, plateSideSegments);

            output.insertConstruct<EdgeSegments>(radialEdges[0].at(4), "",
                                                 param.second.numSegments[0]);
            output.insertConstruct<EdgeSegments>(radialEdges[0].at(5), "",
                                                 param.second.numSegments[0]);
            output.insertConstruct<EdgeSegments>(radialEdges[0].at(6), "",
                                                 param.second.numSegments[0]);
            output.insertConstruct<EdgeSegments>(radialEdges[0].at(7), "",
                                                 param.second.numSegments[0]);
          } else {
            output.insertConstruct<EdgeSegments>(radialEdges[0].at(2), "",
                                                 param.second.numSegments[0]);
            output.insertConstruct<EdgeSegments>(radialEdges[0].at(3), "",
                                                 param.second.numSegments[0]);
          }
        } else {
          if (!ssName.empty()) {
            output.insertConstruct<SideSetEdge>(circumferentialEdges[0].at(2),
                                                ssName);
            output.insertConstruct<SideSetEdge>(
                circumferentialEdges[0].at(plates[0].fuelThickness > 0 ? 13
                                                                       : 7),
                ssName);
          }
        }
      }
      if (param.first == "Fuel" && plates[0].fuelThickness > 0) {
        std::string mat{};
        std::string ssName{};
        getMaterialAndSideset(param.first, materials, sideset, mat, ssName);
        switch (meshingType) {
          case PlateMeshOption::MeshingType::STRUCT:
          case PlateMeshOption::MeshingType::QUAD:
            output.insertConstruct<QuadMeshSurface>(faces[0].at(7), mat);
            break;
          case PlateMeshOption::MeshingType::TRI:
          default: output.insertConstruct<TriMeshSurface>(faces[0].at(7), mat);
        }
        if (meshingType == PlateMeshOption::MeshingType::STRUCT) {
          // Specifying only 1 circum edge even though these faces have 4
          // circum edges total. Leaving the others to be defined by the
          // 'Fuel' mesh options.
          output.insertConstruct<EdgeSegments>(circumferentialEdges[0].at(6),
                                               ssName,
                                               param.second.numSegments[1]);
          output.insertConstruct<EdgeSegments>(circumferentialEdges[0].at(9),
                                               ssName,
                                               param.second.numSegments[1]);

          output.insertConstruct<EdgeSegments>(radialEdges[0].at(8), "",
                                               param.second.numSegments[0]);
          output.insertConstruct<EdgeSegments>(radialEdges[0].at(9), "",
                                               param.second.numSegments[0]);
          output.insertConstruct<EdgeSegments>(radialEdges[0].at(10), "",
                                               param.second.numSegments[0]);
          output.insertConstruct<EdgeSegments>(radialEdges[0].at(11), "",
                                               param.second.numSegments[0]);
        } else {
          if (!ssName.empty()) {
            output.insertConstruct<SideSetEdge>(circumferentialEdges[0].at(6),
                                                ssName);
            output.insertConstruct<SideSetEdge>(circumferentialEdges[0].at(9),
                                                ssName);
          }
        }
      }
    }
  }

  // For plate 2 and more
  auto prevCircEdges = circumferentialEdges.begin();
  auto iterCircEdges =
      plates.size() >= 2 ? prevCircEdges + 1 : circumferentialEdges.end();
  auto iterRadialEdges = radialEdges.begin() + 1;
  auto iterFaces = faces.begin() + 1;
  auto prevPlate = plates.begin();
  auto iterPlate = plates.size() >= 2 ? prevPlate + 1 : plates.end();
  for (; iterPlate != plates.end(); prevCircEdges = iterCircEdges,
                                    ++iterCircEdges, ++iterRadialEdges,
                                    ++iterFaces, ++iterPlate, ++prevPlate) {
    // Current plate meshing params map
    auto &meshingParams = iterPlate->meshType_map;
    auto &materials = iterPlate->materials_map;
    auto &sideset = iterPlate->sideset_map;
    // Previous plate meshing params map
    //auto &prevMeshingParams = prevPlate->meshType_map;

    for (const auto &param : meshingParams) {
      // Get the meshing type
      auto meshingType = param.second.meshingType;
      if (param.first == "End") {
        std::string mat{};
        std::string ssName{};
        getMaterialAndSideset(param.first, materials, sideset, mat, ssName);
        switch (meshingType) {
          case PlateMeshOption::MeshingType::STRUCT:
          case PlateMeshOption::MeshingType::QUAD:
            output.insertConstruct<QuadMeshSurfaceComposite>(iterFaces->at(0),
                                                             mat);
            output.insertConstruct<QuadMeshSurfaceComposite>(iterFaces->at(3),
                                                             mat);
            break;
          case PlateMeshOption::MeshingType::TRI:
          default:
            output.insertConstruct<TriMeshSurface>(iterFaces->at(0), mat);
            output.insertConstruct<TriMeshSurface>(iterFaces->at(3), mat);
        }
        if (meshingType == PlateMeshOption::MeshingType::STRUCT) {
          output.insertConstruct<EdgeSegments>(
              iterCircEdges->at(iterPlate->fuelThickness > 0 ? 6 : 0), ssName,
              param.second.numSegments[1]);
          output.insertConstruct<EdgeSegments>(
              iterCircEdges->at(iterPlate->fuelThickness > 0 ? 10 : 4), ssName,
              param.second.numSegments[1]);
          if (iterPlate->fuelThickness > 0) {
            for (int i = 0; i < 4; ++i) {
              output.insertConstruct<EdgeSegments>(iterRadialEdges->at(i), "",
                                                   param.second.numSegments[0]);
            }
          } else {
            for (int i = 0; i < 6; ++i) {
              if (i == 2 || i == 3)
                continue;
              else
                output.insertConstruct<EdgeSegments>(
                    iterRadialEdges->at(i), "", param.second.numSegments[0]);
            }
          }
        } else {
          if (!ssName.empty()) {
            output.insertConstruct<SideSetEdge>(
                iterCircEdges->at(iterPlate->fuelThickness > 0 ? 6 : 0),
                ssName);
            output.insertConstruct<SideSetEdge>(
                iterCircEdges->at(iterPlate->fuelThickness > 0 ? 10 : 4),
                ssName);
          }
        }
      } else if (param.first == "Side") {
        std::string mat{};
        std::string ssName{};
        getMaterialAndSideset(param.first, materials, sideset, mat, ssName);
        switch (meshingType) {
          case PlateMeshOption::MeshingType::STRUCT:
          case PlateMeshOption::MeshingType::QUAD:
            output.insertConstruct<QuadMeshSurfaceComposite>(iterFaces->at(1),
                                                             mat);
            output.insertConstruct<QuadMeshSurfaceComposite>(iterFaces->at(2),
                                                             mat);
            break;
          case PlateMeshOption::MeshingType::TRI:
          default:
            output.insertConstruct<TriMeshSurface>(iterFaces->at(1), mat);
            output.insertConstruct<TriMeshSurface>(iterFaces->at(2), mat);
        }
        if (meshingType == PlateMeshOption::MeshingType::STRUCT) {
          output.insertConstruct<EdgeSegments>(
              iterCircEdges->at(iterPlate->fuelThickness > 0 ? 7 : 1), ssName,
              param.second.numSegments[1]);
          output.insertConstruct<EdgeSegments>(
              iterCircEdges->at(iterPlate->fuelThickness > 0 ? 9 : 3), ssName,
              param.second.numSegments[1]);

          // Specifying only 2 radial edges even though these faces have 4
          // radial edges total. Leaving the others to be defined by the
          // 'Plate' and 'Fuel' mesh options.
          output.insertConstruct<EdgeSegments>(iterRadialEdges->at(1), "",
                                               param.second.numSegments[0]);
          output.insertConstruct<EdgeSegments>(
              iterRadialEdges->at(iterPlate->fuelThickness > 0 ? 2 : 4), "",
              param.second.numSegments[0]);
        } else {
          if (!ssName.empty()) {
            output.insertConstruct<SideSetEdge>(
                iterCircEdges->at(iterPlate->fuelThickness > 0 ? 7 : 1),
                ssName);
            output.insertConstruct<SideSetEdge>(
                iterCircEdges->at(iterPlate->fuelThickness > 0 ? 9 : 3),
                ssName);
          }
        }
      } else if (param.first == "Plate") {
        std::string mat{};
        std::string ssName{};
        getMaterialAndSideset(param.first, materials, sideset, mat, ssName);
        // Also get mesh params for Fuel
        int fuelCircSegments = 0;
        int plateSideSegments = 1;
        auto srch2 = meshingParams.find("Fuel");
        if (srch2 == meshingParams.end() && iterPlate->fuelThickness > 0)
          std::cerr << "Mesh params for 'Fuel' not found" << std::endl;
        else {
          fuelCircSegments = srch2->second.numSegments[1];
          if (param.second.numSegments[1] > fuelCircSegments) {
            plateSideSegments =
                (param.second.numSegments[1] - fuelCircSegments) / 2;
          }
        }
        switch (meshingType) {
          case PlateMeshOption::MeshingType::STRUCT:
          case PlateMeshOption::MeshingType::QUAD:
            if (iterPlate->fuelThickness > 0) {
              output.insertConstruct<QuadMeshSurfaceComposite>(iterFaces->at(4),
                                                               mat);
              output.insertConstruct<QuadMeshSurfaceComposite>(iterFaces->at(5),
                                                               mat);
              output.insertConstruct<QuadMeshSurfaceComposite>(iterFaces->at(6),
                                                               mat);
              output.insertConstruct<QuadMeshSurfaceComposite>(iterFaces->at(8),
                                                               mat);
            } else
              output.insertConstruct<QuadMeshSurfaceComposite>(iterFaces->at(4),
                                                               mat);
            break;
          case PlateMeshOption::MeshingType::TRI:
          default:
            if (iterPlate->fuelThickness > 0) {
              output.insertConstruct<TriMeshSurface>(iterFaces->at(4), mat);
              output.insertConstruct<TriMeshSurface>(iterFaces->at(5), mat);
              output.insertConstruct<TriMeshSurface>(iterFaces->at(6), mat);
              output.insertConstruct<TriMeshSurface>(iterFaces->at(8), mat);
            } else
              output.insertConstruct<TriMeshSurface>(iterFaces->at(4), mat);
        }
        if (meshingType == PlateMeshOption::MeshingType::STRUCT) {
          // Specifying only 1 circum edge even though these faces have 4
          // circum edges total. Leaving the others to be defined by the
          // 'Fuel' mesh options.
          output.insertConstruct<EdgeSegmentsDistribution>(
              iterCircEdges->at(iterPlate->fuelThickness > 0 ? 8 : 2), ssName,
              param.second.numSegments[1]);

          if (iterPlate->fuelThickness > 0) {
            output.insertConstruct<EdgeSegmentsDistribution>(
                iterCircEdges->at(0), ssName, plateSideSegments);
            output.insertConstruct<EdgeSegmentsDistribution>(
                iterCircEdges->at(2), ssName, plateSideSegments);
            output.insertConstruct<EdgeSegmentsDistribution>(
                iterCircEdges->at(3), ssName, plateSideSegments);
            output.insertConstruct<EdgeSegmentsDistribution>(
                iterCircEdges->at(5), ssName, plateSideSegments);

            for (int i = 4; i < 8; ++i) {
              output.insertConstruct<EdgeSegmentsComposite>(
                  iterRadialEdges->at(i), "", param.second.numSegments[0]);
            }
          } else {
            output.insertConstruct<EdgeSegments>(iterRadialEdges->at(2), "",
                                                 param.second.numSegments[0]);
            output.insertConstruct<EdgeSegments>(iterRadialEdges->at(3), "",
                                                 param.second.numSegments[0]);
          }
        } else {
          if (!ssName.empty()) {
            output.insertConstruct<SideSetEdge>(
                iterCircEdges->at(iterPlate->fuelThickness > 0 ? 8 : 2),
                ssName);
          }
        }
      }
      if (iterPlate->fuelThickness > 0) {
        if (param.first == "Fuel") {
          std::string mat{};
          std::string ssName{};
          getMaterialAndSideset(param.first, materials, sideset, mat, ssName);
          switch (meshingType) {
            case PlateMeshOption::MeshingType::STRUCT:
            case PlateMeshOption::MeshingType::QUAD:
              output.insertConstruct<QuadMeshSurfaceComposite>(iterFaces->at(7),
                                                               mat);
              break;
            case PlateMeshOption::MeshingType::TRI:
            default:
              output.insertConstruct<TriMeshSurface>(iterFaces->at(7), mat);
          }
          if (meshingType == PlateMeshOption::MeshingType::STRUCT) {
            // Specifying only 1 circum edge even though these faces have 4
            // circum edges total. Leaving the others to be defined by the
            // 'Fuel' mesh options.
            output.insertConstruct<EdgeSegmentsComposite>(
                iterCircEdges->at(1), ssName, param.second.numSegments[1]);
            output.insertConstruct<EdgeSegmentsComposite>(
                iterCircEdges->at(4), ssName, param.second.numSegments[1]);

            output.insertConstruct<EdgeSegmentsComposite>(
                iterRadialEdges->at(8), "", param.second.numSegments[0]);
            output.insertConstruct<EdgeSegmentsComposite>(
                iterRadialEdges->at(9), "", param.second.numSegments[0]);
            output.insertConstruct<EdgeSegmentsComposite>(
                iterRadialEdges->at(10), "", param.second.numSegments[0]);
            output.insertConstruct<EdgeSegmentsComposite>(
                iterRadialEdges->at(11), "", param.second.numSegments[0]);
          }
        }
      }
    }
  }

  /*
  STEPControl_Writer writer;
  STEPControl_StepModelType mode = STEPControl_AsIs;
  // mode = STEPControl_BrepWithVoids;
  IFSelect_ReturnStatus stat;
  for (auto &i : vertices) {
    for (auto &vert : i) {
      stat = writer.Transfer(vert, mode);
      if (!stat) {
        std::cerr << " - Error" << std::endl;
      }
    }
  }
  for (auto &i : circumferentialEdges) {
    for (auto &edge : i) {
      std::cout << "writing edge" << std::endl;
      stat = writer.Transfer(edge, mode);
      if (!stat) {
        std::cerr << " - Error" << std::endl;
      }
    }
  }
  for (auto &i : radialEdges) {
    for (auto &edge : i) {
      std::cout << "writing edge" << std::endl;
      stat = writer.Transfer(edge, mode);
      if (!stat) {
        std::cerr << " - Error" << std::endl;
      }
    }
  }
  */

  /*
  STEPControl_Writer writer;
  STEPControl_StepModelType mode = STEPControl_AsIs;
  // mode = STEPControl_BrepWithVoids;
  IFSelect_ReturnStatus stat;
  for (auto &i : faces) {
    for (auto &face : i) {
      std::cout << "writing edge" << std::endl;
      stat = writer.Transfer(face, mode);
      if (!stat) {
        std::cerr << " - Error" << std::endl;
      }
    }
  }
  stat = writer.Write("FuelElement_Test.step");
  if (stat) {
    std::cout << " - File FuelElement_Test.step written successfully."
              << std::endl;
    // return output;
  } else {
    std::cout << " - File  FuelElement_Test.step writing failed." << std::endl;
    // return output;
  }
  */

  return output;
}

NEM::GEO::GeoManager FuelElement::createGeo() const {
  return drawNested(
      startRadius_, angle_, plates_, this->getCenter(),
      [this](const Plate &plate, double currentRadius, int plateIndex) {
        return getPlateVerticesAndCircEdges(
            this->getCenter(), currentRadius, this->angle_, plateIndex,
            plate.fuelThickness, plate.thickness, plate.endWidth,
            plate.sideWidth, plate.fuelOffset);
      });
}

}  // namespace NUCMESH
}  // namespace NEM
