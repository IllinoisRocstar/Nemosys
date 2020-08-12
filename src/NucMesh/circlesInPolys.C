#define _USE_MATH_DEFINES
#include "circlesInPolys.H"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <utility>

#include <gmsh.h>

namespace NEM {
namespace GEO {

/*Creates concentric polygons of n sides.
  Arguments:  number of sides,
              center point,
              radii of points,
              mesh type,
              number of elements,
              region name,
              rotation angle (degrees) */

namespace occ = gmsh::model::occ;
namespace mesh = gmsh::model::mesh;

circlesInPolys::circlesInPolys(int nsides, std::vector<double> cen,
                               std::vector<double> circle_rad,
                               std::vector<double> poly_rad,
                               std::vector<std::string> mtype,
                               std::vector<std::pair<int, int>> el,
                               std::vector<std::string> name, double rot) {
  _nSides = nsides;                       // Number of sides
  _center = std::move(cen);               // Center of circle
  _circle_radii = std::move(circle_rad);  // Radii for concentric circles
  _poly_radii = std::move(poly_rad);      // Radii for concentric polygon
  _meshType = std::move(mtype);           // Mesh type for each layer/circle
  _elems = std::move(el);                 // Number of element for transfinite
  _names = std::move(name);               // Physical group names (Region names)
  _rotation = rot;                        // Rotation angle (degrees)
}

// Draws the polygon
void circlesInPolys::draw() {
  //----------- Create Points ------------//
  //--------------------------------------//
  double incr = 360.0 / _nSides;                  // angle increment
  double defaultRot = -1.0 * (180 - incr) / 2.0;  // default rotation
  int total = _circle_radii.size() + _poly_radii.size();

  // If _center vector has 5 args, use polar position
  if (_center.size() == 5) {
    double angle = _center[4] * M_PI / 180.0;
    _center[0] = _center[0] + _center[3] * std::cos(angle);
    _center[1] = _center[1] + _center[3] * std::sin(angle);
  }

  int p = getMaxID(0) + 1;
  // Create center point
  occ::addPoint(_center[0], _center[1], _center[2], 0, p);

  // Loop through circle radii first and create points
  for (int i = 0; i < _circle_radii.size(); ++i) {
    int k = i * _nSides;
    for (int j = 0; j < _nSides; ++j) {
      double theta = (j * incr + _rotation + defaultRot) * M_PI / 180.0;
      double a = _circle_radii[i] * std::sin(theta);
      double b = _circle_radii[i] * std::cos(theta);
      occ::addPoint(_center[0] + b, _center[1] + a, _center[2], 0,
                    p + j + k + 1);
    }
  }

  int p1 = p + _nSides * _circle_radii.size() + 1;
  // Then loop through polygon radii
  for (int i = 0; i < _poly_radii.size(); ++i) {
    int k = i * _nSides;
    for (int j = 0; j < _nSides; ++j) {
      double theta = (j * incr + _rotation + defaultRot) * M_PI / 180.0;
      double a = _poly_radii[i] * std::sin(theta);
      double b = _poly_radii[i] * std::cos(theta);
      occ::addPoint(_center[0] + b, _center[1] + a, _center[2], 0, p1 + j + k);
    }
  }
  // calling sync to add points
  occ::synchronize();

  //----------- Create Lines/Arcs ------------//
  //------------------------------------------//
  // Get max line ID
  int l = getMaxID(1) + 1;

  // Make cicle arcs
  for (int i = 0; i < _circle_radii.size(); ++i) {
    int k = i * _nSides;
    for (int j = 0; j < _nSides - 1; ++j) {
      occ::addCircleArc(p + j + k + 1, p, p + k + j + 2, l + k + j);
    }
    occ::addCircleArc(p + _nSides + k, p, p + k + 1, l + _nSides - 1 + k);
  }

  // Make poly lines
  int q = _nSides * _circle_radii.size();
  for (int i = 0; i < _poly_radii.size(); ++i) {
    int k = i * _nSides;
    for (int j = 0; j < _nSides - 1; ++j)
      occ::addLine(q + p + j + k + 1, q + p + k + j + 2, q + l + k + j);
    occ::addLine(q + p + _nSides + k, q + p + k + 1, q + l + _nSides - 1 + k);
  }

  // Make radial lines
  for (int i = 0; i < total - 1; ++i) {
    int k = i * _nSides;
    int q = _nSides * (total + i) + 1;
    for (int j = 0; j < _nSides; ++j)
      occ::addLine(p + j + k + 1, p + k + j + _nSides + 1, l + q + j - 1);
  }
  occ::synchronize();

  //----------- Create Lines Loops ------------//
  //-------------------------------------------//

  // Line Loop ID, to be defined...
  int ll;
  int n = _nSides;              // Number of sides
  int nSurf = _nSides * total;  // number of surfaces

  std::vector<int> lineTags;  // vector for line tags
  std::vector<int> loopTags;  // vector for line loop tags

  // Get the first/inner loop
  lineTags.reserve(n);
  for (int i = 0; i < n; ++i) lineTags.push_back(l + i);

  std::reverse(lineTags.begin(), lineTags.end());
  // Automatically set the loop tag
  ll = occ::addCurveLoop(lineTags);
  loopTags.push_back(ll);

  // Get all other loops
  for (int i = 0; i < total - 1; ++i) {
    int k = i * n;
    for (int j = 0; j < n - 1; ++j) {
      lineTags = {l + j + k, l + j + nSurf + k, l + j + n + k,
                  l + j + nSurf + 1 + k};
      ll = occ::addCurveLoop(lineTags);
      loopTags.push_back(ll);
    }
    lineTags = {l + n + k - 1, l + n + nSurf + k - 1, l + n + n + k - 1,
                l + n + nSurf - n + k};
    ll = occ::addCurveLoop(lineTags);
    loopTags.push_back(ll);
  }

  //----------- Create Surfaces ------------//
  //----------------------------------------//

  // Get max surfaces ID
  int s = getMaxID(2) + 1;

  // container for circlesInPoly surface ID's
  _surfaces.clear();

  // Inner surfaces
  int i = 0;
  occ::addPlaneSurface({loopTags[i]}, s + i);
  _surfaces.push_back(s + i);

  // All other surfaces
  nSurf = (total - 1) * _nSides;
  for (int i = 1; i < nSurf + 1; i++) {
    occ::addPlaneSurface({loopTags[i]}, s + i);
    _surfaces.push_back(s + i);
  }

  // std::cout << "Surface" << std::endl;
  //
  // for (int i = 0; i < surfaces.size(); i++)
  //   std::cout << "poly surf " << surfaces[i] << std::endl;

  // Sync to add surfaces
  occ::synchronize();
}

// Applies the mesh type (tri,quad,struct) to surfaces
void circlesInPolys::applyMeshType() {
  // Gather all the lines of the polys for meshing
  std::vector<std::vector<std::pair<int, int>>> allLines;  // lines of circles
  std::vector<std::pair<int, int>> circleLines;

  std::vector<std::pair<int, int>> v = {{2, _surfaces[0]}};

  // Get the boundary (lines) of the _center surface
  gmsh::model::getBoundary(v, circleLines, true, false, false);

  allLines.push_back(circleLines);
  circleLines.clear();

  // Get the lines for all other surfaces
  for (int i = 1; i < _surfaces.size(); i = i + _nSides) {
    std::vector<std::pair<int, int>> v2;
    for (int j = i; j < i + _nSides; ++j) v2.emplace_back(2, _surfaces[j]);

    gmsh::model::getBoundary(v2, circleLines, false, false, false);

    allLines.push_back(circleLines);
    circleLines.clear();
  }

  // Apply _meshType: "Tri"=tris, "Quad"=quads, "Struct"=structured
  for (int i = 0; i < _meshType.size(); i++) {
    int n = _nSides;
    //------------Transfinite Mesh ----------//
    //---------------------------------------//
    if (_meshType[i] == "Struct" || _meshType[i] == "S") {
      // declare containers for meshing lines
      std::vector<int> circum, radial;
      std::vector<std::pair<int, int>> out;

      for (const auto &j : allLines[i]) {
        std::vector<double> pts1, pts2;
        gmsh::model::getBoundary({j}, out, false, false, false);

        gmsh::model::getValue(out[0].first, out[0].second, {}, pts1);
        gmsh::model::getValue(out[1].first, out[1].second, {}, pts2);

        double a = pts1[0] - _center[0];
        double b = pts1[1] - _center[1];
        double c = pts1[2] - _center[2];
        double aa = pts2[0] - _center[0];
        double bb = pts2[1] - _center[1];
        double cc = pts2[2] - _center[2];

        double dist1 = sqrt(pow(a, 2) + pow(b, 2) + pow(c, 2));
        double dist2 = sqrt(pow(aa, 2) + pow(bb, 2) + pow(cc, 2));
        double diff = dist1 - dist2;
        double adiff = fabs(diff);
        double eps = 1e-9;

        if (adiff > eps && adiff < 1e-4) {
          std::cout << "Warning: Mesh type may not be applied correctly."
                    << std::endl;
          std::cout << "         Please use higher precision for shape Center "
                       "and Radii."
                    << std::endl;
        }

        if (adiff < eps)
          circum.push_back(j.second);
        else
          radial.push_back(j.second);
      }

      for (const auto &j : radial)
        mesh::setTransfiniteCurve(j, _elems[i].first + 1);

      for (const auto &j : circum)
        mesh::setTransfiniteCurve(j, _elems[i].second + 1);

      for (int j = n * (i - 1) + 1; j < n * i + 1; ++j) {
        mesh::setTransfiniteSurface(_surfaces[j]);
        mesh::setRecombine(2, _surfaces[j]);
      }

      circum.clear();
      radial.clear();
      out.clear();
    }

    //---------------Quad Mesh---------------//
    //---------------------------------------//
    else if (_meshType[i] == "Quad" || _meshType[i] == "Q") {
      if (i == 0)
        mesh::setRecombine(2, _surfaces[i]);
      else {
        for (int j = n * (i - 1) + 1; j < n * i + 1; ++j)
          mesh::setRecombine(2, _surfaces[j]);
      }
    }
    //---------------Tri Mesh---------------//
    //---------------------------------------//
    else if (_meshType[i] == "Tri" || _meshType[i] == "T") {
      continue;
    } else {
      std::cout << "Mesh Type not recognized. Using Triangles." << std::endl;
    }
  }
}

std::map<int, int> circlesInPolys::getPhysSurf(
    std::map<std::string, int> phystag_map, std::map<int, int> physSurf_map) {
  //------------Physical Groups------------//
  //---------------------------------------//
  int physTag;
  if (!_names.empty()) {
    for (int i = 0; i < _names.size(); ++i) {
      auto it = phystag_map.find(_names[i]);
      if (it != phystag_map.end()) {
        physTag = it->second;
        if (i == 0) {
          // Add new surface id to outSurfaces
          physSurf_map.insert(std::pair<int, int>(_surfaces[i], physTag));
        } else {
          for (int j = (i - 1) * _nSides + 1; j < (i)*_nSides + 1; ++j)
            physSurf_map.insert(std::pair<int, int>(_surfaces[j], physTag));
        }
      } else
        std::cout << "physical tag not in phystag_map" << std::endl;
    }
  }
  return physSurf_map;
  occ::synchronize();
}

// Returns the max ID/Tag of entity with dimension dim
int circlesInPolys::getMaxID(int dim) {
  std::vector<std::pair<int, int>> retTags;
  gmsh::model::getEntities(retTags, dim);
  int max = 0;
  int tag;
  if (!retTags.empty()) {
    for (const auto &retTag : retTags) {
      tag = retTag.second;
      if (tag > max)
        max = tag;
      else
        continue;
    }
    return max;
  } else {
    return 0;
  }
}

}  // namespace GEO
}  // namespace NEM
