#define _USE_MATH_DEFINES
#include "circles.H"

#include <gmsh.h>

#include <cmath>
#include <iostream>
#include <utility>

namespace NEM {
namespace GEO {

/*Creates either a single circle or concentric circles.
  Arguments:  center point,
              radii for each circle,
              mesh type for each circle
              number of elements if structured (radial, circum),
              region name,
              bools to make sidesets */

namespace occ = gmsh::model::occ;
namespace mesh = gmsh::model::mesh;

circles::circles(std::vector<double> cen, std::vector<double> rad,
                 std::vector<std::string> mtype,
                 std::vector<std::pair<int, int>> el,
                 std::vector<std::string> name, std::vector<std::string> ss) {
  _center = std::move(cen);      // Center of circle
  _radii = std::move(rad);       // Radii for concentric circles
  _meshType = std::move(mtype);  // Mesh type for each layer/circle
  _elems = std::move(el);        // Number of element for transfinite
  _names = std::move(name);      // Physical group names (Region names)
  _sideset = std::move(ss);      // Bools to make sideset
}

void circles::draw() {
  //----------- Create Points ------------//
  //--------------------------------------//

  // If _center vector has 5 args, use polar position
  if (_center.size() == 5) {
    double angle = _center[4] * M_PI / 180.0;
    _center[0] = _center[0] + _center[3] * std::cos(angle);
    _center[1] = _center[1] + _center[3] * std::sin(angle);
  }

  // automatically set first point tag
  int p = occ::addPoint(_center[0], _center[1], _center[2], 0, -1);

  // loop through radii to create other points
  int j = 1;
  for (const auto &i : _radii) {
    occ::addPoint(_center[0] - i, _center[1], _center[2], 0, p + j);
    j++;
    occ::addPoint(_center[0] + i, _center[1], _center[2], 0, p + j);
    j++;
  }
  // calling sync to add points
  occ::synchronize();

  //----------- Create Lines/Arcs ------------//
  //------------------------------------------//

  int numLines = 4 * _radii.size() - 2;  // number of lines
  int i = 1;                             // index of Line points
  int q = 1;                             // index of Circle points

  // Get max Line ID
  int l = getMaxID(1) + 1;

  for (int j = 0; j < numLines - 2; j = j + 2) {
    occ::addLine(p + i, p + i + 2, l + j);
    i++;

    if (j < numLines / 2) {
      occ::addCircleArc(p + q, p, p + q + 1, l + j + 1);
      q = q + 2;
    } else {
      q--;
      occ::addCircleArc(p + q, p, p + q - 1, l + j + 1);
      q--;
    }
  }
  j = numLines - 3;
  q--;

  if (numLines < 3) {
    j = 0;
    occ::addCircleArc(p + 1, p, p + 2, l + j);
    occ::addCircleArc(p + 2, p, p + 1, l + j + 1);
  } else {
    occ::addCircleArc(p + q, p, p + q - 1, l + j + 1);
    q = q - 2;
    j++;
    occ::addCircleArc(p + q, p, p + q - 1, l + j + 1);
  }

  // std::cout << "Lines" << std::endl;

  //----------- Create Lines Loops ------------//
  //-------------------------------------------//

  // Line Loop ID identifier, to be defined...
  int ll = -1;
  int nSurf = (_radii.size() - 1) * 2;  // number of surfaces
  std::vector<int> lineTags;

  // If only 2 lines, make a disk
  if (numLines < 3) {
    lineTags = {l + 0, l + 1};
    ll = occ::addCurveLoop(lineTags);
  } else {
    j = 0;
    int h = 3;

    // Top line loops
    for (int i = 0; i < nSurf / 2; i++) {
      lineTags = {l + j, l + h, l + j + 2, l + h - 2};

      if (i == 0) {
        // initialize the line loop tag
        ll = occ::addCurveLoop(lineTags);
        j = j + 4;
        h = h + 2;
      } else {
        occ::addCurveLoop(lineTags, ll + i);
        j = j + 4;
        h = h + 2;
      }
    }

    // Bottom line loops
    j = 0;
    h = numLines - 1;
    for (int i = nSurf / 2; i < nSurf; i++) {
      if (i < nSurf / 2 + 2) {
        // lineTags = {l + j, l + h - 1, l + j + 2, l + h};
        lineTags = {l + j, l + h, l + j + 2, l + h - 1};
        occ::addCurveLoop(lineTags, ll + i);
        j = j + 4;
        h--;
      } else {
        // lineTags = {l + j, l + h - 2, l + j + 2, l + h};
        lineTags = {l + j, l + h, l + j + 2, l + h - 2};
        occ::addCurveLoop(lineTags, ll + i);
        j = j + 4;
        h = h - 2;
      }
    }
    // Inner line loop
    i = nSurf;
    lineTags = {l + 1, l + numLines - 1};
    occ::addCurveLoop(lineTags, ll + i);
  }

  // std::cout << "Line Loops" << std::endl;

  //----------- Create Surfaces ------------//
  //----------------------------------------//

  // Get max surfaces ID
  int s = getMaxID(2) + 1;

  // container for circle surface ID's
  _surfaces.clear();

  // Inner surfaces
  i = 0;
  occ::addPlaneSurface({ll + i + nSurf}, s + i);
  _surfaces.push_back(s + i);

  // All other surfaces
  for (int i = 1; i < nSurf / 2 + 1; i++) {
    occ::addPlaneSurface({ll + i - 1}, s + i);
    _surfaces.push_back(s + i);
    occ::addPlaneSurface({ll + i + nSurf / 2 - 1}, s + i + nSurf / 2);
    _surfaces.push_back(s + i + nSurf / 2);
  }

  // std::cout << "Surface" << std::endl;
  //
  // for (int i = 0; i < surfaces.size(); i++)
  //   std::cout << "circ surf " << surfaces[i] << std::endl;

  // Sync to add surfaces
  occ::synchronize();
}

// Applies the mesh type (tri,quad,struct) to surfaces
void circles::applyMeshType() {
  // Gather all the lines of the circles for meshing
  std::vector<std::vector<std::pair<int, int>>> allLines;  // lines of circles
  allLines = getLines();

  // std::cout << " Gathered lines for circle surfaces" << std::endl;

  // Apply meshType: "Tri"=tris, "Quad"=quads, "Struct"=structured
  for (int i = 0; i < _meshType.size(); i++) {
    if (_meshType[i] == "Struct" || _meshType[i] == "S") {
      //------------Transfinite Mesh ----------//
      //---------------------------------------//
      // declare containers for meshing lines
      std::vector<int> circum, radial;

      for (int j = 0; j < allLines[i].size(); j = j + 2)
        radial.push_back(allLines[i][j].second);
      for (int j = 1; j < allLines[i].size(); j = j + 2)
        circum.push_back(allLines[i][j].second);

      for (const auto &j : radial)
        mesh::setTransfiniteCurve(j, _elems[i].first + 1);

      for (const auto &j : circum)
        mesh::setTransfiniteCurve(j, _elems[i].second + 1);

      int j = i * 2 - 1;
      mesh::setTransfiniteSurface(_surfaces[j]);
      mesh::setTransfiniteSurface(_surfaces[j + 1]);

      mesh::setRecombine(2, _surfaces[j]);
      mesh::setRecombine(2, _surfaces[j + 1]);
    } else if (_meshType[i] == "Quad" || _meshType[i] == "Q") {
      //---------------Quad Mesh---------------//
      //---------------------------------------//
      if (i == 0)
        mesh::setRecombine(2, _surfaces[i]);
      else {
        int j = i * 2 - 1;
        mesh::setRecombine(2, _surfaces[j]);
        mesh::setRecombine(2, _surfaces[j + 1]);
      }
    } else if (_meshType[i] == "Tri" || _meshType[i] == "T") {
      //---------------Tri Mesh---------------//
      //--------------------------------------//
      continue;
    } else {
      std::cout << "Warning: Mesh Type not recognized. Using Triangles."
                << std::endl;
    }
  }
}

std::vector<std::vector<std::pair<int, int>>> circles::getLines() const {
  // Gather all the lines of the circles for meshing
  std::vector<std::vector<std::pair<int, int>>> allLines;  // lines of circles
  std::vector<std::pair<int, int>> circleLines;

  std::vector<std::pair<int, int>> v = {{2, _surfaces[0]}};

  // Get the boundary (lines) of the _center surface
  gmsh::model::getBoundary(v, circleLines, false, false, false);
  allLines.push_back(circleLines);
  circleLines.clear();

  // Get the lines for all other surfaces
  for (int i = 1; i < _surfaces.size(); i = i + 2) {
    std::vector<std::pair<int, int>> v2 = {{2, _surfaces[i]},
                                           {2, _surfaces[i + 1]}};

    gmsh::model::getBoundary(v2, circleLines, false, false, false);
    allLines.push_back(circleLines);
    circleLines.clear();
  }
  return allLines;
}

void circles::getPhysSurf(const std::map<std::string, int> &phystag_map,
                          std::map<int, int> &physSurf_map) const {
  int physTag;
  if (!_names.empty()) {
    for (int i = 0; i < _names.size(); ++i) {
      auto it = phystag_map.find(_names[i]);
      if (it != phystag_map.end()) {
        physTag = it->second;
        if (i == 0)
          physSurf_map.insert(std::pair<int, int>(_surfaces[i], physTag));
        else {
          int j = i * 2 - 1;
          physSurf_map.insert(std::pair<int, int>(_surfaces[j], physTag));
          physSurf_map.insert(std::pair<int, int>(_surfaces[j + 1], physTag));
        }
      }
    }
  }
}

void circles::getPhysLine(const std::map<std::string, int> &sstag_map,
                          std::map<int, int> &physLine_map) const {
  std::vector<std::vector<std::pair<int, int>>> allLines;
  allLines = getLines();

  int ssTag;
  if (!_sideset.empty()) {
    for (int i = 0; i < _sideset.size(); ++i) {
      auto it = sstag_map.find(_sideset[i]);
      if (it != sstag_map.end()) {
        ssTag = it->second;

        // std::cout << "name = " << _sideset[i] << ", ssTag = " << ssTag <<
        // std::endl;

        if (i == 0 && _sideset[i] != "None") {
          physLine_map.insert(
              std::pair<int, int>(allLines[i][0].second, ssTag));
          physLine_map.insert(
              std::pair<int, int>(allLines[i][1].second, ssTag));
          // std::cout << "lines " << allLines[i][0].second << ", " <<
          // allLines[i][1].second << std::endl;
        }
        if (i > 0 && _sideset[i] != "None") {
          physLine_map.insert(
              std::pair<int, int>(allLines[i][1].second, ssTag));
          physLine_map.insert(
              std::pair<int, int>(allLines[i][5].second, ssTag));
          // std::cout << "line " << allLines[i][1].second << std::endl;
          // std::cout << "line " << allLines[i][5].second << std::endl;
        }
      }
    }
  }
}

// Returns the max ID/Tag of entity with dimension dim
int circles::getMaxID(int dim) const {
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
