#include "NucMesh/shape.H"

namespace NEM {
namespace GEO {

void shape::updateSurfaces(const std::vector<std::pair<int, int>> &oldNew_vec) {
  for (auto &&_surface : _surfaces) {
    // get the current surface id
    int currID = _surface;

    // loop through vector of old/new id pairs
    for (const auto &j : oldNew_vec)
      // if current id is equal to old id, change to new id
      if (currID == j.first) _surface = j.second;
  }
}

}  // namespace GEO
}  // namespace NEM
