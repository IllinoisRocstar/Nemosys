#include "gmshParams.H"

#include <iostream>

namespace NEM {

namespace GEN {

gmshParams::gmshParams() {
  // default parameters
  ofname = "";

  minSize = 0.01;

  maxSize = 50.0;

  algo2D = "Frontal";

  algo3D = "Delaunay";

  extSizeFromBoundary = false;

  sizeFromCurvature = false;

  minElePer2Pi = 6;

  optimize = false;

  optimizeThreshold = 0.3;

  mSizeField = false;

  bgField = -1;

  meshExtensions = {".inp", ".unv", ".p3d", ".stl", ".vtk", ".su2"};
}

}  // namespace GEN

}  // namespace NEM
