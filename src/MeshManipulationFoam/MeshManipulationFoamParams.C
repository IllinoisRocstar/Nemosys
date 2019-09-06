#include "MeshManipulationFoamParams.H"


MeshManipulationFoamParams::MeshManipulationFoamParams()
{
  // Default Parameters
  // SurfLambdaMuSmooth
  lambda = 1;
  mu = 1;
  slmsIterations = 50;

  // splitMeshRegions
  _overwriteMsh = true;
  _cellZones = true;


  // mergeMeshes
  _overwriteMergeMsh = true;


  // createPatch
  _overwritecpMsh = true;
}