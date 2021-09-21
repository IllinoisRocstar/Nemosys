#include <boost/filesystem.hpp>
#include <cstdlib>
#include <iostream>
#include <set>
#include <string>

#include "Mesh/geoMeshBase.H"
#include "Mesh/geoMeshFactory.H"
#include "MeshGeneration/blockMeshGen.H"
#include "MeshGeneration/blockMeshParams.H"

// OpenFOAM headers
#include <Time.H>
#include <argList.H>
#include <blockMsh.H>
#include <boundBox.H>
#include <fileName.H>
#include <fvMesh.H>
#include <getDicts.H>
#include <triSurf.H>
#include <triSurfModifier.H>
//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

constexpr double Pi = 3.14159265358979323846;
constexpr double tan30 = 0.577350269189625731;
constexpr double sin45 = 0.707106781186547462;

blockMeshGen::blockMeshGen()  // Default constructor body
{
  // booleans
  params_ = new blockMeshParams();
  defaults = true;

  // initializing the Foam environment
  initialize();
}

// Constructor with user define parameters
blockMeshGen::blockMeshGen(blockMeshParams *params)
    : defaults(false), params_(params) {
  // Initialize foam environment
  initialize();
}

blockMeshGen::~blockMeshGen() {}

void blockMeshGen::initialize() {
  // create dictionaries needed in memory
  bool writeDicts;
  if (params_->isPackMesh)
    writeDicts = true;
  else
    writeDicts = false;
  std::unique_ptr<getDicts> initFoam;
  initFoam = std::unique_ptr<getDicts>(new getDicts());
  controlDict_ = initFoam->createControlDict(writeDicts);
  fvSchemes_ = initFoam->createFvSchemes(writeDicts);
  fvSolution_ = initFoam->createFvSolution(writeDicts);

  // creates dictionaries in constant folder
  createBlockMshDict(writeDicts);

  runTime_ =
      std::unique_ptr<Foam::Time>(new Foam::Time(controlDict_.get(), ".", "."));

  //- 2d cartesian mesher cannot be run in parallel
  Foam::argList::noParallel();
}

int blockMeshGen::createMeshFromSTL(const char *fname) {
  auto mshObj = blockMsh();

  bool writeMsh;
  if (params_->isPackMesh)
    writeMsh = true;
  else
    writeMsh = false;

  fmesh_ = std::unique_ptr<Foam::fvMesh>(
      mshObj.generate(runTime_, blockMshDict_, writeMsh));

  // Create foamGeoMesh
  if (writeMsh) {
    gmData.reset();
    gmData = std::unique_ptr<NEM::MSH::geoMeshBase>(NEM::MSH::Read(".foam"));
  } else {
    auto fgm_ = std::unique_ptr<NEM::MSH::foamGeoMesh>(
        new NEM::MSH::foamGeoMesh(fmesh_.get(), ""));
    gmData = std::unique_ptr<NEM::MSH::geoMeshBase>(
        dynamic_cast<NEM::MSH::geoMeshBase *>(fgm_.get()));
  }

  return 0;
}

void blockMeshGen::createBlockMshDict(const bool &write) {
  // header
  std::string contText =
      "\
/*--------------------------------*- C++ -*----------------------------------*\n\
| =========                   |                                                |\n\
| \\\\      /  F ield         | NEMoSys: blockMesh interface                   |\n\
|  \\\\    /   O peration     |                                                |\n\
|   \\\\  /    A nd           |                                                |\n\
|    \\\\/     M anipulation  |                                                |\n\
\\*---------------------------------------------------------------------------*/\n\
\n\
FoamFile\n\
{\n\
    version   2.0;\n\
    format    ascii;\n\
    class     dictionary;\n\
    location  \"system\";\n\
    object    blockMeshDict;\n\
}\n\n";

  contText =
      contText +
      "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";

  if (auto box = std::dynamic_pointer_cast<bmBox>(params_->shape)) {
    double finalX;
    double finalY;
    double finalZ;

    // Automatically Generates Box Around Packs and Meshes that Hexahedrally.
    if (box->autoGenerate.has_value()) {
      const auto &autoGenerate = box->autoGenerate.value();
      Foam::fileName inFileName((autoGenerate.packFileName));

      Foam::Module::triSurf origSurface(inFileName);
      Foam::Module::triSurfModifier sMod(origSurface);
      Foam::pointField &points = sMod.pointsAccess();

      const Foam::boundBox bb(points);

      Foam::vector negOffset, posOffset;

      negOffset[0] = (autoGenerate.offset[0]);
      negOffset[1] = (autoGenerate.offset[1]);
      negOffset[2] = (autoGenerate.offset[2]);
      posOffset[0] = (autoGenerate.offset[0]);
      posOffset[1] = (autoGenerate.offset[1]);
      posOffset[2] = (autoGenerate.offset[2]);

      const Foam::boundBox newBB(bb.min() - negOffset, bb.max() + posOffset);

      (box->init[0]) = newBB.min()[0];
      (box->init[1]) = newBB.min()[1];
      (box->init[2]) = newBB.min()[2];

      finalX = newBB.max()[0];
      finalY = newBB.max()[1];
      finalZ = newBB.max()[2];

      std::array<double, 3> minPointsBox{newBB.min()[0], newBB.min()[1],
                                         newBB.min()[2]};
      std::array<double, 3> maxPointsBox{newBB.max()[0], newBB.max()[1],
                                         newBB.max()[2]};

      box->coordsBox = std::make_pair(minPointsBox, maxPointsBox);

      if (params_->cellSize.has_value()) {
        double xLength = std::sqrt(((box->init[0]) - (finalX)) *
                                   ((box->init[0]) - (finalX)));
        double yLength = std::sqrt(((box->init[1]) - (finalY)) *
                                   ((box->init[1]) - (finalY)));
        double zLength = std::sqrt(((box->init[2]) - (finalZ)) *
                                   ((box->init[2]) - (finalZ)));

        (params_->nCells[0]) =
            static_cast<int>(std::round(xLength / (params_->cellSize.value())));
        (params_->nCells[1]) =
            static_cast<int>(std::round(yLength / (params_->cellSize.value())));
        (params_->nCells[2]) =
            static_cast<int>(std::round(zLength / (params_->cellSize.value())));
      }
    } else {
      finalX = (box->init[0]) + (box->len[0]);
      finalY = (box->init[1]) + (box->len[1]);
      finalZ = (box->init[2]) + (box->len[2]);

      std::array<double, 3> minPointsBox{box->init[0], box->init[1],
                                         box->init[2]};
      std::array<double, 3> maxPointsBox{finalX, finalY, finalZ};
      box->coordsBox = std::make_pair(minPointsBox, maxPointsBox);
    }

    // Box data
    contText = contText + "convertToMeters " +
               std::to_string(params_->cnvrtToMeters) + ";\n";
    contText = contText + "\nvertices\n";
    contText = contText + "(\n\n";

    contText = contText + "\t(" + std::to_string(box->init[0]) + " " +
               std::to_string(box->init[1]) + " " +
               std::to_string(box->init[2]) + ")\n";

    contText = contText + "\t(" + std::to_string(finalX) + " " +
               std::to_string(box->init[1]) + " " +
               std::to_string(box->init[2]) + ")\n";

    contText = contText + "\t(" + std::to_string(finalX) + " " +
               std::to_string(finalY) + " " + std::to_string(box->init[2]) +
               ")\n";

    contText = contText + "\t(" + std::to_string(box->init[0]) + " " +
               std::to_string(finalY) + " " + std::to_string(box->init[2]) +
               ")\n";

    contText = contText + "\t(" + std::to_string(box->init[0]) + " " +
               std::to_string(box->init[1]) + " " + std::to_string(finalZ) +
               ")\n";

    contText = contText + "\t(" + std::to_string(finalX) + " " +
               std::to_string(box->init[1]) + " " + std::to_string(finalZ) +
               ")\n";

    contText = contText + "\t(" + std::to_string(finalX) + " " +
               std::to_string(finalY) + " " + std::to_string(finalZ) + ")\n";

    contText = contText + "\t(" + std::to_string(box->init[0]) + " " +
               std::to_string(finalY) + " " + std::to_string(finalZ) + ")\n";
    contText = contText + "\n);\n";

    contText = contText + "\nblocks\n(";
    contText = contText + "\n\thex (0 1 2 3 4 5 6 7) (" +
               std::to_string(params_->nCells[0]) + " " +
               std::to_string(params_->nCells[1]) + " " +
               std::to_string(params_->nCells[2]) + ") simpleGrading (" +
               std::to_string(box->smplGrading[0]) + " " +
               std::to_string(box->smplGrading[1]) + " " +
               std::to_string(box->smplGrading[2]) + ")\n\n";
    contText = contText + ");\n";

    contText = contText + "\nedges\n";
    contText = contText + "(\n\n);\n";

    contText = contText + "\npatches";
    contText = contText + "\n(\n";
    contText = contText + "\n\tpatch up\n";
    contText = contText + "\t(\n";
    contText = contText + "\t\t(3 7 6 2)\n";
    contText = contText + "\t)\n";

    contText = contText + "\tpatch down\n";
    contText = contText + "\t(\n";
    contText = contText + "\t\t(0 4 5 1)\n";
    contText = contText + "\t)\n";

    contText = contText + "\tpatch in\n";
    contText = contText + "\t(\n";
    contText = contText + "\t\t(3 7 4 0)\n";
    contText = contText + "\t)\n";

    contText = contText + "\tpatch out\n";
    contText = contText + "\t(\n";
    contText = contText + "\t\t(2 6 5 1)\n";
    contText = contText + "\t)\n";

    contText = contText + "\tpatch front\n";
    contText = contText + "\t(\n";
    contText = contText + "\t\t(0 1 2 3)\n";
    contText = contText + "\t)\n";

    contText = contText + "\tpatch back\n";
    contText = contText + "\t(\n";
    contText = contText + "\t\t(4 5 6 7)\n";
    contText = contText + "\t)\n";

    contText = contText + "\n);\n";

    contText = contText + "\nmergePatchPairs\n(\n);\n";
    contText = contText +
               "\n //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * "
               "* * *//";

  }

  else if (auto sphere = std::dynamic_pointer_cast<bmSphere>(params_->shape)) {
    // Put sphere implementation here
    contText = contText + "\ncx " + std::to_string(sphere->center[0]) + ";";
    contText = contText + "\ncy " + std::to_string(sphere->center[1]) + ";";
    contText = contText + "\ncz " + std::to_string(sphere->center[2]) + ";";
    contText = contText + "\nrad " + std::to_string(sphere->radius) + ";";

    contText = contText + "\ngeometry\n{\n";
    contText = contText + "\tsphere\n\t{\n";
    contText = contText + "\t\ttype searchableSphere;\n";
    contText = contText + "\t\tcentre ($cx $cy $cz);\n";
    contText = contText + "\t\tradius $rad;\n";
    contText = contText + "\t}\n}\n";

    contText = contText + "\nconvertToMeters " +
               std::to_string(params_->cnvrtToMeters) + ";\n";

    double vx = (sphere->center[0]) + tan30 * (sphere->radius);
    double mvx = (sphere->center[0]) - tan30 * (sphere->radius);
    contText = contText + "\nvx " + std::to_string(vx) + ";\n";
    contText = contText + "mvx " + std::to_string(mvx) + ";\n";

    double vy = (sphere->center[1]) + tan30 * (sphere->radius);
    double mvy = (sphere->center[1]) - tan30 * (sphere->radius);
    contText = contText + "\nvy " + std::to_string(vy) + ";\n";
    contText = contText + "mvy " + std::to_string(mvy) + ";\n";

    double vz = (sphere->center[2]) + tan30 * (sphere->radius);
    double mvz = (sphere->center[2]) - tan30 * (sphere->radius);
    contText = contText + "\nvz " + std::to_string(vz) + ";\n";
    contText = contText + "mvz " + std::to_string(mvz) + ";\n";

    double ax = (sphere->center[0]) + sin45 * (sphere->radius);
    double max = (sphere->center[0]) - sin45 * (sphere->radius);
    contText = contText + "\nax " + std::to_string(ax) + ";\n";
    contText = contText + "max " + std::to_string(max) + ";\n";

    double ay = (sphere->center[1]) + sin45 * (sphere->radius);
    double may = (sphere->center[1]) - sin45 * (sphere->radius);
    contText = contText + "\nay " + std::to_string(ay) + ";\n";
    contText = contText + "may " + std::to_string(may) + ";\n";

    double az = (sphere->center[2]) + sin45 * (sphere->radius);
    double naz = (sphere->center[2]) - sin45 * (sphere->radius);
    contText = contText + "\naz " + std::to_string(az) + ";\n";
    contText = contText + "maz " + std::to_string(naz) + ";\n";

    contText = contText + "\nvertices\n";
    contText = contText + "(\n";
    contText = contText + "\t($mvx $mvy $mvz)\n";
    contText = contText + "\t( $vx $mvy $mvz)\n";
    contText = contText + "\t( $vx  $vy $mvz)\n";
    contText = contText + "\t($mvx  $vy $mvz)\n";
    contText = contText + "\t($mvx $mvy  $vz)\n";
    contText = contText + "\t( $vx $mvy  $vz)\n";
    contText = contText + "\t( $vx  $vy  $vz)\n";
    contText = contText + "\t($mvx  $vy  $vz)\n";
    contText = contText + ");\n";

    contText = contText + "\nblocks\n";
    contText = contText + "(\n";
    contText = contText + "\thex (0 1 2 3 4 5 6 7) (" +
               std::to_string(params_->nCells[0]) + " " +
               std::to_string(params_->nCells[1]) + " " +
               std::to_string(params_->nCells[2]) + ") simpleGrading (" +
               std::to_string(sphere->sphrGrading[0]) + " " +
               std::to_string(sphere->sphrGrading[1]) + " " +
               std::to_string(sphere->sphrGrading[2]) + ")\n\n";

    contText = contText + ");\n";

    contText = contText + "\nedges\n";
    contText = contText + "(\n";
    contText = contText + "\tarc 0 1 ($cx $may $maz)\n";
    contText = contText + "\tarc 2 3 ($cx $ay  $maz)\n";
    contText = contText + "\tarc 6 7 ($cx $ay  $az)\n";
    contText = contText + "\tarc 4 5 ($cx $may $az)\n";

    contText = contText + "\tarc 0 3 ($max $cy $maz)\n";
    contText = contText + "\tarc 1 2 ($ax  $cy $maz)\n";
    contText = contText + "\tarc 5 6 ($ax  $cy $az)\n";
    contText = contText + "\tarc 4 7 ($max $cy $az)\n";

    contText = contText + "\tarc 0 4 ($max $may $cz)\n";
    contText = contText + "\tarc 1 5 ($ax  $may $cz)\n";
    contText = contText + "\tarc 2 6 ($ax  $ay  $cz)\n";
    contText = contText + "\tarc 3 7 ($max $ay  $cz)\n";

    contText = contText + ");\n";

    contText = contText + "\nfaces\n(\n";
    contText = contText + "\tproject (0 4 7 3) sphere\n";
    contText = contText + "\tproject (2 6 5 1) sphere\n";
    contText = contText + "\tproject (1 5 4 0) sphere\n";
    contText = contText + "\tproject (3 7 6 2) sphere\n";
    contText = contText + "\tproject (0 3 2 1) sphere\n";
    contText = contText + "\tproject (4 5 6 7) sphere\n";
    contText = contText + ");\n";

    contText = contText + "\nboundary\n(\n";
    contText = contText + "\twalls\n\t{\n";
    contText = contText + "\t\ttype wall;\n";
    contText = contText + "\t\tfaces\n\t\t(\n";
    contText = contText + "\t\t\t(0 4 7 3)\n";
    contText = contText + "\t\t\t(2 6 5 1)\n";
    contText = contText + "\t\t\t(1 5 4 0)\n";
    contText = contText + "\t\t\t(3 7 6 2)\n";
    contText = contText + "\t\t\t(0 3 2 1)\n";
    contText = contText + "\t\t\t(4 5 6 7)\n";
    contText = contText + "\t\t);\n";
    contText = contText + "\t}\n";
    contText = contText + ");\n";
    contText = contText +
               "\n //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * "
               "* * *//";
  } else if (auto cylTaperedCone =
                 std::dynamic_pointer_cast<bmCylTaperedCone>(params_->shape)) {
    contText = contText + "\n\ncx " +
               std::to_string(cylTaperedCone->centerCyl[0]) + ";\n";
    contText =
        contText + "cy " + std::to_string(cylTaperedCone->centerCyl[1]) + ";\n";
    contText =
        contText + "cz " + std::to_string(cylTaperedCone->centerCyl[2]) + ";\n";
    contText =
        contText + "rad1 " + std::to_string(cylTaperedCone->radius1) + ";\n";
    contText = contText + "rad2 " +
               std::to_string(
                   cylTaperedCone->radius2.value_or(cylTaperedCone->radius1)) +
               ";\n";
    contText = contText + "h " + std::to_string(cylTaperedCone->height) + ";\n";
    contText = contText + "cellX " + std::to_string(params_->nCells[0]) + ";\n";
    contText = contText + "cellY " + std::to_string(params_->nCells[1]) + ";\n";
    contText = contText + "cellZ " + std::to_string(params_->nCells[2]) + ";\n";
    contText = contText + "grdX " +
               std::to_string(cylTaperedCone->cylGrading[0]) + ";\n";
    contText = contText + "grdY " +
               std::to_string(cylTaperedCone->cylGrading[1]) + ";\n";
    contText = contText + "grdZ " +
               std::to_string(cylTaperedCone->cylGrading[2]) + ";\n\n";

    double rad1 = cylTaperedCone->radius1;
    double rad2 = cylTaperedCone->radius2.value_or(cylTaperedCone->radius1);
    double X0 = (cylTaperedCone->centerCyl[0]) - sin45 * rad1;
    contText = contText + "\nX0 " + std::to_string(X0) + ";";

    double Y0 = (cylTaperedCone->centerCyl[1]) - sin45 * rad1;
    contText = contText + "\nY0 " + std::to_string(Y0) + ";";

    double X1 = (cylTaperedCone->centerCyl[0]) + sin45 * rad1;
    contText = contText + "\nX1 " + std::to_string(X1) + ";";

    double Y1 = (cylTaperedCone->centerCyl[1]) - sin45 * rad1;
    contText = contText + "\nY1 " + std::to_string(Y1) + ";";

    double X2 = (cylTaperedCone->centerCyl[0]) - sin45 * 0.3 * rad1;
    contText = contText + "\nX2 " + std::to_string(X2) + ";";

    double Y2 = (cylTaperedCone->centerCyl[1]) - sin45 * 0.3 * rad1;
    contText = contText + "\nY2 " + std::to_string(Y2) + ";";

    double X3 = (cylTaperedCone->centerCyl[0]) + sin45 * 0.3 * rad1;
    contText = contText + "\nX3 " + std::to_string(X3) + ";";

    double Y3 = (cylTaperedCone->centerCyl[1]) - sin45 * 0.3 * rad1;
    contText = contText + "\nY3 " + std::to_string(Y3) + ";";

    double X4 = (cylTaperedCone->centerCyl[0]) + sin45 * 0.3 * rad1;
    contText = contText + "\nX4 " + std::to_string(X4) + ";";

    double Y4 = (cylTaperedCone->centerCyl[1]) + sin45 * 0.3 * rad1;
    contText = contText + "\nY4 " + std::to_string(Y4) + ";";

    double X5 = (cylTaperedCone->centerCyl[0]) + sin45 * rad1;
    contText = contText + "\nX5 " + std::to_string(X5) + ";";

    double Y5 = (cylTaperedCone->centerCyl[1]) + sin45 * rad1;
    contText = contText + "\nY5 " + std::to_string(Y5) + ";";

    double X6 = (cylTaperedCone->centerCyl[0]) - sin45 * 0.3 * rad1;
    contText = contText + "\nX6 " + std::to_string(X6) + ";";

    double Y6 = (cylTaperedCone->centerCyl[1]) + sin45 * 0.3 * rad1;
    contText = contText + "\nY6 " + std::to_string(Y6) + ";";

    double X7 = (cylTaperedCone->centerCyl[0]) - sin45 * rad1;
    contText = contText + "\nX7 " + std::to_string(X7) + ";";

    double Y7 = (cylTaperedCone->centerCyl[1]) + sin45 * rad1;
    contText = contText + "\nY7 " + std::to_string(Y7) + ";";

    double X8 = (cylTaperedCone->centerCyl[0]) - sin45 * rad2;
    contText = contText + "\nX8 " + std::to_string(X8) + ";";

    double Y8 = (cylTaperedCone->centerCyl[1]) - sin45 * rad2;
    contText = contText + "\nY8 " + std::to_string(Y8) + ";";

    double X9 = (cylTaperedCone->centerCyl[0]) + sin45 * rad2;
    contText = contText + "\nX9 " + std::to_string(X9) + ";";

    double Y9 = (cylTaperedCone->centerCyl[1]) - sin45 * rad2;
    contText = contText + "\nY9 " + std::to_string(Y9) + ";";

    double X10 = (cylTaperedCone->centerCyl[0]) - sin45 * 0.3 * rad2;
    contText = contText + "\nX10 " + std::to_string(X10) + ";";

    double Y10 = (cylTaperedCone->centerCyl[1]) - sin45 * 0.3 * rad2;
    contText = contText + "\nY10 " + std::to_string(Y10) + ";";

    double X11 = (cylTaperedCone->centerCyl[0]) + sin45 * 0.3 * rad2;
    contText = contText + "\nX11 " + std::to_string(X11) + ";";

    double Y11 = (cylTaperedCone->centerCyl[1]) - sin45 * 0.3 * rad2;
    contText = contText + "\nY11 " + std::to_string(Y11) + ";";

    double X12 = (cylTaperedCone->centerCyl[0]) + sin45 * 0.3 * rad2;
    contText = contText + "\nX12 " + std::to_string(X12) + ";";

    double Y12 = (cylTaperedCone->centerCyl[1]) + sin45 * 0.3 * rad2;
    contText = contText + "\nY12 " + std::to_string(Y12) + ";";

    double X13 = (cylTaperedCone->centerCyl[0]) + sin45 * rad2;
    contText = contText + "\nX13 " + std::to_string(X13) + ";";

    double Y13 = (cylTaperedCone->centerCyl[1]) + sin45 * rad2;
    contText = contText + "\nY13 " + std::to_string(Y13) + ";";

    double X14 = (cylTaperedCone->centerCyl[0]) - sin45 * 0.3 * rad2;
    contText = contText + "\nX14 " + std::to_string(X14) + ";";

    double Y14 = (cylTaperedCone->centerCyl[1]) + sin45 * 0.3 * rad2;
    contText = contText + "\nY14 " + std::to_string(Y14) + ";";

    double X15 = (cylTaperedCone->centerCyl[0]) - sin45 * rad2;
    contText = contText + "\nX15 " + std::to_string(X15) + ";";

    double Y15 = (cylTaperedCone->centerCyl[1]) + sin45 * rad2;
    contText = contText + "\nY15 " + std::to_string(Y15) + ";";

    double arcX1 = cylTaperedCone->centerCyl[0];
    contText = contText + "\narcX1 " + std::to_string(arcX1) + ";";

    double arcY1 = (cylTaperedCone->centerCyl[1] - rad1);
    contText = contText + "\narcY1 " + std::to_string(arcY1) + ";";

    double arcX2 = (cylTaperedCone->centerCyl[0] + rad1);
    contText = contText + "\narcX2 " + std::to_string(arcX2) + ";";

    double arcY2 = cylTaperedCone->centerCyl[1];
    contText = contText + "\narcY2 " + std::to_string(arcY2) + ";";

    double arcX3 = cylTaperedCone->centerCyl[0];
    contText = contText + "\narcX3 " + std::to_string(arcX3) + ";";

    double arcY3 = (cylTaperedCone->centerCyl[1] + rad1);
    contText = contText + "\narcY3 " + std::to_string(arcY3) + ";";

    double arcX4 = (cylTaperedCone->centerCyl[0] - rad1);
    contText = contText + "\narcX4 " + std::to_string(arcX4) + ";";

    double arcY4 = cylTaperedCone->centerCyl[1];
    contText = contText + "\narcY4 " + std::to_string(arcY4) + ";";

    double arcX5 = cylTaperedCone->centerCyl[0];
    contText = contText + "\narcX5 " + std::to_string(arcX5) + ";";

    double arcY5 = (cylTaperedCone->centerCyl[1] - rad2);
    contText = contText + "\narcY5 " + std::to_string(arcY5) + ";";

    double arcX6 = (cylTaperedCone->centerCyl[0] + rad2);
    contText = contText + "\narcX6 " + std::to_string(arcX6) + ";";

    double arcY6 = cylTaperedCone->centerCyl[1];
    contText = contText + "\narcY6 " + std::to_string(arcY6) + ";";

    double arcX7 = cylTaperedCone->centerCyl[0];
    contText = contText + "\narcX7 " + std::to_string(arcX7) + ";";

    double arcY7 = (cylTaperedCone->centerCyl[1] + rad2);
    contText = contText + "\narcY7 " + std::to_string(arcY7) + ";";

    double arcX8 = (cylTaperedCone->centerCyl[0] - rad2);
    contText = contText + "\narcX8 " + std::to_string(arcX8) + ";";

    double arcY8 = cylTaperedCone->centerCyl[1];
    contText = contText + "\narcY8 " + std::to_string(arcY8) + ";";

    double fz = (cylTaperedCone->centerCyl[1] + (cylTaperedCone->height));
    contText = contText + "\nfz " + std::to_string(fz) + ";";

    contText = contText + "\nconvertToMeters " +
               std::to_string(params_->cnvrtToMeters) + ";\n";

    contText = contText + "\nvertices\n(\n\n";
    contText = contText + "\t($X0 $Y0 $cz)\n";
    contText = contText + "\t($X1 $Y1 $cz)\n";
    contText = contText + "\t($X2 $Y2 $cz)\n";
    contText = contText + "\t($X3 $Y3 $cz)\n";
    contText = contText + "\t($X4 $X4 $cz)\n";
    contText = contText + "\t($X5 $Y5 $cz)\n";
    contText = contText + "\t($X6 $Y6 $cz)\n";
    contText = contText + "\t($X7 $Y7 $cz)\n";
    contText = contText + "\t($X8 $Y8 $fz)\n";
    contText = contText + "\t($X9 $Y9 $fz)\n";
    contText = contText + "\t($X10 $Y10 $fz)\n";
    contText = contText + "\t($X11 $Y11 $fz)\n";
    contText = contText + "\t($X12 $X12 $fz)\n";
    contText = contText + "\t($X13 $Y13 $fz)\n";
    contText = contText + "\t($X14 $Y14 $fz)\n";
    contText = contText + "\t($X15 $Y15 $fz)\n";
    contText = contText + ");\n\n";

    contText = contText + "\nblocks\n(\n";
    contText = contText +
               "\thex (2 3 4 6 10 11 12 14) ($cellX $cellY $cellZ) " +
               "simpleGrading ($grdX $grdY $grdZ)\n";
    contText = contText + "\thex (0 1 3 2 8 9 11 10) ($cellX $cellY $cellZ) " +
               "simpleGrading ($grdX $grdY $grdZ)\n";
    contText = contText + "\thex (3 1 5 4 11 9 13 12) ($cellX $cellY $cellZ) " +
               "simpleGrading ($grdX $grdY $grdZ)\n";
    contText = contText +
               "\thex (4 5 7 6 12 13 15 14) ($cellX $cellY $cellZ) " +
               "simpleGrading ($grdX $grdY $grdZ)\n";
    contText = contText + "\thex (6 7 0 2 14 15 8 10) ($cellX $cellY $cellZ) " +
               "simpleGrading ($grdX $grdY $grdZ)\n";
    contText = contText + ");\n";

    contText = contText + "\nedges\n(\n";
    contText = contText + "\tarc 0 1 ($arcX1 $arcY1 $cz)\n";
    contText = contText + "\tarc 1 5 ($arcX2 $arcY2 $cz)\n";
    contText = contText + "\tarc 5 7 ($arcX3 $arcY3 $cz)\n";
    contText = contText + "\tarc 7 0 ($arcX4 $arcY4 $cz)\n";
    contText = contText + "\tarc 8 9 ($arcX5 $arcY5 $fz)\n";
    contText = contText + "\tarc 9 13 ($arcX6 $arcY6 $fz)\n";
    contText = contText + "\tarc 13 15 ($arcX7 $arcY7 $fz)\n";
    contText = contText + "\tarc 15 8 ($arcX8 $arcY8 $fz)\n";
    contText = contText + ");\n\n";

    contText = contText + "boundary\n(\n";
    contText = contText + "\tinlet\n\t{\n";
    contText = contText + "\t\ttype patch;\n";
    contText = contText + "\t\tfaces\n\t\t(\n";
    contText = contText + "\t\t\t(0 1 2 3)\n";
    contText = contText + "\t\t\t(1 3 4 5)\n";
    contText = contText + "\t\t\t(4 5 7 6)\n";
    contText = contText + "\t\t\t(0 2 6 7)\n";
    contText = contText + "\t\t\t(2 3 4 6)\n";
    contText = contText + "\t\t);\n";
    contText = contText + "\t}\n";

    contText = contText + "\toutlet\n\t{\n";
    contText = contText + "\t\ttype patch;\n";
    contText = contText + "\t\tfaces\n\t\t(\n";
    contText = contText + "\t\t\t(8 9 10 11)\n";
    contText = contText + "\t\t\t(11 9 12 13)\n";
    contText = contText + "\t\t\t(12 13 14 15)\n";
    contText = contText + "\t\t\t(14 15 10 8)\n";
    contText = contText + "\t\t\t(10 11 14 12)\n";
    contText = contText + "\t\t);\n";
    contText = contText + "\t}\n";

    contText = contText + "\twalls\n\t{\n";
    contText = contText + "\t\ttype patch;\n";
    contText = contText + "\t\tfaces\n\t\t(\n";
    contText = contText + "\t\t\t(1 5 13 9)\n";
    contText = contText + "\t\t\t(5 7 15 13)\n";
    contText = contText + "\t\t\t(7 0 8 15)\n";
    contText = contText + "\t\t\t(0 1 9 8)\n";
    contText = contText + "\t\t);\n";
    contText = contText + "\t}\n";
    contText = contText + ");\n";

    contText =
        contText +
        "\n //* * * * * * * * * * * * * * * * * * * * * * * * * * * * *//";
  } else {
    std::cerr << "User cannot select multiple geometries in single run!\n"
              << std::endl;
    exit(1);
  }

  // Keep this dictionary in memory
  Foam::dictionary tmptmpDuc =
      new Foam::dictionary(Foam::IStringStream(contText)(), true);
  blockMshDict_ =
      std::unique_ptr<Foam::dictionary>(new Foam::dictionary("blockMeshDict"));
  blockMshDict_->merge(tmptmpDuc);

  // Write in traditional way due to formatting issues
  if (write) {
    // creating a base system directory
    const char dir_path[] = "./system";
    boost::filesystem::path dir(dir_path);
    try {
      boost::filesystem::create_directory(dir);
    } catch (boost::filesystem::filesystem_error &e) {
      std::cerr << "Problem in creating system directory for the cfMesh"
                << "\n";
      std::cerr << e.what() << std::endl;
      throw;
    }
    std::ofstream contDict;
    contDict.open(std::string(dir_path) + "/blockMeshDict");
    contDict << contText;
    contDict.close();
  }
}
