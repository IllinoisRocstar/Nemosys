#define _USE_MATH_DEFINES

#include "meshBase.H"
#include <ANN/ANN.h>
#include <Eigen/Geometry>
#include <chrono>
#include <fstream>
#include <hmxShape.H>
#include <icosidodecahedronShape.H>
#include <iostream>
#include <map>
#include <petnShape.H>
#include <rocPack.H>
#include <rocPackShape.H>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <vtkCell.h>
#include <vtkCell3D.h>
#include <vtkMesh.H>

// GMSH Header
#include <gmsh.h>

namespace NEM {

namespace GEO {

// class constructor
rocPack::rocPack(const std::string &fname, const std::string &outName) {
  std::cout << "rocPack class constructed!" << std::endl;
  InFile = fname;
  OutFile = outName;
}

void rocPack::rocPack2Surf() {
  // Parses output file from Rocpack
  rocParser();

  // Generates geometry from parsed database.
  rocToGeom();

  // Generates surface mesh for created geometry and writes in surface files
  geomToSurf();

  gmsh::finalize();
  std::cout << " - End of process!" << std::endl;
}

void rocPack::rocPack2Periodic3D() {
  periodic3D = true;

  // Parses output file from Rocpack
  rocParser();

  // Generates geometry from parsed database.
  rocToGeom();

  // Generates 3D periodic mesh for created geometry and writes in surface files
  geomToPeriodic3D();

  // Creates cohesive elements if specified
  if (internalCohesiveBool)
    createCohesiveElements(OutFile,OutFile);

  // Writes periodic mesh maps to CSV files
  writePeriodicNodes();

  gmsh::finalize();
  std::cout << " - End of process!" << std::endl;
}

void rocPack::removeBoundaryVolumes() { removeBoundaryPacks = true; }

void rocPack::enableCohesiveElements() { internalCohesiveBool = true; }

void rocPack::setPeriodicGeometry() { enablePeriodicity = true; }

void rocPack::setPeriodicMesh() { periodic3D = true; }

void rocPack::shrinkVolumes(const double percntg) {
  shrinkScale = percntg;
  enableScaling = true;
}

void rocPack::smoothSurfaces(const int smoothingParam) {
  smoothingIter = smoothingParam;
  enableSmoothing = true;
}

void rocPack::setMeshSize(const double size) {
  meshSz = size;
}

void rocPack::enablePhysicalGrps() { enablePhysGrp = true; }

void rocPack::enableTwoPhysGrps() { just2Physgrps = true; }

void rocPack::enableSurfacePatches() { assignSidePatches = true; }

void rocPack::initialize() {
  gmsh::initialize();
  gmsh::model::add("Pack");
  gmsh::option::setNumber("General.NumThreads", 10);
  gmsh::option::setNumber("Geometry.OCCBooleanPreserveNumbering", 1);
  gmsh::option::setNumber("Geometry.OCCParallel", 1);
}

void rocPack::rocParser() {
  // Check if file exist
  std::cout << " - Parsing output file ... " << std::endl;
  std::ifstream rocOut(InFile);
  if (rocOut.is_open()) {
    if (cstmDomain) {
      // Nothing
    } else {
      // Getting box dimensions
      std::vector<std::string> rocTokens;
      std::string line = findWord("boundary");
      rocTokens = Tokenize(line, ' ');
      if (rocTokens[6] != "periodic") {
        std::cerr << "Please select output file with periodic geometries!"
                  << std::endl;
        throw;
      }

      Xdim = std::atof(strToChar(rocTokens[3]));
      Ydim = std::atof(strToChar(rocTokens[4]));
      Zdim = std::atof(strToChar(rocTokens[5]));
      boxPt.push_back((-Xdim / 2) + xUDF);
      boxPt.push_back((-Ydim / 2) + yUDF);
      boxPt.push_back((-Zdim / 2) + zUDF);
    }

    // Stores all lines in file to myLines
    std::string linesOut;             // Temporary string for multiple uses
    std::vector<std::string> myLines; // Stores whole file line by line
    while (std::getline(rocOut, linesOut))
      myLines.push_back(linesOut);

    // Checking for user-specified options
    for (int i = 0; i < myLines.size(); i++) {
      if (myLines[i].find("SetPeriodicity") != std::string::npos)
        if (myLines[i].find("true") != std::string::npos ||
            myLines[i].find("True") != std::string::npos ||
            myLines[i].find("TRUE") != std::string::npos)
          enablePeriodicity = true;

      if (myLines[i].find("RemoveBoundaryPacks") != std::string::npos)
        if (myLines[i].find("true") != std::string::npos ||
            myLines[i].find("True") != std::string::npos ||
            myLines[i].find("TRUE") != std::string::npos)
          removeBoundaryPacks = true;
    }

    if (periodic3D == true && enablePeriodicity == false) {
      std::cerr << "Periodicity needed if using periodic meshing!" << std::endl
                << "Please enable SetPeriodicity option to"
                << " get 3D periodic mesh" << std::endl;
      throw;
    }

    // Checking if crystal shapes are present
    // crystalNames and verts/faces go parallel in terms of indexing
    int nCrystals = 0;
    for (int i = 0; i < myLines.size(); i++) {
      if (myLines[i].find("#import") != std::string::npos ||
          myLines[i].find("import") == std::string::npos) {
        // Nothing
      } else if (myLines[i].find("import") != std::string::npos ||
                 myLines[i].find("#import") == std::string::npos) {
        nCrystals++;
        std::vector<std::string> importStr;
        importStr = Tokenize(myLines[i], "//\"");
        crystalNames.push_back(importStr[2]);
      } else {
        // Nothing
      }
    }

    if (nCrystals > 0) {
      verts.resize(nCrystals);
      faces.resize(nCrystals);
      for (int i = 0; i < nCrystals; i++) {
        rocPackShape *getData = rocPackShape::getShape(crystalNames[i]);
        verts[i] = getData->getVertices();
        faces[i] = getData->getFaces();
      }
      normalizeVerts();
    }

    // Checking if any shapes are present. If present, extract imp data
    for (int i = 0; i < myLines.size(); i++) {
      if (myLines[i].find("shape ") != std::string::npos) {
        std::vector<std::string> baseShapes;
        baseShapes = Tokenize(myLines[i], "         ");
        shapeNames.push_back(baseShapes[3]);
        uniqueNames.push_back(baseShapes[1]);

        if (baseShapes[3] == "cylinder") {
          cylParams.resize(2);
          cylParams[0] = std::atof(strToChar(baseShapes[5]));
          cylParams[1] = std::atof(strToChar(baseShapes[6]));
        }

        if (baseShapes[3] == "ellipsoid") {
          ellipsoidRad.resize(3);
          ellipsoidRad[0] = std::atof(strToChar(baseShapes[5]));
          ellipsoidRad[1] = std::atof(strToChar(baseShapes[6]));
          ellipsoidRad[2] = std::atof(strToChar(baseShapes[7]));
          ellipsoidPresent = true;
        }
      } else {
      }
    }

    if (ellipsoidPresent)
      std::cout << "    !!Warning!!" << std::endl
                << "    This geometry will not be periodic due to presence of "
                << "Ellipsoid shapes!" << std::endl
                << "    Ellipsoid shapes cannot be made periodic as of now!"
                << std::endl;

    // Getting line number where pack data starts
    // Also counting total number of packs
    int numPacks = 0;
    int iter = 0;
    for (int i = 0; i < myLines.size(); i++) {
      if (myLines[i].find("translate") != std::string::npos)
        numPacks++;
      if (myLines[i].find("translate") != std::string::npos && iter == 0)
        iter = i;
      else {
      }
    }
    iter = iter - 1; // Start parsing from "iter"

    // Pack shape data storage initialization
    // Getting pack data
    translateParams.resize(numPacks);
    rotateParams.resize(numPacks);
    scaleOfPack.resize(numPacks);
    nameOfPcks.resize(numPacks);
    int h = 0;
    std::string pckNameStr = " ";
    std::string translateStr = "<,,>";
    std::string rotateStr = "<,,>";
    std::string scaleStr = "       ";

    std::vector<double> udfTranslate;
    udfTranslate.resize(3);
    udfTranslate[0] = xUDF;
    udfTranslate[1] = yUDF;
    udfTranslate[2] = zUDF;

    for (int i = iter; i < myLines.size(); i++) {
      translateParams[h].resize(3);
      rotateParams[h].resize(4);
      std::vector<std::string> pckNmDataTokens;
      pckNmDataTokens.resize(2);
      std::vector<std::string> trnsltDataTokens;
      trnsltDataTokens.resize(4);
      std::vector<std::string> rttDataTokens;
      rttDataTokens.resize(7);
      std::vector<std::string> scaleDataTokens;
      scaleDataTokens.resize(8);

      // Name of pack shape
      pckNmDataTokens = getShapeData(i, pckNameStr, myLines);
      nameOfPcks[h] = pckNmDataTokens[0];

      // Translate Parameters
      trnsltDataTokens = getShapeData(i + 1, translateStr, myLines);

      for (int j = 0; j < 3; j++)
        translateParams[h][j] =
            std::atof(strToChar(trnsltDataTokens[j + 1])) + udfTranslate[j];

      trnsltDataTokens.clear();

      // Rotate Parameters
      rttDataTokens = getShapeData(i + 2, rotateStr, myLines);

      for (int j = 0; j < 4; j++)
        rotateParams[h][j] = std::atof(strToChar(rttDataTokens[j + 1]));

      rttDataTokens.clear();

      // Scale Parameters
      scaleDataTokens = getShapeData(i + 3, scaleStr, myLines);
      scaleOfPack[h] = std::atof(strToChar(scaleDataTokens[1]))*shrinkScale;
      scaleDataTokens.clear();

      i = i + 5;
      h++;
    }
    rocOut.close();
  } else {
    std::cerr << "Cannot open/find " << InFile << " file!" << std::endl;
    throw;
  }
}

// Loops through parsed data and makes gmsh geometries.
void rocPack::rocToGeom() {
  std::cout << " - Creating pack geometries ... " << std::endl;
  // Starting Gmsh commands
  initialize();

  for (int i = 0; i < nameOfPcks.size(); i++) {
    if (nameOfPcks[i] == "sphere") {
      makeSphere(i);
    }

    if (uniqueNames.size() > 0) {
      for (int j = 0; j < uniqueNames.size(); j++) {
        if (nameOfPcks[i] == uniqueNames[j] && shapeNames[j] == "ellipsoid")
          makeEllipsoid(i);
        else if (nameOfPcks[i] == uniqueNames[j] && shapeNames[j] == "cylinder")
          makeCylinder(i);
        else {
        }
      }
    }

    if (crystalNames.size() > 0) {
      for (int k = 0; k < crystalNames.size(); k++)
        if (crystalNames[k] == nameOfPcks[i])
          makeCrystalShape(i, k);
    }
  }

  // Tagging packs intersecting domain boundary
  tagBoundaryPacks();

  if (removeBoundaryPacks == true && ellipsoidPresent == false) {
    gmsh::model::occ::remove(bndryPackTags, true);
    gmsh::model::occ::synchronize();
  } else if (removeBoundaryPacks == true && ellipsoidPresent == true) {
    gmsh::model::geo::remove(bndryPackTags, true);
    gmsh::model::geo::synchronize();
  } else {
    // Nothing
  }

  if (removeBoundaryPacks == true && periodic3D == true) {
    std::vector<std::pair<int, int>> tagsInsidePacks;
    gmsh::model::getEntities(tagsInsidePacks, 3);
    mapPeriodicSurfaces(tagsInsidePacks);
  }

  if (ellipsoidPresent == true || enablePeriodicity == false ||
      cstmDomain == true) {
    // Nothing
  } else {
    std::cout << " - Ensuring periodicity of geometry ..." << std::endl;
    makePeriodic(removeBoundaryPacks);
    gmsh::model::occ::synchronize();
  }
}

void rocPack::geomToSurf() {
  std::cout << " - Writing triangulated surface files ..." << std::endl;

  std::vector<std::pair<int, int>> tagsPts;
  gmsh::model::getEntities(tagsPts, 0);
  gmsh::model::mesh::setSize(tagsPts, meshSz);
  gmsh::option::setNumber("Mesh.Algorithm", meshingAlgorithm);
  gmsh::model::mesh::generate(2);
  gmsh::model::mesh::refine();
  gmsh::model::mesh::refine();
  gmsh::model::mesh::removeDuplicateNodes();

  if (OutFile.find(".stl") != std::string::npos)
    geomToSTL(OutFile);
  else if (OutFile.find(".vtu") != std::string::npos ||
           OutFile.find(".vtk") != std::string::npos) {
    size_t pos = std::string::npos;
    while ((pos = OutFile.find(".vtu")) != std::string::npos) {
      OutFile.erase(pos, 4);
    }
    while ((pos = OutFile.find(".vtk")) != std::string::npos) {
      OutFile.erase(pos, 4);
    }
    geomToVTK(OutFile + ".msh");
  } else if (OutFile.find(".msh") != std::string::npos) {
    geomToMsh(OutFile);
  } else {
    gmsh::write(OutFile);
  }

  if (defOutputs) {
    size_t lastindex = OutFile.find_last_of(".");
    std::string rawname = OutFile.substr(0, lastindex);
    std::string stlFile = rawname + ".stl";
    std::string vtkFile = rawname + ".msh";
    std::string mshFile = rawname + ".msh";
    geomToSTL(stlFile);
    geomToMsh(mshFile);
    geomToVTK(vtkFile);
  }
}

void rocPack::geomToPeriodic3D() {
  // Figure out how to save geometry
  std::cout << " - Writing Periodic mesh files ..." << std::endl;

  std::vector<std::pair<int, int>> tagsPts;
  gmsh::model::getEntities(tagsPts, 0);
  gmsh::model::mesh::setSize(tagsPts, meshSz);

  if ((enablePhysGrp) || (just2Physgrps) || (physGrpPerShape))
    gmsh::option::setNumber("Mesh.StlOneSolidPerSurface", 2);

  gmsh::option::setNumber("Mesh.Algorithm3D", meshingAlgorithm);
  gmsh::option::setNumber("Mesh.SaveGroupsOfNodes", 1);
  gmsh::model::mesh::generate(3);
  gmsh::model::mesh::removeDuplicateNodes();

  if (OutFile.find(".stl") != std::string::npos)
    geomToSTL(OutFile);
  else if (OutFile.find(".vtu") != std::string::npos ||
           OutFile.find(".vtk") != std::string::npos) {
    size_t pos = std::string::npos;
    while ((pos = OutFile.find(".vtu")) != std::string::npos) {
      OutFile.erase(pos, 4);
    }
    while ((pos = OutFile.find(".vtk")) != std::string::npos) {
      OutFile.erase(pos, 4);
    }
    geomToVTK(OutFile + ".msh");
  } else if (OutFile.find(".msh") != std::string::npos) {
    geomToMsh(OutFile);
  } else {
    gmsh::write(OutFile);
  }

  if (defOutputs) {
    size_t lastindex = OutFile.find_last_of(".");
    std::string rawname = OutFile.substr(0, lastindex);
    std::string stlFile = rawname + ".stl";
    std::string vtkFile = rawname + ".msh";
    std::string mshFile = rawname + ".msh";
    geomToSTL(stlFile);
    geomToMsh(mshFile);
    geomToVTK(vtkFile);
  }
}

// Private methods
void rocPack::geomToSTL(const std::string &writeFile) {
  // Writes created periodic geometries into STL file
  gmsh::write(writeFile);
}

void rocPack::geomToVTK(const std::string &writeFile) {
  // Writes created periodic geometries into VTU file with metadata
  size_t lastindex = writeFile.find_last_of(".");
  std::string rawname = writeFile.substr(0, lastindex);
  rawname = rawname + "_oldMSH.msh";
  gmsh::option::setNumber("Mesh.MshFileVersion", 2.2);
  geomToMsh(rawname);
  meshBase *mb = meshBase::exportGmshToVtk(rawname);
  mb->write(rawname);

  if (mb)
    delete mb;
}

void rocPack::geomToMsh(const std::string &writeFile) {
  // Writes created periodic geometries into .msh file
  gmsh::write(writeFile);
}

void rocPack::makePeriodic(const bool rmbPacks) {
  if (rmbPacks == false) {

    // To get processor clocktime
    std::clock_t start;
    double duration;
    start = std::clock();

    std::vector<int> xTranslate;
    std::vector<int> yTranslate;
    std::vector<int> zTranslate;

    xTranslate.push_back(Xdim);   // 0 -> X+
    yTranslate.push_back(0);      // 0 -> X+
    zTranslate.push_back(0);      // 0 -> X+
    xTranslate.push_back(Xdim);   //// 1 -> X+ Z+
    yTranslate.push_back(0);      //// 1 -> X+ Z+
    zTranslate.push_back(Zdim);   //// 1 -> X+ Z+
    xTranslate.push_back(Xdim);   // 2 -> X+ Z-
    yTranslate.push_back(0);      // 2 -> X+ Z-
    zTranslate.push_back(-Zdim);  // 2 -> X+ Z-

    xTranslate.push_back(-Xdim);  // 3 -> X-
    yTranslate.push_back(0);      // 3 -> X-
    zTranslate.push_back(0);      // 3 -> X-
    xTranslate.push_back(-Xdim);  //// 4 -> X- Z+
    yTranslate.push_back(0);      //// 4 -> X- Z+
    zTranslate.push_back(Zdim);   //// 4 -> X- Z+
    xTranslate.push_back(-Xdim);  // 5 -> X- Z-
    yTranslate.push_back(0);      // 5 -> X- Z-
    zTranslate.push_back(-Zdim);  // 5 -> X- Z-

    xTranslate.push_back(0);      // 6 -> Y+
    yTranslate.push_back(Ydim);   // 6 -> Y+
    zTranslate.push_back(0);      // 6 -> Y+
    xTranslate.push_back(0);      //// 7 -> Y+ Z+
    yTranslate.push_back(Ydim);   //// 7 -> Y+ Z+
    zTranslate.push_back(Zdim);   //// 7 -> Y+ Z+
    xTranslate.push_back(0);      // 8 -> Y+ Z-
    yTranslate.push_back(Ydim);   // 8 -> Y+ Z-
    zTranslate.push_back(-Zdim);  // 8 -> Y+ Z-

    xTranslate.push_back(0);      // 9 -> Y-
    yTranslate.push_back(-Ydim);  // 9 -> Y-
    zTranslate.push_back(0);      // 9 -> Y-
    xTranslate.push_back(0);      //// 10 -> Y- Z+
    yTranslate.push_back(-Ydim);  //// 10 -> Y- Z+
    zTranslate.push_back(Zdim);   //// 10 -> Y- Z+
    xTranslate.push_back(0);      // 11 -> Y- Z-
    yTranslate.push_back(-Ydim);  // 11 -> Y- Z-
    zTranslate.push_back(-Zdim);  // 11 -> Y- Z-

    xTranslate.push_back(Xdim);   // 12 -> X+ Y+
    yTranslate.push_back(Ydim);   // 12 -> X+ Y+
    zTranslate.push_back(0);      // 12 -> X+ Y+
    xTranslate.push_back(Xdim);   //// 13 -> X+ Y+ Z+
    yTranslate.push_back(Ydim);   //// 13 -> X+ Y+ Z+
    zTranslate.push_back(Zdim);   //// 13 -> X+ Y+ Z+
    xTranslate.push_back(Xdim);   // 14 -> X+ Y+ Z-
    yTranslate.push_back(Ydim);   // 14 -> X+ Y+ Z-
    zTranslate.push_back(-Zdim);  // 14 -> X+ Y+ Z-

    xTranslate.push_back(Xdim);   // 15 -> X+ Y-
    yTranslate.push_back(-Ydim);  // 15 -> X+ Y-
    zTranslate.push_back(0);      // 15 -> X+ Y-
    xTranslate.push_back(Xdim);   //// 16 -> X+ Y- Z+
    yTranslate.push_back(-Ydim);  //// 16 -> X+ Y- Z+
    zTranslate.push_back(Zdim);   //// 16 -> X+ Y- Z+
    xTranslate.push_back(Xdim);   // 17 -> X+ Y- Z-
    yTranslate.push_back(-Ydim);  // 17 -> X+ Y- Z-
    zTranslate.push_back(-Zdim);  // 17 -> X+ Y- Z-

    xTranslate.push_back(-Xdim);  // 18 -> X- Y+
    yTranslate.push_back(Ydim);   // 18 -> X- Y+
    zTranslate.push_back(0);      // 18 -> X- Y+
    xTranslate.push_back(-Xdim);  //// 19 -> X- Y+ Z+
    yTranslate.push_back(Ydim);   //// 19 -> X- Y+ Z+
    zTranslate.push_back(Zdim);   //// 19 -> X- Y+ Z+
    xTranslate.push_back(-Xdim);  // 20 -> X- Y+ Z-
    yTranslate.push_back(Ydim);   // 20 -> X- Y+ Z-
    zTranslate.push_back(-Zdim);  // 20 -> X- Y+ Z-

    xTranslate.push_back(-Xdim);  // 21 -> X- Y-
    yTranslate.push_back(-Ydim);  // 21 -> X- Y-
    zTranslate.push_back(0);      // 21 -> X- Y-
    xTranslate.push_back(-Xdim);  //// 22 -> X- Y- Z+
    yTranslate.push_back(-Ydim);  //// 22 -> X- Y- Z+
    zTranslate.push_back(Zdim);   //// 22 -> X- Y- Z+
    xTranslate.push_back(-Xdim);  // 23 -> X- Y- Z-
    yTranslate.push_back(-Ydim);  // 23 -> X- Y- Z-
    zTranslate.push_back(-Zdim);  // 23 -> X- Y- Z-

    xTranslate.push_back(0);      // 24 -> Z+
    yTranslate.push_back(0);      // 24 -> Z+
    zTranslate.push_back(Zdim);   // 24 -> Z+
    xTranslate.push_back(0);      //// 25 -> Z-
    yTranslate.push_back(0);      //// 25 -> Z-
    zTranslate.push_back(-Zdim);  //// 25 -> Z-

    /*// Need some implementation of picking up boundary shapes using get Entities
    // in bounding box command and translate only those. Aim to get points in a
    // bounding box and iteratively find parent entity.

    // If there are less than 25 total volumes, just use the whole tagged
    // boundary volumes, otherwise use the selected volumes from each side
    // (left, right, up, etc..) 

    double tol = 0.5;
    double edgeTol = 0.001;
    // double Xdim2 = 0;
    // double Ydim2 = 0;
    // double Zdim2 = 0;
    
    // left
    std::vector<std::pair<int, int>> entityLeft;
    gmsh::model::getEntitiesInBoundingBox(boxPt[0] - tol, boxPt[1] - tol, 
                                          boxPt[2] - tol, boxPt[0] + edgeTol,
                                          boxPt[1] + Ydim + tol,
                                          Zdim + boxPt[2] + tol, entityLeft,0);

    // Right
    std::vector<std::pair<int, int>> entityRight;
    gmsh::model::getEntitiesInBoundingBox(Xdim + boxPt[0] - edgeTol, boxPt[1] - tol, 
                                          boxPt[2] - tol, Xdim + boxPt[0] + tol,
                                          Ydim + boxPt[1] + tol,
                                          Zdim + boxPt[2] + tol, entityRight,0);

    // Up
    std::vector<std::pair<int, int>> entityUp;
    gmsh::model::getEntitiesInBoundingBox(boxPt[0] - tol, Ydim + boxPt[1] - edgeTol,
                                          boxPt[2] - tol, Xdim + boxPt[0] + tol,
                                          Ydim + boxPt[1] + tol,
                                          Zdim + boxPt[2] + tol, entityUp, 0);

    // Down
    std::vector<std::pair<int, int>> entityDown;
    gmsh::model::getEntitiesInBoundingBox(boxPt[0] - tol, boxPt[1] - tol, 
                                          boxPt[2] - tol, Xdim + boxPt[0] + tol,
                                          boxPt[1] + edgeTol, Zdim + boxPt[2] + tol, 
                                          entityDown, 0);

    // Back
    std::vector<std::pair<int, int>> entityBack;
    gmsh::model::getEntitiesInBoundingBox(boxPt[0] - tol, boxPt[1] - tol,
                                          boxPt[2] - tol, Xdim + boxPt[0] + tol,
                                          Ydim + boxPt[1] + tol, boxPt[2] + edgeTol,
                                          entityBack, 0);

    //Front
    std::vector<std::pair<int, int>> entityFront;
    gmsh::model::getEntitiesInBoundingBox(boxPt[0] - tol, boxPt[1] - tol,
                                          Zdim + boxPt[2] - edgeTol,
                                          Xdim + boxPt[0] + tol,
                                          Ydim + boxPt[1] + tol, 
                                          Zdim + boxPt[2] + tol,
                                          entityFront, 0);

    // Now start getting parent entities for captured points in other vectors
    std::vector<std::pair<int,int>> leftShapes;
    std::vector<std::pair<int,int>> rightShapes;
    std::vector<std::pair<int,int>> UpShapes;
    std::vector<std::pair<int,int>> DownShapes;
    std::vector<std::pair<int,int>> backShapes;
    std::vector<std::pair<int,int>> frontShapes;

    for (int i=0; i<bndryPackTags.size(); i++) {
      std::vector<std::pair<int,int>> tmpTagVec;
      std::vector<std::pair<int,int>> outPoints;
      tmpTagVec.push_back(std::make_pair(bndryPackTags[i].first,
                                         bndryPackTags[i].second));
      gmsh::model::getBoundary(tmpTagVec,outPoints,true,true,true);

      int matchFound = 0;

      // Left
      for (int j=0; j<outPoints.size(); j++) {
        for (int k=0; k<entityLeft.size(); k++) {
          if (outPoints[j].second == entityLeft[k].second) {
            matchFound++;
            break;
          }
        }
        if (matchFound > 0)
          break;
      }
      if (matchFound > 0) {
        leftShapes.push_back(std::make_pair(3,bndryPackTags[i].second));
      }


      // Right
      matchFound = 0;
      for (int j=0; j<outPoints.size(); j++) {
        for (int k=0; k<entityRight.size(); k++) {
          if (outPoints[j].second == entityRight[k].second) {
            matchFound++;
            break;
          }
        }

        if (matchFound > 0)
          break;
      }

      if (matchFound > 0) {
        rightShapes.push_back(std::make_pair(3,bndryPackTags[i].second));
      }


      // Up
      matchFound = 0;
      for (int j=0; j<outPoints.size(); j++) {
        for (int k=0; k<entityUp.size(); k++) {
          if (outPoints[j].second == entityUp[k].second) {
            matchFound++;
            break;
          }
        }

        if (matchFound > 0)
          break;
      }

      if (matchFound > 0) {
        UpShapes.push_back(std::make_pair(3,bndryPackTags[i].second));
      }


      // Down
      matchFound = 0;
      for (int j=0; j<outPoints.size(); j++) {
        for (int k=0; k<entityDown.size(); k++) {
          if (outPoints[j].second == entityDown[k].second) {
            matchFound++;
            break;
          }
        }

        if (matchFound > 0)
          break;
      }

      if (matchFound > 0) {
        DownShapes.push_back(std::make_pair(3,bndryPackTags[i].second));
      }


      // Back
      matchFound = 0;
      for (int j=0; j<outPoints.size(); j++) {
        for (int k=0; k<entityBack.size(); k++) {
          if (outPoints[j].second == entityBack[k].second) {
            matchFound++;
            break;
          }
        }

        if (matchFound > 0)
          break;
      }

      if (matchFound > 0) {
        backShapes.push_back(std::make_pair(3,bndryPackTags[i].second));
      }


      // Front
      matchFound = 0;
      for (int j=0; j<outPoints.size(); j++) {
        for (int k=0; k<entityFront.size(); k++) {
          if (outPoints[j].second == entityFront[k].second) {
            matchFound++;
            break;
          }
        }

        if (matchFound > 0)
          break;
      }

      if (matchFound > 0) {
        frontShapes.push_back(std::make_pair(3,bndryPackTags[i].second));
      }
    }

   std::ofstream fileInspect;
   fileInspect.open("Inspect.csv");

   fileInspect << "Left";
   for (int i=0; i<leftShapes.size(); i++)
     fileInspect << "," << leftShapes[i];
   fileInspect << std::endl;

   fileInspect << "Right";
   for (int i=0; i<rightShapes.size(); i++)
     fileInspect << "," << rightShapes[i];
   fileInspect << std::endl;

   fileInspect << "Up";
   for (int i=0; i<UpShapes.size(); i++)
     fileInspect << "," << UpShapes[i];
   fileInspect << std::endl;

   fileInspect << "Down";
   for (int i=0; i<DownShapes.size(); i++)
     fileInspect << "," << DownShapes[i];
   fileInspect << std::endl;

   fileInspect << "Back";
   for (int i=0; i<backShapes.size(); i++)
     fileInspect << "," << backShapes[i];
   fileInspect << std::endl;

   fileInspect << "Front";
   for (int i=0; i<frontShapes.size(); i++)
     fileInspect << "," << frontShapes[i];
   fileInspect << std::endl;

    fileInspect.close();

    //gmsh::model::occ::addBox(boxPt[0], boxPt[1], boxPt[2], Xdim, Ydim, Zdim);
    //gmsh::model::occ::synchronize();
    //gmsh::fltk::run();

    // Start Translating
    std::vector<int> indTra;
    indTra.resize(6);

    indTra.push_back(0); indTra.push_back(3); indTra.push_back(6);
    indTra.push_back(9); indTra.push_back(24); indTra.push_back(25);

    // Translation is not working, Need to think something.

    if (leftShapes.size() > 0) {
      // Left side (X+, X-, Y+, Y-, Z+, Z-)
      for (int h=0; h<indTra.size(); h++) {
        std::vector<std::pair<int, int>> tagsCopy;
        gmsh::model::occ::copy(leftShapes, tagsCopy);
        gmsh::model::occ::translate(tagsCopy, xTranslate[indTra[h]], yTranslate[indTra[h]],
                                    zTranslate[indTra[h]]);
      }
    }

    if (rightShapes.size() > 0) {
      // Right Side (X+, X-, Y+, Y-, Z+, Z-)
      for (int h=0; h<indTra.size(); h++) {
        std::vector<std::pair<int, int>> tagsCopy;
        gmsh::model::occ::copy(rightShapes, tagsCopy);
        gmsh::model::occ::translate(tagsCopy, xTranslate[indTra[h]], yTranslate[indTra[h]],
                                    zTranslate[indTra[h]]);
      }
    }

    if (UpShapes.size() > 0) {
      // Up (X+, X-, Y+, Y-, Z+, Z-)
      for (int h=0; h<indTra.size(); h++) {
        std::vector<std::pair<int, int>> tagsCopy;
        gmsh::model::occ::copy(UpShapes, tagsCopy);
        gmsh::model::occ::translate(tagsCopy, xTranslate[indTra[h]], yTranslate[indTra[h]],
                                    zTranslate[indTra[h]]);
      }
    }

    if (DownShapes.size() > 0) {
      // Down (X+, X-, Y+, Y-, Z+, Z-)
      for (int h=0; h<indTra.size(); h++) {
        std::vector<std::pair<int, int>> tagsCopy;
        gmsh::model::occ::copy(DownShapes, tagsCopy);
        gmsh::model::occ::translate(tagsCopy, xTranslate[indTra[h]], yTranslate[indTra[h]],
                                    zTranslate[indTra[h]]);
      }
    }

    if (backShapes.size() > 0) {
      // Back (X+, X-, Y+, Y-, Z+, Z-)
      for (int h=0; h<indTra.size(); h++) {
        std::vector<std::pair<int, int>> tagsCopy;
        gmsh::model::occ::copy(backShapes, tagsCopy);
        gmsh::model::occ::translate(tagsCopy, xTranslate[indTra[h]], yTranslate[indTra[h]],
                                    zTranslate[indTra[h]]);
      }
    }

    if (frontShapes.size() > 0) {
      // Front (X+, X-, Y+, Y-, Z+, Z-)
      for (int h=0; h<indTra.size(); h++) {
        std::vector<std::pair<int, int>> tagsCopy;
        gmsh::model::occ::copy(frontShapes, tagsCopy);
        gmsh::model::occ::translate(tagsCopy, xTranslate[indTra[h]], yTranslate[indTra[h]],
                                    zTranslate[indTra[h]]);
      }
    }
    gmsh::model::occ::synchronize();*/

    // Translating shapes
    std::vector<int> indexVols;
    std::vector<int> indexShapes;
    if (physGrpPerShape) {
      for (int g=0; g<bndryPackTags.size(); g++) {
        std::map<int, int>::iterator itrMp = storeShapeNames.begin();
        itrMp = storeShapeNames.find(bndryPackTags[g].second);

        if (itrMp != storeShapeNames.end()) {
          indexShapes.push_back(itrMp->second);
          indexVols.push_back(g);
        }
      }
    }

    std::vector<int> tmpVec;
    if ((enablePhysGrp) || (internalCohesiveBool)) {
      linkMultiPhysGrps.resize(26*bndryPackTags.size()+nameOfPcks.size());
      for (int i=0; i<nameOfPcks.size(); i++) {
        tmpVec.push_back(i+1);
        linkMultiPhysGrps[i] = i+1;
      }
    }

    int ptg = 1;
    std::cerr << "    Progress --> [0%";

    for (int h = 0; h < 26; h++) {
      std::vector<std::pair<int, int>> tagsCopy;
      gmsh::model::occ::copy(bndryPackTags, tagsCopy);
      gmsh::model::occ::translate(tagsCopy, xTranslate[h], yTranslate[h],
                                  zTranslate[h]);

      // Links shape type with volume number
      if (physGrpPerShape)
        for (int g=0; g<indexVols.size(); g++)
          storeShapeNames[tagsCopy[indexVols[g]].second] = indexShapes[g];

      if ((enablePhysGrp) || (internalCohesiveBool))
        for (int g=0; g<bndryPackTags.size(); g++)
          linkMultiPhysGrps[tagsCopy[g].second-1] = tmpVec[bndryPackTags[g].second-1];

      if (h % 3 == 0) {
        std::cerr.precision(3);
        std::cerr << ".." << 10.7142857 * (ptg) << "%";
        ptg++;
      }
    }
    gmsh::model::occ::synchronize();

    // Containers for boolean operation
    std::vector<std::pair<int, int>> tagsPacks;
    std::vector<std::pair<int, int>> outBoolean;
    std::vector<std::vector<std::pair<int, int>>> outBoolMap;
    gmsh::model::getEntities(tagsPacks, 3);
    std::vector<std::pair<int, int>> tagsBox;
    int tmpTag = gmsh::model::occ::addBox(boxPt[0], boxPt[1], boxPt[2], Xdim,
                                          Ydim, Zdim);
    tagsBox.push_back(std::make_pair(3, tmpTag));

    // Boolean Intersection
    gmsh::model::occ::intersect(tagsPacks, tagsBox, outBoolean, outBoolMap);
    std::cout << "..100%]" << std::endl;

    gmsh::model::occ::synchronize();

    if (physGrpPerShape) {
      for (int mp=0; mp<outBoolean.size(); mp++) {
        std::map<int, int>::iterator itrMp = storeShapeNames.begin();

        itrMp = storeShapeNames.find(outBoolean[mp].second);

        if (itrMp != storeShapeNames.end()) {
          if (itrMp->second == 0)
            spherePhysicalGroup.push_back(itrMp->first);
          if (itrMp->second == 1)
            ellipsoidPhysicalGroup.push_back(itrMp->first);
          if (itrMp->second == 2)
            cylindersPhysicalGroup.push_back(itrMp->first);
          if (itrMp->second == 3)
            hmxPhysicalGroup.push_back(itrMp->first);
          if (itrMp->second == 4)
            petnPhysicalGroup.push_back(itrMp->first);
          if (itrMp->second == 5)
            icosidodecahedronPhysicalGroup.push_back(itrMp->first);
        }
      }
    }

    if ((enablePhysGrp) || (internalCohesiveBool)) {
      // storeMultiPhysGrps
      for (int k=0; k<tmpVec.size(); k++) {
        std::vector<int> oneVols;
        for (int j=0; j<outBoolean.size(); j++) {
          if (tmpVec[k] == linkMultiPhysGrps[outBoolean[j].second-1]) {
            oneVols.push_back(outBoolean[j].second);
          }
        }
        storeMultiPhysGrps.push_back(oneVols);
      }
    }

    if (periodic3D == true) {
      std::cout << " - Mapping Periodic Surfaces" << std::endl;
      mapPeriodicSurfaces(outBoolean);
    }
  } else {
    // Nothing
  }
}

// Finds particular word in file and returns that line
std::string rocPack::findWord(const std::string &word) {
  std::ifstream input(InFile);
  std::string line;
  while (std::getline(input, line))
    if (line.find(word) == 0)
      return line;
  return "";
}

// Tokenized string in a string vector
std::vector<std::string> rocPack::Tokenize(const std::string &lineIn,
                                           const char &delim) {
  std::vector<std::string> rocTokens;

  std::string intermediate;
  std::stringstream checkStr(lineIn);

  // Tokenizing
  while (std::getline(checkStr, intermediate, delim)) {
    rocTokens.push_back(intermediate);
  }
  return rocTokens;
}

std::vector<std::string> rocPack::Tokenize(const std::string &lineIn,
                                           const std::string &delims) {
  std::vector<std::string> rocTokens;
  std::vector<char> delimStr;
  delimStr.reserve(delims.size());

  for (int i = 0; i < delims.size(); i++)
    delimStr.push_back(delims[i]);

  std::string intermediate;
  std::stringstream checkStr(lineIn);

  // Tokenizing
  int iter = 0;
  while (std::getline(checkStr, intermediate, delimStr[iter])) {
    rocTokens.push_back(intermediate);
    iter++;
  }
  return rocTokens;
}

std::vector<std::string>
rocPack::getShapeData(const int &iter, const std::string &a,
                      const std::vector<std::string> &L) {
  std::vector<std::string> LTokens = Tokenize(L[iter], a);
  return LTokens;
}

void rocPack::makeSphere(const int &n) {
  double sphereRad = 1 * scaleOfPack[n];
  int newVol =
      gmsh::model::occ::addSphere(translateParams[n][0], translateParams[n][1],
                                  translateParams[n][2], sphereRad);

  if (physGrpPerShape)
    storeShapeNames[newVol] = 0;

  gmsh::model::occ::synchronize();
}

void rocPack::makeEllipsoid(const int &n) {
  std::vector<int> pts;
  std::vector<int> arcs;
  std::vector<int> loops;
  std::vector<int> surfs;
  std::vector<int> tmpCurveLoop;
  std::vector<std::vector<double>> elVerts;
  std::vector<std::vector<double>> rotatedVerts;
  elVerts.resize(7);
  rotatedVerts.resize(7);
  pts.resize(7);
  arcs.resize(12);
  loops.resize(8);
  surfs.resize(8);

  // Scaling vertices
  elVerts[0].resize(3);
  elVerts[0][0] = 0;
  elVerts[0][1] = 0;
  elVerts[0][2] = 0;
  elVerts[1].resize(3);
  elVerts[1][0] = scaleOfPack[n] * ellipsoidRad[0];
  elVerts[1][1] = 0;
  elVerts[1][2] = 0;
  elVerts[2].resize(3);
  elVerts[2][0] = 0;
  elVerts[2][1] = scaleOfPack[n] * ellipsoidRad[1];
  elVerts[2][2] = 0;
  elVerts[3].resize(3);
  elVerts[3][0] = 0;
  elVerts[3][1] = 0;
  elVerts[3][2] = scaleOfPack[n] * ellipsoidRad[2];
  elVerts[4].resize(3);
  elVerts[4][0] = -scaleOfPack[n] * ellipsoidRad[0];
  elVerts[4][1] = 0;
  elVerts[4][2] = 0;
  elVerts[5].resize(3);
  elVerts[5][0] = 0;
  elVerts[5][1] = -scaleOfPack[n] * ellipsoidRad[1];
  elVerts[5][2] = 0;
  elVerts[6].resize(3);
  elVerts[6][0] = 0;
  elVerts[6][1] = 0;
  elVerts[6][2] = -scaleOfPack[n] * ellipsoidRad[2];

  // Rotating vertices
  rocQuaternion q = toQuaternion(rotateParams[n]);
  for (int i = 0; i < 7; i++)
    rotatedVerts[i] = rotateByQuaternion(q, elVerts[i]);

  // Adding translated points

  for (int i = 0; i < rotatedVerts.size(); i++)
    pts[i] =
        gmsh::model::geo::addPoint(rotatedVerts[i][0] + translateParams[n][0],
                                   rotatedVerts[i][1] + translateParams[n][1],
                                   rotatedVerts[i][2] + translateParams[n][2]);

  // Try circle arc here
  arcs[0] = gmsh::model::geo::addEllipseArc(pts[1], pts[0], pts[6], pts[6]);
  arcs[1] = gmsh::model::geo::addEllipseArc(pts[6], pts[0], pts[4], pts[4]);
  arcs[2] = gmsh::model::geo::addEllipseArc(pts[4], pts[0], pts[3], pts[3]);
  arcs[3] = gmsh::model::geo::addEllipseArc(pts[3], pts[0], pts[1], pts[1]);
  arcs[4] = gmsh::model::geo::addEllipseArc(pts[1], pts[0], pts[2], pts[2]);
  arcs[5] = gmsh::model::geo::addEllipseArc(pts[2], pts[0], pts[4], pts[4]);
  arcs[6] = gmsh::model::geo::addEllipseArc(pts[4], pts[0], pts[5], pts[5]);
  arcs[7] = gmsh::model::geo::addEllipseArc(pts[5], pts[0], pts[1], pts[1]);
  arcs[8] = gmsh::model::geo::addEllipseArc(pts[6], pts[0], pts[2], pts[2]);
  arcs[9] = gmsh::model::geo::addEllipseArc(pts[2], pts[0], pts[3], pts[3]);
  arcs[10] = gmsh::model::geo::addEllipseArc(pts[3], pts[0], pts[5], pts[5]);
  arcs[11] = gmsh::model::geo::addEllipseArc(pts[5], pts[0], pts[6], pts[6]);

  // Adding line loops
  tmpCurveLoop.push_back(arcs[4]);
  tmpCurveLoop.push_back(arcs[9]);
  tmpCurveLoop.push_back(arcs[3]);
  loops[0] = gmsh::model::geo::addCurveLoop(tmpCurveLoop);
  tmpCurveLoop.clear();
  tmpCurveLoop.push_back(arcs[8]);
  tmpCurveLoop.push_back(-arcs[4]);
  tmpCurveLoop.push_back(arcs[0]);
  loops[1] = gmsh::model::geo::addCurveLoop(tmpCurveLoop);
  tmpCurveLoop.clear();
  tmpCurveLoop.push_back(-arcs[9]);
  tmpCurveLoop.push_back(arcs[5]);
  tmpCurveLoop.push_back(arcs[2]);
  loops[2] = gmsh::model::geo::addCurveLoop(tmpCurveLoop);
  tmpCurveLoop.clear();
  tmpCurveLoop.push_back(-arcs[5]);
  tmpCurveLoop.push_back(-arcs[8]);
  tmpCurveLoop.push_back(arcs[1]);
  loops[3] = gmsh::model::geo::addCurveLoop(tmpCurveLoop);
  tmpCurveLoop.clear();
  tmpCurveLoop.push_back(arcs[7]);
  tmpCurveLoop.push_back(-arcs[3]);
  tmpCurveLoop.push_back(arcs[10]);
  loops[4] = gmsh::model::geo::addCurveLoop(tmpCurveLoop);
  tmpCurveLoop.clear();
  tmpCurveLoop.push_back(arcs[11]);
  tmpCurveLoop.push_back(-arcs[7]);
  tmpCurveLoop.push_back(-arcs[0]);
  loops[5] = gmsh::model::geo::addCurveLoop(tmpCurveLoop);
  tmpCurveLoop.clear();
  tmpCurveLoop.push_back(-arcs[10]);
  tmpCurveLoop.push_back(-arcs[2]);
  tmpCurveLoop.push_back(arcs[6]);
  loops[6] = gmsh::model::geo::addCurveLoop(tmpCurveLoop);
  tmpCurveLoop.clear();
  tmpCurveLoop.push_back(-arcs[1]);
  tmpCurveLoop.push_back(-arcs[6]);
  tmpCurveLoop.push_back(-arcs[11]);
  loops[7] = gmsh::model::geo::addCurveLoop(tmpCurveLoop);
  tmpCurveLoop.clear();

  // Adding compound surfaces
  tmpCurveLoop.push_back(loops[0]);
  surfs[0] = gmsh::model::geo::addSurfaceFilling(tmpCurveLoop);
  tmpCurveLoop.clear();
  tmpCurveLoop.push_back(loops[1]);
  surfs[1] = gmsh::model::geo::addSurfaceFilling(tmpCurveLoop);
  tmpCurveLoop.clear();
  tmpCurveLoop.push_back(loops[2]);
  surfs[2] = gmsh::model::geo::addSurfaceFilling(tmpCurveLoop);
  tmpCurveLoop.clear();
  tmpCurveLoop.push_back(loops[3]);
  surfs[3] = gmsh::model::geo::addSurfaceFilling(tmpCurveLoop);
  tmpCurveLoop.clear();
  tmpCurveLoop.push_back(loops[4]);
  surfs[4] = gmsh::model::geo::addSurfaceFilling(tmpCurveLoop);
  tmpCurveLoop.clear();
  tmpCurveLoop.push_back(loops[5]);
  surfs[5] = gmsh::model::geo::addSurfaceFilling(tmpCurveLoop);
  tmpCurveLoop.clear();
  tmpCurveLoop.push_back(loops[6]);
  surfs[6] = gmsh::model::geo::addSurfaceFilling(tmpCurveLoop);
  tmpCurveLoop.clear();
  tmpCurveLoop.push_back(loops[7]);
  surfs[7] = gmsh::model::geo::addSurfaceFilling(tmpCurveLoop);
  tmpCurveLoop.clear();

  // Surface loop
  tmpCurveLoop.push_back(surfs[0]);
  tmpCurveLoop.push_back(surfs[1]);
  tmpCurveLoop.push_back(surfs[2]);
  tmpCurveLoop.push_back(surfs[3]);
  tmpCurveLoop.push_back(surfs[4]);
  tmpCurveLoop.push_back(surfs[5]);
  tmpCurveLoop.push_back(surfs[6]);
  tmpCurveLoop.push_back(surfs[7]);
  int surfLoop = gmsh::model::geo::addSurfaceLoop(tmpCurveLoop);
  tmpCurveLoop.clear();

  // Volume
  tmpCurveLoop.push_back(surfLoop);
  int newVol = gmsh::model::geo::addVolume(tmpCurveLoop);
  gmsh::model::geo::synchronize();

  if (physGrpPerShape)
      storeShapeNames[newVol] = 1;

  gmsh::model::geo::synchronize();
}

void rocPack::makeCylinder(const int &n) {
  // "cylParams" is vector with 2 values in shape definition
  // cylParams[0] = Radius
  // cylParams[1] = Height

  // Scaling everything on thr fly
  // Cylinder Center
  std::vector<double> cylCenter = std::vector<double>(3);
  cylCenter[0] = 0;
  cylCenter[1] = 0 - (cylParams[1] / 2) * scaleOfPack[n];
  cylCenter[2] = 0;

  // Cylinder Translation Axis
  std::vector<double> translation_axis = std::vector<double>(3);
  translation_axis[0] = 0;
  translation_axis[1] = (cylParams[1]) * scaleOfPack[n];
  translation_axis[2] = 0;

  // Cylinder Vertices
  std::vector<std::vector<double>> cylVertices;
  cylVertices.resize(4);
  cylVertices[0].resize(3);
  cylVertices[1].resize(3);
  cylVertices[2].resize(3);
  cylVertices[3].resize(3);

  cylVertices[0][0] = 0;
  cylVertices[0][1] = 0 - (cylParams[1] / 2) * scaleOfPack[n];
  cylVertices[0][2] = 0 + cylParams[0] * scaleOfPack[n];
  cylVertices[1][0] = 0 + cylParams[0] * scaleOfPack[n];
  cylVertices[1][1] = 0 - (cylParams[1] / 2) * scaleOfPack[n];
  cylVertices[1][2] = 0;
  cylVertices[2][0] = 0;
  cylVertices[2][1] = 0 - (cylParams[1] / 2) * scaleOfPack[n];
  cylVertices[2][2] = 0 - cylParams[0] * scaleOfPack[n];
  cylVertices[3][0] = 0 - cylParams[0] * scaleOfPack[n];
  cylVertices[3][1] = 0 - (cylParams[1] / 2) * scaleOfPack[n];
  cylVertices[3][2] = 0;

  // Rotate By Quaternion
  rocQuaternion q = toQuaternion(rotateParams[n]);

  std::vector<double> rotatedCylCenter;
  std::vector<double> rotatedTranslationAxis;
  std::vector<std::vector<double>> rotatedCylVertices;
  rotatedCylVertices.resize(4);

  // Center
  rotatedCylCenter = rotateByQuaternion(q, cylCenter);

  // Translation_Axis
  rotatedTranslationAxis = rotateByQuaternion(q, translation_axis);

  // Vertices
  for (int i = 0; i < 4; i++)
    rotatedCylVertices[i] = rotateByQuaternion(q, cylVertices[i]);

  // Start adding geometry
  std::vector<int> pts = std::vector<int>(5);

  pts[0] =
      gmsh::model::occ::addPoint(rotatedCylCenter[0] + translateParams[n][0],
                                 rotatedCylCenter[1] + translateParams[n][1],
                                 rotatedCylCenter[2] + translateParams[n][2]);

  for (int i = 0; i < rotatedCylVertices.size(); i++)
    pts[i + 1] = gmsh::model::occ::addPoint(
        rotatedCylVertices[i][0] + translateParams[n][0],
        rotatedCylVertices[i][1] + translateParams[n][1],
        rotatedCylVertices[i][2] + translateParams[n][2]);

  std::vector<int> arcs = std::vector<int>(4);
  arcs[0] = gmsh::model::occ::addCircleArc(pts[1], pts[0], pts[2]);
  arcs[1] = gmsh::model::occ::addCircleArc(pts[2], pts[0], pts[3]);
  arcs[2] = gmsh::model::occ::addCircleArc(pts[3], pts[0], pts[4]);
  arcs[3] = gmsh::model::occ::addCircleArc(pts[4], pts[0], pts[1]);

  std::vector<int> curves = std::vector<int>(4);
  curves[0] = arcs[0];
  curves[1] = arcs[1];
  curves[2] = arcs[2];
  curves[3] = arcs[3];
  int loop = gmsh::model::occ::addCurveLoop(curves);

  std::vector<int> loops = std::vector<int>(1);
  loops[0] = loop;
  int surf = gmsh::model::occ::addPlaneSurface(loops);

  gmsh::model::occ::synchronize();
  std::vector<std::pair<int, int>> tagsSurfs;
  std::vector<std::pair<int, int>> outVols;
  tagsSurfs.push_back(std::make_pair(2, surf));

  gmsh::model::occ::extrude(tagsSurfs, rotatedTranslationAxis[0],
                            rotatedTranslationAxis[1],
                            rotatedTranslationAxis[2], outVols);

  for (int i=0; i<outVols.size(); i++)
    if (outVols[i].first == 3)
      if (physGrpPerShape)
        storeShapeNames[outVols[i].second] = 2;

  gmsh::model::occ::synchronize();
}

void rocPack::makeCrystalShape(const int &n, const int &index) {
  std::vector<std::vector<double>> scaledVerts;
  scaledVerts.resize(verts[index].size());

  // Scaling Vertices
  for (int i = 0; i < verts[index].size(); i++) {
    scaledVerts[i].resize(3);
    for (int j = 0; j < verts[index][i].size(); j++)
      scaledVerts[i][j] = verts[index][i][j] * scaleOfPack[n];
  }

  // Rotating Vertices using Quaternions
  rocQuaternion q = toQuaternion(rotateParams[n]);
  std::vector<std::vector<double>> rotatedVerts;
  rotatedVerts.resize(verts[index].size());

  for (int i = 0; i < verts[index].size(); i++)
    rotatedVerts[i] = rotateByQuaternion(q, scaledVerts[i]);

  // Translating Vertices
  for (int i = 0; i < rotatedVerts.size(); i++)
    for (int j = 0; j < rotatedVerts[i].size(); j++)
      rotatedVerts[i][j] = rotatedVerts[i][j] + translateParams[n][j];

  // Creating geometry
  std::vector<int> pts;
  pts.resize(rotatedVerts.size());
  for (int i = 0; i < rotatedVerts.size(); i++) {
    pts[i] = gmsh::model::occ::addPoint(rotatedVerts[i][0], rotatedVerts[i][1],
                                        rotatedVerts[i][2]);
  }

  std::vector<std::vector<int>> lines;
  lines.resize(faces[index].size());
  std::vector<int> surfs;
  surfs.resize(faces[index].size());
  std::vector<int> lineLoops;
  lineLoops.resize(faces[index].size());

  std::map<std::pair<int, int>, int> lineMap;
  // Adding lines and surfaces
  for (int i = 0; i < faces[index].size(); i++) {
    lines[i].resize(faces[index][i].size());
    for (int j = 0; j < faces[index][i].size(); j++) {
      // Finding for existing line
      std::map<std::pair<int, int>, int>::iterator it;
      int a = pts[faces[index][i][j]];
      int b = pts[faces[index][i][(j + 1) % faces[index][i].size()]];

      it = lineMap.find(std::make_pair(a, b));

      if (it == lineMap.end())
        it = lineMap.find(std::make_pair(b, a));

      if (it == lineMap.end()) {
        lines[i][j] = gmsh::model::occ::addLine(a, b);
        lineMap.insert(std::make_pair(std::make_pair(a, b), lines[i][j]));
      } else
        lines[i][j] = it->second;
    }

    // Adding Line Loop
    lineLoops[i] = gmsh::model::occ::addCurveLoop(lines[i]);
    // Adding Surface
    std::vector<int> tagsSurfFill = std::vector<int>(1);
    tagsSurfFill[0] = lineLoops[i];
    if (nameOfPcks[n] == "icosidodecahedron")
      surfs[i] = gmsh::model::occ::addSurfaceFilling(lineLoops[i]);
    else
      surfs[i] = gmsh::model::occ::addPlaneSurface(tagsSurfFill);
  }

  // Adding Surface Loop
  int surfLoop = gmsh::model::occ::addSurfaceLoop(surfs);

  // Adding Volumes
  std::vector<int> surfLoopTags = std::vector<int>(1);
  surfLoopTags[0] = surfLoop;
  int newVol = gmsh::model::occ::addVolume(surfLoopTags);


  if (physGrpPerShape) {
    if (nameOfPcks[n] == "hmx")
      storeShapeNames[newVol] = 3;

    if (nameOfPcks[n] == "petn")
      storeShapeNames[newVol] = 4;

    if (nameOfPcks[n] == "icosidodecahedron")
      storeShapeNames[newVol] = 5;
  }

  gmsh::model::occ::synchronize();
}

void rocPack::normalizeVerts() {
  // Normalize
  for (int k = 0; k < verts.size(); k++) {
    double maxNorm = 0.0;
    for (int i = 0; i < verts[k].size(); i++) {
      double norm = sqrt(pow(verts[k][i][0], 2) + pow(verts[k][i][1], 2) +
                         pow(verts[k][i][2], 2));
      if (norm > maxNorm)
        maxNorm = norm;
    }

    for (int i = 0; i < verts[k].size(); i++) {
      verts[k][i][0] = (verts[k][i][0] / maxNorm);
      verts[k][i][1] = (verts[k][i][1] / maxNorm);
      verts[k][i][2] = (verts[k][i][2] / maxNorm);
    }
  }
}

void rocPack::tagBoundaryPacks() {
  // Add method for removing boundary packs
  std::vector<std::pair<int, int>> tagsWithinBoundary;
  std::vector<std::pair<int, int>> tagsAll;

  std::vector<int> insidePacks;
  std::vector<int> allPacks;
  std::vector<int> bndryPackVols;

  gmsh::model::getEntitiesInBoundingBox(boxPt[0], boxPt[1], boxPt[2],
                                        boxPt[0] + Xdim, boxPt[1] + Ydim,
                                        boxPt[2] + Zdim, tagsWithinBoundary, 3);

  gmsh::model::getEntities(tagsAll, 3);

  for (auto iter : tagsWithinBoundary)
    insidePacks.push_back(iter.second);

  for (auto iter2 : tagsAll)
    allPacks.push_back(iter2.second);

  for (int i = 0; i < allPacks.size(); i++) {
    int ht = 0;
    for (int j = 0; j < insidePacks.size(); j++)
      if (allPacks[i] == insidePacks[j]) {
        ht++;
        break;
      }

    if (ht == 0)
      bndryPackVols.push_back(allPacks[i]);
  }

  if ((tagsWithinBoundary.size() == 0 && sntChk) && (removeBoundaryPacks)) {
    if (just2Physgrps) {
      std::cerr << "There are no volumes left inside the Box!" << std::endl;
      std::cerr << "Try increasing the domain size or pack density!"
                << std::endl;
      std::cerr << "Or disable physical group options" << std::endl;
      throw;
    }
  }

  // Converting volumes into pair
  for (int i = 0; i < bndryPackVols.size(); i++)
    bndryPackTags.push_back(std::make_pair(3, bndryPackVols[i]));
}

char *rocPack::strToChar(const std::string &strng) {
  char *tab = new char[strng.length() + 1];
  std::strcpy(tab, strng.c_str());

  return tab;
}

std::vector<double> rocPack::rotateByQuaternion(const rocQuaternion &q,
                                                const std::vector<double> &v) {
  Eigen::VectorXd InputVec(3);
  Eigen::VectorXd RotatedVec(3);
  std::vector<double> returnVec = std::vector<double>(3);

  InputVec[0] = v[0];
  InputVec[1] = v[1];
  InputVec[2] = v[2];

  Eigen::Quaternion<double> Quaternion(q.w, q.x, q.y, q.z);

  RotatedVec = Quaternion._transformVector(InputVec);

  returnVec[0] = RotatedVec[0];
  returnVec[1] = RotatedVec[1];
  returnVec[2] = RotatedVec[2];

  return returnVec;
}

rocQuaternion rocPack::toQuaternion(const std::vector<double> &r) {
  rocQuaternion q;                            // Declaring Quaternion
  double angle = r[3] * 3.141592653589 / 180; // Declaring Angle in Radians

  q.x = r[0] * sin(angle / 2);
  q.y = r[1] * sin(angle / 2);
  q.z = r[2] * sin(angle / 2);
  q.w = cos(angle / 2);

  return q;
}

void rocPack::mapPeriodicSurfaces(
    const std::vector<std::pair<int, int>> &prevTags) {

  // Boolean Fragment
  gmsh::model::occ::synchronize();
  std::vector<std::pair<int, int>> tagsBox2;
  int tmpTag2 =
      gmsh::model::occ::addBox(boxPt[0], boxPt[1], boxPt[2], Xdim, Ydim, Zdim);
  tagsBox2.push_back(std::make_pair(3, tmpTag2));
  std::vector<std::pair<int, int>> outBoolean2;
  std::vector<std::vector<std::pair<int, int>>> outBoolMap2;
  gmsh::model::occ::fragment(tagsBox2, prevTags, outBoolean2, outBoolMap2);
  gmsh::model::occ::synchronize();

  if (assignSidePatches) {
    // Get boundaries of box
    std::vector<std::pair<int, int>> tagstmpVOL;
    tagstmpVOL.push_back(std::make_pair(3, tmpTag2));
    std::vector<std::pair<int, int>> tagsAllsurfs;

    gmsh::model::getBoundary(tagstmpVOL, tagsAllsurfs);

    int totalSurfsinBox = tagsAllsurfs.size() - 1;

    // Last 6 surfaces
    std::vector<int> surfTagUp;
    surfTagUp.push_back(tagsAllsurfs[totalSurfsinBox - 2].second);
    int Up = gmsh::model::addPhysicalGroup(2, surfTagUp);
    gmsh::model::setPhysicalName(2, Up, "Up");

    std::vector<int> surfTagDown;
    surfTagDown.push_back(tagsAllsurfs[totalSurfsinBox - 4].second);
    int Down = gmsh::model::addPhysicalGroup(2, surfTagDown);
    gmsh::model::setPhysicalName(2, Down, "Down");

    std::vector<int> surfTagLeft;
    surfTagLeft.push_back(tagsAllsurfs[totalSurfsinBox - 5].second);
    int Left = gmsh::model::addPhysicalGroup(2, surfTagLeft);
    gmsh::model::setPhysicalName(2, Left, "Left");

    std::vector<int> surfTagRight;
    surfTagRight.push_back(tagsAllsurfs[totalSurfsinBox].second);
    int Right = gmsh::model::addPhysicalGroup(2, surfTagRight);
    gmsh::model::setPhysicalName(2, Right, "Right");

    std::vector<int> surfTagFront;
    surfTagFront.push_back(tagsAllsurfs[totalSurfsinBox - 3].second);
    int Front = gmsh::model::addPhysicalGroup(2, surfTagFront);
    gmsh::model::setPhysicalName(2, Front, "Front");

    std::vector<int> surfTagBack;
    surfTagBack.push_back(tagsAllsurfs[totalSurfsinBox - 1].second);
    int Back = gmsh::model::addPhysicalGroup(2, surfTagBack);
    gmsh::model::setPhysicalName(2, Back, "Back");
  }

  if (((enablePhysGrp) && (just2Physgrps)) || ((enablePhysGrp) && (physGrpPerShape))
     || ((just2Physgrps) && (physGrpPerShape))) {
    std::cerr << "Please select only one option for physical group"
              << std::endl;
    throw;
  }

  if (enablePhysGrp) {
    // Creating Physical Groups
    // Same volume on both periodic boundaries needs to be same physical group
    std::vector<int> vecTags;
    vecTags.push_back(tmpTag2);
    surroundingGrp = gmsh::model::addPhysicalGroup(3, vecTags);
    gmsh::model::setPhysicalName(3, surroundingGrp, "Surrounding");

    for (int i=0; i<storeMultiPhysGrps.size(); i++) {
      std::vector<int> newTags;
      std::string name = "Vol" + std::to_string(i);

      for (int j=0; j<storeMultiPhysGrps[i].size(); j++)
        newTags.push_back(storeMultiPhysGrps[i][j]);

      int newGrp = gmsh::model::addPhysicalGroup(3, newTags);
      multiGrpIndices.push_back(newGrp);
      gmsh::model::setPhysicalName(3, newGrp, name);
    }
  } else if (just2Physgrps) {
    // Creating Physical Groups
    std::vector<int> vecTags;
    vecTags.push_back(tmpTag2);
    surroundingGrp = gmsh::model::addPhysicalGroup(3, vecTags);
    gmsh::model::setPhysicalName(3, surroundingGrp, "Surrounding");

    std::vector<std::pair<int, int>> tagsAll;
    gmsh::model::getEntities(tagsAll, 3);
    std::vector<int> newTags;

    for (int i = 0; i < tagsAll.size(); i++) {
      if (tagsAll[i].second == tmpTag2) {
      } else {
        newTags.push_back(tagsAll[i].second);
      }
    }

    packGrp = gmsh::model::addPhysicalGroup(3, newTags);
    gmsh::model::setPhysicalName(3, packGrp, "MicroStructures");
  } else if (physGrpPerShape) {
    // Defining Physical Groups Per Shape
    std::vector<int> vecTags;
    vecTags.push_back(tmpTag2);
    surroundingGrp = gmsh::model::addPhysicalGroup(3, vecTags);
    gmsh::model::setPhysicalName(3, surroundingGrp, "Surrounding");

    // Spheres
    if (spherePhysicalGroup.size() > 0){
      int newGrp = gmsh::model::addPhysicalGroup(3, spherePhysicalGroup);
      multiGrpIndices.push_back(newGrp);
      gmsh::model::setPhysicalName(3, newGrp, "Spheres");
    }

    // Cylinders
    if (cylindersPhysicalGroup.size() > 0){
      int newGrp = gmsh::model::addPhysicalGroup(3, cylindersPhysicalGroup);
      multiGrpIndices.push_back(newGrp);
      gmsh::model::setPhysicalName(3, newGrp, "Cylinders");
    }

    // Ellipsoids
    if (ellipsoidPhysicalGroup.size() > 0){
      int newGrp = gmsh::model::addPhysicalGroup(3, ellipsoidPhysicalGroup);
      multiGrpIndices.push_back(newGrp);
      gmsh::model::setPhysicalName(3, newGrp, "Ellipsoids");
    }

    // PETN
    if (petnPhysicalGroup.size() > 0){
      int newGrp = gmsh::model::addPhysicalGroup(3, petnPhysicalGroup);
      multiGrpIndices.push_back(newGrp);
      gmsh::model::setPhysicalName(3, newGrp, "PETN");
    }

    // Icosidodecahedron
    if (icosidodecahedronPhysicalGroup.size() > 0){
      int newGrp = gmsh::model::addPhysicalGroup(3, icosidodecahedronPhysicalGroup);
      multiGrpIndices.push_back(newGrp);
      gmsh::model::setPhysicalName(3, newGrp, "Icosidodecahedron");
    }

    // HMX
    if (hmxPhysicalGroup.size() > 0){
      int newGrp = gmsh::model::addPhysicalGroup(3, hmxPhysicalGroup);
      multiGrpIndices.push_back(newGrp);
      gmsh::model::setPhysicalName(3, newGrp, "HMX");
    }

  } else {
    surroundingGrp = tmpTag2;

    std::vector<std::pair<int, int>> tagsAll;
    gmsh::model::getEntities(tagsAll, 3);

    for (int i = 0; i < tagsAll.size(); i++) {
      if (tagsAll[i].second == tmpTag2) {
      } else {
        multiGrpIndices.push_back(tagsAll[i].second);
      }
    }
  }

  // Finding surfaces at each side
  double tol = 0.001;

  // Left
  std::vector<std::pair<int, int>> entityLeft;
  gmsh::model::getEntitiesInBoundingBox(
      boxPt[0] - tol, boxPt[1] - tol, boxPt[2] - tol, boxPt[0] + tol,
      boxPt[1] + Ydim + tol, Zdim + boxPt[2] + tol, entityLeft, 2);

  // Right
  std::vector<std::pair<int, int>> entityRight;
  gmsh::model::getEntitiesInBoundingBox(Xdim + boxPt[0] - tol, boxPt[1] - tol,
                                        boxPt[2] - tol, Xdim + boxPt[0] + tol,
                                        Ydim + boxPt[1] + tol,
                                        Zdim + boxPt[2] + tol, entityRight, 2);

  // Up
  std::vector<std::pair<int, int>> entityUp;
  gmsh::model::getEntitiesInBoundingBox(boxPt[0] - tol, Ydim + boxPt[1] - tol,
                                        boxPt[2] - tol, Xdim + boxPt[0] + tol,
                                        Ydim + boxPt[1] + tol,
                                        Zdim + boxPt[2] + tol, entityUp, 2);

  // Down
  std::vector<std::pair<int, int>> entityDown;
  gmsh::model::getEntitiesInBoundingBox(
      boxPt[0] - tol, boxPt[1] - tol, boxPt[2] - tol, Xdim + boxPt[0] + tol,
      boxPt[1] + tol, Zdim + boxPt[2] + tol, entityDown, 2);

  // Back
  std::vector<std::pair<int, int>> entityBack;
  gmsh::model::getEntitiesInBoundingBox(
      boxPt[0] - tol, boxPt[1] - tol, boxPt[2] - tol, Xdim + boxPt[0] + tol,
      Ydim + boxPt[1] + tol, boxPt[2] + tol, entityBack, 2);

  // Front
  std::vector<std::pair<int, int>> entityFront;
  gmsh::model::getEntitiesInBoundingBox(
      boxPt[0] - tol, boxPt[1] - tol, Zdim + boxPt[2] - tol,
      Xdim + boxPt[0] + tol, Ydim + boxPt[1] + tol, Zdim + boxPt[2] + tol,
      entityFront, 2);

  // Getting points of individual surfaces
  std::vector<std::vector<std::pair<int, int>>> vertsLeft =
      getAllPoints(entityLeft);
  std::vector<std::vector<std::pair<int, int>>> vertsRight =
      getAllPoints(entityRight);
  std::vector<std::vector<std::pair<int, int>>> vertsUp =
      getAllPoints(entityUp);
  std::vector<std::vector<std::pair<int, int>>> vertsDown =
      getAllPoints(entityDown);
  std::vector<std::vector<std::pair<int, int>>> vertsFront =
      getAllPoints(entityFront);
  std::vector<std::vector<std::pair<int, int>>> vertsBack =
      getAllPoints(entityBack);

  // Getting Periodic Surfaces
  std::vector<std::pair<int, int>> periodicSurfsX;
  std::vector<std::pair<int, int>> periodicSurfsY;
  std::vector<std::pair<int, int>> periodicSurfsZ;

  periodicSurfsX = getPeriodicSurfs(vertsLeft, vertsRight, 0, Xdim);
  periodicSurfsY = getPeriodicSurfs(vertsDown, vertsUp, 1, Ydim);
  periodicSurfsZ = getPeriodicSurfs(vertsBack, vertsFront, 2, Zdim);

  // Enforcing Periodicity in GMSH model
  for (int i = 0; i < periodicSurfsX.size(); i++) {
    slaveX.push_back(periodicSurfsX[i].first);
    MasterX.push_back(periodicSurfsX[i].second);
  }

  for (int i = 0; i < periodicSurfsY.size(); i++) {
    slaveY.push_back(periodicSurfsY[i].first);
    MasterY.push_back(periodicSurfsY[i].second);
  }

  for (int i = 0; i < periodicSurfsZ.size(); i++) {
    slaveZ.push_back(periodicSurfsZ[i].first);
    MasterZ.push_back(periodicSurfsZ[i].second);
  }

  std::vector<std::vector<double>> affineTransform;
  affineTransform.resize(3);
  affineTransform[0].resize(16); // X
  affineTransform[1].resize(16); // Y
  affineTransform[2].resize(16); // Z
  int xT = -1 * Xdim;
  int yT = -1 * Ydim;
  int zT = -1 * Zdim;

  affineTransform[0][0] = 1;
  affineTransform[0][1] = 0;
  affineTransform[0][2] = 0;
  affineTransform[0][3] = xT;
  affineTransform[0][4] = 0;
  affineTransform[0][5] = 1;
  affineTransform[0][6] = 0;
  affineTransform[0][7] = 0;
  affineTransform[0][8] = 0;
  affineTransform[0][9] = 0;
  affineTransform[0][10] = 1;
  affineTransform[0][11] = 0;
  affineTransform[0][12] = 0;
  affineTransform[0][13] = 0;
  affineTransform[0][14] = 0;
  affineTransform[0][15] = 1;

  affineTransform[1][0] = 1;
  affineTransform[1][1] = 0;
  affineTransform[1][2] = 0;
  affineTransform[1][3] = 0;
  affineTransform[1][4] = 0;
  affineTransform[1][5] = 1;
  affineTransform[1][6] = 0;
  affineTransform[1][7] = yT;
  affineTransform[1][8] = 0;
  affineTransform[1][9] = 0;
  affineTransform[1][10] = 1;
  affineTransform[1][11] = 0;
  affineTransform[1][12] = 0;
  affineTransform[1][13] = 0;
  affineTransform[1][14] = 0;
  affineTransform[1][15] = 1;

  affineTransform[2][0] = 1;
  affineTransform[2][1] = 0;
  affineTransform[2][2] = 0;
  affineTransform[2][3] = 0;
  affineTransform[2][4] = 0;
  affineTransform[2][5] = 1;
  affineTransform[2][6] = 0;
  affineTransform[2][7] = 0;
  affineTransform[2][8] = 0;
  affineTransform[2][9] = 0;
  affineTransform[2][10] = 1;
  affineTransform[2][11] = zT;
  affineTransform[2][12] = 0;
  affineTransform[2][13] = 0;
  affineTransform[2][14] = 0;
  affineTransform[2][15] = 1;

  gmsh::model::mesh::setPeriodic(2, slaveX, MasterX, affineTransform[0]);
  gmsh::model::mesh::setPeriodic(2, slaveY, MasterY, affineTransform[1]);
  gmsh::model::mesh::setPeriodic(2, slaveZ, MasterZ, affineTransform[2]);
  gmsh::model::occ::synchronize();
}

std::vector<std::vector<std::pair<int, int>>> 
rocPack::getAllPoints(std::vector<std::pair<int, int>> surfaces) {
  std::vector<std::vector<std::pair<int, int>>> verts;
  verts.resize(surfaces.size());
  for (int j = 0; j < surfaces.size(); j++) {
    std::vector<std::pair<int, int>> inputTags;
    std::vector<std::pair<int, int>> outTags;
    inputTags.push_back(std::make_pair(2, surfaces[j].second));

    gmsh::model::getBoundary(inputTags, outTags, true, true, true);
    verts[j].resize(outTags.size());
    for (int i = 0; i < outTags.size(); i++)
      verts[j][i] = std::make_pair(surfaces[j].second, outTags[i].second);
  }

  return verts;
}

std::vector<std::pair<int, int>> rocPack::getPeriodicSurfs(
    const std::vector<std::vector<std::pair<int, int>>> &vertsOneSide,
    const std::vector<std::vector<std::pair<int, int>>> &vertsOtherSide,
    const int &indexTranslate, const double &amountTranslate) {
  // Comparing vertices and mapping periodic surfaces
  std::vector<double> paramCoords;
  std::vector<std::vector<double>> pointsOneSide;
  std::vector<std::pair<int, int>> periodicSurfs;

  for (int i = 0; i < vertsOneSide.size(); i++) {
    pointsOneSide.resize(vertsOneSide[i].size());

    for (int s = 0; s < vertsOneSide[i].size(); s++) {
      gmsh::model::getValue(0, vertsOneSide[i][s].second, paramCoords,
                            pointsOneSide[s]);
      pointsOneSide[s][indexTranslate] =
          pointsOneSide[s][indexTranslate] + amountTranslate;
    }

    for (int j = 0; j < vertsOtherSide.size(); j++) {
      int know = 0;
      for (int l = 0; l < vertsOtherSide[j].size(); l++) {
        std::vector<double> pointsOtherSide;
        gmsh::model::getValue(0, vertsOtherSide[j][l].second, paramCoords,
                              pointsOtherSide);

        for (int g = 0; g < pointsOneSide.size(); g++) {
          double a = pointsOtherSide[0] - pointsOneSide[g][0];
          double b = pointsOtherSide[1] - pointsOneSide[g][1];
          double c = pointsOtherSide[2] - pointsOneSide[g][2];

          if (std::sqrt(a * a) < 1e-10 &&
              pointsOneSide.size() == vertsOtherSide[j].size())
            if (std::sqrt(b * b) < 1e-10 &&
                pointsOneSide.size() == vertsOtherSide[j].size())
              if (std::sqrt(c * c) < 1e-10 &&
                  pointsOneSide.size() == vertsOtherSide[j].size())
                know++;
        }
      }
      if (know == pointsOneSide.size())
        periodicSurfs.push_back(std::make_pair(vertsOneSide[i][0].first,
                                               vertsOtherSide[j][0].first));
    }
  }
  return periodicSurfs;
}

void rocPack::writePeriodicNodes() {

  assignPeriodicEqNodes();

  // Change this method to write periodic equations.
  std::cout << " - Writing periodic equation file"
            << std::endl;

  // Getting all volumes
  std::vector<std::pair<int, int>> allVolumes;
  gmsh::model::getEntities(allVolumes, 3);

  // Creating pairs of <Volume,Surface> for node mapping
  std::vector<std::pair<int, int>> volSurfLinksX;
  std::vector<std::pair<int, int>> volSurfLinksY;
  std::vector<std::pair<int, int>> volSurfLinksZ;

  // A map for linking master surfaces with corrosponding volumes
  // <Surface, Volume>
  std::map<int, int> masterLinksX;
  std::map<int, int> masterLinksY;
  std::map<int, int> masterLinksZ;

  // Linking surfaces with volumes
  for (int i = 0; i < allVolumes.size(); i++) {
    std::vector<std::pair<int, int>> inputTags;
    std::vector<std::pair<int, int>> outTags;
    inputTags.push_back(std::make_pair(3, allVolumes[i].second));

    gmsh::model::getBoundary(inputTags, outTags, true, true, false);

    for (int j = 0; j < slaveX.size(); j++)
      for (int jj = 0; jj < outTags.size(); jj++)
        if (outTags[jj].first == 2 && outTags[jj].second == slaveX[j])
          volSurfLinksX.push_back(
              std::make_pair(allVolumes[i].second, slaveX[j]));

    for (int k = 0; k < slaveY.size(); k++)
      for (int kk = 0; kk < outTags.size(); kk++)
        if (outTags[kk].first == 2 && outTags[kk].second == slaveY[k])
          volSurfLinksY.push_back(
              std::make_pair(allVolumes[i].second, slaveY[k]));

    for (int l = 0; l < slaveZ.size(); l++)
      for (int ll = 0; ll < outTags.size(); ll++)
        if (outTags[ll].first == 2 && outTags[ll].second == slaveZ[l])
          volSurfLinksZ.push_back(
              std::make_pair(allVolumes[i].second, slaveZ[l]));

    for (int j = 0; j < MasterX.size(); j++)
      for (int jj = 0; jj < outTags.size(); jj++)
        if (outTags[jj].first == 2 && outTags[jj].second == MasterX[j])
          masterLinksX.insert(std::make_pair(MasterX[j], allVolumes[i].second));

    for (int k = 0; k < MasterY.size(); k++)
      for (int kk = 0; kk < outTags.size(); kk++)
        if (outTags[kk].first == 2 && outTags[kk].second == MasterY[k])
          masterLinksY.insert(std::make_pair(MasterY[k], allVolumes[i].second));

    for (int l = 0; l < MasterZ.size(); l++)
      for (int ll = 0; ll < outTags.size(); ll++)
        if (outTags[ll].first == 2 && outTags[ll].second == MasterZ[l])
          masterLinksZ.insert(std::make_pair(MasterZ[l], allVolumes[i].second));
  }

  // Gathering node data

  // Writing periodic.equ file (TODO: Do not include n0,nx,ny,nz as pairs in equations)
  std::ofstream periodicEquation;
  periodicEquation.open("periodic.equ");
  periodicEquation << "**set definitions" << std::endl;
  periodicEquation << "*nset, nset=n0" << std::endl;
  periodicEquation << eqRefNodes[0] << std::endl;
  periodicEquation << "*nset, nset=nx" << std::endl;
  periodicEquation << eqRefNodes[1] << std::endl;
  periodicEquation << "*nset, nset=ny" << std::endl;
  periodicEquation << eqRefNodes[2] << std::endl;
  periodicEquation << "*nset, nset=nz" << std::endl;
  periodicEquation << eqRefNodes[3] << std::endl;
  periodicEquation << "*equation" << std::endl;

  // X Direction (Master(Right) -> Slave(Left))
  int forSurfTestX = round(randomSurfTest * volSurfLinksX.size() / 100);
  int forSurfTestY = round(randomSurfTest * volSurfLinksY.size() / 100);
  int forSurfTestZ = round(randomSurfTest * volSurfLinksZ.size() / 100);
  for (int srf = 0; srf < volSurfLinksX.size(); srf++) {
    int masterSurfTag;
    std::vector<std::size_t> slaveNodes;
    std::vector<std::size_t> masterNodes;
    std::vector<double> affineTransformData;
    gmsh::model::mesh::getPeriodicNodes(2, volSurfLinksX[srf].second,
                                        masterSurfTag, slaveNodes, masterNodes,
                                        affineTransformData);

    int volSlave = volSurfLinksX[srf].first;
    int surfSlave = volSurfLinksX[srf].second;

    auto mapIterator = masterLinksX.find(masterSurfTag);
    int volMaster = mapIterator->second;
    int surfMaster = masterSurfTag;

    if (srf == forSurfTestX && masterNodes.size() != 0) {
      int indReturn = round(randomNodeX * masterNodes.size() / 100);
      std::vector<double> Mcoords;
      std::vector<double> Scoords;
      std::vector<double> paramCoords;

      gmsh::model::mesh::getNode(masterNodes[indReturn], Mcoords, paramCoords);
      gmsh::model::mesh::getNode(slaveNodes[indReturn], Scoords, paramCoords);

      if ((std::sqrt(Mcoords[0] * Mcoords[0]) -
           std::sqrt(Scoords[0] * Scoords[0])) < 1e-05)
        if ((std::sqrt(Mcoords[1] * Mcoords[1]) -
             std::sqrt(Scoords[1] * Scoords[1])) < 1e-05)
          if ((std::sqrt(Mcoords[2] * Mcoords[2]) -
               std::sqrt(Scoords[2] * Scoords[2])) < 1e-05)
            matchingCoordsX = true;
    }

    // Check slave nodes for any interface nodes and if found, add duplicate node
    // spawned from that interface node and link against a duplicate node spawned
    // from its opposite interface node.
    std::vector<int> newMasterNodes;
    std::vector<int> newSlaveNodes;
    if (internalCohesiveBool) {
      for (int i=0; i<masterNodes.size(); i++) {
        if (slaveInterfaceId[masterNodes[i]-1] == 1) {
          newMasterNodes.push_back(ptsReplacer[masterNodes[i]-1]);
          newSlaveNodes.push_back(ptsReplacer[slaveNodes[i]-1]);
        }
      }
      for (int i=0; i<newMasterNodes.size(); i++) {
        masterNodes.push_back(newMasterNodes[i]);
        slaveNodes.push_back(newSlaveNodes[i]);
      }
    }

    for (int i = 0; i < masterNodes.size(); i++) {

      if ((slaveNodes[i] != eqRefNodes[0]) && (slaveNodes[i] != eqRefNodes[1])
      &&  (slaveNodes[i] != eqRefNodes[2]) && (slaveNodes[i] != eqRefNodes[3])) {
        periodicEquation << "3" << std::endl;
        periodicEquation << slaveNodes[i] << ",1,-1," << masterNodes[i] << ",1,1,"
        << eqRefNodes[1] << ",1,1" << std::endl;

        periodicEquation << "3" << std::endl;
        periodicEquation << slaveNodes[i] << ",2,-1," << masterNodes[i] << ",2,1,"
        << eqRefNodes[1] << ",2,1" << std::endl;

        periodicEquation << "3" << std::endl;
        periodicEquation << slaveNodes[i] << ",3,-1," << masterNodes[i] << ",3,1,"
        << eqRefNodes[1] << ",3,1" << std::endl;
      } 
    }
  }

  // Y Direction (Master(Up) -> Slave(Down))
  for (int srf = 0; srf < volSurfLinksY.size(); srf++) {
    int masterSurfTag;
    std::vector<std::size_t> slaveNodes;
    std::vector<std::size_t> masterNodes;
    std::vector<double> affineTransformData;
    gmsh::model::mesh::getPeriodicNodes(2, volSurfLinksY[srf].second,
                                        masterSurfTag, slaveNodes, masterNodes,
                                        affineTransformData);

    int volSlave = volSurfLinksY[srf].first;
    int surfSlave = volSurfLinksY[srf].second;

    auto mapIterator = masterLinksY.find(masterSurfTag);
    int volMaster = mapIterator->second;
    int surfMaster = masterSurfTag;

    if (srf == forSurfTestY && masterNodes.size() != 0) {
      int indReturn = round(randomNodeY * masterNodes.size() / 100);
      std::vector<double> Mcoords;
      std::vector<double> Scoords;
      std::vector<double> paramCoords;

      gmsh::model::mesh::getNode(masterNodes[indReturn], Mcoords, paramCoords);
      gmsh::model::mesh::getNode(slaveNodes[indReturn], Scoords, paramCoords);

      if ((std::sqrt(Mcoords[0] * Mcoords[0]) -
           std::sqrt(Scoords[0] * Scoords[0])) < 1e-05)
        if ((std::sqrt(Mcoords[1] * Mcoords[1]) -
             std::sqrt(Scoords[1] * Scoords[1])) < 1e-05)
          if ((std::sqrt(Mcoords[2] * Mcoords[2]) -
               std::sqrt(Scoords[2] * Scoords[2])) < 1e-05)
            matchingCoordsY = true;
    }

    // Check slave nodes for any interface nodes and if found, add duplicate node
    // spawned from that interface node and link against a duplicate node spawned
    // from its opposite interface node.
    std::vector<int> newMasterNodes;
    std::vector<int> newSlaveNodes;
    if (internalCohesiveBool) {
      for (int i=0; i<masterNodes.size(); i++) {
        if (slaveInterfaceId[masterNodes[i]-1] == 1) {
          newMasterNodes.push_back(ptsReplacer[masterNodes[i]-1]);
          newSlaveNodes.push_back(ptsReplacer[slaveNodes[i]-1]);
        }
      }

      for (int i=0; i<newMasterNodes.size(); i++) {
        masterNodes.push_back(newMasterNodes[i]);
        slaveNodes.push_back(newSlaveNodes[i]);
      }
    }

    for (int i = 0; i < masterNodes.size(); i++) {
      if ((slaveNodes[i] != eqRefNodes[0]) && (slaveNodes[i] != eqRefNodes[1])
      &&  (slaveNodes[i] != eqRefNodes[2]) && (slaveNodes[i] != eqRefNodes[3])) {
        periodicEquation << "3" << std::endl;
        periodicEquation << slaveNodes[i] << ",1,-1," << masterNodes[i] << ",1,1,"
        << eqRefNodes[2] << ",1,1" << std::endl;

        periodicEquation << "3" << std::endl;
        periodicEquation << slaveNodes[i] << ",2,-1," << masterNodes[i] << ",2,1,"
        << eqRefNodes[2] << ",2,1" << std::endl;

        periodicEquation << "3" << std::endl;
        periodicEquation << slaveNodes[i] << ",3,-1," << masterNodes[i] << ",3,1,"
        << eqRefNodes[2] << ",3,1" << std::endl;
      }
    }
  }

  // Z Direction (Master(Front) -> Slave(Back))
  for (int srf = 0; srf < volSurfLinksZ.size(); srf++) {
    int masterSurfTag;
    std::vector<std::size_t> slaveNodes;
    std::vector<std::size_t> masterNodes;
    std::vector<double> affineTransformData;
    gmsh::model::mesh::getPeriodicNodes(2, volSurfLinksZ[srf].second,
                                        masterSurfTag, slaveNodes, masterNodes,
                                        affineTransformData);

    int volSlave = volSurfLinksZ[srf].first;
    int surfSlave = volSurfLinksZ[srf].second;

    auto mapIterator = masterLinksZ.find(masterSurfTag);
    int volMaster = mapIterator->second;
    int surfMaster = masterSurfTag;

    if (srf == forSurfTestZ && masterNodes.size() != 0) {
      int indReturn = round(randomNodeZ * masterNodes.size() / 100);
      std::vector<double> Mcoords;
      std::vector<double> Scoords;
      std::vector<double> paramCoords;

      gmsh::model::mesh::getNode(masterNodes[indReturn], Mcoords, paramCoords);
      gmsh::model::mesh::getNode(slaveNodes[indReturn], Scoords, paramCoords);

      if ((std::sqrt(Mcoords[0] * Mcoords[0]) -
           std::sqrt(Scoords[0] * Scoords[0])) < 1e-05)
        if ((std::sqrt(Mcoords[1] * Mcoords[1]) -
             std::sqrt(Scoords[1] * Scoords[1])) < 1e-05)
          if ((std::sqrt(Mcoords[2] * Mcoords[2]) -
               std::sqrt(Scoords[2] * Scoords[2])) < 1e-05)
            matchingCoordsZ = true;
    }

    // Check slave nodes for any interface nodes and if found, add duplicate node
    // spawned from that interface node and link against a duplicate node spawned
    // from its opposite interface node.
    std::vector<int> newMasterNodes;
    std::vector<int> newSlaveNodes;
    if (internalCohesiveBool) {
      for (int i=0; i<masterNodes.size(); i++) {
        if (slaveInterfaceId[masterNodes[i]-1] == 1) {
          newMasterNodes.push_back(ptsReplacer[masterNodes[i]-1]);
          newSlaveNodes.push_back(ptsReplacer[slaveNodes[i]-1]);
        }
      }

      for (int i=0; i<newMasterNodes.size(); i++) {
        masterNodes.push_back(newMasterNodes[i]);
        slaveNodes.push_back(newSlaveNodes[i]);
      }
    }

    for (int i = 0; i < masterNodes.size(); i++) {
      if ((slaveNodes[i] != eqRefNodes[0]) && (slaveNodes[i] != eqRefNodes[1])
      &&  (slaveNodes[i] != eqRefNodes[2]) && (slaveNodes[i] != eqRefNodes[3])) {
        periodicEquation << "3" << std::endl;
        periodicEquation << slaveNodes[i] << ",1,-1," << masterNodes[i] << ",1,1,"
        << eqRefNodes[3] << ",1,1" << std::endl;

        periodicEquation << "3" << std::endl;
        periodicEquation << slaveNodes[i] << ",2,-1," << masterNodes[i] << ",2,1,"
        << eqRefNodes[3] << ",2,1" << std::endl;

        periodicEquation << "3" << std::endl;
        periodicEquation << slaveNodes[i] << ",3,-1," << masterNodes[i] << ",3,1,"
        << eqRefNodes[3] << ",3,1" << std::endl;
      }
    }
  }
  periodicEquation.close();
}

void rocPack::setNodeLocations(const int &x, const int &y, const int &z) {
  randomNodeX = x;
  randomNodeY = y;
  randomNodeZ = z;
}

void rocPack::setRandomSurface(const int &surf) { randomSurfTest = surf; }

bool rocPack::getTestResult() {
  if ((matchingCoordsX) && (matchingCoordsY) && (matchingCoordsZ))
    return true;
  else
    return false;
}

void rocPack::scaleVols(const int &vol, const int &index) {

  std::vector<std::pair<int, int>> tagsIndividual;
  tagsIndividual.push_back(std::make_pair(3, vol));
  gmsh::model::occ::dilate(tagsIndividual, translateParams[index][0],
                           translateParams[index][1], translateParams[index][2],
                           shrinkScale, shrinkScale, shrinkScale);
}

void rocPack::performSmoothing() {
  std::vector<std::pair<int, int>> tagsSurfs;
  gmsh::model::getEntities(tagsSurfs, 2);

  for (int i = 0; i < tagsSurfs.size(); i++) {
    gmsh::model::mesh::setSmoothing(tagsSurfs[i].first, tagsSurfs[i].second,
                                    smoothingIter);
  }
  gmsh::model::occ::synchronize();
}

void rocPack::createCohesiveElements(const std::string &filename,
                                     const std::string &outname) {

  // Pre-processing to get some important data
  if ((sntChk) && (just2Physgrps)) {
    gmsh::model::mesh::getNodesForPhysicalGroup(3, surroundingGrp, surrNodeTags,
                                                surrCoords);
    gmsh::model::mesh::getNodesForPhysicalGroup(3, packGrp, geomsNodeTags,
                                                geomsCoords);
  } else if ((sntChk) && (enablePhysGrp)) {
    gmsh::model::mesh::getNodesForPhysicalGroup(3, surroundingGrp, surrNodeTags,
                                                surrCoords);

    for (int i=0; i<multiGrpIndices.size(); i++) {
      std::vector<std::size_t> tmpNodeTags;
      gmsh::model::mesh::getNodesForPhysicalGroup(3, multiGrpIndices[i], tmpNodeTags,
                                                  geomsCoords);
      for (int j=0; j<tmpNodeTags.size(); j++)
        geomsNodeTags.push_back(tmpNodeTags[j]);
    }
  } else if ((sntChk) && (physGrpPerShape)) {
    gmsh::model::mesh::getNodesForPhysicalGroup(3, surroundingGrp, surrNodeTags,
                                                surrCoords);

    for (int i=0; i<multiGrpIndices.size(); i++) {
      std::vector<std::size_t> tmpNodeTags;
      gmsh::model::mesh::getNodesForPhysicalGroup(3, multiGrpIndices[i], tmpNodeTags,
                                                  geomsCoords);

      for (int j=0; j<tmpNodeTags.size(); j++)
        geomsNodeTags.push_back(tmpNodeTags[j]);
    }
  } else if (((sntChk) && !(enablePhysGrp)) && ((sntChk) && !(physGrpPerShape)) 
         && ((sntChk) && !(just2Physgrps))) {
    std::vector<double> noUse;
    gmsh::model::mesh::getNodes(surrNodeTags,surrCoords,noUse,3,surroundingGrp,true,false);

    for (int i=0; i<multiGrpIndices.size(); i++) {
      std::vector<std::size_t> tmpNodes;
      gmsh::model::mesh::getNodes(tmpNodes,geomsCoords,noUse,3,multiGrpIndices[i],true,false);

      for (int j=0; j<tmpNodes.size(); j++)
        geomsNodeTags.push_back(tmpNodes[j]);
    }

    for (int k=0; k<multiGrpIndices.size(); k++) {
      // Getting element tags in surrounding -> surrElementIds
      std::vector<int> elementTypes;
      std::vector<std::vector<std::size_t>> elementTags;
      std::vector<std::vector<std::size_t>> nodeTags;

      gmsh::model::mesh::getElements(elementTypes,elementTags,nodeTags,3,multiGrpIndices[k]);

      for (int i=0; i<elementTags.size(); i++)
        for (int j=0; j<elementTags[i].size(); j++)
        geomElementIds.push_back(elementTags[0][j]);
    }
  } else {
    // Nothing
  }

  // Reading in mesh file
  size_t lastindex = filename.find_last_of(".");
  std::string rawname = filename.substr(0, lastindex);
  rawname = rawname + "_oldMSH.msh";
  meshBase *mb = meshBase::Create(rawname);

  // Changing output name for cohesive elements
  size_t lastindex2 = outname.find_last_of(".");
  std::string rawname2 = outname.substr(0, lastindex2);
  rawname2 = rawname2 + "_withCohesive.vtu";

  // Creating vtk data set for mesh
  vtkSmartPointer<vtkDataSet> dataSetSurr;
  dataSetSurr = mb->getDataSet();

  int numCells = dataSetSurr->GetNumberOfCells();
  int numPoints = dataSetSurr->GetNumberOfPoints();

  // Material assignment vector here
  for (int i=0; i<numPoints; i++)
    ptsCohesiveGrp.push_back(-1);

  // storeMultiPhysGrps is vector of vectors, storing volumes in each group
  int countParts = 100;
  for (int i=0; i<storeMultiPhysGrps.size(); i++) {
    std::vector<std::size_t> storeNodes;
    for (int j=0; j<storeMultiPhysGrps[i].size(); j++) {
      std::vector<std::size_t> getNodeTags;
      std::vector<double> notUseful;
      std::vector<double> notUseful2;
      gmsh::model::mesh::getNodes(getNodeTags,notUseful,notUseful2,3,
                                  storeMultiPhysGrps[i][j],true,false);

      for (int k=0; k<getNodeTags.size(); k++)
        storeNodes.push_back(getNodeTags[k]);
    }

    for (int j=0; j<storeNodes.size(); j++) {
      //if (ptsCohesiveGrp[storeNodes[j]-1] != -1) {
      //  std::cerr << "Exception at material assignement" << std::endl;
      //  throw;
      //} else {
        ptsCohesiveGrp[storeNodes[j]-1] = countParts;
      //}
    }
    countParts++;
  }

  // Fatching physical groups information
  std::vector<double> grpIds(mb->getNumberOfCells(), surroundingGrp);
  if ((just2Physgrps) || (enablePhysGrp) || (physGrpPerShape)) {
    mb->getCellDataArray("PhysGrpId", grpIds);
  } else {
    for (int i=0; i<geomElementIds.size(); i++)
      grpIds[i] = 2;
  }

  std::vector<int> newPtIds;

  // Transfer old dataset to new one with new points
  vtkSmartPointer<vtkUnstructuredGrid> dataSetCohesive =
      vtkSmartPointer<vtkUnstructuredGrid>::New();

  vtkSmartPointer<vtkUnstructuredGrid> dataSetViz =
      vtkSmartPointer<vtkUnstructuredGrid>::New();

  vtkSmartPointer<vtkPoints> pointsViz = vtkSmartPointer<vtkPoints>::New();

  // Provides means of identifying if the node at given index is new duplicate.
  std::vector<int> ptsInterfaceID;
  // Replaces duplicate nodes at their indexes with their surrounding counterparts
  std::vector<int> traceBackToSurr;
  // new vector
  std::vector<int> identifyInterfaceNodes;

  for (int i = 0; i < numPoints; i++) {
    std::vector<double> getPt = std::vector<double>(3);
    dataSetSurr->GetPoint(i, &getPt[0]);
    pointsViz->InsertNextPoint(getPt[0], getPt[1], getPt[2]);
    identifyInterfaceNodes.push_back(0);
    ptsReplacer.push_back(i);
    ptsInterfaceID.push_back(0);
    traceBackToSurr.push_back(i);
    slaveInterfaceId.push_back(0);
  }

  std::clock_t start;
  double duration;
  start = std::clock();

  // A faster method for finding interface nodes (0.031592 Seconds vs 813 Seconds for 61782 nodes)
  for (int i=0; i<surrNodeTags.size(); i++)
    identifyInterfaceNodes[surrNodeTags[i]-1] = 1;

  for (int i=0; i<geomsNodeTags.size(); i++) {
    if (identifyInterfaceNodes[geomsNodeTags[i]-1] == 1) {
      interfaceNodes.push_back(geomsNodeTags[i]-1);
      std::vector<double> getPt = std::vector<double>(3);
      dataSetSurr->GetPoint(geomsNodeTags[i]-1, &getPt[0]);
      int newP = pointsViz->InsertNextPoint(getPt[0], getPt[1], getPt[2]);
      newPtIds.push_back(newP);
      ptsReplacer[geomsNodeTags[i]-1] = newP;
      traceBackToSurr.push_back(geomsNodeTags[i]-1);
      slaveInterfaceId[geomsNodeTags[i]-1] = 1;
    }
  }

  duration = (std::clock() - start ) / (double) CLOCKS_PER_SEC;
  std::cout << "Process Finished!" << std::endl;
  std::cout << "Total " << interfaceNodes.size() 
            << " Nodes duplicated in " << duration 
            << " seconds!" << std::endl;

  dataSetCohesive->SetPoints(pointsViz);
  dataSetViz->SetPoints(pointsViz);

  for (int i = 0; i < interfaceNodes.size(); i++)
    ptsInterfaceID.push_back(1); // Updating points identification vector

  // Once we have interface nodes and duplicated the nodes, Surrounding side
  // of surface cells needs to be replaced by new nodes.
  // Find all cells containing interface nodes and rule out those in physical 
  // group surrounding.
  std::vector<int> cellIdentification;

  // Initialize with zero
  for (int i=0; i<numCells; i++)
    cellIdentification.push_back(0);

  // Finding all cells in microstructure group that atleast contains one interface
  // node.
  std::vector<int> pickCells;
  for (int i=0; i<interfaceNodes.size(); i++) {
    vtkIdList *cells = vtkIdList::New();
    dataSetSurr->GetPointCells(interfaceNodes[i], cells);
    for (int j=0; j<cells->GetNumberOfIds(); j++) {
      if (grpIds[cells->GetId(j)] != surroundingGrp) {
        cellIdentification[cells->GetId(j)] = 1;
        pickCells.push_back(cells->GetId(j));
      }
    }
  }
  sort(pickCells.begin(), pickCells.end());
  pickCells.erase(unique(pickCells.begin(),pickCells.end()),pickCells.end());

  // Now selected cells contain anywhere between 1 to 3 intefaceNodes that will
  // need replacement.
  // Make a mesh and keep everything same but replace interface
  // nodes with duplicate nodes when the selected cells are encountered.
  for (int i = 0; i < numCells; i++) {
    vtkCell *cell;
    vtkPoints *fp;
    if (cellIdentification[i] == 1) {
      cell = dataSetSurr->GetCell(i);
      vtkIdList *pts = cell->GetPointIds();
      std::vector<int> ptForCells;
      ptForCells.resize(4);
      for (int j = 0; j < 4; j++)
        ptForCells[j] = ptsReplacer[pts->GetId(j)];
      createVtkCell(dataSetCohesive, 10, ptForCells);
    } else {
      cell = dataSetSurr->GetCell(i);
      vtkIdList *pts = cell->GetPointIds();
      std::vector<int> ptForCells;
      ptForCells.resize(4);
      for (int j = 0; j < 4; j++)
        ptForCells[j] = pts->GetId(j);
      createVtkCell(dataSetCohesive, 10, ptForCells);
    }
  }

  // Mesh separated successfully!
  // Finding out interface faces within picked cells to reveal node order.
  std::vector<int> seqNodeMap;
  std::vector<int> cellFinal;
  std::vector<int> cellsNotPicked;
  int numOfCohElm = 0;

  for (int i=0; i<pickCells.size(); i++) {
    vtkCell *cell;
    cell = dataSetCohesive->GetCell(pickCells[i]);
    vtkIdList *pts = cell->GetPointIds();
    vtkCell3D *cell3d = static_cast<vtkCell3D *>(cell);

    int atleastOne = 0;
    // Looping through all faces of current cell
    for (int j = 0; j < 4; j++) {
      std::vector<int> tmpSeqStore;
      int *ptFaces = nullptr;
      cell3d->GetFacePoints(j, ptFaces);
      for (int k=0; k<3; k++)
        if (ptsInterfaceID[pts->GetId(ptFaces[k])] == 1)
          tmpSeqStore.push_back(pts->GetId(ptFaces[k]));

      if (tmpSeqStore.size() == 3)
        for (int n = 0; n < 3; n++)
          seqNodeMap.push_back(tmpSeqStore[n]);
    }
  }

  // Update the material group array here.
  std::vector<double> groupId;
  for (int i=0; i<grpIds.size(); i++)
    groupId.push_back(grpIds[i]);

  // Create cohesive elements
  int add = 0;
  for (int i = 0; i < (seqNodeMap.size()/3); i++) {
    std::vector<int> ptForCells;
    ptForCells.resize(6);
    for (int j = 0; j < 3 ; j++) {
      ptForCells[j] = seqNodeMap[j+add];
      ptForCells[j+3] = traceBackToSurr[seqNodeMap[j+add]];
    }

    // Search for periodic cohesive elements (existing on periodic boundaries)
    std::vector<double> getPt1 = std::vector<double>(3);
    std::vector<double> getPt2 = std::vector<double>(3);
    std::vector<double> getPt3 = std::vector<double>(3);
    dataSetCohesive->GetPoint(ptForCells[0], &getPt1[0]);
    dataSetCohesive->GetPoint(ptForCells[1], &getPt2[0]);
    dataSetCohesive->GetPoint(ptForCells[2], &getPt3[0]);
    if (existsOnPeriodicBoundary(getPt1,getPt2,getPt3)) {
      // Do not create that element
    } else {
      if (ptsCohesiveGrp[ptForCells[3]] == ptsCohesiveGrp[ptForCells[4]])
        if (ptsCohesiveGrp[ptForCells[4]] == ptsCohesiveGrp[ptForCells[5]])
          if (ptsCohesiveGrp[ptForCells[3]] == ptsCohesiveGrp[ptForCells[5]])
            groupId.push_back(ptsCohesiveGrp[ptForCells[3]]);

      createVtkCell(dataSetViz, 13, ptForCells);
      createVtkCell(dataSetCohesive, 13, ptForCells);
    }
    add = add + 3;
  }

  /*// Write cohesive elements in gmsh mesh
  std::vector<int> elementTypes;
  elementTypes.push_back(6);
  std::vector<std::vector<std::size_t>> elementTags;
  elementTags.resize(1);
  std::vector<std::vector<std::size_t>> nodeTags;
  nodeTags.resize(1);

  int add2 = 0;
  for (int i = 0; i < (seqNodeMap.size()/3); i++) {
   std::vector<int> ptForCells;
   ptForCells.resize(6);
   for (int j = 0; j < 3 ; j++) {
     ptForCells[j] = seqNodeMap[j+add2];
     ptForCells[j+3] = traceBackToSurr[seqNodeMap[j+add2]];
   }

   elementTags[0].push_back(numCells+i+1);

   for (int j=0; j<6; j++)
    nodeTags[0].push_back(ptForCells[j]+1);
   add2 = add2 + 3;
  }

  gmsh::model::mesh::addElements(3,surroundingGrp,elementTypes,elementTags,nodeTags);

  gmsh::write("GMSH_WITH_COHESIVE.msh");*/

  if (groupId.size() != dataSetCohesive->GetNumberOfCells()) {
    std::cerr << "Exception at group id vector size" << std::endl;
    throw;
  }

  vtkMesh *vm = new vtkMesh(dataSetViz, "CohesiveElements.vtu");
  vm->report();
  vm->write();

  vtkMesh *vm2 = new vtkMesh(dataSetCohesive, rawname2);
  vm2->setCellDataArray("PhysGrpId",groupId);
  vm2->report();
  vm2->write();
}

void rocPack::createVtkCell(vtkSmartPointer<vtkUnstructuredGrid> dataSet,
                            const int cellType, std::vector<int> &vrtIds) {
  vtkSmartPointer<vtkIdList> vtkCellIds = vtkSmartPointer<vtkIdList>::New();
  vtkCellIds->SetNumberOfIds(vrtIds.size());
  for (auto pit = vrtIds.begin(); pit != vrtIds.end(); pit++)
    vtkCellIds->SetId(pit - vrtIds.begin(), *pit);
  dataSet->InsertNextCell(cellType, vtkCellIds);
}

void rocPack::translateAll(const double &X, const double &Y, const double &Z) {
  xUDF = X;
  yUDF = Y;
  zUDF = Z;
}

void rocPack::setCustomDomain(const std::vector<double> &domainBounds) {
  cstmDomain = true;

  boxPt.clear();
  boxPt.resize(3);

  boxPt[0] = domainBounds[0] + xUDF;
  boxPt[1] = domainBounds[1] + yUDF;
  boxPt[2] = domainBounds[2] + zUDF;

  Xdim = domainBounds[3];
  Ydim = domainBounds[4];
  Zdim = domainBounds[5];
}

void rocPack::setMeshingAlgorithm(const int &mshAlg) {
  meshingAlgorithm = mshAlg;
}

void rocPack::enableDefOuts() { defOutputs = true; }

void rocPack::sanityCheckOn() { sntChk = true; }

void rocPack::enablePhysicalGroupsPerShape() { physGrpPerShape = true; }

void rocPack::assignPeriodicEqNodes() {
  std::size_t cellTag;
  int cellType = 0;
  std::vector<std::size_t> nodeIds;
  double u,v,w;
  int elemDim = -1;

  gmsh::model::mesh::getElementByCoordinates(boxPt[0],boxPt[1],boxPt[2],cellTag,cellType,
                                              nodeIds,u,v,w,elemDim,true);

  for (int i=0; i<nodeIds.size(); i++) {
    std::vector<double> coords;
    std::vector<double> paramCoords;

    gmsh::model::mesh::getNode(nodeIds[i],coords,paramCoords);
    
    if (coords[0] == boxPt[0])
      if (coords[1] == boxPt[1])
        if (coords[2] == boxPt[2])
          eqRefNodes[0] = nodeIds[i];
  }

  nodeIds.clear();

  gmsh::model::mesh::getElementByCoordinates(boxPt[0]+Xdim,boxPt[1],boxPt[2],cellTag,cellType,
                                              nodeIds,u,v,w,elemDim,true);

  for (int i=0; i<nodeIds.size(); i++) {
    std::vector<double> coords;
    std::vector<double> paramCoords;

    gmsh::model::mesh::getNode(nodeIds[i],coords,paramCoords);
    
    if (coords[0] == boxPt[0]+Xdim) {
      if (coords[1] == boxPt[1]) {
        if (coords[2] == boxPt[2]) {
          eqRefNodes[1] = nodeIds[i];
          break;
        }
      }
    }
  }
  nodeIds.clear();

  gmsh::model::mesh::getElementByCoordinates(boxPt[0],boxPt[1]+Ydim,boxPt[2],cellTag,cellType,
                                              nodeIds,u,v,w,elemDim,true);
  for (int i=0; i<nodeIds.size(); i++) {
    std::vector<double> coords;
    std::vector<double> paramCoords;

    gmsh::model::mesh::getNode(nodeIds[i],coords,paramCoords);
    
    if (coords[0] == boxPt[0]) {
      if (coords[1] == boxPt[1]+Ydim) {
        if (coords[2] == boxPt[2]) {
          eqRefNodes[2] = nodeIds[i];
          break;
        }
      }
    }
  }
  nodeIds.clear();

  gmsh::model::mesh::getElementByCoordinates(boxPt[1],boxPt[2],boxPt[2]+Zdim,cellTag,cellType,
                                              nodeIds,u,v,w,elemDim,true);
  for (int i=0; i<nodeIds.size(); i++) {
    std::vector<double> coords;
    std::vector<double> paramCoords;

    gmsh::model::mesh::getNode(nodeIds[i],coords,paramCoords);
    
    if (coords[0] == boxPt[0]) {
      if (coords[1] == boxPt[1]) {
        if (coords[2] == boxPt[2]+Zdim) {
          eqRefNodes[3] = nodeIds[i];
          break;
        }
      }
    }
  }
  nodeIds.clear();
}

bool rocPack::existsOnPeriodicBoundary(const std::vector<double> &getPt1,
                                       const std::vector<double> &getPt2,
                                       const std::vector<double> &getPt3) {
  // Check if all 3 points exist on any of these 6 surfaces
  // 1. X = boxPt[0]
  if ((getPt1[0] == boxPt[0]) && (getPt2[0] == boxPt[0]) && (getPt3[0] == boxPt[0]))
    return true;
  // 2. X = boxPt[0] + Xdim
  else if ((getPt1[0] == boxPt[0] + Xdim) && (getPt2[0] == boxPt[0] + Xdim) && (getPt3[0] == boxPt[0] + Xdim))
    return true;
  // 3. Y = boxPt[1]
  else if ((getPt1[1] == boxPt[1]) && (getPt2[1] == boxPt[1]) && (getPt3[1] == boxPt[1]))
    return true;
  // 4. Y = boxPt[1] + Ydim
  else if ((getPt1[1] == boxPt[1] + Ydim) && (getPt2[1] == boxPt[1] + Ydim) && (getPt3[1] == boxPt[1] + Ydim))
    return true;
  // 5. Z = boxPt[2]
  else if ((getPt1[2] == boxPt[2]) && (getPt2[2] == boxPt[2]) && (getPt3[2] == boxPt[2]))
    return true;
  // 6. Z = boxPt[2] + Zdim
  else if ((getPt1[2] == boxPt[2] + Zdim) && (getPt2[2] == boxPt[2] + Zdim) && (getPt3[2] == boxPt[2] + Zdim))
    return true;
  else
    return false;
}

} // namespace GEO

} // namespace NEM
