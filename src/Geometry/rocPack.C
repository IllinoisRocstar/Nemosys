#define _USE_MATH_DEFINES

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <rocPack.H>
#include <Eigen/Geometry>
#include <rocPackShape.H>
#include <hmxShape.H>
#include <petnShape.H>
#include <icosidodecahedronShape.H>

// GMSH Header
#include <gmsh.h>

namespace NEM {

namespace GEO {

// class constructor
rocPack::rocPack(const std::string &fname, const std::string &outName)
{
  std::cout << "rocPack class constructed!" << std::endl;

  InFile = fname;
  OutFile = outName;
}

void rocPack::rocPack2Surf()
{
  // Parses output file from Rocpack
  rocParser();

  // Generates geometry from parsed database and writes in STL/VTK files.
  rocToGeom();
}

void rocPack::initialize()
{
  gmsh::initialize();
  gmsh::model::add("Packs");
  gmsh::option::setNumber("General.NumThreads", 10);
  gmsh::option::setNumber("Geometry.OCCBooleanPreserveNumbering", 1);
  gmsh::option::setNumber("Geometry.OCCParallel", 1);
  gmsh::option::setNumber("Mesh.FlexibleTransfinite", 0);
  gmsh::option::setNumber("Mesh.MshFileVersion", 2.2);
}

void rocPack::rocParser()
{
  // Check if file exist
  std::cout << " - Parsing output file ... " << std::endl;
  std::ifstream rocOut(InFile);
  if (rocOut.is_open())
  {
    // Getting box dimensions
    std::string line = findWord("boundary");
    std::vector<std::string> rocTokens;

    rocTokens = Tokenize(line, ' ');

    Xdim = std::atof(strToChar(rocTokens[3]));
    Ydim = std::atof(strToChar(rocTokens[4]));
    Zdim = std::atof(strToChar(rocTokens[5]));

    if (rocTokens[6] != "periodic")
    {
      std::cerr << "Please select output file with periodic geometries!"
                << std::endl;
      throw;
    }

    boxPt.push_back(-Xdim/2);
    boxPt.push_back(-Ydim/2);
    boxPt.push_back(-Zdim/2);

    // Stores all lines in file to myLines
    std::string linesOut; // Temporary string for multiple uses
    std::vector<std::string> myLines; // Stores whole file line by line
    while (std::getline(rocOut, linesOut))
      myLines.push_back(linesOut);


    // Checking for user-specified options
    for (int i=0; i<myLines.size(); i++)
    {
      if (myLines[i].find("NoPeriodicity") != std::string::npos)
        if (myLines[i].find("true") != std::string::npos ||
            myLines[i].find("True") != std::string::npos ||
            myLines[i].find("TRUE") != std::string::npos)
          noPeriodicity = true;

      if (myLines[i].find("RemoveBoundaryPacks") != std::string::npos)
        if (myLines[i].find("true") != std::string::npos ||
            myLines[i].find("True") != std::string::npos ||
            myLines[i].find("TRUE") != std::string::npos)
          removeBoundaryPacks = true;
    }

    // Checking if crystal shapes are present
    // crystalNames and verts/faces go parallel in terms of indexing
    int nCrystals = 0;
    for (int i=0; i<myLines.size(); i++)
    {
      if (myLines[i].find("#import") != std::string::npos ||
          myLines[i].find("import") == std::string::npos)
      {
        // Nothing
      }
      else if (myLines[i].find("import") != std::string::npos ||
               myLines[i].find("#import") == std::string::npos)
      {
        nCrystals++;
        std::vector<std::string> importStr;
        importStr = Tokenize(myLines[i], "//\"");
        crystalNames.push_back(importStr[2]);
      }
      else
      {
        
      }
    }

    if (nCrystals > 0)
    {
      verts.resize(nCrystals);
      faces.resize(nCrystals);
      
      for (int i=0; i<nCrystals; i++)
      {
        rocPackShape *getData = rocPackShape::getShape(crystalNames[i]);
        verts[i] = getData->getVertices();
        faces[i] = getData->getFaces();
      }

      normalizeVerts();
    }

    // Checking if any shapes are present. If present, extract imp data
    for (int i=0; i<myLines.size(); i++)
    {
      if (myLines[i].find("shape ") != std::string::npos)
      {
        std::vector<std::string> baseShapes;
        baseShapes = Tokenize(myLines[i],"         ");
        shapeNames.push_back(baseShapes[3]);
        uniqueNames.push_back(baseShapes[1]);

        if (baseShapes[3] == "cylinder")
        {
          cylParams.resize(2);
          cylParams[0] = std::atof(strToChar(baseShapes[5]));
          cylParams[1] = std::atof(strToChar(baseShapes[6]));
        }

        if (baseShapes[3] == "ellipsoid")
        {
          ellipsoidRad.resize(3);
          ellipsoidRad[0] = std::atof(strToChar(baseShapes[5]));
          ellipsoidRad[1] = std::atof(strToChar(baseShapes[6]));
          ellipsoidRad[2] = std::atof(strToChar(baseShapes[7]));
          ellipsoidPresent = true;
        }
      }
      else
      {}
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
    for (int i=0; i<myLines.size(); i++)
    {
      if (myLines[i].find("translate") != std::string::npos)
        numPacks++;
      if (myLines[i].find("translate") != std::string::npos && iter == 0)
        iter = i;
      else
        {}
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

    for (int i=iter; i<myLines.size(); i++)
    {
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
      trnsltDataTokens = getShapeData(i+1, translateStr, myLines);

      for (int j=0; j<3; j++)
        translateParams[h][j] = std::atof(strToChar(trnsltDataTokens[j+1]));

      trnsltDataTokens.clear();

      // Rotate Parameters
      rttDataTokens = getShapeData(i+2, rotateStr, myLines);

      for (int j=0; j<4; j++)
        rotateParams[h][j] = std::atof(strToChar(rttDataTokens[j+1]));

      rttDataTokens.clear();

      // Scale Parameters
      scaleDataTokens = getShapeData(i+3, scaleStr, myLines);
      scaleOfPack[h] = std::atof(strToChar(scaleDataTokens[1]));
      scaleDataTokens.clear();

      i = i + 5;
      h++;
    }   
  rocOut.close();
  }
  else
  {
    std::cerr << "Cannot open/find " << InFile << " file!" << std::endl;
    throw;
  }
}

// Loops through parsed data and makes gmsh geometries.
void rocPack::rocToGeom()
{
  std::cout << " - Creating pack geometries ... " << std::endl;
  // Starting Gmsh commands
  initialize();

  for (int i=0; i<nameOfPcks.size(); i++)
  {
    if (nameOfPcks[i] == "sphere")
      makeSphere(i);
    
    if (uniqueNames.size() > 0)
    {
      for (int j=0; j<uniqueNames.size(); j++)
      {
        if (nameOfPcks[i] == uniqueNames[j] && shapeNames[j] == "ellipsoid")
          makeEllipsoid(i);
        else if (nameOfPcks[i] == uniqueNames[j] && shapeNames[j] == "cylinder")
          makeCylinder(i);
        else {}
      }
    }

    if (crystalNames.size() > 0)
    {
      for (int k=0; k<crystalNames.size(); k++)
        if (crystalNames[k] == nameOfPcks[i])
          makeCrystalShape(i,k);
    }
  }

  // Tagging packs intersecting domain boundary
  tagBoundaryPacks();

  if (removeBoundaryPacks == true)
  {
    gmsh::model::occ::remove(bndryPackTags,true);
    gmsh::model::occ::synchronize();
  }

  if (ellipsoidPresent==true || noPeriodicity == true)
  {
    // Nothing
  }
  else
  {
    std::cout << " - Ensuring periodicity of geometry ..." << std::endl;
    makePeriodic(removeBoundaryPacks);
    gmsh::model::occ::synchronize();
  }

  std::cout << " - Writing triangulated surface files ..." << std::endl;
  gmsh::model::mesh::generate(2);
  gmsh::model::mesh::refine();
  gmsh::model::mesh::refine();
  gmsh::model::mesh::removeDuplicateNodes();
  gmsh::option::setNumber("Mesh.StlOneSolidPerSurface", 0);

  if (OutFile.find(".stl") != std::string::npos)
    geomToSTL(OutFile);
  else if (OutFile.find(".vtk") != std::string::npos)
    geomToVTK(OutFile);
  else if (OutFile.find(".msh") != std::string::npos)
    geomToMsh(OutFile);
  else {
    std::string stlFile = OutFile + ".stl";
    std::string vtkFile = OutFile + ".vtk";
    std::string mshFile = OutFile + ".msh";
    geomToSTL(stlFile); geomToVTK(vtkFile); geomToMsh(mshFile);
  }
  gmsh::finalize();
  std::cout << " - End of process!" << std::endl;

}

// Private methods
void rocPack::geomToSTL(const std::string &writeFile)
{
  // Writes created periodic geometries into STL file
  gmsh::write(writeFile);
}

void rocPack::geomToVTK(const std::string &writeFile)
{
  // Writes created periodic geometries into VTK file
  gmsh::write(writeFile);
}

void rocPack::geomToMsh(const std::string &writeFile)
{
  // Writes created periodic geometries into .msh file
  gmsh::write(writeFile);
}

void rocPack::makePeriodic(const bool rmbPacks)
{
  if (rmbPacks == false)
  {
    std::vector<int> xTranslate;
    std::vector<int> yTranslate;
    std::vector<int> zTranslate;

    xTranslate.push_back(1); yTranslate.push_back(0); zTranslate.push_back(0);
    xTranslate.push_back(1); yTranslate.push_back(0); zTranslate.push_back(1);
    xTranslate.push_back(1); yTranslate.push_back(0); zTranslate.push_back(-1);

    xTranslate.push_back(-1); yTranslate.push_back(0); zTranslate.push_back(0);
    xTranslate.push_back(-1); yTranslate.push_back(0); zTranslate.push_back(1);
    xTranslate.push_back(-1); yTranslate.push_back(0); zTranslate.push_back(-1);

    xTranslate.push_back(0); yTranslate.push_back(1); zTranslate.push_back(0);
    xTranslate.push_back(0); yTranslate.push_back(1); zTranslate.push_back(1);
    xTranslate.push_back(0); yTranslate.push_back(1); zTranslate.push_back(-1);

    xTranslate.push_back(0); yTranslate.push_back(-1); zTranslate.push_back(0);
    xTranslate.push_back(0); yTranslate.push_back(-1); zTranslate.push_back(1);
    xTranslate.push_back(0); yTranslate.push_back(-1); zTranslate.push_back(-1);

    xTranslate.push_back(1); yTranslate.push_back(1); zTranslate.push_back(0);
    xTranslate.push_back(1); yTranslate.push_back(1); zTranslate.push_back(1);
    xTranslate.push_back(1); yTranslate.push_back(1); zTranslate.push_back(-1);

    xTranslate.push_back(1); yTranslate.push_back(-1); zTranslate.push_back(0);
    xTranslate.push_back(1); yTranslate.push_back(-1); zTranslate.push_back(1);
    xTranslate.push_back(1); yTranslate.push_back(-1); zTranslate.push_back(-1);

    xTranslate.push_back(-1); yTranslate.push_back(1); zTranslate.push_back(0);
    xTranslate.push_back(-1); yTranslate.push_back(1); zTranslate.push_back(1);
    xTranslate.push_back(-1); yTranslate.push_back(1); zTranslate.push_back(-1);

    xTranslate.push_back(-1); yTranslate.push_back(-1); zTranslate.push_back(0);
    xTranslate.push_back(-1); yTranslate.push_back(-1); zTranslate.push_back(1);
    xTranslate.push_back(-1); yTranslate.push_back(-1); zTranslate.push_back(-1);

    xTranslate.push_back(0); yTranslate.push_back(0); zTranslate.push_back(1);
    xTranslate.push_back(0); yTranslate.push_back(0); zTranslate.push_back(-1);

    // Translating shapes
    int ptg = 1;
    std::cerr << "    Progress --> [0%";


    for (int h=0; h<26; h++)
    {
      std::vector<std::pair<int,int>> tagsCopy;
      gmsh::model::occ::copy(bndryPackTags, tagsCopy);
      gmsh::model::occ::translate(tagsCopy, xTranslate[h], 
                                  yTranslate[h], zTranslate[h]);

      if (h%3 == 0)
      {
        std::cerr.precision(3);
        std::cerr << ".." << 10.7142857*(ptg) << "%";
        ptg++;
      }
    }

    gmsh::model::occ::synchronize();

    // Containers for boolean operation
    std::vector<std::pair<int,int>> tagsPacks;
    std::vector<std::pair<int,int>> outBoolean;
    std::vector<std::vector<std::pair<int, int>>> outBoolMap;
    gmsh::model::getEntities(tagsPacks, 3);
    std::vector<std::pair<int,int>> tagsBox;
    int tmpTag =
          gmsh::model::occ::addBox(boxPt[0],boxPt[1],boxPt[2],Xdim,Ydim,Zdim);
    tagsBox.push_back(std::make_pair(3,tmpTag));

    // Boolean Intersection
    gmsh::model::occ::intersect(tagsPacks, tagsBox, outBoolean, outBoolMap);
    std::cout << "..100%]" << std::endl;
  }
  else
  {
    // Nothing
  }
}

// Finds particular word in file and returns that line
std::string rocPack::findWord(const std::string &word)
{
  std::ifstream input(InFile);
  std::string line;
  while (std::getline(input, line))
    if (line.find(word) == 0)
      return line;
  return "";
}


// Tokenized string in a string vector
std::vector<std::string> rocPack::Tokenize(const std::string &lineIn, 
                                           const char &delim )
{
  std::vector<std::string> rocTokens;

  std::string intermediate;
  std::stringstream checkStr(lineIn);  
      
  // Tokenizing
  while(std::getline(checkStr, intermediate, delim))
  { 
    rocTokens.push_back(intermediate);
  } 

  return rocTokens;
}

std::vector<std::string> rocPack::Tokenize(const std::string &lineIn, 
                                           const std::string &delims)
{
  std::vector<std::string> rocTokens;
  std::vector<char> delimStr;
  delimStr.reserve(delims.size());

  for (int i=0; i<delims.size(); i++)
    delimStr.push_back(delims[i]);

  std::string intermediate;
  std::stringstream checkStr(lineIn);  
      
  // Tokenizing
  int iter = 0;
  while(std::getline(checkStr, intermediate, delimStr[iter]))
  { 
    rocTokens.push_back(intermediate);
    iter++;
  }

  return rocTokens;
}


std::vector<std::string> rocPack::getShapeData(const int &iter,
                                               const std::string &a,
                                               const std::vector<std::string> &L)
{
  std::vector<std::string> LTokens = Tokenize(L[iter], a);

  return LTokens;
}


void rocPack::makeSphere(const int &n)
{

  double sphereRad = 1*scaleOfPack[n];

  gmsh::model::occ::addSphere(translateParams[n][0],
                              translateParams[n][1],
                              translateParams[n][2],
                              sphereRad);

  gmsh::model::occ::synchronize();

}


void rocPack::makeEllipsoid(const int &n)
{
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
  elVerts[1][0] = scaleOfPack[n]*ellipsoidRad[0];
  elVerts[1][1] = 0;
  elVerts[1][2] = 0;
  elVerts[2].resize(3);
  elVerts[2][0] = 0;
  elVerts[2][1] = scaleOfPack[n]*ellipsoidRad[1];
  elVerts[2][2] = 0;
  elVerts[3].resize(3);
  elVerts[3][0] = 0;
  elVerts[3][1] = 0;
  elVerts[3][2] = scaleOfPack[n]*ellipsoidRad[2];
  elVerts[4].resize(3);
  elVerts[4][0] = -scaleOfPack[n]*ellipsoidRad[0];
  elVerts[4][1] = 0;
  elVerts[4][2] = 0;
  elVerts[5].resize(3);
  elVerts[5][0] = 0;
  elVerts[5][1] = -scaleOfPack[n]*ellipsoidRad[1];
  elVerts[5][2] = 0;
  elVerts[6].resize(3);
  elVerts[6][0] = 0;
  elVerts[6][1] = 0;
  elVerts[6][2] = -scaleOfPack[n]*ellipsoidRad[2];


  // Rotating vertices
  rocQuaternion q = toQuaternion(rotateParams[n]);
  for (int i=0; i<7; i++)
    rotatedVerts[i] = rotateByQuaternion(q,elVerts[i]);


  // Adding translated points

  for (int i=0; i<rotatedVerts.size(); i++)
    pts[i] = gmsh::model::geo::addPoint(rotatedVerts[i][0]+translateParams[n][0],
                                      rotatedVerts[i][1]+translateParams[n][1],
                                      rotatedVerts[i][2]+translateParams[n][2]);


  // Try circle arc here
  arcs[0] = gmsh::model::geo::addEllipseArc(pts[1],pts[0],pts[6],pts[6]);
  arcs[1] = gmsh::model::geo::addEllipseArc(pts[6],pts[0],pts[4],pts[4]);
  arcs[2] = gmsh::model::geo::addEllipseArc(pts[4],pts[0],pts[3],pts[3]);
  arcs[3] = gmsh::model::geo::addEllipseArc(pts[3],pts[0],pts[1],pts[1]);
  arcs[4] = gmsh::model::geo::addEllipseArc(pts[1],pts[0],pts[2],pts[2]);
  arcs[5] = gmsh::model::geo::addEllipseArc(pts[2],pts[0],pts[4],pts[4]);
  arcs[6] = gmsh::model::geo::addEllipseArc(pts[4],pts[0],pts[5],pts[5]);
  arcs[7] = gmsh::model::geo::addEllipseArc(pts[5],pts[0],pts[1],pts[1]);
  arcs[8] = gmsh::model::geo::addEllipseArc(pts[6],pts[0],pts[2],pts[2]);
  arcs[9] = gmsh::model::geo::addEllipseArc(pts[2],pts[0],pts[3],pts[3]);
  arcs[10] = gmsh::model::geo::addEllipseArc(pts[3],pts[0],pts[5],pts[5]);
  arcs[11] = gmsh::model::geo::addEllipseArc(pts[5],pts[0],pts[6],pts[6]);

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
  tmpCurveLoop.push_back(surfs[0]); tmpCurveLoop.push_back(surfs[1]);
  tmpCurveLoop.push_back(surfs[2]); tmpCurveLoop.push_back(surfs[3]);
  tmpCurveLoop.push_back(surfs[4]); tmpCurveLoop.push_back(surfs[5]);
  tmpCurveLoop.push_back(surfs[6]); tmpCurveLoop.push_back(surfs[7]);
  int surfLoop = gmsh::model::geo::addSurfaceLoop(tmpCurveLoop);
  tmpCurveLoop.clear();

  // Volume
  tmpCurveLoop.push_back(surfLoop);
  gmsh::model::geo::addVolume(tmpCurveLoop);
  gmsh::model::geo::synchronize();

}


void rocPack::makeCylinder(const int &n)
{
  // "cylParams" is vector with 2 values in shape definition
  // cylParams[0] = Radius
  // cylParams[1] = Height

  // Scaling everything on thr fly
  // Cylinder Center
  std::vector<double> cylCenter = std::vector<double>(3);
  cylCenter[0] = 0;
  cylCenter[1] = 0 - (cylParams[1]/2)*scaleOfPack[n];
  cylCenter[2] = 0;

  // Cylinder Translation Axis
  std::vector<double> translation_axis = std::vector<double>(3);
  translation_axis[0] = 0;
  translation_axis[1] = (cylParams[1])*scaleOfPack[n];
  translation_axis[2] = 0;

  //Cylinder Vertices
  std::vector<std::vector<double>> cylVertices;
  cylVertices.resize(4);
  cylVertices[0].resize(3); cylVertices[1].resize(3);
  cylVertices[2].resize(3); cylVertices[3].resize(3);

  cylVertices[0][0] = 0;
  cylVertices[0][1] = 0 - (cylParams[1]/2)*scaleOfPack[n];
  cylVertices[0][2] = 0 + cylParams[0]*scaleOfPack[n];
  cylVertices[1][0] = 0 + cylParams[0]*scaleOfPack[n];
  cylVertices[1][1] = 0 - (cylParams[1]/2)*scaleOfPack[n];
  cylVertices[1][2] = 0;
  cylVertices[2][0] = 0;
  cylVertices[2][1] = 0 - (cylParams[1]/2)*scaleOfPack[n];
  cylVertices[2][2] = 0 - cylParams[0]*scaleOfPack[n];
  cylVertices[3][0] = 0 - cylParams[0]*scaleOfPack[n];
  cylVertices[3][1] = 0 - (cylParams[1]/2)*scaleOfPack[n];
  cylVertices[3][2] = 0;

  // Rotate By Quaternion
  rocQuaternion q = toQuaternion(rotateParams[n]);

  std::vector<double> rotatedCylCenter;
  std::vector<double> rotatedTranslationAxis;
  std::vector<std::vector<double>> rotatedCylVertices;
  rotatedCylVertices.resize(4);

  // Center
  rotatedCylCenter = rotateByQuaternion(q,cylCenter);

  // Translation_Axis
  rotatedTranslationAxis = rotateByQuaternion(q,translation_axis);

  // Vertices
  for (int i=0; i<4; i++)
    rotatedCylVertices[i] = rotateByQuaternion(q,cylVertices[i]);

  // Start adding geometry
  std::vector<int> pts = std::vector<int>(5);

  pts[0] = gmsh::model::occ::addPoint(rotatedCylCenter[0]
                                      +translateParams[n][0],
                                      rotatedCylCenter[1]
                                      +translateParams[n][1],
                                      rotatedCylCenter[2]
                                      +translateParams[n][2]);

  for (int i=0; i<rotatedCylVertices.size(); i++)
    pts[i+1] = gmsh::model::occ::addPoint(rotatedCylVertices[i][0]
                                        +translateParams[n][0],
                                        rotatedCylVertices[i][1]
                                        +translateParams[n][1],
                                        rotatedCylVertices[i][2]
                                        +translateParams[n][2]);

  std::vector<int> arcs = std::vector<int>(4);
  arcs[0] = gmsh::model::occ::addCircleArc(pts[1],pts[0],pts[2]);
  arcs[1] = gmsh::model::occ::addCircleArc(pts[2],pts[0],pts[3]);
  arcs[2] = gmsh::model::occ::addCircleArc(pts[3],pts[0],pts[4]);
  arcs[3] = gmsh::model::occ::addCircleArc(pts[4],pts[0],pts[1]);

  std::vector<int> curves = std::vector<int>(4);
  curves[0] = arcs[0]; curves[1] = arcs[1]; curves[2] = arcs[2]; 
  curves[3] = arcs[3];
  int loop = gmsh::model::occ::addCurveLoop(curves);

  std::vector<int> loops = std::vector<int>(1);
  loops[0] = loop;
  int surf = gmsh::model::occ::addPlaneSurface(loops);

  gmsh::model::occ::synchronize();
  std::vector<std::pair<int,int>> tagsSurfs;
  std::vector<std::pair<int,int>> outVols;
  tagsSurfs.push_back(std::make_pair(2,surf));

  gmsh::model::occ::extrude(tagsSurfs,rotatedTranslationAxis[0],
                            rotatedTranslationAxis[1],
                            rotatedTranslationAxis[2],
                            outVols);

  gmsh::model::occ::synchronize();

}

void rocPack::makeCrystalShape(const int &n, const int &index)
{
  std::vector<std::vector<double>> scaledVerts;
  scaledVerts.resize(verts[index].size());

  // Scaling Vertices
  for (int i=0; i<verts[index].size(); i++)
  {
    scaledVerts[i].resize(3);
    for (int j=0; j<verts[index][i].size(); j++)
      scaledVerts[i][j] = verts[index][i][j]*scaleOfPack[n];
  }

  // Rotating Vertices using Quaternions
  rocQuaternion q = toQuaternion(rotateParams[n]);
  std::vector<std::vector<double>> rotatedVerts;
  rotatedVerts.resize(verts[index].size());

  for (int i=0; i<verts[index].size(); i++)
    rotatedVerts[i] = rotateByQuaternion(q,scaledVerts[i]);

  // Translating Vertices
  for (int i=0; i<rotatedVerts.size(); i++)
    for (int j=0; j<rotatedVerts[i].size(); j++)
      rotatedVerts[i][j] = rotatedVerts[i][j] + translateParams[n][j];

  // Creating geometry
  std::vector<int> pts;
  pts.resize(rotatedVerts.size());
  for (int i=0; i<rotatedVerts.size(); i++)
  {
    pts[i] = gmsh::model::occ::addPoint(rotatedVerts[i][0],
                                        rotatedVerts[i][1],
                                        rotatedVerts[i][2]);
  }

  std::vector<std::vector<int>> lines;
  lines.resize(faces[index].size());
  std::vector<int> surfs;
  surfs.resize(faces[index].size());
  std::vector<int> lineLoops;
  lineLoops.resize(faces[index].size());

  std::map<std::pair<int,int>,int> lineMap;
  // Adding lines and surfaces
  for (int i=0; i<faces[index].size(); i++)
  {
    lines[i].resize(faces[index][i].size());
    for (int j=0; j<faces[index][i].size(); j++)
    {
      // Finding for existing line
      std::map<std::pair<int,int>,int>::iterator it;
      int a = pts[faces[index][i][j]];
      int b = pts[faces[index][i][(j+1)%faces[index][i].size()]];

      it = lineMap.find(std::make_pair(a,b));

      if (it == lineMap.end())
        it = lineMap.find(std::make_pair(b,a));

      if (it == lineMap.end())
      {
        lines[i][j] = gmsh::model::occ::addLine(a,b);
        lineMap.insert(std::make_pair(std::make_pair(a,b),lines[i][j]));
      }
      else
        lines[i][j] =  it->second;
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
  gmsh::model::occ::addVolume(surfLoopTags);

  gmsh::model::occ::synchronize();

}


void rocPack::normalizeVerts()
{
  // Normalize
  for (int k=0; k<verts.size(); k++)
  {
    double maxNorm = 0.0;
    for (int i=0; i<verts[k].size(); i++)
    {
      double norm = 
      sqrt(pow(verts[k][i][0],2) + pow(verts[k][i][1],2) + pow(verts[k][i][2],2));

      if (norm > maxNorm)
        maxNorm = norm;
    }

    for (int i=0; i<verts[k].size(); i++)
    {
      verts[k][i][0] = (verts[k][i][0]/maxNorm);
      verts[k][i][1] = (verts[k][i][1]/maxNorm);
      verts[k][i][2] = (verts[k][i][2]/maxNorm);
    }
  }
}


void rocPack::tagBoundaryPacks()
{
  // Add method for removing boundary packs
  std::vector<std::pair<int,int>> tagsWithinBoundary;
  std::vector<std::pair<int,int>> tagsAll;
  
  std::vector<int> insidePacks;
  std::vector<int> allPacks;
  std::vector<int> bndryPackVols;

  gmsh::model::getEntitiesInBoundingBox(boxPt[0],boxPt[1],boxPt[2],
                                        boxPt[0]+Xdim,boxPt[1]+Ydim,
                                        boxPt[2]+Zdim,tagsWithinBoundary,3);

  gmsh::model::getEntities(tagsAll,3);

  for (auto iter : tagsWithinBoundary)
    insidePacks.push_back(iter.second);

  for (auto iter2 : tagsAll)
    allPacks.push_back(iter2.second);

  for (int i=0; i<allPacks.size(); i++)
  {
    int ht = 0;
    for (int j=0; j<insidePacks.size(); j++)
      if (allPacks[i] == insidePacks[j])
      {
        ht++;
        break;
      }

    if (ht == 0)
      bndryPackVols.push_back(allPacks[i]);
  }

  // Converting volumes into pair
  for (int i=0; i<bndryPackVols.size(); i++)
  {
    bndryPackTags.push_back(std::make_pair(3,bndryPackVols[i]));
  }

}


char* rocPack::strToChar(const std::string &strng)
{
  char * tab = new char [strng.length()+1];
  std::strcpy(tab, strng.c_str());

  return tab;
}


std::vector<double> rocPack::rotateByQuaternion(const rocQuaternion &q,
                                                const std::vector<double> &v)
{
  Eigen::VectorXd InputVec(3);
  Eigen::VectorXd RotatedVec(3);
  std::vector<double> returnVec = std::vector<double>(3);

  InputVec[0] = v[0]; InputVec[1] = v[1]; InputVec[2] = v[2];

  Eigen::Quaternion<double> Quaternion(q.w,q.x,q.y,q.z);

  RotatedVec = Quaternion._transformVector(InputVec);

  returnVec[0] = RotatedVec[0];
  returnVec[1] = RotatedVec[1];
  returnVec[2] = RotatedVec[2];

  return returnVec;
}

rocQuaternion rocPack::toQuaternion(const std::vector<double> &r)
{
  rocQuaternion q; // Declaring Quaternion
  double angle = r[3]*3.141592653589/180; // Declaring Angle in Radians

  q.x = r[0] * sin(angle/2);
  q.y = r[1] * sin(angle/2);
  q.z = r[2] * sin(angle/2);
  q.w = cos(angle/2);

  return q;
}


}

}