#include <iostream>
#include <string>
#include "snappymeshGen.H"
#include "snappymeshParams.H"
#include <boost/filesystem.hpp>
#include "meshGen.H"
#include <vtkUnstructuredGrid.h>
#include <set>

// openfoam headers
#include "fvCFD.H"
#include "fvMesh.H"
#include "vtkTopo.H"
#include "fileName.H"
#include "argList.H"
#include "Time.H"
#include "decompositionMethod.H"
#include "noDecomp.H"
#include "fvMeshDistribute.H"
#include "wallPolyPatch.H"
#include "refinementParameters.H"
#include "vtkSetWriter.H"
#include "faceSet.H"
#include "motionSmoother.H"
#include "polyTopoChange.H"
#include "cellModeller.H"
#include "uindirectPrimitivePatch.H"
#include "surfZoneIdentifierList.H"
#include "UnsortedMeshedSurface.H"
#include "MeshedSurface.H"
#include "globalIndex.H"
#include "IOmanip.H"
#include "fvMeshTools.H"
#include "meshRefinement.H"

// snappyHexMesh Headers
#include "snappyRefineDriver.H"
#include "snappySnapDriver.H"
#include "snappyLayerDriver.H"
#include "searchableSurfaces.H"
#include "refinementSurfaces.H"
#include "refinementFeatures.H"
#include "shellSurfaces.H"
#include "snapParameters.H"
#include "layerParameters.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// This class has all the methods from snappyHexMesh. Currently, Surface
// Simplyfy method is not enabled.

snappymeshGen::snappymeshGen() // default constructor definition
{
  // Default Meshing Parameters
  _params = new snappymeshParams();
  _params->_withCastMesh = true;
  _params->_withSnap = true;
  _params->_withLayers = false;
  _params->geomFileName = "";
  _params->geomPatchName = "";
  _params->surfRefPatch = "";
  _params->maxLCells = 2000000;
  _params->maxGCells = 4000000;
  _params->minRefCells = 0;
  _params->cellsBetnLvls = 3;
  _params->refSurfLvlMin = 0;
  _params->refSurfLvlMax = 0;
  _params->featAngle = 60;
  _params->locMeshX = 0;
  _params->locMeshY = 0;
  _params->locMeshZ = 0;
  _params->_alwFreeZone = true;
  _params->snapSmthPatch = 4;
  _params->snapTol = 0.5;
  _params->solveSnapIter = 200;
  _params->relaxSnapIter = 6;
  _params->_relSize = true;
  _params->expRatio = 1.3;
  _params->finLThick = 1;
  _params->minThick = 0.1;
  _params->nGrow = 0;
  _params->lyrFeatAngle = 30;
  _params->lyrRelaxIter = 3;
  _params->lyrSmthSurfNorm = 1;
  _params->lyrSmthNorm = 3;
  _params->lyrSmthThick = 2;
  _params->lyrMaxFcTR = 0.5;
  _params->lyrMaxThickTMR = 1;
  _params->lyrMinMedAngl = 90;
  _params->lyrBuffrCells = 0;
  _params->lyrIter = 50;
  _params->qcMaxNOrtho = 65;
  _params->qcMaxBndrySkew = 20;
  _params->qcMaxIntSkew = 4;
  _params->qcMaxConc = 80;
  _params->qcMinVol = 1e-13;
  _params->qcMinTetQ = 1e-15;
  _params->qcMinArea = -1;
  _params->qcMinTwist = 0.02;
  _params->qcMinFaceW = 0.05;
  _params->qcMinVolRto = 0.01;
  _params->qcMinDet = 0.001;
  _params->qcMinTrTwist = -1;
  _params->qcSmthScale = 5;
  _params->qcErrRedctn = 0.75;

  // Initialization Tasks
  initialize();

}

// constructor definition with parameters input
snappymeshGen::snappymeshGen(snappymeshParams* params):_params(params)
{
  
  initialize(); // Foam Initializatiob

}


snappymeshGen::~snappymeshGen() // destructor definition
{
  if (defaults)
    delete _params;
}


void snappymeshGen::initialize()
{
  caseName = "_snappymesh";

  // create dictionaries needed
  createControlDict();
  createSnappyDict();
  createfvSchemesDict();
  createfvSolutionDict();
}


void snappymeshGen::createControlDict()
{
  // creating a base system directory
  const char dir_path[] = "./system";
  boost::filesystem::path dir(dir_path);
  try
  {
    boost::filesystem::create_directory(dir);
  }
  catch (boost::filesystem::filesystem_error &e)
  {
    std::cerr << "Problem in creating system directory for the cfMesh" << "\n";
    std::cerr << e.what() << std::endl;
    throw;
  }

  std::ofstream contDict;
  contDict.open(std::string(dir_path)+"/controlDict");
  std::string contText=
    "\
/*--------------------------------*- C++ -*----------------------------------*\n\
| =========                 |                                                |\n\
| \\\\      /  F ield         | NEMoSys: cfMesh interface                      |\n\
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
    object    controlDict;\n\
}\n\n\
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n\
deltaT  1;\n\n\
startTime   0;\n\n\
writeInterval   1;\n\n\
// ***************************************************************** //";
    contDict << contText;
    contDict.close();
}

void snappymeshGen::createfvSchemesDict()
{
  // creating a base system directory
  const char dir_path[] = "./system";
  boost::filesystem::path dir(dir_path);
  try
  {
    boost::filesystem::create_directory(dir);
  }
  catch (boost::filesystem::filesystem_error &e)
  {
    std::cerr << "Problem in creating system directory for the cfMesh" << "\n";
    std::cerr << e.what() << std::endl;
    throw;
  }

  std::ofstream contDict;
  contDict.open(std::string(dir_path)+"/fvSchemes");
  std::string contText=
    "\
/*--------------------------------*- C++ -*----------------------------------*\n\
| =========                 |                                                |\n\
| \\\\      /  F ield         | NEMoSys: cfMesh interface                      |\n\
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
    object    fvSchemes;\n\
}\n\n\
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n\
gradSchemes\n\
{\n\
    default         Gauss linear;\n\
    grad(p)         Gauss linear;\n\
}\n\
\n\
divSchemes\n\
{\n\
    default         none;\n\
    div(phi,U)      Gauss linear;\n\
}\n\
\n\
laplacianSchemes\n\
{\n\
    default         none;\n\
    laplacian(nu,U) Gauss linear corrected;\n\
    laplacian((1|A(U)),p) Gauss linear corrected;\n\
}\n\
// ******************************************************************* //";
    contDict << contText;
    contDict.close();
}

void snappymeshGen::createfvSolutionDict()
{
  // creating a base system directory
  const char dir_path[] = "./system";
  boost::filesystem::path dir(dir_path);
  try
  {
    boost::filesystem::create_directory(dir);
  }
  catch (boost::filesystem::filesystem_error &e)
  {
    std::cerr << "Problem in creating system directory for the cfMesh" << "\n";
    std::cerr << e.what() << std::endl;
    throw;
  }

  std::ofstream contDict;
  contDict.open(std::string(dir_path)+"/fvSolution");
  std::string contText=
    "\
/*--------------------------------*- C++ -*----------------------------------*\n\
| =========                 |                                                |\n\
| \\\\      /  F ield         | NEMoSys: cfMesh interface                      |\n\
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
    object    fvSolution;\n\
}\n\n\
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n\
\n\
// ********************************************************************** //";
  contDict << contText;
  contDict.close();
}

void snappymeshGen::createSnappyDict()
{
  // creating snappyHexMeshDict

  // creating a base system directory
  const char dir_path[] = "./system";
  boost::filesystem::path dir(dir_path);
  try
  {
    boost::filesystem::create_directory(dir);
  }
  catch (boost::filesystem::filesystem_error &e)
  {
    std::cerr << "Problem in creating system directory for the snappyHexMesh" 
              << "\n";
    std::cerr << e.what() << std::endl;
    throw;
  }

    // creating mesh dictionary file
  std::ofstream contDict;
  contDict.open(std::string(dir_path)+"/snappyHexMeshDict");
  // header
  std::string contText=
    "\
/*--------------------------------*- C++ -*----------------------------------*\n\
| =========                 |                                                |\n\
| \\\\      /  F ield         | NEMoSys: snappyHexMesh interface               |\n\
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
    object    snappyHexMeshDict;\n\
}\n\n";

  contText = contText + 
  "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";

  // initial booleans
  if ((_params->_withCastMesh)){
    contText = contText + "\ncastellatedMesh  true;\n";
  }
  else{
    contText = contText + "\ncastellatedMesh  false;\n";
  }

  if ((_params->_withSnap)){
    contText = contText + "\nsnap       true;\n";
  }else{
    contText = contText + "\nsnap       false;\n";
  }

  if ((_params->_withLayers)){
    contText = contText + "\naddLayers      true;\n";
  }
  else{
    contText = contText + "\naddLayers      false;\n";
  }

  // Geometry declaration
  contText = contText + "\ngeometry\n{\n";
  contText = contText + "\t" + (_params->geomFileName);
  contText = contText + "\n\t{\n\t\ttype\ttriSurfaceMesh;\n";
  contText = contText + "\t\tname\t" + (_params->geomPatchName) + ";\n\t}\n";

  if ((_params->_withGeomRefReg))
  {
    int i=1;
    for (auto pt=(_params->geomRefs).begin(); 
                          pt!=(_params->geomRefs).end(); pt++)
    {
      contText = contText + "\t" + (pt->patchNm);
      contText = contText + "\n\t{";
      contText = contText + "\n\t\ttype " + (pt->searchableName) + ";\n";
      contText = contText + "\t\t" + (pt->shapeParameters1) + ";\n";

      if ((pt->searchableName) == "searchableCylinder")
      {
        contText = contText + "\t\t" + (pt->shapeParameters2) + ";\n";
      }

      contText = contText + "\t\tradius " + std::to_string(pt->rad) + ";\n";
      contText = contText + "\n\t}\n";
    }
  }

  contText = contText + "};\n";

  // Castellated Mesh Controls
  contText = contText + "\n\ncastellatedMeshControls\n{\n";
  contText = contText + "\tmaxLocalCells\t"
          + std::to_string(_params->maxLCells) + ";\n";
  contText = contText + "\tmaxGlobalCells\t"
          + std::to_string(_params->maxGCells) + ";\n";
  contText = contText + "\tnCellsBetweenLevels\t"
          + std::to_string(_params->cellsBetnLvls) + ";\n";
  contText = contText + "\tminRefinementCells\t"
          + std::to_string(_params->minRefCells) + ";\n";
  contText = contText + "\n\tfeatures\n\t(\n\t);\n";
  contText = contText + "\n\trefinementSurfaces\n\t{\n";
  contText = contText + "\t\t" + (_params->surfRefPatch) + "\n\t\t{";
  contText = contText + "\n\t\t\tlevel ("
          + std::to_string(_params->refSurfLvlMin) + " "
          + std::to_string(_params->refSurfLvlMax) + ");\n";

  if ((_params->_withCellZones))
  {
    contText = contText + "\t\t\tfaceZone\t"
          + (_params->surfRefPatch) + ";\n";
    contText = contText + "\t\t\tcellZone\t"
          + (_params->surfRefPatch) +";\n";
    contText = contText + "\t\t\tcellZoneInside\tinside;\n";
  }

  contText = contText + "\t\t}\n\t}\n";
  contText = contText + "\tresolveFeatureAngle\t"
          + std::to_string(_params->featAngle) + ";\n";
    

  if ((_params->_withSurfRefReg))
  {
    contText = contText + "\n\trefinementSurfaces\n";
    contText = contText + "\t{\n";

    for (auto pt=(_params->surfRefs).begin();
                    pt!=(_params->surfRefs).end(); pt++)
    {
      contText = contText + "\n\t\t" + (pt->refPatchNm);
      contText = contText + "\n\t\t{";
      contText = contText + "\t\t\tlevels (("
              + std::to_string(pt->minLvl)
              + " " + std::to_string(pt->maxLvl) + "));\n";
      contText = contText + "\t\t}";
      }

      contText = contText + "\n\t}\n";
  }


  contText = contText + "\n\trefinementRegions\n\t{\n";


  if ((_params->_withGeomRefReg))
  {
    for (auto pt=(_params->geomRefs).begin();
                      pt!=(_params->geomRefs).end(); pt++)
    {
      contText = contText + "\n\t\t" + (pt->patchNm);
      contText = contText + "\n\t\t{";
      contText = contText + "\n\t\t\tmode " + (pt->mode) + ";\n";
      contText = contText + "\t\t\tlevels ((" + std::to_string(pt->minLvl)
                            + " " + std::to_string( pt->maxLvl) + "));\n";
      contText = contText + "\t\t}\n";
    }
  }

  contText = contText + "\t}";

  contText = contText + "\n\tlocationInMesh\t(" 
          + std::to_string(_params->locMeshX)
          + " " + std::to_string(_params->locMeshY) + " " 
          + std::to_string(_params->locMeshZ) + ");\n";
  if ((_params->_alwFreeZone) == 1)
    contText = contText + "\tallowFreeStandingZoneFaces\ttrue;\n";
  if ((_params->_alwFreeZone) == 0)
    contText = contText + "\tallowFreeStandingZoneFaces\tfalse;\n";
  contText = contText + "}\n";


  // Snap Controls
  contText = contText + "\n\nsnapControls\n{\n";
  contText = contText + "\tnSmoothPatch\t"
        + std::to_string(_params->snapSmthPatch) + ";\n";
  contText = contText + "\ttolerance\t"
        + std::to_string(_params->snapTol) + ";\n";
  contText = contText + "\tnSolveIter\t"
        + std::to_string(_params->solveSnapIter) + ";\n";
  contText = contText + "\tnRelaxIter\t"
        + std::to_string(_params->relaxSnapIter) + ";\n";
  contText = contText + "}\n";


  // Layer Controls
  contText = contText + "\n\naddLayersControls\n{\n";
  if ((_params->_relSize) == 1)
    contText = contText + "\trelativeSizes\ttrue;\n";
  if ((_params->_relSize) == 0)
    contText = contText + "\trelativeSizes\tfalse;\n";
  contText = contText + "\tlayers\n\t{\n\t}\n";
  contText = contText + "\texpansionRatio\t"
          + std::to_string(_params->expRatio) + ";\n";
  contText = contText + "\tfinalLayerThickness\t"
          + std::to_string(_params->finLThick) + ";\n";
  contText = contText + "\tminThickness\t"
          + std::to_string(_params->minThick) + ";\n";
  contText = contText + "\tnGrow\t" + std::to_string(_params->nGrow) + ";\n";
  contText = contText + "\tfeatureAngle\t"
          + std::to_string(_params->lyrFeatAngle) + ";\n";
  contText = contText + "\tnRelaxIter\t"
          + std::to_string(_params->lyrRelaxIter) + ";\n";
  contText = contText + "\tnSmoothSurfaceNormals\t"
          + std::to_string(_params->lyrSmthSurfNorm) + ";\n";
  contText = contText + "\tnSmoothNormals\t"
          + std::to_string(_params->lyrSmthNorm) + ";\n";
  contText = contText + "\tnSmoothThickness\t"
          + std::to_string(_params->lyrSmthThick) + ";\n";
  contText = contText + "\tmaxFaceThicknessRatio\t"
          + std::to_string(_params->lyrMaxFcTR) + ";\n";
  contText = contText + "\tmaxThicknessToMedialRatio\t"
          + std::to_string(_params->lyrMaxThickTMR) + ";\n";
  contText = contText + "\tminMedialAxisAngle\t"
          + std::to_string(_params->lyrMinMedAngl) + ";\n";
  contText = contText + "\tnBufferCellsNoExtrude\t"
          + std::to_string(_params->lyrBuffrCells) + ";\n";
  contText = contText + "\tnLayerIter\t"
          + std::to_string(_params->lyrIter) + ";\n";
  contText = contText + "}\n";


  // Mesh Quality Controls
  contText = contText + "\n\nmeshQualityControls\n{\n";
  contText = contText + "\tmaxNonOrtho\t" 
          + std::to_string(_params->qcMaxNOrtho) + ";\n";
  contText = contText + "\tmaxBoundarySkewness\t" 
          + std::to_string(_params->qcMaxBndrySkew) + ";\n";
  contText = contText + "\tmaxInternalSkewness\t" 
          + std::to_string(_params->qcMaxIntSkew) + ";\n";
  contText = contText + "\tmaxConcave\t" 
          + std::to_string(_params->qcMaxConc) + ";\n";
  contText = contText + "\tminVol\t" 
          + std::to_string(_params->qcMinVol) + ";\n";
  contText = contText + "\tminTetQuality\t" 
          + std::to_string(_params->qcMinTetQ) + ";\n";
  contText = contText + "\tminArea\t" 
          + std::to_string(_params->qcMinArea) + ";\n";
  contText = contText + "\tminTwist\t" 
          + std::to_string(_params->qcMinTwist) + ";\n";
  contText = contText + "\tminFaceWeight\t" 
          + std::to_string(_params->qcMinFaceW) + ";\n";
  contText = contText + "\tminVolRatio\t" 
          + std::to_string(_params->qcMinVolRto) + ";\n";
  contText = contText + "\tminDeterminant\t" 
          + std::to_string(_params->qcMinDet) + ";\n";
  contText = contText + "\tminTriangleTwist\t" 
          + std::to_string(_params->qcMinTrTwist) + ";\n";
  contText = contText + "\tnSmoothScale\t" 
          + std::to_string(_params->qcSmthScale) + ";\n";
  contText = contText + "\terrorReduction\t" 
          + std::to_string(_params->qcErrRedctn) + ";\n";
  contText = contText + "}\n";

  contText = contText + "\nmergeTolerance\t"
          + std::to_string(_params->mergeTol) + ";\n";

  contText = contText +
  "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //";
  contDict << contText;
  contDict.close();

}

// Convert size (as fraction of defaultCellSize) to refinement level
Foam::label snappymeshGen::sizeCoeffToRefinement
(
    const Foam::scalar level0Coeff, // ratio of hex cell size 
                                    // v.s. defaultCellSize
    const Foam::scalar sizeCoeff
)
{
     using namespace Foam;
     return round(::log(level0Coeff/sizeCoeff)/::log(2));
}


Foam::autoPtr<refinementSurfaces> snappymeshGen::createRefinementSurfaces
(
    const Foam::searchableSurfaces& allGeometry,
    const Foam::dictionary& surfacesDict,
    const Foam::dictionary& shapeControlDict,
    const Foam::label gapLevelIncrement,
    const Foam::scalar level0Coeff
)
{
  using namespace Foam;

  autoPtr<refinementSurfaces> surfacePtr;

  // Count number of surfaces.
  label surfI = 0;
  forAll(allGeometry.names(), geomI)
  {
    const word& geomName = allGeometry.names()[geomI];

    if (surfacesDict.found(geomName))
    {
      surfI++;
    }
  }

  labelList surfaces(surfI);
  wordList names(surfI);
  PtrList<surfaceZonesInfo> surfZones(surfI);

  labelList regionOffset(surfI);

  labelList globalMinLevel(surfI, 0);
  labelList globalMaxLevel(surfI, 0);
  labelList globalLevelIncr(surfI, 0);
  PtrList<dictionary> globalPatchInfo(surfI);
  List<Map<label>> regionMinLevel(surfI);
  List<Map<label>> regionMaxLevel(surfI);
  List<Map<label>> regionLevelIncr(surfI);
  List<Map<scalar>> regionAngle(surfI);
  List<Map<autoPtr<dictionary>>> regionPatchInfo(surfI);

  HashSet<word> unmatchedKeys(surfacesDict.toc());

  surfI = 0;
  forAll(allGeometry.names(), geomI)
  {
    const word& geomName = allGeometry.names()[geomI];

    const entry* ePtr = surfacesDict.lookupEntryPtr(geomName, false, true);

    if (ePtr)
    {
      const dictionary& shapeDict = ePtr->dict();
      unmatchedKeys.erase(ePtr->keyword());

      names[surfI] = geomName;
      surfaces[surfI] = geomI;

      const searchableSurface& surface = allGeometry[geomI];

      // Find the index in shapeControlDict
      // Invert surfaceCellSize to get the refinementLevel
      // OF-4.x does not have "optionalSubDict Method"

      const word scsFuncName =
          shapeDict.lookup("surfaceCellSizeFunction");
      const dictionary& scsDict =
          //shapeDict.optionalSubDict(scsFuncName + "Coeffs");
          shapeDict.subDict(scsFuncName + "Coeffs");

      const scalar surfaceCellSize =
          readScalar(scsDict.lookup("surfaceCellSizeCoeff"));

      const label refLevel = sizeCoeffToRefinement
      (
          level0Coeff,
          surfaceCellSize
      );

      globalMinLevel[surfI] = refLevel;
      globalMaxLevel[surfI] = refLevel;
      globalLevelIncr[surfI] = gapLevelIncrement;

      // Surface zones
      surfZones.set(surfI, new surfaceZonesInfo(surface, shapeDict));


      // Global perpendicular angle
      if (shapeDict.found("patchInfo"))
      {
        globalPatchInfo.set
        (
          surfI,
          shapeDict.subDict("patchInfo").clone()
        );
      }


      // Per region override of patchInfo

      if (shapeDict.found("regions"))
      {
        const dictionary& regionsDict = shapeDict.subDict("regions");
        const wordList& regionNames =
            allGeometry[surfaces[surfI]].regions();

        forAll(regionNames, regionI)
        {
          if (regionsDict.found(regionNames[regionI]))
          {
            // Get the dictionary for region
            const dictionary& regionDict = regionsDict.subDict
            (
              regionNames[regionI]
            );

            if (regionDict.found("patchInfo"))
            {
              regionPatchInfo[surfI].insert
              (
                regionI,
                regionDict.subDict("patchInfo").clone()
              );
            }
          }
        }
      }

            // Per region override of cellSize
      if (shapeDict.found("regions"))
      {
        const dictionary& shapeControlRegionsDict =
            shapeDict.subDict("regions");
        const wordList& regionNames =
            allGeometry[surfaces[surfI]].regions();

        forAll(regionNames, regionI)
        {
          if (shapeControlRegionsDict.found(regionNames[regionI]))
          {
            const dictionary& shapeControlRegionDict =
              shapeControlRegionsDict.subDict
              (
                regionNames[regionI]
              );

            const word scsFuncName =
              shapeControlRegionDict.lookup
              (
                "surfaceCellSizeFunction"
              );
            const dictionary& scsDict =
              shapeControlRegionDict.subDict
              (
                scsFuncName + "Coeffs"
              );

            const scalar surfaceCellSize =
              readScalar
              (
                scsDict.lookup("surfaceCellSizeCoeff")
              );

            const label refLevel = sizeCoeffToRefinement
              (
                level0Coeff,
                surfaceCellSize
              );

            regionMinLevel[surfI].insert(regionI, refLevel);
            regionMaxLevel[surfI].insert(regionI, refLevel);
            regionLevelIncr[surfI].insert(regionI, 0);
          }
        }
      }

      surfI++;
    }
  }

    // Calculate local to global region offset
    label nRegions = 0;

    forAll(surfaces, surfI)
    {
        regionOffset[surfI] = nRegions;
        nRegions += allGeometry[surfaces[surfI]].regions().size();
    }

    // Rework surface specific information into information per global region
    labelList minLevel(nRegions, 0);
    labelList maxLevel(nRegions, 0);
    labelList gapLevel(nRegions, -1);
    PtrList<dictionary> patchInfo(nRegions);

    forAll(globalMinLevel, surfI)
    {
        label nRegions = allGeometry[surfaces[surfI]].regions().size();

        // Initialise to global (i.e. per surface)
        for (label i = 0; i < nRegions; i++)
        {
            label globalRegionI = regionOffset[surfI] + i;
            minLevel[globalRegionI] = globalMinLevel[surfI];
            maxLevel[globalRegionI] = globalMaxLevel[surfI];
            gapLevel[globalRegionI] =
                maxLevel[globalRegionI]
              + globalLevelIncr[surfI];

            if (globalPatchInfo.set(surfI))
            {
                patchInfo.set
                (
                    globalRegionI,
                    globalPatchInfo[surfI].clone()
                );
            }
        }

        // Overwrite with region specific information
        forAllConstIter(Map<label>, regionMinLevel[surfI], iter)
        {
            label globalRegionI = regionOffset[surfI] + iter.key();

            minLevel[globalRegionI] = iter();
            maxLevel[globalRegionI] = regionMaxLevel[surfI][iter.key()];
            gapLevel[globalRegionI] =
                maxLevel[globalRegionI]
              + regionLevelIncr[surfI][iter.key()];
        }

        const Map<autoPtr<dictionary>>& localInfo = regionPatchInfo[surfI];
        forAllConstIter(Map<autoPtr<dictionary>>, localInfo, iter)
        {
            label globalRegionI = regionOffset[surfI] + iter.key();
            patchInfo.set(globalRegionI, iter()().clone());
        }
    }

    surfacePtr.set
    (
        new refinementSurfaces
        (
            allGeometry,
            surfaces,
            names,
            surfZones,
            regionOffset,
            minLevel,
            maxLevel,
            gapLevel,
            scalarField(nRegions, -GREAT),  //perpendicularAngle,
            patchInfo
        )
    );


    const refinementSurfaces& rf = surfacePtr();

    // Determine maximum region name length
    label maxLen = 0;
    forAll(rf.surfaces(), surfI)
    {
        label geomI = rf.surfaces()[surfI];
        const wordList& regionNames = allGeometry.regionNames()[geomI];
        forAll(regionNames, regionI)
        {
            maxLen = Foam::max(maxLen, label(regionNames[regionI].size()));
        }
    }


    Info<< Foam::setw(maxLen) << "Region"
        << Foam::setw(10) << "Min Level"
        << Foam::setw(10) << "Max Level"
        << Foam::setw(10) << "Gap Level" << nl
        << Foam::setw(maxLen) << "------"
        << Foam::setw(10) << "---------"
        << Foam::setw(10) << "---------"
        << Foam::setw(10) << "---------" << endl;

    forAll(rf.surfaces(), surfI)
    {
        label geomI = rf.surfaces()[surfI];

        Info<< rf.names()[surfI] << ':' << nl;

        const wordList& regionNames = allGeometry.regionNames()[geomI];

        forAll(regionNames, regionI)
        {
            label globalI = rf.globalRegion(surfI, regionI);

            Info<< Foam::setw(maxLen) << regionNames[regionI]
                << Foam::setw(10) << rf.minLevel()[globalI]
                << Foam::setw(10) << rf.maxLevel()[globalI]
                << Foam::setw(10) << rf.gapLevel()[globalI] << endl;
        }
    }


    return surfacePtr;
}



void snappymeshGen::writeMesh
(
    const Foam::string& msg,
    const Foam::meshRefinement& meshRefiner,
    const Foam::meshRefinement::debugType debugLevel,
    const Foam::meshRefinement::writeType writeLevel
)
{
    using namespace Foam;
    const fvMesh& mesh = meshRefiner.mesh();

    meshRefiner.printMeshInfo(debugLevel, msg);
    Info<< "Writing mesh to time " << meshRefiner.timeName() << endl;

    meshRefiner.write
    (
        debugLevel,
        meshRefinement::writeType(writeLevel | meshRefinement::WRITEMESH),
        mesh.time().path()/meshRefiner.timeName()
    );
    Info<< "Wrote mesh in = "
        << mesh.time().cpuTimeIncrement() << " s." << endl;
}


void snappymeshGen::removeZeroSizedPatches(Foam::fvMesh& mesh)
{
    using namespace Foam;
    // Remove any zero-sized ones. Assumes
    // - processor patches are already only there if needed
    // - all other patches are available on all processors
    // - but coupled ones might still be needed, even if zero-size
    //   (e.g. processorCyclic)
    // See also logic in createPatch.
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    labelList oldToNew(pbm.size(), -1);
    label newPatchi = 0;
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];

        if (!isA<processorPolyPatch>(pp))
        {
            if
            (
                isA<coupledPolyPatch>(pp)
             || returnReduce(pp.size(), sumOp<label>())
            )
            {
                // Coupled (and unknown size) or uncoupled and used
                oldToNew[patchi] = newPatchi++;
            }
        }
    }

    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];

        if (isA<processorPolyPatch>(pp))
        {
            oldToNew[patchi] = newPatchi++;
        }
    }


    const label nKeepPatches = newPatchi;

    // Shuffle unused ones to end
    if (nKeepPatches != pbm.size())
    {
        Info<< endl
            << "Removing zero-sized patches:" << endl << incrIndent;

        forAll(oldToNew, patchi)
        {
            if (oldToNew[patchi] == -1)
            {
                Info<< indent << pbm[patchi].name()
                    << " type " << pbm[patchi].type()
                    << " at position " << patchi << endl;
                oldToNew[patchi] = newPatchi++;
            }
        }
        Info<< decrIndent;

        fvMeshTools::reorderPatches(mesh, oldToNew, nKeepPatches, true);
        Info<< endl;
    }
}


Foam::scalar snappymeshGen::getMergeDistance(const Foam::polyMesh& mesh, 
                                            const Foam::scalar mergeTol)
{
    using namespace Foam;
    const boundBox& meshBb = mesh.bounds();
    scalar mergeDist = mergeTol * meshBb.mag();

    Info<< nl
        << "Overall mesh bounding box  : " << meshBb << nl
        << "Relative tolerance         : " << mergeTol << nl
        << "Absolute matching distance : " << mergeDist << nl
        << endl;

    // check writing tolerance
    if (mesh.time().writeFormat() == IOstream::ASCII)
    {
        const scalar writeTol = std::pow
        (
            scalar(10.0),
            -scalar(IOstream::defaultPrecision())
        );

        if (mergeTol < writeTol)
        {
            FatalErrorInFunction
                << "Your current settings specify ASCII writing with "
                << IOstream::defaultPrecision() << " digits precision." << nl
                << "Your merging tolerance (" << mergeTol
                << ") is finer than this." << nl
                << "Change to binary writeFormat, "
                << "or increase the writePrecision" << endl
                << "or adjust the merge tolerance (mergeTol)."
                << exit(FatalError);
        }
    }

    return mergeDist;
}


void snappymeshGen::extractSurface
(
    const Foam::polyMesh& mesh,
    const Foam::Time& runTime,
    const Foam::labelHashSet& includePatches,
    const Foam::fileName& outFileName
)
{
    using namespace Foam;
    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

    // Collect sizes. Hash on names to handle local-only patches (e.g.
    //  processor patches)
    HashTable<label> patchSize(1000);
    label nFaces = 0;
    forAllConstIter(labelHashSet, includePatches, iter)
    {
        const polyPatch& pp = bMesh[iter.key()];
        patchSize.insert(pp.name(), pp.size());
        nFaces += pp.size();
    }
    Pstream::mapCombineGather(patchSize, plusEqOp<label>());


    // Allocate zone/patch for all patches
    HashTable<label> compactZoneID(1000);
    forAllConstIter(HashTable<label>, patchSize, iter)
    {
        label sz = compactZoneID.size();
        compactZoneID.insert(iter.key(), sz);
    }
    Pstream::mapCombineScatter(compactZoneID);


    // Rework HashTable into labelList just for speed of conversion
    labelList patchToCompactZone(bMesh.size(), -1);
    forAllConstIter(HashTable<label>, compactZoneID, iter)
    {
        label patchi = bMesh.findPatchID(iter.key());
        if (patchi != -1)
        {
            patchToCompactZone[patchi] = iter();
        }
    }

    // Collect faces on zones
    DynamicList<label> faceLabels(nFaces);
    DynamicList<label> compactZones(nFaces);
    forAllConstIter(labelHashSet, includePatches, iter)
    {
        const polyPatch& pp = bMesh[iter.key()];
        forAll(pp, i)
        {
            faceLabels.append(pp.start()+i);
            compactZones.append(patchToCompactZone[pp.index()]);
        }
    }

    // Addressing engine for all faces
    uindirectPrimitivePatch allBoundary
    (
        UIndirectList<face>(mesh.faces(), faceLabels),
        mesh.points()
    );


    // Find correspondence to master points
    labelList pointToGlobal;
    labelList uniqueMeshPoints;
    autoPtr<globalIndex> globalNumbers = mesh.globalData().mergePoints
    (
        allBoundary.meshPoints(),
        allBoundary.meshPointMap(),
        pointToGlobal,
        uniqueMeshPoints
    );

    // Gather all unique points on master
    List<pointField> gatheredPoints(Pstream::nProcs());
    gatheredPoints[Pstream::myProcNo()] = pointField
    (
        mesh.points(),
        uniqueMeshPoints
    );
    Pstream::gatherList(gatheredPoints);

    // Gather all faces
    List<faceList> gatheredFaces(Pstream::nProcs());
    gatheredFaces[Pstream::myProcNo()] = allBoundary.localFaces();
    forAll(gatheredFaces[Pstream::myProcNo()], i)
    {
        inplaceRenumber(pointToGlobal, gatheredFaces[Pstream::myProcNo()][i]);
    }
    Pstream::gatherList(gatheredFaces);

    // Gather all ZoneIDs
    List<labelList> gatheredZones(Pstream::nProcs());

   #ifdef HAVE_OF4
    gatheredZones[Pstream::myProcNo()] = compactZones.xfer();
   #endif
   #ifdef HAVE_OF5
    gatheredZones[Pstream::myProcNo()] = compactZones.xfer();
   #endif
   #ifdef HAVE_OF6
    gatheredZones[Pstream::myProcNo()] = compactZones.xfer();
   #endif
   #ifdef HAVE_OF7
    gatheredZones[Pstream::myProcNo()] = move(compactZones);
   #endif
    Pstream::gatherList(gatheredZones);

    // On master combine all points, faces, zones
    if (Pstream::master())
    {
        pointField allPoints = ListListOps::combine<pointField>
        (
            gatheredPoints,
            accessOp<pointField>()
        );
        gatheredPoints.clear();

        faceList allFaces = ListListOps::combine<faceList>
        (
            gatheredFaces,
            accessOp<faceList>()
        );
        gatheredFaces.clear();

        labelList allZones = ListListOps::combine<labelList>
        (
            gatheredZones,
            accessOp<labelList>()
        );
        gatheredZones.clear();


        // Zones
        surfZoneIdentifierList surfZones(compactZoneID.size());
        forAllConstIter(HashTable<label>, compactZoneID, iter)
        {
            surfZones[iter()] = surfZoneIdentifier(iter.key(), iter());
            Info<< "surfZone " << iter()  <<  " : " << surfZones[iter()].name()
                << endl;
        }

      #ifdef HAVE_OF4
        UnsortedMeshedSurface<face> unsortedFace
        (
            xferMove(allPoints),
            xferMove(allFaces),
            xferMove(allZones),
            xferMove(surfZones)
        );
      #endif
      #ifdef HAVE_OF5
        UnsortedMeshedSurface<face> unsortedFace
        (
            xferMove(allPoints),
            xferMove(allFaces),
            xferMove(allZones),
            xferMove(surfZones)
        );
      #endif
      #ifdef HAVE_OF6
        UnsortedMeshedSurface<face> unsortedFace
        (
            xferMove(allPoints),
            xferMove(allFaces),
            xferMove(allZones),
            xferMove(surfZones)
        );
      #endif
      #ifdef HAVE_OF7
        UnsortedMeshedSurface<face> unsortedFace
        (
            clone(allPoints),
            clone(allFaces),
            clone(allZones),
            clone(surfZones)
        );
      #endif


        MeshedSurface<face> sortedFace(unsortedFace);

        fileName globalCasePath
        (
            runTime.processorCase()
          ? runTime.path()/".."/outFileName
          : runTime.path()/outFileName
        );
        globalCasePath.clean();

        Info<< "Writing merged surface to " << globalCasePath << endl;

        sortedFace.write(globalCasePath);
    }
}


int snappymeshGen::createMeshFromSTL(const char* fname)
{
    using namespace Foam;

    int argc = 1;
    char** argv = new char*[2];
    argv[0] = new char[100];
    strcpy(argv[0], "NONE");
    Foam::argList args(argc, argv);
    Foam::Info<< "Create time\n" << Foam::endl;
    Foam::argList::noParallel();

    Time runTime
    (
        Time::controlDictName,
        "",
        ""
    );

    const bool overwrite = true;

    autoPtr<fvMesh> meshPtr;

    {
        Foam::Info
            << "Create mesh for time = "
            << runTime.timeName() << Foam::nl << Foam::endl;

        meshPtr.set
        (
            new fvMesh
            (
                Foam::IOobject
                (
                    Foam::fvMesh::defaultRegion,
                    runTime.timeName(),
                    runTime,
                    Foam::IOobject::MUST_READ
                )
            )
        );
    }

    fvMesh& mesh = meshPtr();

    Info<< "Read mesh in = "
        << runTime.cpuTimeIncrement() << " s" << endl;

    // Check patches and faceZones are synchronised
    mesh.boundaryMesh().checkParallelSync(true);
    meshRefinement::checkCoupledFaceZones(mesh);


    // Read meshing dictionary
    const word dictName("snappyHexMeshDict");
    #include "setSystemMeshDictionaryIO.H"
    const IOdictionary meshDict(dictIO);


    // all surface geometry
    const dictionary& geometryDict = meshDict.subDict("geometry");

    // refinement parameters
    const dictionary& refineDict = meshDict.subDict("castellatedMeshControls");

    // mesh motion and mesh quality parameters
    const dictionary& motionDict = meshDict.subDict("meshQualityControls");

    // snap-to-surface parameters
    const dictionary& snapDict = meshDict.subDict("snapControls");

    // layer addition parameters
    const dictionary& layerDict = meshDict.subDict("addLayersControls");

    // absolute merge distance
    const scalar mergeDist = getMergeDistance
    (
        mesh,
        readScalar(meshDict.lookup("mergeTolerance"))
    );

    const Switch keepPatches(meshDict.lookupOrDefault("keepPatches", false));



    // Read decomposePar dictionary
    dictionary decomposeDict;
    {
        if (Pstream::parRun())
        {
            decomposeDict = IOdictionary
            (
                IOobject
                (
                    "decomposeParDict",
                    runTime.system(),
                    mesh,
                    IOobject::MUST_READ_IF_MODIFIED,
                    IOobject::NO_WRITE
                )
            );
        }
        else
        {
            decomposeDict.add("method", "none");
            decomposeDict.add("numberOfSubdomains", 1);
        }
    }


    // Debug
    // ~~~~~

    // Set debug level
    meshRefinement::debugType debugLevel = meshRefinement::debugType
    (
        meshDict.lookupOrDefault<label>
        (
            "debug",
            0
        )
    );
    {
        wordList flags;
        if (meshDict.readIfPresent("debugFlags", flags))
        {
            debugLevel = meshRefinement::debugType
            (
                meshRefinement::readFlags
                (
                    meshRefinement::IOdebugTypeNames,
                    flags
                )
            );
        }
    }
    if (debugLevel > 0)
    {
        meshRefinement::debug   = debugLevel;
        snappyRefineDriver::debug = debugLevel;
        snappySnapDriver::debug   = debugLevel;
        snappyLayerDriver::debug  = debugLevel;
    }

    // Set file writing level
    {
        wordList flags;
        if (meshDict.readIfPresent("writeFlags", flags))
        {
            meshRefinement::writeLevel
            (
                meshRefinement::writeType
                (
                    meshRefinement::readFlags
                    (
                        meshRefinement::IOwriteTypeNames,
                        flags
                    )
                )
            );
        }
    }

    // Set output level
    {
        wordList flags;
        if (meshDict.readIfPresent("outputFlags", flags))
        {
            meshRefinement::outputLevel
            (
                meshRefinement::outputType
                (
                    meshRefinement::readFlags
                    (
                        meshRefinement::IOoutputTypeNames,
                        flags
                    )
                )
            );
        }
    }


    // Read geometry
    // ~~~~~~~~~~~~~

    searchableSurfaces allGeometry
    (
        IOobject
        (
            "abc",                      // dummy name
            mesh.time().constant(),     // instance
            //mesh.time().findInstance("triSurface", word::null),// instance
            "triSurface",               // local
            mesh.time(),                // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        geometryDict,
        meshDict.lookupOrDefault("singleRegionName", true)
    );


    // Read refinement surfaces
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    autoPtr<refinementSurfaces> surfacesPtr;

    Info<< "Reading refinement surfaces." << endl;

    /*if (surfaceSimplify)
    {
        IOdictionary foamyHexMeshDict
        (
           IOobject
           (
                "foamyHexMeshDict",
                runTime.system(),
                runTime,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
           )
        );

        const dictionary& conformationDict =
            foamyHexMeshDict.subDict("surfaceConformation").subDict
            (
                "geometryToConformTo"
            );

        const dictionary& motionDict =
            foamyHexMeshDict.subDict("motionControl");

        const dictionary& shapeControlDict =
            motionDict.subDict("shapeControlFunctions");

        // Calculate current ratio of hex cells v.s. wanted cell size
        const scalar defaultCellSize =
            readScalar(motionDict.lookup("defaultCellSize"));

        const scalar initialCellSize = ::pow(meshPtr().V()[0], 1.0/3.0);

        //Info<< "Wanted cell size  = " << defaultCellSize << endl;
        //Info<< "Current cell size = " << initialCellSize << endl;
        //Info<< "Fraction          = " << initialCellSize/defaultCellSize
        //    << endl;

        surfacesPtr =
            createRefinementSurfaces
            (
                allGeometry,
                conformationDict,
                shapeControlDict,
                refineDict.lookupOrDefault("gapLevelIncrement", 0),
                initialCellSize/defaultCellSize
            );
    }*/
    //else
    //{
        surfacesPtr.set
        (
            new refinementSurfaces
            (
                allGeometry,
                refineDict.subDict("refinementSurfaces"),
                refineDict.lookupOrDefault("gapLevelIncrement", 0)
            )
        );

        Info<< "Read refinement surfaces in = "
            << mesh.time().cpuTimeIncrement() << " s" << nl << endl;
    //}

    refinementSurfaces& surfaces = surfacesPtr();


    // Checking only?

    /*if (checkGeometry)
    {
        // Extract patchInfo
        List<wordList> patchTypes(allGeometry.size());

        const PtrList<dictionary>& patchInfo = surfaces.patchInfo();
        const labelList& surfaceGeometry = surfaces.surfaces();
        forAll(surfaceGeometry, surfI)
        {
            label geomI = surfaceGeometry[surfI];
            const wordList& regNames = allGeometry.regionNames()[geomI];

            patchTypes[geomI].setSize(regNames.size());
            forAll(regNames, regionI)
            {
                label globalRegionI = surfaces.globalRegion(surfI, regionI);

                if (patchInfo.set(globalRegionI))
                {
                    patchTypes[geomI][regionI] =
                        word(patchInfo[globalRegionI].lookup("type"));
                }
                else
                {
                    patchTypes[geomI][regionI] = wallPolyPatch::typeName;
                }
            }
        }

        // Write some stats
        allGeometry.writeStats(patchTypes, Info);
        // Check topology
        allGeometry.checkTopology(true);
        // Check geometry
        allGeometry.checkGeometry
        (
            100.0,      // max size ratio
            1e-9,       // intersection tolerance
            autoPtr<writer<scalar>>(new vtkSetWriter<scalar>()),
            0.01,       // min triangle quality
            true
        );

        return 0;
    }*/



    // Read refinement shells
    // ~~~~~~~~~~~~~~~~~~~~~~

    Info<< "Reading refinement shells." << endl;
    shellSurfaces shells
    (
        allGeometry,
        refineDict.subDict("refinementRegions")
    );
    Info<< "Read refinement shells in = "
        << mesh.time().cpuTimeIncrement() << " s" << nl << endl;


    Info<< "Setting refinement level of surface to be consistent"
        << " with shells." << endl;
    surfaces.setMinLevelFields(shells);
    Info<< "Checked shell refinement in = "
        << mesh.time().cpuTimeIncrement() << " s" << nl << endl;


    // Read feature meshes
    // ~~~~~~~~~~~~~~~~~~~

    Info<< "Reading features." << endl;
    refinementFeatures features
    (
        mesh,
        refineDict.lookup("features")
    );
    Info<< "Read features in = "
        << mesh.time().cpuTimeIncrement() << " s" << nl << endl;



    // Refinement engine
    // ~~~~~~~~~~~~~~~~~

    Info<< nl
        << "Determining initial surface intersections" << nl
        << "-----------------------------------------" << nl
        << endl;

    // Main refinement engine
    meshRefinement meshRefiner
    (
        mesh,
        mergeDist,          // tolerance used in sorting coordinates
        overwrite,          // overwrite mesh files?
        surfaces,           // for surface intersection refinement
        features,           // for feature edges/point based refinement
        shells              // for volume (inside/outside) refinement
    );
    Info<< "Calculated surface intersections in = "
        << mesh.time().cpuTimeIncrement() << " s" << nl << endl;

    // Some stats
    meshRefiner.printMeshInfo(debugLevel, "Initial mesh");

    meshRefiner.write
    (
        meshRefinement::debugType(debugLevel&meshRefinement::OBJINTERSECTIONS),
        meshRefinement::writeType(0),
        mesh.time().path()/meshRefiner.timeName()
    );


    // Add all the surface regions as patches
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //- Global surface region to patch (non faceZone surface) or patches
    //  (faceZone surfaces)
    labelList globalToMasterPatch;
    labelList globalToSlavePatch;
    {
      Info<< nl
          << "Adding patches for surface regions" << nl
          << "----------------------------------" << nl
          << endl;

      // From global region number to mesh patch.
      globalToMasterPatch.setSize(surfaces.nRegions(), -1);
      globalToSlavePatch.setSize(surfaces.nRegions(), -1);

      Info<< setf(ios_base::left)
          << Foam::setw(6) << "Patch"
          << Foam::setw(20) << "Type"
          << Foam::setw(30) << "Region" << nl
          << Foam::setw(6) << "-----"
          << Foam::setw(20) << "----"
          << Foam::setw(30) << "------" << endl;

      const labelList& surfaceGeometry = surfaces.surfaces();
      const PtrList<dictionary>& surfacePatchInfo = surfaces.patchInfo();

      forAll(surfaceGeometry, surfI)
      {
        label geomI = surfaceGeometry[surfI];

        const wordList& regNames = allGeometry.regionNames()[geomI];

        Info<< surfaces.names()[surfI] << ':' << nl << nl;

        if (surfaces.surfZones()[surfI].faceZoneName().empty())
        {
          // 'Normal' surface
          forAll(regNames, i)
          {
            label globalRegionI = surfaces.globalRegion(surfI, i);

            label patchi;

            if (surfacePatchInfo.set(globalRegionI))
            {
              patchi = meshRefiner.addMeshedPatch
              (
                regNames[i],
                surfacePatchInfo[globalRegionI]
              );
            }
            else
            {
                dictionary patchInfo;
                patchInfo.set("type", wallPolyPatch::typeName);

                patchi = meshRefiner.addMeshedPatch
                (
                    regNames[i],
                    patchInfo
                );
            }

            Info<< setf(ios_base::left)
                << Foam::setw(6) << patchi
                << Foam::setw(20) << mesh.boundaryMesh()[patchi].type()
                << Foam::setw(30) << regNames[i] << nl;

                    globalToMasterPatch[globalRegionI] = patchi;
                    globalToSlavePatch[globalRegionI] = patchi;
                }
            }
            else
            {
                // Zoned surface
              forAll(regNames, i)
              {
                label globalRegionI = surfaces.globalRegion(surfI, i);

                // Add master side patch
                {
                  label patchi;

                  if (surfacePatchInfo.set(globalRegionI))
                  {
                      patchi = meshRefiner.addMeshedPatch
                      (
                          regNames[i],
                          surfacePatchInfo[globalRegionI]
                      );
                  }
                  else
                  {
                    dictionary patchInfo;
                    patchInfo.set("type", wallPolyPatch::typeName);

                    patchi = meshRefiner.addMeshedPatch
                    (
                        regNames[i],
                        patchInfo
                    );
                  }

                    Info<< setf(ios_base::left)
                        << Foam::setw(6) << patchi
                        << Foam::setw(20) << mesh.boundaryMesh()[patchi].type()
                        << Foam::setw(30) << regNames[i] << nl;

                    globalToMasterPatch[globalRegionI] = patchi;
                }
                // Add slave side patch
                {
                    const word slaveName = regNames[i] + "_slave";
                    label patchi;

                    if (surfacePatchInfo.set(globalRegionI))
                    {
                        patchi = meshRefiner.addMeshedPatch
                        (
                            slaveName,
                            surfacePatchInfo[globalRegionI]
                        );
                    }
                    else
                    {
                        dictionary patchInfo;
                        patchInfo.set("type", wallPolyPatch::typeName);

                        patchi = meshRefiner.addMeshedPatch
                        (
                            slaveName,
                            patchInfo
                        );
                    }

                    Info<< setf(ios_base::left)
                        << Foam::setw(6) << patchi
                        << Foam::setw(20) << mesh.boundaryMesh()[patchi].type()
                        << Foam::setw(30) << slaveName << nl;

                    globalToSlavePatch[globalRegionI] = patchi;
                }
              }
            }

            Info<< nl;
        }
        Info<< "Added patches in = "
            << mesh.time().cpuTimeIncrement() << " s" << nl << endl;
    }


    // Parallel
    // ~~~~~~~~

    // Decomposition
    autoPtr<decompositionMethod> decomposerPtr
    (
        decompositionMethod::New
        (
            decomposeDict
        )
    );
    decompositionMethod& decomposer = decomposerPtr();

    if (Pstream::parRun() && !decomposer.parallelAware())
    {
        FatalErrorInFunction
            << "You have selected decomposition method "
            << decomposer.typeName
            << " which is not parallel aware." << endl
            << "Please select one that is (hierarchical, ptscotch)"
            << exit(FatalError);
    }

    // Mesh distribution engine (uses tolerance to reconstruct meshes)
    fvMeshDistribute distributor(mesh, mergeDist);





    // Now do the real work -refinement -snapping -layers
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const Switch wantRefine(meshDict.lookup("castellatedMesh"));
    const Switch wantSnap(meshDict.lookup("snap"));
    const Switch wantLayers(meshDict.lookup("addLayers"));

    // Refinement parameters
    const refinementParameters refineParams(refineDict);

    // Snap parameters
    const snapParameters snapParams(snapDict);

    // Layer addition parameters
    const layerParameters layerParams(layerDict, mesh.boundaryMesh());


    if (wantRefine)
    {
        cpuTime timer;

        snappyRefineDriver refineDriver
        (
            meshRefiner,
            decomposer,
            distributor,
            globalToMasterPatch,
            globalToSlavePatch
        );


        if (!overwrite && !debugLevel)
        {
            const_cast<Time&>(mesh.time())++;
        }


        refineDriver.doRefine
        (
            refineDict,
            refineParams,
            snapParams,
            refineParams.handleSnapProblems(),
            motionDict
        );


        if (!keepPatches && !wantSnap && !wantLayers)
        {
            removeZeroSizedPatches(mesh);
        }

        writeMesh
        (
            "Refined mesh",
            meshRefiner,
            debugLevel,
            meshRefinement::writeLevel()
        );

        Info<< "Mesh refined in = "
            << timer.cpuTimeIncrement() << " s." << endl;
    }

    if (wantSnap)
    {
        cpuTime timer;

        snappySnapDriver snapDriver
        (
            meshRefiner,
            globalToMasterPatch,
            globalToSlavePatch
        );

        if (!overwrite && !debugLevel)
        {
            const_cast<Time&>(mesh.time())++;
        }

        // Use the resolveFeatureAngle from the refinement parameters
        scalar curvature = refineParams.curvature();
        scalar planarAngle = refineParams.planarAngle();

        snapDriver.doSnap
        (
            snapDict,
            motionDict,
            curvature,
            planarAngle,
            snapParams
        );

        if (!keepPatches && !wantLayers)
        {
            removeZeroSizedPatches(mesh);
        }

        writeMesh
        (
            "Snapped mesh",
            meshRefiner,
            debugLevel,
            meshRefinement::writeLevel()
        );

        Info<< "Mesh snapped in = "
            << timer.cpuTimeIncrement() << " s." << endl;
    }

    if (wantLayers)
    {
        cpuTime timer;

        snappyLayerDriver layerDriver
        (
            meshRefiner,
            globalToMasterPatch,
            globalToSlavePatch
        );

        // Use the maxLocalCells from the refinement parameters
        bool preBalance = returnReduce
        (
            (mesh.nCells() >= refineParams.maxLocalCells()),
            orOp<bool>()
        );


        if (!overwrite && !debugLevel)
        {
            const_cast<Time&>(mesh.time())++;
        }

        layerDriver.doLayers
        (
            layerDict,
            motionDict,
            layerParams,
            preBalance,
            decomposer,
            distributor
        );

        if (!keepPatches)
        {
            removeZeroSizedPatches(mesh);
        }

        writeMesh
        (
            "Layer mesh",
            meshRefiner,
            debugLevel,
            meshRefinement::writeLevel()
        );

        Info<< "Layers added in = "
            << timer.cpuTimeIncrement() << " s." << endl;
    }


    {
        // Check final mesh
        Info<< "Checking final mesh ..." << endl;
        faceSet wrongFaces(mesh, "wrongFaces", mesh.nFaces()/100);
        motionSmoother::checkMesh(false, mesh, motionDict, wrongFaces);
        const label nErrors = returnReduce
        (
            wrongFaces.size(),
            sumOp<label>()
        );

        if (nErrors > 0)
        {
            Info<< "Finished meshing with " << nErrors << " illegal faces"
                << " (concave, zero area or negative cell pyramid volume)"
                << endl;
            wrongFaces.write();
        }
        else
        {
            Info<< "Finished meshing without any errors" << endl;
        }
    }

    // Comment out surface simplify if not going to use
    /*if (surfaceSimplify)
    {
        const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

        labelHashSet includePatches(bMesh.size());

        if (args.optionFound("patches"))
        {
            includePatches = bMesh.patchSet
            (
                wordReList(args.optionLookup("patches")())
            );
        }
        else
        {
            forAll(bMesh, patchi)
            {
                const polyPatch& patch = bMesh[patchi];

                if (!isA<processorPolyPatch>(patch))
                {
                    includePatches.insert(patchi);
                }
            }
        }

        // Comment out args if needed
        fileName outFileName
        (
            args.optionLookupOrDefault<fileName>
            (
                "outFile",
                "constant/triSurface/simplifiedSurface.stl"
            )
        );

        extractSurface
        (
            mesh,
            runTime,
            includePatches,
            outFileName
        );

        pointIOField cellCentres
        (
            IOobject
            (
                "internalCellCentres",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh.cellCentres()
        );

        cellCentres.write();
    }*/


    Info<< "Finished meshing in = "
        << runTime.elapsedCpuTime() << " s." << endl;

    Info<< "End\n" << endl;

    //return 0;

    readSnappyFoamMesh();

}


void snappymeshGen::readSnappyFoamMesh()
{
    Foam::Info<< "Create time\n" << Foam::endl;
    Foam::argList::noParallel();

    Time runTime
    (
        Time::controlDictName,
        "",
        ""
    );
    // reading mesh database and converting
    Foam::word regionName;
    regionName = Foam::fvMesh::defaultRegion;
    Foam::Info
        << "Create mesh for time = "
        << runTime.timeName() << Foam::nl << Foam::endl;

    _fmesh = new Foam::fvMesh 
    (
        Foam::IOobject
        (
            regionName,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );

    genMshDB();
}

void snappymeshGen::genMshDB()
{
    // Obtaining the mesh data and generating requested output file
    std::cout << "Number of points " << _fmesh->nPoints() << std::endl;
    std::cout << "Number of cells "<< _fmesh->nCells() << std::endl;

    // declare vtk dataset
    vtkSmartPointer<vtkUnstructuredGrid> dataSet_tmp 
        = vtkSmartPointer<vtkUnstructuredGrid>::New();

    // decomposition
    Foam::vtkTopo::decomposePoly = false;

    // creating equivalent vtk topology from fvMesh
    // by default polyhedral cells will be decomposed to 
    // tets and pyramids. Additional points will be added
    // to underlying fvMesh.
    std::cout << "Performing topological decomposition.\n";
    Foam::vtkTopo topo(*_fmesh);

    // point coordinates
    Foam::pointField pf = _fmesh->points();
    vtkSmartPointer<vtkPoints> points 
        = vtkSmartPointer<vtkPoints>::New(); 
    for (int ipt=0; ipt<_fmesh->nPoints(); ipt++)
        points->InsertNextPoint(
                pf[ipt].x(), 
                pf[ipt].y(), 
                pf[ipt].z() 
                );
    dataSet_tmp->SetPoints(points);

    // cell types
    std::vector<int> pntIds;
    int nCelPnts = 0;
    for (int icl=0; icl<topo.vertLabels().size(); icl++)
    {
      if ( topo.cellTypes()[icl] != VTK_POLYHEDRON )
      {
        nCelPnts = topo.vertLabels()[icl].size();
        pntIds.resize(nCelPnts, -1);
        for (int ip=0; ip< nCelPnts; ip++)
            pntIds[ip] = topo.vertLabels()[icl][ip];
        createVtkCell(dataSet_tmp, topo.cellTypes()[icl], pntIds);
      }
      else
      {
        // polyhedral cells treated differently in vtk
        // faces should be defined for them
        int nFace = topo.vertLabels()[icl][0]; 
        vtkSmartPointer<vtkCellArray> faces =
            vtkSmartPointer<vtkCellArray>::New();
        std::vector<vtkIdType> faceIds;
        std::set<vtkIdType> pntIds;
        int dataId = 1;
        for (int iFace=0; iFace<nFace; iFace++)
        {
          faceIds.clear();
          int nFaceId = topo.vertLabels()[icl][dataId++];
          int pntId;
          for (int ifid=0; ifid<nFaceId; ifid++)
          {
            pntId = topo.vertLabels()[icl][dataId++];
            pntIds.insert(pntId);
            faceIds.push_back(pntId);
          }
          faces->InsertNextCell(nFaceId, &faceIds[0]);
          }
          std::vector<vtkIdType> pntIdsVec(pntIds.begin(), pntIds.end());
          dataSet_tmp->InsertNextCell(VTK_POLYHEDRON, pntIdsVec.size(),
                      &pntIdsVec[0], nFace, faces->GetPointer());
      }
  }
  dataSet = dataSet_tmp;
}


void snappymeshGen::createVtkCell(vtkSmartPointer<vtkUnstructuredGrid> dataSet,
                              const int cellType, 
                              std::vector<int>& vrtIds)
{
    vtkSmartPointer<vtkIdList> vtkCellIds = vtkSmartPointer<vtkIdList>::New();
    vtkCellIds->SetNumberOfIds(vrtIds.size());
    for (auto pit = vrtIds.begin(); pit!= vrtIds.end(); pit++)
        vtkCellIds->SetId(pit-vrtIds.begin(), *pit);
    dataSet->InsertNextCell(cellType,vtkCellIds);
}