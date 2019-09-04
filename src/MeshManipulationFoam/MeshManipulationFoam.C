#include <iostream>
#include <string>
#include <tuple>
#include <vector>
#include <vtkUnstructuredGrid.h>
#include <set>
#include <unordered_map>
#include "meshBase.H"
#include "foamMesh.H"
#include "MeshManipulationFoam.H"
#include "MeshManipulationFoamParams.H"
#include <vtkCellData.h>
#include <vtkFieldData.h>
#include <vtkCell.h>
#include <vtkMesh.H>
#include <vtkCellType.h>
#include <vtkPointData.h>
#include <vtkPolyDataNormals.h>
#include <vtkGeometryFilter.h>
#include <AuxiliaryFunctions.H>
#include <boost/filesystem.hpp>
#include <vtkCell3D.h>
#include <vtkDataSetTriangleFilter.h>
#include <vtkAppendFilter.h>
#include <vtkDataSet.h>
#include <vtkBooleanOperationPolyDataFilter.h>
#include <vtkSTLReader.h>
#include <vtkSTLWriter.h>
#include <vtkTriangleFilter.h>

// New

#include <vtkActor.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>


// openfoam headers
#include "fvCFD.H"
#include "fvMesh.H"
#include "fileName.H"
#include "vtkTopo.H"

//SurfLambdaMuSmooth
#include "argList.H"
#include "boundBox.H"
#include "edgeMesh.H"
#include "matchPoints.H"
#include "MeshedSurfaces.H"

//splitMeshByRegions
#include "SortableList.H"
#include "regionSplit.H"
#include "fvMeshSubset.H"
#include "IOobjectList.H"
#include "volFields.H"
#include "faceSet.H"
#include "cellSet.H"
#include "polyTopoChange.H"
#include "removeCells.H"
#include "EdgeMap.H"
#include "syncTools.H"
#include "ReadFields.H"
#include "mappedWallPolyPatch.H"
#include "fvMeshTools.H"
#include "zeroGradientFvPatchFields.H"

//mergeMeshes
#include "Time.H"
#include "mergePolyMesh.H"

//createPatch
#include "cyclicPolyPatch.H"
#include "syncTools.H"
#include "polyMesh.H"
#include "SortableList.H"
#include "OFstream.H"
#include "meshTools.H"
#include "IOPtrList.H"
#include "polyModifyFace.H"
#include "wordReList.H"
#include "IOdictionary.H"

//splitMeshByTopology
#include "triSurface.H"

// foamToSurface
#include "timeSelector.H"
#include "polyMesh.H"
#include "IOdictionary.H"
#include "polyPatch.H"

// third party
#include <ANN/ANN.h>

// TODO
// 1. Two of these methods (mergeMesh and CreatePatch) are little bit pack
//    mesh specific but can be modified to accept very broad user
//    arguments and perform mesh manipulations. - Akash
// 2. Extend cohesive element methods for tetrahedrons.

//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

MeshManipulationFoam::~MeshManipulationFoam()
{
  //class destructor
}

/*
  SurfaceLambdaMuSmooth takes input surface file and user defined values of
  lambda, mu, and surface smoothing iterations to perform laplacian smoothing
  on surface of geometry. This utility accepts following extensions as input:
  .ofs, .obj, .inp, .stl, .tri, .off, .stlb, .nas, .bdf, .gts, .vtk, and .ac.
  It will write output file at userdefined path
*/
void MeshManipulationFoam::surfLambdaMuSmooth()
{
  using namespace Foam;

  int argc = 1;
  char** argv = new char*[2];
  argv[0] = new char[100];
  strcpy(argv[0], "NONE");
  Foam::argList args(argc, argv);
  Foam::Info<< "Create time\n" << Foam::endl;
  Foam::Time runTime(Foam::Time::controlDictName, args);
  Foam::argList::noParallel();

  const fileName surfFileName = (_mshMnipPrms->slmssurfaceFile);
  const fileName outFileName = (_mshMnipPrms->slmsoutputFile);
  const scalar lambda = (_mshMnipPrms->lambda);
  const scalar mu = (_mshMnipPrms->mu);
  const label  iters = (_mshMnipPrms->slmsIterations);

  if (lambda < 0 || lambda > 1)
  {
      FatalErrorInFunction
          << lambda << endl
          << "0: no change   1: move vertices to average of neighbours"
          << exit(FatalError);
  }
  if (mu < 0 || mu > 1)
  {
      FatalErrorInFunction
          << mu << endl
          << "0: no change   1: move vertices to average of neighbours"
          << exit(FatalError);
  }
  Info<< "lambda      : " << lambda << nl
      << "mu          : " << mu << nl
      << "Iters       : " << iters << nl
      << "Reading surface from " << surfFileName << " ..." << endl;

  meshedSurface surf1(surfFileName);

  Info<< "Faces       : " << surf1.size() << nl
      << "Vertices    : " << surf1.nPoints() << nl
      << "Bounding Box: " << boundBox(surf1.localPoints()) << endl;

  PackedBoolList fixedPoints(surf1.localPoints().size(), false);

  if (_mshMnipPrms->_addFeatureFile)
  {
    const fileName featureFileName(_mshMnipPrms->_addFeatureFile);
    Info<< "Reading features from " << featureFileName << " ..." << endl;

    edgeMesh feMesh(featureFileName);

    getFixedPoints(feMesh, surf1.localPoints(), fixedPoints);

    Info<< "Number of fixed points on surface = " << fixedPoints.count()
        << endl;
  }

  pointField newPoints(surf1.localPoints());

  for (label iter = 0; iter < iters; iter++)
  {
    // Lambda
    {
      pointField newLocalPoints
      (
          (1 - lambda)*surf1.localPoints()
        + lambda*avg(surf1, fixedPoints)
      );

      pointField newPoints(surf1.points());
      UIndirectList<point>(newPoints, surf1.meshPoints()) =
      newLocalPoints;

      surf1.movePoints(newPoints);
    }

    // Mu
    if (mu != 0)
    {
      pointField newLocalPoints
      (
        (1 + mu)*surf1.localPoints()
        - mu*avg(surf1, fixedPoints)
      );

      pointField newPoints(surf1.points());
      UIndirectList<point>(newPoints, surf1.meshPoints()) =
      newLocalPoints;

      surf1.movePoints(newPoints);
    }
  }

  Info<< "Writing surface to " << outFileName << " ..." << endl;
  surf1.write(outFileName);

  Info<< "End\n" << endl;
}


/*
  SplitMeshRegions utility splits mesh into multiple regions. It takes cell-face
  -cell walk through the mesh, separates different cellZones into domains, and
  writes them into constant directory. Currently it supports splitting of mesh
  using cellZones. This function returns an integer to be used in mergeMeshes
  utility. SplitMeshRegions utility assigns random numbering to split domains
  from 0 to number of domains and skips a number automatically. This number is 
  returned and used in mergeMeshes function so that it knows which directory is
  missing from constant folder.
*/
std::pair<int,int> MeshManipulationFoam::splitMshRegions()
{
  int impVar;
  using namespace Foam;

  int argc = 1;
  char** argv = new char*[2];
  argv[0] = new char[100];
  strcpy(argv[0], "NONE");
  Foam::argList args(argc, argv);
  Foam::Info<< "Create time\n" << Foam::endl;
  Foam::Time runTime(Foam::Time::controlDictName, args);
  Foam::argList::noParallel();

  #include "createNamedMesh.H"
  const word oldInstance = mesh.pointsInstance();
  word blockedFacesName;

  const bool makeCellZones    = false;
  const bool largestOnly      = false;
  const bool insidePoint      = false;
  const bool useCellZones     = (_mshMnipPrms->_cellZones);
  const bool useCellZonesOnly = false;
  const bool useCellZonesFile = false;
  const bool overwrite        = (_mshMnipPrms->_overwriteMsh);
  const bool detectOnly       = false;
  const bool sloppyCellZones  = false;
  const bool useFaceZones     = false;
  const bool prefixRegion     = false;

  if
  (
    (useCellZonesOnly || useCellZonesFile)
   && (useCellZones || blockedFacesName.size())
  )
  {
    FatalErrorInFunction
      << "You cannot specify both -cellZonesOnly or -cellZonesFileOnly"
      << " (which specify complete split)"
      << " in combination with -blockedFaces or -cellZones"
      << " (which imply a split based on topology)"
      << exit(FatalError);
  }

  if (useFaceZones)
  {
    Info<< "Using current faceZones to divide inter-region interfaces"
        << " into multiple patches."
        << nl << endl;
  }
  else
  {
    Info<< "Creating single patch per inter-region interface."
        << nl << endl;
  }

  if (insidePoint && largestOnly)
  {
    FatalErrorInFunction
      << "You cannot specify both -largestOnly"
      << " (keep region with most cells)"
      << " and -insidePoint (keep region containing point)"
      << exit(FatalError);
  }

  const cellZoneMesh& cellZones = mesh.cellZones();

  // Existing zoneID
  labelList zoneID(mesh.nCells(), -1);
  // Neighbour zoneID.
  labelList neiZoneID(mesh.nFaces()-mesh.nInternalFaces());
  getZoneID(mesh, cellZones, zoneID, neiZoneID);

  // Determine per cell the region it belongs to
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // cellRegion is the labelList with the region per cell.
  labelList cellRegion;
  // Region per zone
  labelList regionToZone;
  // Name of region
  wordList regionNames;
  // Zone to region
  labelList zoneToRegion;

  label nCellRegions = 0;
  if (useCellZonesOnly)
  {
    Info<< "Using current cellZones to split mesh into regions."
        << " This requires all"
        << " cells to be in one and only one cellZone." << nl << endl;

    label unzonedCelli = findIndex(zoneID, -1);
    if (unzonedCelli != -1)
    {
      FatalErrorInFunction
        << "For the cellZonesOnly option all cells "
        << "have to be in a cellZone." << endl
        << "Cell " << unzonedCelli
        << " at" << mesh.cellCentres()[unzonedCelli]
        << " is not in a cellZone. There might be more unzoned cells."
        << exit(FatalError);
    }

    cellRegion = zoneID;
    nCellRegions = gMax(cellRegion)+1;
    regionToZone.setSize(nCellRegions);
    regionNames.setSize(nCellRegions);
    zoneToRegion.setSize(cellZones.size(), -1);

    for (label regionI = 0; regionI < nCellRegions; regionI++)
    {
        regionToZone[regionI] = regionI;
        zoneToRegion[regionI] = regionI;
        regionNames[regionI] = cellZones[regionI].name();
    }
  }

  // Will be implemented in future
  else if (useCellZonesFile) 
  {
    const word zoneFile = args.option("cellZonesFileOnly");
    Info<< "Reading split from cellZones file " << zoneFile << endl
        << "This requires all"
        << " cells to be in one and only one cellZone." << nl << endl;

    cellZoneMesh newCellZones
    (
      IOobject
      (
        zoneFile,
        mesh.facesInstance(),
        polyMesh::meshSubDir,
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
      ),
      mesh
    );

    labelList newZoneID(mesh.nCells(), -1);
    labelList newNeiZoneID(mesh.nFaces()-mesh.nInternalFaces());
    getZoneID(mesh, newCellZones, newZoneID, newNeiZoneID);

    label unzonedCelli = findIndex(newZoneID, -1);
    if (unzonedCelli != -1)
    {
      FatalErrorInFunction
          << "For the cellZonesFileOnly option all cells "
          << "have to be in a cellZone." << endl
          << "Cell " << unzonedCelli
          << " at" << mesh.cellCentres()[unzonedCelli]
          << " is not in a cellZone. There might be more unzoned cells."
          << exit(FatalError);
    }
    cellRegion = newZoneID;
    nCellRegions = gMax(cellRegion)+1;
    zoneToRegion.setSize(newCellZones.size(), -1);
    regionToZone.setSize(nCellRegions);
    regionNames.setSize(nCellRegions);

    for (label regionI = 0; regionI < nCellRegions; regionI++)
    {
        regionToZone[regionI] = regionI;
        zoneToRegion[regionI] = regionI;
        regionNames[regionI] = newCellZones[regionI].name();
    }
  }
  else
  {
    // Determine connected regions
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Mark additional faces that are blocked
    boolList blockedFace;

    // Read from faceSet
    if (blockedFacesName.size())
    {
      faceSet blockedFaceSet(mesh, blockedFacesName);
      Info<< "Read "
          << returnReduce(blockedFaceSet.size(), sumOp<label>())
          << " blocked faces from set " << blockedFacesName << nl << endl;

      blockedFace.setSize(mesh.nFaces(), false);

      forAllConstIter(faceSet, blockedFaceSet, iter)
      {
        blockedFace[iter.key()] = true;
      }
    }

    // Imply from differing cellZones
    if (useCellZones)
    {
      blockedFace.setSize(mesh.nFaces(), false);

      for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
      {
        label own = mesh.faceOwner()[facei];
        label nei = mesh.faceNeighbour()[facei];

        if (zoneID[own] != zoneID[nei])
        {
          blockedFace[facei] = true;
        }
      }

      // Different cellZones on either side of processor patch.
      forAll(neiZoneID, i)
      {
        label facei = i+mesh.nInternalFaces();

        if (zoneID[mesh.faceOwner()[facei]] != neiZoneID[i])
        {
          blockedFace[facei] = true;
        }
      }
    }

    // Do a topological walk to determine regions
    regionSplit regions(mesh, blockedFace);
    nCellRegions = regions.nRegions();
    cellRegion.transfer(regions);

    // Make up region names. If possible match them to existing zones.
    impVar = matchRegions
    (
      sloppyCellZones,
      mesh,
      nCellRegions,
      cellRegion,

      regionToZone,
      regionNames,
      zoneToRegion
    );

    // Override any default region names if single region selected
    if (largestOnly || insidePoint)
    {
      forAll(regionToZone, regionI)
      {
        if (regionToZone[regionI] == -1)
        {
          if (overwrite)
          {
            regionNames[regionI] = polyMesh::defaultRegion;
          }
          else if (insidePoint)
          {
            regionNames[regionI] = "insidePoint";
          }
          else if (largestOnly)
          {
            regionNames[regionI] = "largestOnly";
          }
        }
      }
    }
  }

  Info<< endl << "Number of regions:" << nCellRegions << nl << endl;


  // Write decomposition to file
  writeCellToRegion(mesh, cellRegion);


  // Sizes per region
  // ~~~~~~~~~~~~~~~~

  labelList regionSizes(nCellRegions, 0);

  forAll(cellRegion, celli)
  {
    regionSizes[cellRegion[celli]]++;
  }
  forAll(regionSizes, regionI)
  {
    reduce(regionSizes[regionI], sumOp<label>());
  }

  Info<< "Region\tCells" << nl
      << "------\t-----" << endl;

  forAll(regionSizes, regionI)
  {
    Info<< regionI << '\t' << regionSizes[regionI] << nl;
  }
  Info<< endl;


  // Print region to zone
  Info<< "Region\tZone\tName" << nl
      << "------\t----\t----" << endl;
  forAll(regionToZone, regionI)
  {
    Info<< regionI << '\t' << regionToZone[regionI] << '\t'
        << regionNames[regionI] << nl;
  }
  Info<< endl;


  // Since we're going to mess with patches and zones make sure all
  // is synchronised
  mesh.boundaryMesh().checkParallelSync(true);
  mesh.faceZones().checkParallelSync(true);


  // Interfaces
  // ----------
  // per interface:
  // - the two regions on either side
  // - the name
  // - the (global) size
  edgeList interfaces;
  List<Pair<word>> interfaceNames;
  labelList interfaceSizes;
  // per face the interface
  labelList faceToInterface;

  getInterfaceSizes
  (
    mesh,
    useFaceZones,
    cellRegion,
    regionNames,

    interfaces,
    interfaceNames,
    interfaceSizes,
    faceToInterface
  );

  Info<< "Sizes of interfaces between regions:" << nl << nl
      << "Interface\tRegion\tRegion\tFaces" << nl
      << "---------\t------\t------\t-----" << endl;

  forAll(interfaces, interI)
  {
    const edge& e = interfaces[interI];

    Info<< interI
        << "\t\t" << e[0] << '\t' << e[1]
        << '\t' << interfaceSizes[interI] << nl;
  }
  Info<< endl;

  if (detectOnly)
  {
    exit(1);//return 0;
  }

  // Read objects in time directory
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IOobjectList objects(mesh, runTime.timeName());

  // Read vol fields.

  PtrList<volScalarField> vsFlds;
  ReadFields(mesh, objects, vsFlds);

  PtrList<volVectorField> vvFlds;
  ReadFields(mesh, objects, vvFlds);

  PtrList<volSphericalTensorField> vstFlds;
  ReadFields(mesh, objects, vstFlds);

  PtrList<volSymmTensorField> vsymtFlds;
  ReadFields(mesh, objects, vsymtFlds);

  PtrList<volTensorField> vtFlds;
  ReadFields(mesh, objects, vtFlds);

  // Read surface fields.

  PtrList<surfaceScalarField> ssFlds;
  ReadFields(mesh, objects, ssFlds);

  PtrList<surfaceVectorField> svFlds;
  ReadFields(mesh, objects, svFlds);

  PtrList<surfaceSphericalTensorField> sstFlds;
  ReadFields(mesh, objects, sstFlds);

  PtrList<surfaceSymmTensorField> ssymtFlds;
  ReadFields(mesh, objects, ssymtFlds);

  PtrList<surfaceTensorField> stFlds;
  ReadFields(mesh, objects, stFlds);

  Info<< endl;


  // Remove any demand-driven fields ('S', 'V' etc)
  mesh.clearOut();


  if (nCellRegions == 1)
  {
    Info<< "Only one region. Doing nothing." << endl;
  }
  else if (makeCellZones)
  {
    Info<< "Putting cells into cellZones instead of splitting mesh."
        << endl;

    // Check if region overlaps with existing zone. If so keep.

    for (label regionI = 0; regionI < nCellRegions; regionI++)
    {
      label zoneI = regionToZone[regionI];

      if (zoneI != -1)
      {
        Info<< "    Region " << regionI << " : corresponds to existing"
            << " cellZone "
            << zoneI << ' ' << cellZones[zoneI].name() << endl;
      }
      else
      {
        // Create new cellZone.
        labelList regionCells = findIndices(cellRegion, regionI);

        word zoneName = "region" + Foam::name(regionI);

        zoneI = cellZones.findZoneID(zoneName);

        if (zoneI == -1)
        {
          zoneI = cellZones.size();
          mesh.cellZones().setSize(zoneI+1);
          mesh.cellZones().set
          (
            zoneI,
            new cellZone
            (
              zoneName,           //name
              regionCells,        //addressing
              zoneI,              //index
              cellZones           //cellZoneMesh
            )
          );
        }
        else
        {
          mesh.cellZones()[zoneI].clearAddressing();
          mesh.cellZones()[zoneI] = regionCells;
        }
        Info<< "    Region " << regionI << " : created new cellZone "
            << zoneI << ' ' << cellZones[zoneI].name() << endl;
      }
    }

    mesh.cellZones().writeOpt() = IOobject::AUTO_WRITE;

    if (!overwrite)
    {
      runTime++;
      mesh.setInstance(runTime.timeName());
    }
    else
    {
      mesh.setInstance(oldInstance);
    }

    Info<< "Writing cellZones as new mesh to time " << runTime.timeName()
      << nl << endl;

    mesh.write();


    // Write cellSets for convenience
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Info<< "Writing cellSets corresponding to cellZones." << nl << endl;

    forAll(cellZones, zoneI)
    {
      const cellZone& cz = cellZones[zoneI];

      cellSet(mesh, cz.name(), cz).write();
    }
  }
  else
  {
    // Add patches for interfaces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Add all possible patches. Empty ones get filtered later on.
    Info<< nl << "Adding patches" << nl << endl;

    labelList interfacePatches
    (
      addRegionPatches
      (
        mesh,
        regionNames,
        interfaces,
        interfaceNames
      )
    );


    if (!overwrite)
    {
      runTime++;
    }


    // Create regions
    // ~~~~~~~~~~~~~~

    if (insidePoint)
    {
      // Will be implemented in future
      const point insidePoint = args.optionRead<point>("insidePoint");

      label regionI = -1;

      (void)mesh.tetBasePtIs();

      label celli = mesh.findCell(insidePoint);

      Info<< nl << "Found point " << insidePoint << " in cell " << celli
                << endl;

      if (celli != -1)
      {
        regionI = cellRegion[celli];
      }

      reduce(regionI, maxOp<label>());

      Info<< nl
          << "Subsetting region " << regionI
          << " containing point " << insidePoint << endl;

      if (regionI == -1)
      {
        FatalErrorInFunction
          << "Point " << insidePoint
          << " is not inside the mesh." << nl
          << "Bounding box of the mesh:" << mesh.bounds()
          << exit(FatalError);
      }

      createAndWriteRegion
      (
        mesh,
        cellRegion,
        regionNames,
        prefixRegion,
        faceToInterface,
        interfacePatches,
        regionI,
        (overwrite ? oldInstance : runTime.timeName())
      );
    }
    else if (largestOnly)
    {
      label regionI = findMax(regionSizes);

      Info<< nl
          << "Subsetting region " << regionI
          << " of size " << regionSizes[regionI]
          << " as named region " << regionNames[regionI] << endl;

      createAndWriteRegion
      (
        mesh,
        cellRegion,
        regionNames,
        prefixRegion,
        faceToInterface,
        interfacePatches,
        regionI,
        (overwrite ? oldInstance : runTime.timeName())
      );
    }
    else
    {
      // Split all
      for (label regionI = 0; regionI < nCellRegions; regionI++)
      {
        Info<< nl
            << "Region " << regionI << nl
            << "-------- " << endl;

        createAndWriteRegion
        (
          mesh,
          cellRegion,
          regionNames,
          prefixRegion,
          faceToInterface,
          interfacePatches,
          regionI,
          (overwrite ? oldInstance : runTime.timeName())
        );
      }
    }
  }

    Info<< "End\n" << endl;
    
return std::make_pair(impVar,nCellRegions);

}

/*
MergeMeshes utility takes master region and slave region names/paths from user
and merges the slave regions to master region one by one.
*/
void MeshManipulationFoam::mergeMeshes(int dirStat, int nDomains)
{
  using namespace Foam;

  int argc = 1;
  char** argv = new char*[2];
  argv[0] = new char[100];
  strcpy(argv[0], "NONE");
  Foam::argList args(argc, argv);
  Foam::Info<< "Create time\n" << Foam::endl;
  Foam::Time runTime(Foam::Time::controlDictName, args);
  Foam::argList::noParallel();

  int ndoms;

  if ((_mshMnipPrms->numDomains) == -1)
  {
    ndoms = nDomains;
  }
  else{
    ndoms = (_mshMnipPrms->numDomains);
  }

  std::vector<std::string> addCases;


  /*if ((_mshMnipPrms->numDomains) == 2)
  {
        addCases.push_back(_mshMnipPrms->addCase);
  }
  else
  {
    for (int i=2; i<(ndoms+1); i++)
    {
      if (i == dirStat)
      {
          i++;
      }
      addCases.push_back("domain" + (std::to_string(i)));
    }
    
    addCases.push_back(_mshMnipPrms->addCase);
  }*/


  // Collecting all domain names for automatic merging of all mesh regions.
  // Directory number obtained from splitMeshRegions is used here. Also
  // number of packs are passed through this function for use in loop. 

  if (ndoms == 2)
  {
    addCases.push_back(_mshMnipPrms->addCase);
  }
  else
  {
    for (int i=2; i<(ndoms+1); i++)
    {
      if (i == dirStat)
      {
          i++;
      }
      addCases.push_back("domain" + (std::to_string(i)));
    }
    
    addCases.push_back(_mshMnipPrms->addCase);
  }

// Main for loop. It loops through all the slave regions and adds them to
// master region one by one. At the end, master region directory will have
// all the slave regions added. Keep overwrite boolean true for this. 
for (int j=0; j<(ndoms-1); j++)
{
  
  const bool overwrite = (_mshMnipPrms->_overwriteMergeMsh);
  word masterRegion = polyMesh::defaultRegion;

  fileName masterCase = (_mshMnipPrms->masterCasePath);
  if (ndoms == 2)
  {
    if (dirStat == 1)
      masterRegion = "domain2";
    else
      masterRegion = (_mshMnipPrms->masterCase);
  }
  else
    masterRegion = (_mshMnipPrms->masterCase);

    fileName addCase = (_mshMnipPrms->addCasePath);
    word addRegion = polyMesh::defaultRegion;
    addRegion = addCases[j];

    getRootCase(masterCase);
    getRootCase(addCase);

    Info<< "Master:      " << masterCase << "  region " << masterRegion << nl
        << "mesh to add: " << addCase    << "  region " << addRegion << endl;

    const fileName masterCasePath = masterCase.path();
    const fileName masterCaseName = masterCase.name();

    Time runTimeMaster
    (
        Time::controlDictName,
        masterCasePath,
        masterCaseName
    );
    runTimeMaster.functionObjects().off();

    const fileName addCasePath = addCase.path();
    const fileName addCaseName = addCase.name();

    Time runTimeToAdd
    (
        Time::controlDictName,
        addCasePath,
        addCaseName
    );
    runTimeToAdd.functionObjects().off();

    Info<< "Reading master mesh for time = " << runTimeMaster.timeName() << nl;

    Info<< "Create mesh\n" << endl;
    Foam::mergePolyMesh masterMesh
    (
        IOobject
        (
            masterRegion,
            runTimeMaster.timeName(),
            runTimeMaster
        )
    );
    const word oldInstance = masterMesh.pointsInstance();


    Info<< "Reading mesh to add for time = " << runTimeToAdd.timeName() << nl;

    Info<< "Create mesh\n" << endl;
    polyMesh meshToAdd
    (
        IOobject
        (
            addRegion,
            runTimeToAdd.timeName(),
            runTimeToAdd
        )
    );

    if (!(_mshMnipPrms->_overwriteMergeMsh))
    {
        runTimeMaster++;
    }

    Info<< "Writing combined mesh to " << runTimeMaster.timeName() << endl;

    masterMesh.addMesh(meshToAdd);
    masterMesh.merge();

    if ((_mshMnipPrms->_overwriteMergeMsh))
    {
        masterMesh.setInstance(oldInstance);
    }

    masterMesh.write();

}

    Info<< "End\n" << endl;

}

/*
CreatePatchDict function creates dictionary file for createPatch utility and 
puts them in defined path (i.e system/domainX). Currently this function has 
customizations for Pack Mesh Generation service. More general methods could be
added in future.
*/
void MeshManipulationFoam::createPatchDict(int dirStat)
{
  // creating a base system directory
    const char dir_path[] = "./system/domain0";
  boost::filesystem::path dir(dir_path);
    try
    {
      boost::filesystem::create_directory(dir);
    }
    catch (boost::filesystem::filesystem_error &e)
    {
    std::cerr << "Problem in creating system directory for the snappyHexMesh" << "\n";
        std::cerr << e.what() << std::endl;
        throw;
  }

    // creating mesh dictionary file
    std::ofstream contDict;
    contDict.open(std::string(dir_path)+"/createPatchDict");
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
    object    createPatchDict;\n\
}\n\n";

    contText = contText + 
    "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";

    contText = contText + "\n\npointSync true;\n";
    contText = contText + "\n\npatches\n";
    contText = contText + "(\n";
    contText = contText + "\t{\n";
    contText = contText + "\t\tname " + (_mshMnipPrms->surroundingName)
                        + ";\n";
    contText = contText + "\n\t\tpatchInfo\n";
    contText = contText + "\t\t{\n";
    contText = contText + "\t\t\ttype " + (_mshMnipPrms->srrndngPatchType)
                        + ";\n";
    contText = contText + "\t\t}\n";
    contText = contText + "\n\t\tconstructFrom patches;\n";
    contText = contText + "\n\t\tpatches (\"domain0_to_domain.*\");";
    contText = contText + "\n\t}\n";
    contText = contText + ");\n";

    contText = contText + 
    "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//\n\n";

    contDict << contText;
    contDict.close();


  if (dirStat == 1)
  {
    const char dir_path2[] = "./system/domain2";
    boost::filesystem::path dir2(dir_path2);
  try
  {
    boost::filesystem::create_directory(dir2);
  }
  catch (boost::filesystem::filesystem_error &e)
  {
    std::cerr << "Problem in creating system directory for createPatch" << "\n";
    std::cerr << e.what() << std::endl;
    throw;
  }
    // creating mesh dictionary file
    std::ofstream contDict2;
    contDict2.open(std::string(dir_path2)+"/createPatchDict");
    // header
    std::string contText2=
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
    object    createPatchDict;\n\
}\n\n";

    contText2 = contText2 + 
    "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";

    contText2 = contText2 + "\n\npointSync true;\n";
    contText2 = contText2 + "\n\npatches\n";
    contText2 = contText2 + "(\n";
    contText2 = contText2 + "\t{\n";
    contText2 = contText2 + "\t\tname " + (_mshMnipPrms->packsName) + ";\n";
    contText2 = contText2 + "\n\t\tpatchInfo\n";
    contText2 = contText2 + "\t\t{\n";
    contText2 = contText2 + "\t\t\ttype " + (_mshMnipPrms->packsPatchType)
                          + ";\n";
    contText2 = contText2 + "\t\t}\n";
    contText2 = contText2 + "\n\t\tconstructFrom patches;\n";
    contText2 = contText2 + "\n\t\tpatches (\"domain.*_to_domain0\");";
    contText2 = contText2 + "\n\t}\n";
    contText2 = contText2 + ");\n";

    contText2 = contText2 + 
    "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";

    contDict2 << contText2;
    contDict2.close();
  }
  else
  {
    const char dir_path2[] = "./system/domain1";
    boost::filesystem::path dir2(dir_path2);
  try
  {
    boost::filesystem::create_directory(dir2);
  }
  catch (boost::filesystem::filesystem_error &e)
  {
    std::cerr << "Problem in creating system directory for createPatch" << "\n";
    std::cerr << e.what() << std::endl;
    throw;
  }
    // creating mesh dictionary file
    std::ofstream contDict2;
    contDict2.open(std::string(dir_path2)+"/createPatchDict");
    // header
    std::string contText2=
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
    object    createPatchDict;\n\
}\n\n";

    contText2 = contText2 + 
    "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";

    contText2 = contText2 + "\n\npointSync true;\n";
    contText2 = contText2 + "\n\npatches\n";
    contText2 = contText2 + "(\n";
    contText2 = contText2 + "\t{\n";
    contText2 = contText2 + "\t\tname " + (_mshMnipPrms->packsName) + ";\n";
    contText2 = contText2 + "\n\t\tpatchInfo\n";
    contText2 = contText2 + "\t\t{\n";
    contText2 = contText2 + "\t\t\ttype " + (_mshMnipPrms->packsPatchType)
                          + ";\n";
    contText2 = contText2 + "\t\t}\n";
    contText2 = contText2 + "\n\t\tconstructFrom patches;\n";
    contText2 = contText2 + "\n\t\tpatches (\"domain.*_to_domain0\");";
    contText2 = contText2 + "\n\t}\n";
    contText2 = contText2 + ");\n";

    contText2 = contText2 + 
    "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";

    contDict2 << contText2;
    contDict2.close();
  }

  
}

/*
CreatePatch utility merges multiple patches of certain mesh in a single patch.
It reads createPatchDict from system/mesh_reg_name/createPatchDict. Currently
the utility implementation is wrapped with for loop to execute it twice for
Pack Mesh Generation service.
*/
void MeshManipulationFoam::createPatch(int dirStat)
{
  using namespace Foam;

  int argc = 1;
  char** argv = new char*[2];
  argv[0] = new char[100];
  strcpy(argv[0], "NONE");
  Foam::argList args(argc, argv);
  Foam::Info<< "Create time\n" << Foam::endl;
  Foam::argList::noParallel();

  createPatchDict(dirStat);

  Time runTime
  (
      Time::controlDictName,
      "",
      ""
  );
  

  std::string firstPtch = _mshMnipPrms->pathSurrounding;
  std::string secondPtch;
  
  if (dirStat == 1)
    secondPtch = "domain2";
  else
    secondPtch = _mshMnipPrms->pathPacks;
  
  std::string regnName;

// This loop wraps around utility execution to merge patches of all packs and 
// surounding regions into two distinct patches.
for (int i=0; i<2; i++)
{
  if (i == 0)
  {
    regnName = firstPtch;
  }
  if (i == 1)
  {
    regnName = secondPtch;
  }
  const bool overwrite = (_mshMnipPrms->_overwritecpMsh);

    Foam::word meshRegionName = regnName;

    Foam::word regionName;
    regionName = regnName;
    Foam::Info
        << "Create polyMesh for time = "
        << runTime.timeName() << Foam::nl << Foam::endl;

    Foam::polyMesh mesh
  (
      Foam::IOobject
      (
          regionName,
          runTime.timeName(),
          runTime,
          Foam::IOobject::MUST_READ
      )
  );

    const word oldInstance = mesh.pointsInstance();

    const word dictName("createPatchDict");

  IOobject dictIO
  (
      dictName,
      runTime.system(),
      mesh,
      IOobject::MUST_READ_IF_MODIFIED,
      IOobject::NO_WRITE
  );

  Info<< "Reading " << dictName << nl << endl;

  IOdictionary dict(dictIO);

  // Whether to synchronise points
  const Switch pointSync(dict.lookup("pointSync"));


  const polyBoundaryMesh& patches = mesh.boundaryMesh();

  // If running parallel check same patches everywhere
  patches.checkParallelSync(true);


  dumpCyclicMatch("initial_", mesh);

  // Read patch construct info from dictionary
  PtrList<dictionary> patchSources(dict.lookup("patches"));

  HashSet<word> addedPatchNames;
  forAll(patchSources, addedI)
  {
    const dictionary& dict = patchSources[addedI];
    addedPatchNames.insert(dict.lookup("name"));
  }


  // 1. Add all new patches
  // ~~~~~~~~~~~~~~~~~~~~~~

  if (patchSources.size())
  {
    // Old and new patches.
    DynamicList<polyPatch*> allPatches(patches.size()+patchSources.size());

    label startFacei = mesh.nInternalFaces();

    // Copy old patches.
    forAll(patches, patchi)
    {
      const polyPatch& pp = patches[patchi];

      if (!isA<processorPolyPatch>(pp))
      {
        allPatches.append
        (
          pp.clone
          (
            patches,
            patchi,
            pp.size(),
            startFacei
          ).ptr()
        );
        startFacei += pp.size();
      }
    }

    forAll(patchSources, addedI)
    {
      const dictionary& dict = patchSources[addedI];

      word patchName(dict.lookup("name"));

      label destPatchi = patches.findPatchID(patchName);

      if (destPatchi == -1)
      {
        dictionary patchDict(dict.subDict("patchInfo"));

        destPatchi = allPatches.size();

        Info<< "Adding new patch " << patchName
            << " as patch " << destPatchi
            << " from " << patchDict << endl;

        patchDict.set("nFaces", 0);
        patchDict.set("startFace", startFacei);

        // Add an empty patch.
        allPatches.append
        (
          polyPatch::New
          (
            patchName,
            patchDict,
            destPatchi,
            patches
          ).ptr()
        );
      }
      else
      {
        Info<< "Patch '" << patchName << "' already exists.  Only "
            << "moving patch faces - type will remain the same" << endl;
      }
    }

    // Copy old patches.
    forAll(patches, patchi)
    {
      const polyPatch& pp = patches[patchi];

      if (isA<processorPolyPatch>(pp))
      {
        allPatches.append
        (
          pp.clone
          (
            patches,
            patchi,
            pp.size(),
            startFacei
          ).ptr()
        );
        startFacei += pp.size();
      }
    }

    allPatches.shrink();
    mesh.removeBoundary();
    mesh.addPatches(allPatches);
    Info<< endl;
  }



  // 2. Repatch faces
  // ~~~~~~~~~~~~~~~~

  polyTopoChange meshMod(mesh);


  forAll(patchSources, addedI)
  {
    const dictionary& dict = patchSources[addedI];

    const word patchName(dict.lookup("name"));
    label destPatchi = patches.findPatchID(patchName);

    if (destPatchi == -1)
    {
      FatalErrorInFunction
        << "patch " << patchName << " not added. Problem."
        << abort(FatalError);
    }

    const word sourceType(dict.lookup("constructFrom"));

    if (sourceType == "patches")
    {
      labelHashSet patchSources
      (
        patches.patchSet
        (
          wordReList(dict.lookup("patches"))
        )
      );

      // Repatch faces of the patches.
      forAllConstIter(labelHashSet, patchSources, iter)
      {
        const polyPatch& pp = patches[iter.key()];

        Info<< "Moving faces from patch " << pp.name()
            << " to patch " << destPatchi << endl;

        forAll(pp, i)
        {
          changePatchID
          (
            mesh,
            pp.start() + i,
            destPatchi,
            meshMod
          );
        }
      }
    }
    else if (sourceType == "set")
    {
      const word setName(dict.lookup("set"));

      faceSet faces(mesh, setName);

      Info<< "Read " << returnReduce(faces.size(), sumOp<label>())
          << " faces from faceSet " << faces.name() << endl;

      // Sort (since faceSet contains faces in arbitrary order)
      labelList faceLabels(faces.toc());

      SortableList<label> patchFaces(faceLabels);

      forAll(patchFaces, i)
      {
        label facei = patchFaces[i];

        if (mesh.isInternalFace(facei))
        {
          FatalErrorInFunction
            << "Face " << facei << " specified in set "
            << faces.name()
            << " is not an external face of the mesh." << endl
            << "This application can only repatch existing boundary"
            << " faces." << exit(FatalError);
        }

        changePatchID
        (
          mesh,
          facei,
          destPatchi,
          meshMod
        );
      }
    }
    else
    {
      FatalErrorInFunction
        << "Invalid source type " << sourceType << endl
        << "Valid source types are 'patches' 'set'" << exit(FatalError);
    }
  }
  Info<< endl;


  // Change mesh, use inflation to reforce calculation of transformation
  // tensors.
  Info<< "Doing topology modification to order faces." << nl << endl;
  autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, true);
  mesh.movePoints(map().preMotionPoints());

  dumpCyclicMatch("coupled_", mesh);

    // Synchronise points.
  if (!pointSync)
  {
      Info<< "Not synchronising points." << nl << endl;
  }
  else
  {
    Info<< "Synchronising points." << nl << endl;

    // This is a bit tricky. Both normal and position might be out and
    // current separation also includes the normal
    // ( separation_ = (nf&(Cr - Cf))*nf ).

    // For cyclic patches:
    // - for separated ones use user specified offset vector

    forAll(mesh.boundaryMesh(), patchi)
    {
      const polyPatch& pp = mesh.boundaryMesh()[patchi];

      if (pp.size() && isA<coupledPolyPatch>(pp))
      {
        const coupledPolyPatch& cpp =
          refCast<const coupledPolyPatch>(pp);

        if (cpp.separated())
        {
          Info<< "On coupled patch " << pp.name()
              << " separation[0] was "
              << cpp.separation()[0] << endl;

          if (isA<cyclicPolyPatch>(pp) && pp.size())
          {
            const cyclicPolyPatch& cycpp =
              refCast<const cyclicPolyPatch>(pp);

            if (cycpp.transform() == cyclicPolyPatch::TRANSLATIONAL)
            {
              // Force to wanted separation
              Info<< "On cyclic translation patch " << pp.name()
                << " forcing uniform separation of "
                << cycpp.separationVector() << endl;
              const_cast<vectorField&>(cpp.separation()) =
                pointField(1, cycpp.separationVector());
            }
            else
            {
              const cyclicPolyPatch& nbr = cycpp.neighbPatch();
              const_cast<vectorField&>(cpp.separation()) =
              pointField
              (
                1,
                nbr[0].centre(mesh.points())
                  - cycpp[0].centre(mesh.points())
              );
            }
          }
          Info<< "On coupled patch " << pp.name()
              << " forcing uniform separation of "
              << cpp.separation() << endl;
        }
        else if (!cpp.parallel())
        {
          Info<< "On coupled patch " << pp.name()
              << " forcing uniform rotation of "
              << cpp.forwardT()[0] << endl;

          const_cast<tensorField&>
            (
              cpp.forwardT()
            ).setSize(1);
            const_cast<tensorField&>
            (
              cpp.reverseT()
            ).setSize(1);

            Info<< "On coupled patch " << pp.name()
                << " forcing uniform rotation of "
                << cpp.forwardT() << endl;
          }
        }
    }

      Info<< "Synchronising points." << endl;

      pointField newPoints(mesh.points());

      syncPoints
      (
        mesh,
        newPoints,
        minMagSqrEqOp<vector>(),
        point(GREAT, GREAT, GREAT)
      );

      scalarField diff(mag(newPoints-mesh.points()));
      Info<< "Points changed by average:" << gAverage(diff)
          << " max:" << gMax(diff) << nl << endl;

      mesh.movePoints(newPoints);
  }

    // 3. Remove zeros-sized patches
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Info<< "Removing patches with no faces in them." << nl<< endl;
    filterPatches(mesh, addedPatchNames);


    dumpCyclicMatch("final_", mesh);


    // Set the precision of the points data to 10
    IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));


    if (!overwrite)
    {
        runTime++;
    }
    else
    {
        mesh.setInstance(oldInstance);
    }

    // Write resulting mesh
    Info<< "Writing repatched mesh to " << runTime.timeName() << nl << endl;
    mesh.write();
}
  
    Info<< "End\n" << endl;

}

/*
This utility converts OpenFoam mesh into STL file. Filename is user-input.
It can take paths for output file location and name.
*/
void MeshManipulationFoam::foamToSurface()
{
  using namespace Foam;
    
    int argc = 1;
    char** argv = new char*[2];
    argv[0] = new char[100];
    strcpy(argv[0], "NONE");
    Foam::argList args(argc, argv);
    Foam::Info<< "Create time\n" << Foam::endl;
    Foam::Time runTime(Foam::Time::controlDictName, args);
    Foam::argList::noParallel();
  
  fileName exportName = (_mshMnipPrms->outSurfName);
  
  scalar scaleFactor = 0;
    const bool doTriangulate = true;

    fileName exportBase = exportName.lessExt();
    word exportExt = exportName.ext();

    Time runTimeMaster
    (
        Time::controlDictName,
        ".",
        "."
    );
  
  polyMesh mesh
    (
        IOobject
        (
            Foam::polyMesh::defaultRegion,
            runTimeMaster.timeName(),
            runTimeMaster,
            Foam::IOobject::MUST_READ
        )
    );
  
  polyMesh::readUpdateState state = mesh.readUpdate();
  exportName = exportBase + "." + exportExt;
  
  meshedSurface surf(mesh.boundaryMesh());
    surf.scalePoints(scaleFactor);

    Info<< "writing " << exportName;
  
    if (doTriangulate)
    {
        Info<< " triangulated";
        surf.triangulate();
    }

    if (scaleFactor <= 0)
    {
        Info<< " without scaling" << endl;
    }
    else
    {
        Info<< " with scaling " << scaleFactor << endl;
    }
  
  // Add an error handler for directory check before writing
    surf.write(exportName);
  //surf.write("./constant/triSurface/InputPacks.stl");

    Info<< nl << endl;
  
  Info<< "End\n" << endl;
  
}


/*
surfaceSplitByTopology is a surface file manipulation utility which aims at
splitting multiple disconnected regions in geometry into separate surfaces and
outputs a single surface file containing all regions as well as separate STL
files for all regions. User inputs are input file name and output file name. 
*/
int MeshManipulationFoam::surfSpltByTopology()
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

    fileName surfFileName(_mshMnipPrms->surfFile);
    Info<< "Reading surface from " << surfFileName << endl;

    fileName outFileName(_mshMnipPrms->outSurfFile);
    fileName outFileBaseName = outFileName.lessExt();
    word outExtension = outFileName.ext();

    // Load surface
    triSurface surf(surfFileName);

    bool anyZoneRemoved = false;

    label iterationNo = 0;
    label iterationLimit = 10;

    Info<< "Splitting off baffle parts " << endl;

    do
    {
        anyZoneRemoved = false;

        labelList faceZone;

        const labelListList& edFaces = surf.edgeFaces();
        const labelListList& faceEds = surf.faceEdges();

        boolList multipleEdges(edFaces.size(), false);

        forAll(multipleEdges, i)
        {
            if (edFaces[i].size() > 2)
            {
                multipleEdges[i] = true;
            }
        }

        label nZones = surf.markZones(multipleEdges, faceZone);

        if (nZones < 2)
        {
            break;
        }

        boolList nonBaffle(faceZone.size(), true);
        boolList baffle(faceZone.size(), true);
        labelList pointMap;
        labelList faceMap;


        for (label z = 0; z < nZones; z++)
        {
            bool keepZone = true;

            forAll(faceZone, f)
            {
                if (faceZone[f] == z)
                {
                    forAll(faceEds[f], fe)
                    {
                        if (edFaces[faceEds[f][fe]].size() < 2)
                        {
                            keepZone = false;

                            anyZoneRemoved = true;

                            break;
                        }
                    }
                }

                if (!keepZone)
                {
                    break;
                }
            }

            forAll(faceZone, f)
            {
                if (faceZone[f] == z)
                {
                    nonBaffle[f] = keepZone;
                    baffle[f] = !keepZone;
                }
            }
        }

        Info<< "    Iteration " << iterationNo << endl;

        triSurface baffleSurf = surf.subsetMesh(baffle, pointMap, faceMap);

        if (baffleSurf.size())
        {
            fileName bafflePartFileName =
            outFileBaseName
          + "_bafflePart_"
          + name(iterationNo)
          + "." + outExtension;

            Info<< "    Writing baffle part to " << bafflePartFileName << endl;

            baffleSurf.write(bafflePartFileName);
        }

        surf = surf.subsetMesh(nonBaffle, pointMap, faceMap);

        if (iterationNo == iterationLimit)
        {
            WarningInFunction
            << "Iteration limit of " << iterationLimit << "reached" << endl;
        }

        iterationNo++;

    } while (anyZoneRemoved && iterationNo < iterationLimit);

    Info<< "Writing new surface to " << outFileName << endl;

    surf.write(outFileName);

    labelList faceZone;

    const labelListList& edFaces = surf.edgeFaces();

    boolList multipleEdges(edFaces.size(), false);

    forAll(multipleEdges, i)
    {
        if (edFaces[i].size() > 2)
        {
            multipleEdges[i] = true;
        }
    }

    label nZones = surf.markZones(multipleEdges, faceZone);

    int nofPacks = nZones;

    // nZones is number of pack domains. it will be great to get that out for mergeMeshes.
    // Also, surface renumbering is needed as output STL

    /*Info<< "Splitting remaining multiply connected parts" << endl;

    for (label z = 0; z < nZones; z++)
    {

        boolList include(faceZone.size(), false);
        labelList pointMap;
        labelList faceMap;

        forAll(faceZone, f)
        {
            if (faceZone[f] == z)
            {
                include[f] = true;
            }
        }

        triSurface zoneSurf = surf.subsetMesh(include, pointMap, faceMap);


        fileName remainingPartFileName =
            outFileBaseName
          + "_multiplePart_"
          + name(z)
          + "." + outExtension;

        Info<< "    Writing mulitple part "
            << z << " to " << remainingPartFileName << endl;

        zoneSurf.write(remainingPartFileName);
    }*/

    Info << "End\n" << endl;

    return nofPacks;
}

void MeshManipulationFoam::createControlDict()
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


// Feature in development currently

// Adds connectivity information at shared interfaces.
void MeshManipulationFoam::addCohesiveElements(double tol, const std::string outName)
{
  // Initializing FOAM
  int argc = 1;
  char** argv = new char*[2];
  argv[0] = new char[100];
  strcpy(argv[0], "NONE");
  Foam::argList args(argc, argv);
  Foam::argList::noParallel();

  Time runTime
  (
      Time::controlDictName,
      ".",
      "."
  );

  Foam::word regionName = Foam::polyMesh::defaultRegion;

  auto _fmeshSurr = new Foam::polyMesh
  (
    Foam::IOobject
    (
      "domain0",
      runTime.timeName(),
      runTime,
      Foam::IOobject::MUST_READ
    )
  );

  // New query points (From Packs)
  auto _fmeshPacks = new Foam::polyMesh
  (
    Foam::IOobject
    (
      "domain1",
      runTime.timeName(),
      runTime,
      Foam::IOobject::MUST_READ
    )
  );

  Foam::pointField packsPF = _fmeshPacks->points();

  // declare vtk datasets
  vtkSmartPointer<vtkUnstructuredGrid> dataSetSurr
      = vtkSmartPointer<vtkUnstructuredGrid>::New();


  // creating equivalent vtk topology from fvMesh
  // by default polyhedral cells will be decomposed to 
  // tets and pyramids. Additional points will be added
  // to underlying fvMesh.
  std::cout << "Performing topological decomposition.\n";
  Foam::vtkTopo topo1(*_fmeshSurr);
  Foam::vtkTopo topo2(*_fmeshPacks);

  // point data for surrounding
  Foam::pointField pfSurr = _fmeshSurr->points();
  vtkSmartPointer<vtkPoints> points 
      = vtkSmartPointer<vtkPoints>::New(); 
  for (int ipt=0; ipt<_fmeshSurr->nPoints(); ipt++)
      points->InsertNextPoint(
              pfSurr[ipt].x(), 
              pfSurr[ipt].y(), 
              pfSurr[ipt].z() 
              );

  // point data for packs
  int impPts = (_fmeshSurr->nPoints());
  for (int ipt=0; ipt<_fmeshPacks->nPoints(); ipt++)
      points->InsertNextPoint(
              packsPF[ipt].x(), 
              packsPF[ipt].y(), 
              packsPF[ipt].z() 
              );
  

  dataSetSurr->SetPoints(points);

  // cell data for surrounding
  std::vector<int> pntIds;
  int nCelPnts = 0;
  for (int icl=0; icl<topo1.vertLabels().size(); icl++)
  {
      nCelPnts = topo1.vertLabels()[icl].size();
      pntIds.resize(nCelPnts, -1);
      for (int ip=0; ip< nCelPnts; ip++)
      {
          pntIds[ip] = topo1.vertLabels()[icl][ip];
      }
      createVtkCell(dataSetSurr, topo1.cellTypes()[icl], pntIds);
  }


  // cell data for packs
  int nCelPnts2;
  for (int icl=0; icl<topo2.vertLabels().size(); icl++)
  {
      std::vector<int> pntIds2;
      nCelPnts2 = topo2.vertLabels()[icl].size();
      pntIds2.resize(nCelPnts2, -1);
      for (int ip=0; ip< nCelPnts2; ip++)
      {
        pntIds2[ip] = impPts + topo2.vertLabels()[icl][ip];
      }
      createVtkCell(dataSetSurr, topo2.cellTypes()[icl], pntIds2);
  }

  int numCells = dataSetSurr->GetNumberOfCells();
  int numPoints = dataSetSurr->GetNumberOfPoints();


  int nDim = 3;
  ANNpointArray pntCrd;
  pntCrd = annAllocPts(numPoints, nDim);
  for (int iPnt=0; iPnt < numPoints; iPnt++)
  {
    std::vector<double> getPt = std::vector<double>(3);
    dataSetSurr->GetPoint(iPnt, &getPt[0]);
    pntCrd[iPnt][0] = getPt[0];
    pntCrd[iPnt][1] = getPt[1];
    pntCrd[iPnt][2] = getPt[2];
  }

  ANNkd_tree *kdTree = new ANNkd_tree(pntCrd, numPoints, nDim);

  // Finding duplicate nodes
  double rad = 1e-05;
  std::map<int, int> dupNdeMap;
  ANNpoint qryPnt; // Query point.
  int nNib = 2; // Number of neighbours to return (including query pt. itself).
  ANNidxArray nnIdx = new ANNidx[nNib]; // Nearest neighbour ID array.
  ANNdistArray dists = new ANNdist[nNib]; // Neighbour distance array.
  qryPnt = annAllocPt(nDim); // Initializing query point

  for (int iNde = 0; iNde < impPts; iNde++) {
    std::vector<double> getPt = std::vector<double>(3);
    dataSetSurr->GetPoint(iNde, &getPt[0]);
    qryPnt[0] = getPt[0];
    qryPnt[1] = getPt[1];
    qryPnt[2] = getPt[2];
    kdTree->annkFRSearch(qryPnt, rad, nNib, nnIdx, dists, 0);
    if (dists[1] <= tol)
    {
      if ((nnIdx[1]) >= (_fmeshSurr->nPoints()))
        dupNdeMap[nnIdx[0]] = nnIdx[1]; // Format Surrounding::Pack
      else
        dupNdeMap[nnIdx[1]] = nnIdx[0]; // Format Surrounding::Pack
    }
  }

  int numDupPnts = dupNdeMap.size();
  std::cout << "Found " << numDupPnts << " duplicate nodes.\n";


  // Getting cells in packs
  std::map<int, int>::iterator it = dupNdeMap.begin();
  std::vector<int> cellID;
  for (int i=0; i<(dupNdeMap.size()); i++)
  {
    int nCells, cellNum;
    vtkIdList *cells = vtkIdList::New();

    // get cells the point belongs to
    dataSetSurr->GetPointCells( it->second, cells );
    nCells = cells->GetNumberOfIds();

    for( cellNum=0; cellNum<nCells; cellNum++ )
    {
      cellID.push_back(cells->GetId( cellNum ));  // get cell id from the list
    }
    it++;
  }

  sort( cellID.begin(), cellID.end() );
  cellID.erase( unique( cellID.begin(), cellID.end() ), cellID.end() );

  int know = cellID.size();

  //std::cout << "Total Cells Found Are " << know << std::endl;


  // Determines which cells are useful
  std::vector<int> newCellIds;
  std::vector<int> debugCellIds;

  for (int i=0; i<know; i++)
  {
    // Getting all cell defining points
    std::vector<int> ptIdzz(8);
    int cntr = 0;
    vtkCell *cell;
    vtkPoints *fp;
    cell = dataSetSurr->GetCell( cellID[i] ); 

    vtkIdList *pts = cell->GetPointIds();

    for (int j=0; j<8; j++)
    {
      ptIdzz[j] = pts->GetId(j);
    }

    // Checking if any cell contains any 4 or more of duplicate nodes.
    std::map<int, int>::iterator it2 = dupNdeMap.begin();
    while(it2 != dupNdeMap.end())
    {
      for (int k=0; k<8; k++)
      {
        if ((it2->second) == ptIdzz[k])
        {
          cntr++;
        }
        else
        {
          // Nothing
        }
      }
      it2++;
    }

    if (cntr >= 4)
    {
      newCellIds.push_back(cellID[i]);
    }
    else
    {
      debugCellIds.push_back(cellID[i]);
    }

  }


  int size2 = newCellIds.size();

  // Getting useful faces from cells
  std::vector<int> globalPtIds;
  std::vector<int> surroundingArray;
  std::vector<int> packArray;
  for (int i=0; i<size2; i++)
  {
    vtkCell *cell;
    cell = dataSetSurr->GetCell(newCellIds[i]);
    vtkIdList *pts = cell->GetPointIds();
    
    for (int j=0; j<6; j++)
    {
      int isFour = 0;
      int* ptFaces = nullptr;
      vtkCell3D* cell3d = static_cast<vtkCell3D*>( cell );
      cell3d->GetFacePoints(j,ptFaces);
      std::vector<int> keysDupMap = std::vector<int>(4);
      std::map<int, int>::iterator it3 = dupNdeMap.begin();
      while (it3 != dupNdeMap.end())
      {
        for (int h=0; h<4; h++)
        {
          if ((it3->second) == (pts->GetId(ptFaces[h])))
          {
            isFour++;
            keysDupMap[h] = it3->first;
          }
          else
          {
            // Nothing
          }
        }

        it3++;
      }

      if (isFour == 4)
      {
        for (int k=0; k<4; k++)
        {
          globalPtIds.push_back(keysDupMap[k]);
        }
      }
    }
  }

  std::map<int, int>::iterator it5;

  // Creating node map with sequence.
  std::unordered_multimap<int, int> cohesiveMap;

  for (int i=0; i<globalPtIds.size(); i++)
  {
    it5 = dupNdeMap.find(globalPtIds[i]);

    if (it5 == dupNdeMap.end())
    {
      // Nothing
    }
    else
    {
      surroundingArray.push_back(it5->first);
      packArray.push_back(it5->second);
    }
  }

  //std::cout << "Size = " << packArray.size() << std::endl;
  int newCells = packArray.size()/4;
  int startNum = dataSetSurr->GetNumberOfCells();

  // Finally, creating VTK cells
  int it7 = 0;
  std::vector<int> pntCohesiveIds;
  int nCelCohesivePnts = 0;
  for (int icl=0; icl<newCells; icl++)
  {
      nCelCohesivePnts = 4;
      pntCohesiveIds.resize((nCelCohesivePnts*2), -1);
      for (int ip=0; ip< nCelCohesivePnts; ip++)
      {
          pntCohesiveIds[ip] = packArray[it7];
          pntCohesiveIds[ip+4] = surroundingArray[it7];
          it7++;
      }
      createVtkCell(dataSetSurr, 12, pntCohesiveIds);
  }

  int endNum = dataSetSurr->GetNumberOfCells();

  // ****************************** //
  // Zero volume cells visualization method for debugging purpose.
  std::vector<int> cohCellIds;
  for (int i=startNum; i<endNum; i++)
  {
    cohCellIds.push_back(i);
  }

  int size3 = cohCellIds.size();
  // Takes cell Ids as input and writes VTK mesh
  vtkSmartPointer<vtkUnstructuredGrid> visDataSet
      = vtkSmartPointer<vtkUnstructuredGrid>::New();

  std::vector<int> ptIdzz2;
  for (int i=0; i<size3; i++){
    
    vtkCell *cell;
    vtkPoints *fp;
    cell = dataSetSurr->GetCell( cohCellIds[i] );

    vtkIdList *pts = cell->GetPointIds();

    for (int j=0; j<8; j++)
    {
      ptIdzz2.push_back(pts->GetId(j));
    }

  }

  int lvar = 0;

  vtkSmartPointer<vtkPoints> pointsViz 
        = vtkSmartPointer<vtkPoints>::New();

  for (int i=0; i<(size3*8); i++)
  {
    std::vector<double> getPt = std::vector<double>(3);
    dataSetSurr->GetPoint(ptIdzz2[i], &getPt[0]);
    pointsViz->InsertNextPoint(getPt[0],getPt[1],getPt[2]);
  }

  visDataSet->SetPoints(pointsViz);

  

  for (int i=0; i<size3; i++)
  {
    std::vector<int> ptnewIdz;
    for (int j=0; j<8; j++)
    {
      ptnewIdz.push_back(lvar);
      lvar++;
    }
    createVtkCell(visDataSet, 12, ptnewIdz);
  }


  vtkMesh* vm = new vtkMesh(visDataSet,"CohesiveElements.vtu");
  vm->report();
  vm->write();

  // ********************************************* //
  vtkMesh* vm2 = new vtkMesh(dataSetSurr,outName);
  vm2->report();
  vm2->write();

  if (vm)
    delete vm;
  if (vm2)
    delete vm2;

  bool tetra = false;
  std::string ofnameTet = "Tetra" + outName;

  if (tetra == true)
  {
    // create meshBase object
    std::shared_ptr<meshBase> myMesh = meshBase::CreateShared(outName);

    // Converts hex mesh to tet mesh and writes in VTU file.
    myMesh->convertHexToTetVTK(dataSetSurr);
    myMesh->report();
    myMesh->write(ofnameTet);
  }
  // ******************************************************************************* //

}

void MeshManipulationFoam::addArtificialThicknessElements
(
  double tol,
  const std::string outName,
  double thickness
)
{
  // Initializing FOAM
  int argc = 1;
  char** argv = new char*[2];
  argv[0] = new char[100];
  strcpy(argv[0], "NONE");
  Foam::argList args(argc, argv);
  Foam::argList::noParallel();

  Time runTime
  (
      Time::controlDictName,
      ".",
      "."
  );

  Foam::word regionName = Foam::polyMesh::defaultRegion;

  auto _fmeshSurr = new Foam::polyMesh
  (
    Foam::IOobject
    (
      "domain0",
      runTime.timeName(),
      runTime,
      Foam::IOobject::MUST_READ
    )
  );

  // New query points (From Packs)
  auto _fmeshPacks = new Foam::polyMesh
  (
    Foam::IOobject
    (
      "domain1",
      runTime.timeName(),
      runTime,
      Foam::IOobject::MUST_READ
    )
  );

  Foam::pointField packsPF = _fmeshPacks->points();

  // declare vtk datasets
  vtkSmartPointer<vtkUnstructuredGrid> dataSetSurr
      = vtkSmartPointer<vtkUnstructuredGrid>::New();


  // creating equivalent vtk topology from fvMesh
  // by default polyhedral cells will be decomposed to 
  // tets and pyramids. Additional points will be added
  // to underlying fvMesh.
  std::cout << "Performing topological decomposition.\n";
  Foam::vtkTopo topo1(*_fmeshSurr);
  Foam::vtkTopo topo2(*_fmeshPacks);

  // point data for surrounding
  Foam::pointField pfSurr = _fmeshSurr->points();
  vtkSmartPointer<vtkPoints> points 
      = vtkSmartPointer<vtkPoints>::New(); 
  for (int ipt=0; ipt<_fmeshSurr->nPoints(); ipt++)
      points->InsertNextPoint(
              pfSurr[ipt].x(), 
              pfSurr[ipt].y(), 
              pfSurr[ipt].z() 
              );

  // point data for packs
  int impPts = (_fmeshSurr->nPoints());
  for (int ipt=0; ipt<_fmeshPacks->nPoints(); ipt++)
      points->InsertNextPoint(
              packsPF[ipt].x(), 
              packsPF[ipt].y(), 
              packsPF[ipt].z() 
              );
  

  dataSetSurr->SetPoints(points);

  // cell data for surrounding
  std::vector<int> pntIds;
  int nCelPnts = 0;
  for (int icl=0; icl<topo1.vertLabels().size(); icl++)
  {
      nCelPnts = topo1.vertLabels()[icl].size();
      pntIds.resize(nCelPnts, -1);
      for (int ip=0; ip< nCelPnts; ip++)
      {
          pntIds[ip] = topo1.vertLabels()[icl][ip];
      }
      createVtkCell(dataSetSurr, topo1.cellTypes()[icl], pntIds);
  }


  // cell data for packs
  int nCelPnts2;
  for (int icl=0; icl<topo2.vertLabels().size(); icl++)
  {
      std::vector<int> pntIds2;
      nCelPnts2 = topo2.vertLabels()[icl].size();
      pntIds2.resize(nCelPnts2, -1);
      for (int ip=0; ip< nCelPnts2; ip++)
      {
        pntIds2[ip] = impPts + topo2.vertLabels()[icl][ip];
      }
      createVtkCell(dataSetSurr, topo2.cellTypes()[icl], pntIds2);
  }

  int numCells = dataSetSurr->GetNumberOfCells();
  int numPoints = dataSetSurr->GetNumberOfPoints();

  int nDim = 3;
  ANNpointArray pntCrd;
  pntCrd = annAllocPts(numPoints, nDim);
  for (int iPnt=0; iPnt < numPoints; iPnt++)
  {
    std::vector<double> getPt = std::vector<double>(3);
    dataSetSurr->GetPoint(iPnt, &getPt[0]);
    pntCrd[iPnt][0] = getPt[0];
    pntCrd[iPnt][1] = getPt[1];
    pntCrd[iPnt][2] = getPt[2];
  }

  ANNkd_tree *kdTree = new ANNkd_tree(pntCrd, numPoints, nDim);

  // Finding duplicate nodes
  double rad = 1e-05;
  std::map<int, int> dupNdeMap;
  ANNpoint qryPnt; // Query point.
  int nNib = 2; // Number of neighbours to return (including query pt. itself).
  ANNidxArray nnIdx = new ANNidx[nNib]; // Nearest neighbour ID array.
  ANNdistArray dists = new ANNdist[nNib]; // Neighbour distance array.
  qryPnt = annAllocPt(nDim); // Initializing query point

  for (int iNde = 0; iNde < impPts; iNde++) {
    std::vector<double> getPt = std::vector<double>(3);
    dataSetSurr->GetPoint(iNde, &getPt[0]);
    qryPnt[0] = getPt[0];
    qryPnt[1] = getPt[1];
    qryPnt[2] = getPt[2];
    kdTree->annkFRSearch(qryPnt, rad, nNib, nnIdx, dists, 0);
    if (dists[1] <= tol){
      if ((nnIdx[1]) >= (_fmeshSurr->nPoints()))
      {
        dupNdeMap[nnIdx[0]] = nnIdx[1]; // Format Surrounding::Pack 
      }
      else{
        dupNdeMap[nnIdx[1]] = nnIdx[0]; // Format Surrounding::Pack
      }
      
    }
  }

  int numDupPnts = dupNdeMap.size();
  std::cout << "Found " << numDupPnts << " duplicate nodes.\n";


  // Getting cells in packs
  std::map<int, int>::iterator it = dupNdeMap.begin();
  std::vector<int> cellID;
  for (int i=0; i<(dupNdeMap.size()); i++)
  {
    int nCells, cellNum;
    vtkIdList *cells = vtkIdList::New();

    // get cells the point belongs to
    dataSetSurr->GetPointCells( it->second, cells );
    nCells = cells->GetNumberOfIds();

    for( cellNum=0; cellNum<nCells; cellNum++ )
    {
      cellID.push_back(cells->GetId( cellNum ));  // get cell id from the list
    }
    it++;
  }

  sort( cellID.begin(), cellID.end() );
  cellID.erase( unique( cellID.begin(), cellID.end() ), cellID.end() );

  int know = cellID.size();

  // Determines which cells are useful
  std::vector<int> newCellIds;
  std::vector<int> debugCellIds;

  for (int i=0; i<know; i++)
  {
    // Getting all cell defining points
    std::vector<int> ptIdzz(8);
    int cntr = 0;
    vtkCell *cell;
    vtkPoints *fp;
    cell = dataSetSurr->GetCell( cellID[i] ); 

    vtkIdList *pts = cell->GetPointIds();

    for (int j=0; j<8; j++)
    {
      ptIdzz[j] = pts->GetId(j);
    }

    // Checking if any cell contains any 4 or more of duplicate nodes.
    std::map<int, int>::iterator it2 = dupNdeMap.begin();
    while(it2 != dupNdeMap.end())
    {
      for (int k=0; k<8; k++)
      {
        if ((it2->second) == ptIdzz[k])
        {
          cntr++;
        }
        else
        {
          // Nothing
        }
      }
      it2++;
    }

    if (cntr >= 4)
    {
      newCellIds.push_back(cellID[i]);
    }
    else
    {
      debugCellIds.push_back(cellID[i]);
    }

  }


  int size2 = newCellIds.size();

  // Getting useful faces from cells
  std::vector<int> globalPtIds;
  std::vector<int> surroundingArray;
  std::vector<int> packArray;
  for (int i=0; i<size2; i++)
  {
    vtkCell *cell;
    cell = dataSetSurr->GetCell(newCellIds[i]);
    vtkIdList *pts = cell->GetPointIds();
    
    for (int j=0; j<6; j++)
    {
      int isFour = 0;
      int* ptFaces = nullptr;
      vtkCell3D* cell3d = static_cast<vtkCell3D*>( cell );
      cell3d->GetFacePoints(j,ptFaces);
      std::vector<int> keysDupMap = std::vector<int>(4);
      std::map<int, int>::iterator it3 = dupNdeMap.begin();
      while (it3 != dupNdeMap.end())
      {
        for (int h=0; h<4; h++)
        {
          if ((it3->second) == (pts->GetId(ptFaces[h])))
          {
            isFour++;
            keysDupMap[h] = it3->first;
          }
          else
          {
            // Nothing
          }
        }

        it3++;
      }

      if (isFour == 4)
      {
        for (int k=0; k<4; k++)
        {
          globalPtIds.push_back(keysDupMap[k]);
        }
      }
    }
  }

  std::map<int, int>::iterator it5;

  // Creating node map with sequence.
  std::unordered_multimap<int, int> cohesiveMap;

  for (int i=0; i<globalPtIds.size(); i++)
  {
    it5 = dupNdeMap.find(globalPtIds[i]);

    if (it5 == dupNdeMap.end())
    {
      // Nothing
    }
    else
    {
      surroundingArray.push_back(it5->first);
      packArray.push_back(it5->second);
    }
  }

  std::cout << "End of previous method" << std::endl;

  // Converting vtkUnstructuredData to vtkPolyData
  vtkSmartPointer<vtkGeometryFilter> geometryFilter =  
  vtkSmartPointer<vtkGeometryFilter>::New();
  geometryFilter->SetInputData(dataSetSurr);
  geometryFilter->Update();

  vtkSmartPointer<vtkPolyData> polydata = 
  vtkSmartPointer<vtkPolyData>::New(); 
  polydata = geometryFilter->GetOutput();

  vtkSmartPointer<vtkPolyDataNormals> polyDataNormal = 
  vtkSmartPointer<vtkPolyDataNormals>::New();

  polyDataNormal->SetInputData(polydata);

  // Setting points normals calculation to ON with other options
  polyDataNormal->ComputePointNormalsOn();
  polyDataNormal->ComputeCellNormalsOff();
  //polyDataNormal->SetNonManifoldTraversal(1);
  //polyDataNormal->SetAutoOrientNormals(0);
  polyDataNormal->SetSplitting(0);
  //polyDataNormal->SetFeatureAngle(0.1);
  polyDataNormal->SetConsistency(1);
  //polyDataNormal->SetFlipNormals(0);
  polyDataNormal->Update();

  // Extracting point normals from polyDataNormal
  polydata = polyDataNormal->GetOutput();

  // Count points
  vtkIdType numPointsNew = polydata->GetNumberOfPoints();
  std::cout << "There are " << numPointsNew << " points." << std::endl;
  
  vtkDataArray* normalsGeneric = polydata->GetPointData()->GetNormals();

  // Number of normals should match number of points
  std::cout << "There are " << normalsGeneric->GetNumberOfTuples()
            << " normals in normalsGeneric" << std::endl;


  // Changing point coordinates for adding artificial thickness;

  // Surrounding cells
  for (int i=0; i<surroundingArray.size(); i++)
  {
    double normalOfPt[3];
    normalsGeneric->GetTuple(surroundingArray[i], normalOfPt);

    std::vector<double> getPt = std::vector<double>(3);
    dataSetSurr->GetPoint(surroundingArray[i], &getPt[0]);

    // TODO : Need to figure out thickness formula considering pack and
    // surrounding mesh cell size and unit normal scale.
    double newCoords[3];
    newCoords[0] = getPt[0] + 0.001*normalOfPt[0];
    newCoords[1] = getPt[1] + 0.001*normalOfPt[1];
    newCoords[2] = getPt[2] + 0.001*normalOfPt[2];

    dataSetSurr->GetPoints()->SetPoint(surroundingArray[i], newCoords);
  }


  // Pack cells
  for (int i=0; i<packArray.size(); i++)
  {
    double normalOfPt[3];
    normalsGeneric->GetTuple(packArray[i], normalOfPt);

    std::vector<double> getPt = std::vector<double>(3);
    dataSetSurr->GetPoint(packArray[i], &getPt[0]);

    double newCoords[3];
    newCoords[0] = (getPt[0] + 0.001*normalOfPt[0]);
    newCoords[1] = (getPt[1] + 0.001*normalOfPt[1]);
    newCoords[2] = (getPt[2] + 0.001*normalOfPt[2]);

    dataSetSurr->GetPoints()->SetPoint(packArray[i], newCoords);
  }

  // Making cells now and plotting them
  int newCells = packArray.size()/4;
  int startNum = dataSetSurr->GetNumberOfCells();

  // Finally, creating VTK cells
  int it7 = 0;
  std::vector<int> pntCohesiveIds;
  int nCelCohesivePnts = 0;
  for (int icl=0; icl<newCells; icl++)
  {
      nCelCohesivePnts = 4;
      pntCohesiveIds.resize((nCelCohesivePnts*2), -1);
      for (int ip=0; ip< nCelCohesivePnts; ip++)
      {
          pntCohesiveIds[ip] = packArray[it7];
          pntCohesiveIds[ip+4] = surroundingArray[it7];
          it7++;
      }
      createVtkCell(dataSetSurr, 12, pntCohesiveIds);
  }

  int endNum = dataSetSurr->GetNumberOfCells();

  // ****************************** //
  // Zero volume cells visualization method for debugging purpose.
  std::vector<int> cohCellIds;
  for (int i=startNum; i<endNum; i++)
  {
    cohCellIds.push_back(i);
  }

  int size3 = cohCellIds.size();
  // Takes cell Ids as input and writes VTK mesh
  vtkSmartPointer<vtkUnstructuredGrid> visDataSet
      = vtkSmartPointer<vtkUnstructuredGrid>::New();

  std::vector<int> ptIdzz2;
  for (int i=0; i<size3; i++){
    
    vtkCell *cell;
    vtkPoints *fp;
    cell = dataSetSurr->GetCell( cohCellIds[i] );

    vtkIdList *pts = cell->GetPointIds();

    for (int j=0; j<8; j++)
    {
      ptIdzz2.push_back(pts->GetId(j));
    }

  }

  int lvar = 0;

  vtkSmartPointer<vtkPoints> pointsViz 
        = vtkSmartPointer<vtkPoints>::New();

  for (int i=0; i<(size3*8); i++)
  {
    std::vector<double> getPt = std::vector<double>(3);
    dataSetSurr->GetPoint(ptIdzz2[i], &getPt[0]);
    pointsViz->InsertNextPoint(getPt[0],getPt[1],getPt[2]);
  }

  visDataSet->SetPoints(pointsViz);

  

  for (int i=0; i<size3; i++)
  {
    std::vector<int> ptnewIdz;
    for (int j=0; j<8; j++)
    {
      ptnewIdz.push_back(lvar);
      lvar++;
    }
    createVtkCell(visDataSet, 12, ptnewIdz);
  }

  vtkMesh* vm = new vtkMesh(visDataSet,"ArtificialElements.vtu");
  vm->report();
  vm->write();

  if (vm)
    delete vm;

  vtkMesh* vm2 = new vtkMesh(dataSetSurr,outName);
  vm2->report();
  vm2->write();

  if (vm2)
    delete vm2;

  bool tetra = false;
  std::string ofnameTet = "Tetra" + outName;

  if (tetra == true)
  {
    // create meshBase object
    std::shared_ptr<meshBase> myMesh = meshBase::CreateShared(outName);

    // Converts hex mesh to tet mesh and writes in VTU file.
    myMesh->convertHexToTetVTK(dataSetSurr);
    myMesh->report();
    myMesh->write(ofnameTet);
  }
  

}

// Creates a VTK cell in databse using provided points and cell type.
void MeshManipulationFoam::createVtkCell(
        vtkSmartPointer<vtkUnstructuredGrid> dataSet,
        const int cellType, std::vector<int>& vrtIds)
{
    vtkSmartPointer<vtkIdList> vtkCellIds = vtkSmartPointer<vtkIdList>::New();
    vtkCellIds->SetNumberOfIds(vrtIds.size());
    for (auto pit = vrtIds.begin(); pit!= vrtIds.end(); pit++)
        vtkCellIds->SetId(pit-vrtIds.begin(), *pit);
    dataSet->InsertNextCell(cellType,vtkCellIds);
}


// Generates Periodic Mesh
void MeshManipulationFoam::periodicMeshMapper(std::string patch1, 
  std::string patch2)
{
  // Initializing FOAM
  int argc = 1;
  char** argv = new char*[2];
  argv[0] = new char[100];
  strcpy(argv[0], "NONE");
  Foam::argList args(argc, argv);
  Foam::argList::noParallel();

  Time runTime
  (
      Time::controlDictName,
      ".",
      "."
  );

  vtkSmartPointer<vtkPoints> points 
      = vtkSmartPointer<vtkPoints>::New(); 

  Foam::word regionName = Foam::polyMesh::defaultRegion;

  auto _fmeshPacks = new Foam::polyMesh
  (
    Foam::IOobject
    (
      regionName,
      runTime.timeName(),
      runTime,
      Foam::IOobject::MUST_READ
    )
  );

  // Creating patch objects for periodic boundary patches
  int patchIDIn = _fmeshPacks->boundaryMesh().findPatchID(patch1);
  int patchIDOut = _fmeshPacks->boundaryMesh().findPatchID(patch2);

  const Foam::polyPatch& InPolyPatch = _fmeshPacks->boundaryMesh()[patchIDIn];
  const Foam::polyPatch& OutPolyPatch = _fmeshPacks->boundaryMesh()[patchIDOut];

  // Point IDs from periodic patches.
  Foam::labelList inPts(InPolyPatch.meshPoints());
  Foam::labelList outPts(OutPolyPatch.meshPoints());

  if (inPts.size() != outPts.size())
  {
    std::cerr << "Selected boundaries are not periodic!" 
              << "\nExiting!" << std::endl;
    throw;
  }

  Foam::pointField packsPF = _fmeshPacks->points();

  // Creating map of IDs on both patches
  std::map<int, int> periodicMap;

  for (int i=0; i<inPts.size(); i++)
  {
    periodicMap[inPts[i]] = outPts[i];
  }

  ofstream myfile;
  myfile.open ("PeriodicMap.csv");
  for (int i=0; i<inPts.size(); i++)
  {
    myfile << inPts[i] << "," << outPts[i] << std::endl;
  }
  myfile.close();


  // Method test using VTK unstructured writer
  // Very useful for extracting face node orders. Might use this in
  // cohesive elements. Also, Foam::facePointPatch::pointNormals()
  // can be useful for artificial thickness elements
  vtkSmartPointer<vtkUnstructuredGrid> dataSetTest
      = vtkSmartPointer<vtkUnstructuredGrid>::New();
  vtkSmartPointer<vtkPoints> pointsMesh
      = vtkSmartPointer<vtkPoints>::New();

  for (int i=0; i<packsPF.size(); i++)
  {
    pointsMesh->InsertNextPoint(packsPF[i].x(),packsPF[i].y(),packsPF[i].z());
  }

  dataSetTest->SetPoints(pointsMesh);

  forAll (_fmeshPacks->boundaryMesh()[patchIDIn],facei) 
  {
    const Foam::label& faceID = _fmeshPacks->boundaryMesh()[patchIDIn].start() + facei;
    std::vector<int> pntIds;
    pntIds.resize(4, -1);
    int g = 0;
    forAll (_fmeshPacks->faces()[faceID], nodei)
    {
      const Foam::label& nodeID = _fmeshPacks->faces()[faceID][nodei];
      pntIds[g] = nodeID;
      g++;
    }
    createVtkCell(dataSetTest,9,pntIds);
  }

  forAll (_fmeshPacks->boundaryMesh()[patchIDOut],facei) 
  {
    const Foam::label& faceID = _fmeshPacks->boundaryMesh()[patchIDOut].start() + facei;
    std::vector<int> pntIds;
    pntIds.resize(4, -1);
    int g = 0;
    forAll (_fmeshPacks->faces()[faceID], nodei)
    {
      const Foam::label& nodeID = _fmeshPacks->faces()[faceID][nodei];
      pntIds[g] = nodeID;
      g++;
    }
    createVtkCell(dataSetTest,9,pntIds);
  }

  vtkMesh* vm = new vtkMesh(dataSetTest,"PeriodicMesh.vtu");
  vm->report();
  vm->write();

  if (vm)
    delete vm;

}


// All the private functions are defined here

//SurfaceLambdaMuSmooth
Foam::tmp<Foam::pointField> MeshManipulationFoam::avg
(
  const Foam::meshedSurface& s,
  const Foam::PackedBoolList& fixedPoints
)
{
  using namespace Foam;
  const labelListList& pointEdges = s.pointEdges();

    tmp<pointField> tavg(new pointField(s.nPoints(), Zero));
    pointField& avg = tavg.ref();

    forAll(pointEdges, vertI)
    {
        vector& avgPos = avg[vertI];

        if (fixedPoints[vertI])
        {
            avgPos = s.localPoints()[vertI];
        }
        else
        {
            const labelList& pEdges = pointEdges[vertI];

            forAll(pEdges, myEdgeI)
            {
                const edge& e = s.edges()[pEdges[myEdgeI]];

                label otherVertI = e.otherVertex(vertI);

                avgPos += s.localPoints()[otherVertI];
            }

            avgPos /= pEdges.size();
        }
    }

    return tavg;
}

//SurfaceLambdaMuSmooth
void MeshManipulationFoam::getFixedPoints
(
    const Foam::edgeMesh& feMesh,
    const Foam::pointField& points,
    Foam::PackedBoolList& fixedPoints
)
{
  using namespace Foam;
    scalarList matchDistance(feMesh.points().size(), 1e-1);
    labelList from0To1;

    bool matchedAll = matchPoints
    (
        feMesh.points(),
        points,
        matchDistance,
        false,
        from0To1
    );

    if (!matchedAll)
    {
        WarningInFunction
            << "Did not match all feature points to points on the surface"
            << endl;
    }

    forAll(from0To1, fpI)
    {
        if (from0To1[fpI] != -1)
        {
            fixedPoints[from0To1[fpI]] = true;
        }
    }
}


//splitMeshByRegions
void MeshManipulationFoam::renamePatches
(
    Foam::fvMesh& mesh,
    const Foam::word& prefix,
    const Foam::labelList& patchesToRename
)
{
  using namespace Foam;
    polyBoundaryMesh& polyPatches =
        const_cast<polyBoundaryMesh&>(mesh.boundaryMesh());
    forAll(patchesToRename, i)
    {
        label patchi = patchesToRename[i];
        polyPatch& pp = polyPatches[patchi];

        if (isA<coupledPolyPatch>(pp))
        {
            WarningInFunction
                << "Encountered coupled patch " << pp.name()
                << ". Will only rename the patch itself,"
                << " not any referred patches."
                << " This might have to be done by hand."
                << endl;
        }

        pp.name() = prefix + '_' + pp.name();
    }
    // Recalculate any demand driven data (e.g. group to name lookup)
    polyPatches.updateMesh();
}

//splitMeshByRegions
template<class GeoField>
void MeshManipulationFoam::subsetVolFields
(
    const Foam::fvMesh& mesh,
    const Foam::fvMesh& subMesh,
    const Foam::labelList& cellMap,
    const Foam::labelList& faceMap,
    const Foam::labelHashSet& addedPatches
)
{
  using namespace Foam;
    const labelList patchMap(identity(mesh.boundaryMesh().size()));

    HashTable<const GeoField*> fields
    (
        mesh.objectRegistry::lookupClass<GeoField>()
    );
    forAllConstIter(typename HashTable<const GeoField*>, fields, iter)
    {
        const GeoField& fld = *iter();

        Info<< "Mapping field " << fld.name() << endl;

        tmp<GeoField> tSubFld
        (
            fvMeshSubset::interpolate
            (
                fld,
                subMesh,
                patchMap,
                cellMap,
                faceMap
            )
        );

        // Hack: set value to 0 for introduced patches (since don't
        //       get initialised.
        forAll(tSubFld().boundaryField(), patchi)
        {
            if (addedPatches.found(patchi))
            {
                tSubFld.ref().boundaryFieldRef()[patchi] ==
                    typename GeoField::value_type(Zero);
            }
        }

        // Store on subMesh
        GeoField* subFld = tSubFld.ptr();
        subFld->rename(fld.name());
        subFld->writeOpt() = IOobject::AUTO_WRITE;
        subFld->store();
    }
}

//splitMeshByRegions
template<class GeoField>
void MeshManipulationFoam::subsetSurfaceFields
(
    const Foam::fvMesh& mesh,
    const Foam::fvMesh& subMesh,
    const Foam::labelList& cellMap,
    const Foam::labelList& faceMap,
    const Foam::labelHashSet& addedPatches
)
{
  using namespace Foam;
    const labelList patchMap(identity(mesh.boundaryMesh().size()));

    HashTable<const GeoField*> fields
    (
        mesh.objectRegistry::lookupClass<GeoField>()
    );
    forAllConstIter(typename HashTable<const GeoField*>, fields, iter)
    {
        const GeoField& fld = *iter();

        Info<< "Mapping field " << fld.name() << endl;

        tmp<GeoField> tSubFld
        (
            fvMeshSubset::interpolate
            (
                fld,
                subMesh,
                patchMap,
                cellMap,
                faceMap
            )
        );

        // Hack: set value to 0 for introduced patches (since don't
        //       get initialised.
        forAll(tSubFld().boundaryField(), patchi)
        {
            if (addedPatches.found(patchi))
            {
                tSubFld.ref().boundaryFieldRef()[patchi] ==
                    typename GeoField::value_type(Zero);
            }
        }

        // Store on subMesh
        GeoField* subFld = tSubFld.ptr();
        subFld->rename(fld.name());
        subFld->writeOpt() = IOobject::AUTO_WRITE;
        subFld->store();
    }
}


//splitMeshByRegions
Foam::labelList MeshManipulationFoam::getNonRegionCells
(
  const Foam::labelList& cellRegion,
  const Foam::label regionI
)
{
  using namespace Foam;
    DynamicList<label> nonRegionCells(cellRegion.size());
    forAll(cellRegion, celli)
    {
        if (cellRegion[celli] != regionI)
        {
            nonRegionCells.append(celli);
        }
    }
    return nonRegionCells.shrink();
}

//splitMeshByRegions
void MeshManipulationFoam::addToInterface
(
    const Foam::polyMesh& mesh,
    const Foam::label zoneID,
    const Foam::label ownRegion,
    const Foam::label neiRegion,
    Foam::EdgeMap<Foam::Map<Foam::label>>& regionsToSize
)
{
  using namespace Foam;
    edge interface
    (
        min(ownRegion, neiRegion),
        max(ownRegion, neiRegion)
    );

    EdgeMap<Map<label>>::iterator iter = regionsToSize.find
    (
        interface
    );

    if (iter != regionsToSize.end())
    {
        // Check if zone present
        Map<label>::iterator zoneFnd = iter().find(zoneID);
        if (zoneFnd != iter().end())
        {
            // Found zone. Increment count.
            zoneFnd()++;
        }
        else
        {
            // New or no zone. Insert with count 1.
            iter().insert(zoneID, 1);
        }
    }
    else
    {
        // Create new interface of size 1.
        Map<label> zoneToSize;
        zoneToSize.insert(zoneID, 1);
        regionsToSize.insert(interface, zoneToSize);
    }
}

//splitMeshByRegions
void MeshManipulationFoam::getInterfaceSizes
(
    const Foam::polyMesh& mesh,
    const bool useFaceZones,
    const Foam::labelList& cellRegion,
    const Foam::wordList& regionNames,

    Foam::edgeList& interfaces,
    Foam::List<Foam::Pair<Foam::word>>& interfaceNames,
    Foam::labelList& interfaceSizes,
    Foam::labelList& faceToInterface
)
{
  using namespace Foam;
    // From region-region to faceZone (or -1) to number of faces.

    EdgeMap<Map<label>> regionsToSize;


    // Internal faces
    // ~~~~~~~~~~~~~~

    forAll(mesh.faceNeighbour(), facei)
    {
        label ownRegion = cellRegion[mesh.faceOwner()[facei]];
        label neiRegion = cellRegion[mesh.faceNeighbour()[facei]];

        if (ownRegion != neiRegion)
        {
            addToInterface
            (
                mesh,
                (useFaceZones ? mesh.faceZones().whichZone(facei) : -1),
                ownRegion,
                neiRegion,
                regionsToSize
            );
        }
    }

    // Boundary faces
    // ~~~~~~~~~~~~~~

    // Neighbour cellRegion.
    labelList coupledRegion(mesh.nFaces()-mesh.nInternalFaces());

    forAll(coupledRegion, i)
    {
        label celli = mesh.faceOwner()[i+mesh.nInternalFaces()];
        coupledRegion[i] = cellRegion[celli];
    }
    syncTools::swapBoundaryFaceList(mesh, coupledRegion);

    forAll(coupledRegion, i)
    {
        label facei = i+mesh.nInternalFaces();
        label ownRegion = cellRegion[mesh.faceOwner()[facei]];
        label neiRegion = coupledRegion[i];

        if (ownRegion != neiRegion)
        {
            addToInterface
            (
                mesh,
                (useFaceZones ? mesh.faceZones().whichZone(facei) : -1),
                ownRegion,
                neiRegion,
                regionsToSize
            );
        }
    }


    if (Pstream::parRun())
    {
        if (Pstream::master())
        {
            // Receive and add to my sizes
            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave();
                slave++
            )
            {
                IPstream fromSlave(Pstream::commsTypes::blocking, slave);

                EdgeMap<Map<label>> slaveSizes(fromSlave);

                forAllConstIter(EdgeMap<Map<label>>, slaveSizes, slaveIter)
                {
                    EdgeMap<Map<label>>::iterator masterIter =
                        regionsToSize.find(slaveIter.key());

                    if (masterIter != regionsToSize.end())
                    {
                        // Same inter-region
                        const Map<label>& slaveInfo = slaveIter();
                        Map<label>& masterInfo = masterIter();

                        forAllConstIter(Map<label>, slaveInfo, iter)
                        {
                            label zoneID = iter.key();
                            label slaveSize = iter();

                            Map<label>::iterator zoneFnd = masterInfo.find
                            (
                                zoneID
                            );
                            if (zoneFnd != masterInfo.end())
                            {
                                zoneFnd() += slaveSize;
                            }
                            else
                            {
                                masterInfo.insert(zoneID, slaveSize);
                            }
                        }
                    }
                    else
                    {
                        regionsToSize.insert(slaveIter.key(), slaveIter());
                    }
                }
            }
        }
        else
        {
            // Send to master
            {
                OPstream toMaster
                (
                    Pstream::commsTypes::blocking,
                    Pstream::masterNo()
                );
                toMaster << regionsToSize;
            }
        }
    }

    // Rework

    Pstream::scatter(regionsToSize);



    // Now we have the global sizes of all inter-regions.
    // Invert this on master and distribute.
    label nInterfaces = 0;
    forAllConstIter(EdgeMap<Map<label>>, regionsToSize, iter)
    {
        const Map<label>& info = iter();
        nInterfaces += info.size();
    }

    interfaces.setSize(nInterfaces);
    interfaceNames.setSize(nInterfaces);
    interfaceSizes.setSize(nInterfaces);
    EdgeMap<Map<label>> regionsToInterface(nInterfaces);

    nInterfaces = 0;
    forAllConstIter(EdgeMap<Map<label>>, regionsToSize, iter)
    {
        const edge& e = iter.key();
        const word& name0 = regionNames[e[0]];
        const word& name1 = regionNames[e[1]];

        const Map<label>& info = iter();
        forAllConstIter(Map<label>, info, infoIter)
        {
            interfaces[nInterfaces] = iter.key();
            label zoneID = infoIter.key();
            if (zoneID == -1)
            {
                interfaceNames[nInterfaces] = Pair<word>
                (
                    name0 + "_to_" + name1,
                    name1 + "_to_" + name0
                );
            }
            else
            {
                const word& zoneName = mesh.faceZones()[zoneID].name();
                interfaceNames[nInterfaces] = Pair<word>
                (
                    zoneName + "_" + name0 + "_to_" + name1,
                    zoneName + "_" + name1 + "_to_" + name0
                );
            }
            interfaceSizes[nInterfaces] = infoIter();

            if (regionsToInterface.found(e))
            {
                regionsToInterface[e].insert(zoneID, nInterfaces);
            }
            else
            {
                Map<label> zoneAndInterface;
                zoneAndInterface.insert(zoneID, nInterfaces);
                regionsToInterface.insert(e, zoneAndInterface);
            }
            nInterfaces++;
        }
    }


    // Now all processor have consistent interface information

    Pstream::scatter(interfaces);
    Pstream::scatter(interfaceNames);
    Pstream::scatter(interfaceSizes);
    Pstream::scatter(regionsToInterface);

    // Mark all inter-region faces.
    faceToInterface.setSize(mesh.nFaces(), -1);

    forAll(mesh.faceNeighbour(), facei)
    {
        label ownRegion = cellRegion[mesh.faceOwner()[facei]];
        label neiRegion = cellRegion[mesh.faceNeighbour()[facei]];

        if (ownRegion != neiRegion)
        {
            label zoneID = -1;
            if (useFaceZones)
            {
                zoneID = mesh.faceZones().whichZone(facei);
            }

            edge interface
            (
                min(ownRegion, neiRegion),
                max(ownRegion, neiRegion)
            );

            faceToInterface[facei] = regionsToInterface[interface][zoneID];
        }
    }
    forAll(coupledRegion, i)
    {
        label facei = i+mesh.nInternalFaces();
        label ownRegion = cellRegion[mesh.faceOwner()[facei]];
        label neiRegion = coupledRegion[i];

        if (ownRegion != neiRegion)
        {
            label zoneID = -1;
            if (useFaceZones)
            {
                zoneID = mesh.faceZones().whichZone(facei);
            }

            edge interface
            (
                min(ownRegion, neiRegion),
                max(ownRegion, neiRegion)
            );

            faceToInterface[facei] = regionsToInterface[interface][zoneID];
        }
    }
}


//splitMeshByRegions
Foam::autoPtr<Foam::mapPolyMesh> MeshManipulationFoam::createRegionMesh
(
    const Foam::fvMesh& mesh,
    // Region info
    const Foam::labelList& cellRegion,
    const Foam::label regionI,
    const Foam::word& regionName,
    // Interface info
    const Foam::labelList& interfacePatches,
    const Foam::labelList& faceToInterface,

    Foam::autoPtr<Foam::fvMesh>& newMesh
)
{
  using namespace Foam;
    // Create dummy system/fv*
    {
        IOobject io
        (
            "fvSchemes",
            mesh.time().system(),
            regionName,
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        );

        Info<< "Testing:" << io.objectPath() << endl;

#ifdef HAVE_OF5

        if (!io.typeHeaderOk<IOdictionary>(true))
        // if (!exists(io.objectPath()))
        {
            Info<< "Writing dummy " << regionName/io.name() << endl;
            dictionary dummyDict;
            dictionary divDict;
            dummyDict.add("divSchemes", divDict);
            dictionary gradDict;
            dummyDict.add("gradSchemes", gradDict);
            dictionary laplDict;
            dummyDict.add("laplacianSchemes", laplDict);

            IOdictionary(io, dummyDict).regIOobject::write();
        }
    }
    {
        IOobject io
        (
            "fvSolution",
            mesh.time().system(),
            regionName,
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        );

        if (!io.typeHeaderOk<IOdictionary>(true))
        //if (!exists(io.objectPath()))
        {
            Info<< "Writing dummy " << regionName/io.name() << endl;
            dictionary dummyDict;
            IOdictionary(io, dummyDict).regIOobject::write();
        }
    }


#endif

#ifdef HAVE_OF4

    
    if (!io.headerOk())
        {
            Info<< "Writing dummy " << regionName/io.name() << endl;
            dictionary dummyDict;
            dictionary divDict;
            dummyDict.add("divSchemes", divDict);
            dictionary gradDict;
            dummyDict.add("gradSchemes", gradDict);
            dictionary laplDict;
            dummyDict.add("laplacianSchemes", laplDict);

            IOdictionary(io, dummyDict).regIOobject::write();
        }
    }
    {
        IOobject io
        (
            "fvSolution",
            mesh.time().system(),
            regionName,
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        );

        
        if (!io.headerOk())
        {
            Info<< "Writing dummy " << regionName/io.name() << endl;
            dictionary dummyDict;
            IOdictionary(io, dummyDict).regIOobject::write();
        }
    }

#endif

#ifdef HAVE_OF6

        if (!io.typeHeaderOk<IOdictionary>(true))
        // if (!exists(io.objectPath()))
        {
            Info<< "Writing dummy " << regionName/io.name() << endl;
            dictionary dummyDict;
            dictionary divDict;
            dummyDict.add("divSchemes", divDict);
            dictionary gradDict;
            dummyDict.add("gradSchemes", gradDict);
            dictionary laplDict;
            dummyDict.add("laplacianSchemes", laplDict);

            IOdictionary(io, dummyDict).regIOobject::write();
        }
    }
    {
        IOobject io
        (
            "fvSolution",
            mesh.time().system(),
            regionName,
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        );

        if (!io.typeHeaderOk<IOdictionary>(true))
        //if (!exists(io.objectPath()))
        {
            Info<< "Writing dummy " << regionName/io.name() << endl;
            dictionary dummyDict;
            IOdictionary(io, dummyDict).regIOobject::write();
        }
    }


#endif

#ifdef HAVE_OF7

        if (!io.typeHeaderOk<IOdictionary>(true))
        // if (!exists(io.objectPath()))
        {
            Info<< "Writing dummy " << regionName/io.name() << endl;
            dictionary dummyDict;
            dictionary divDict;
            dummyDict.add("divSchemes", divDict);
            dictionary gradDict;
            dummyDict.add("gradSchemes", gradDict);
            dictionary laplDict;
            dummyDict.add("laplacianSchemes", laplDict);

            IOdictionary(io, dummyDict).regIOobject::write();
        }
    }
    {
        IOobject io
        (
            "fvSolution",
            mesh.time().system(),
            regionName,
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        );

        if (!io.typeHeaderOk<IOdictionary>(true))
        //if (!exists(io.objectPath()))
        {
            Info<< "Writing dummy " << regionName/io.name() << endl;
            dictionary dummyDict;
            IOdictionary(io, dummyDict).regIOobject::write();
        }
    }


#endif

    // Neighbour cellRegion.
    labelList coupledRegion(mesh.nFaces()-mesh.nInternalFaces());

    forAll(coupledRegion, i)
    {
        label celli = mesh.faceOwner()[i+mesh.nInternalFaces()];
        coupledRegion[i] = cellRegion[celli];
    }
    syncTools::swapBoundaryFaceList(mesh, coupledRegion);


    // Topology change container. Start off from existing mesh.
    polyTopoChange meshMod(mesh);

    // Cell remover engine
    removeCells cellRemover(mesh);

    // Select all but region cells
    labelList cellsToRemove(getNonRegionCells(cellRegion, regionI));

    // Find out which faces will get exposed. Note that this
    // gets faces in mesh face order. So both regions will get same
    // face in same order (important!)
    labelList exposedFaces = cellRemover.getExposedFaces(cellsToRemove);

    labelList exposedPatchIDs(exposedFaces.size());
    forAll(exposedFaces, i)
    {
        label facei = exposedFaces[i];
        label interfacei = faceToInterface[facei];

        label ownRegion = cellRegion[mesh.faceOwner()[facei]];
        label neiRegion = -1;

        if (mesh.isInternalFace(facei))
        {
            neiRegion = cellRegion[mesh.faceNeighbour()[facei]];
        }
        else
        {
            neiRegion = coupledRegion[facei-mesh.nInternalFaces()];
        }


        // Check which side is being kept - determines which of the two
        // patches will be used.

        label otherRegion = -1;

        if (ownRegion == regionI && neiRegion != regionI)
        {
            otherRegion = neiRegion;
        }
        else if (ownRegion != regionI && neiRegion == regionI)
        {
            otherRegion = ownRegion;
        }
        else
        {
            FatalErrorInFunction
                << "Exposed face:" << facei
                << " fc:" << mesh.faceCentres()[facei]
                << " has owner region " << ownRegion
                << " and neighbour region " << neiRegion
                << " when handling region:" << regionI
                << exit(FatalError);
        }

        // Find the patch.
        if (regionI < otherRegion)
        {
            exposedPatchIDs[i] = interfacePatches[interfacei];
        }
        else
        {
            exposedPatchIDs[i] = interfacePatches[interfacei]+1;
        }
    }

    // Remove faces
    cellRemover.setRefinement
    (
        cellsToRemove,
        exposedFaces,
        exposedPatchIDs,
        meshMod
    );

    autoPtr<mapPolyMesh> map = meshMod.makeMesh
    (
        newMesh,
        IOobject
        (
            regionName,
            mesh.time().timeName(),
            mesh.time(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    return map;
}

//splitMeshRegions
//void MeshManipulationFoam::createAndWriteRegion
void MeshManipulationFoam::createAndWriteRegion
(
    const Foam::fvMesh& mesh,
    const Foam::labelList& cellRegion,
    const Foam::wordList& regionNames,
    const bool prefixRegion,
    const Foam::labelList& faceToInterface,
    const Foam::labelList& interfacePatches,
    const Foam::label regionI,
    const Foam::word& newMeshInstance
)
{
  using namespace Foam;
    Info<< "Creating mesh for region " << regionI
        << ' ' << regionNames[regionI] << endl;

    autoPtr<fvMesh> newMesh;
    autoPtr<mapPolyMesh> map = createRegionMesh
    (
        mesh,
        cellRegion,
        regionI,
        regionNames[regionI],
        interfacePatches,
        faceToInterface,
        newMesh
    );


    // Make map of all added patches
    labelHashSet addedPatches(2*interfacePatches.size());
    forAll(interfacePatches, interfacei)
    {
        addedPatches.insert(interfacePatches[interfacei]);
        addedPatches.insert(interfacePatches[interfacei]+1);
    }


    Info<< "Mapping fields" << endl;

    // Map existing fields
    newMesh().updateMesh(map());

    // Add subsetted fields
    subsetVolFields<volScalarField>
    (
        mesh,
        newMesh(),
        map().cellMap(),
        map().faceMap(),
        addedPatches
    );
    subsetVolFields<volVectorField>
    (
        mesh,
        newMesh(),
        map().cellMap(),
        map().faceMap(),
        addedPatches
    );
    subsetVolFields<volSphericalTensorField>
    (
        mesh,
        newMesh(),
        map().cellMap(),
        map().faceMap(),
        addedPatches
    );
    subsetVolFields<volSymmTensorField>
    (
        mesh,
        newMesh(),
        map().cellMap(),
        map().faceMap(),
        addedPatches
    );
    subsetVolFields<volTensorField>
    (
        mesh,
        newMesh(),
        map().cellMap(),
        map().faceMap(),
        addedPatches
    );

    subsetSurfaceFields<surfaceScalarField>
    (
        mesh,
        newMesh(),
        map().cellMap(),
        map().faceMap(),
        addedPatches
    );
    subsetSurfaceFields<surfaceVectorField>
    (
        mesh,
        newMesh(),
        map().cellMap(),
        map().faceMap(),
        addedPatches
    );
    subsetSurfaceFields<surfaceSphericalTensorField>
    (
        mesh,
        newMesh(),
        map().cellMap(),
        map().faceMap(),
        addedPatches
    );
    subsetSurfaceFields<surfaceSymmTensorField>
    (
        mesh,
        newMesh(),
        map().cellMap(),
        map().faceMap(),
        addedPatches
    );
    subsetSurfaceFields<surfaceTensorField>
    (
        mesh,
        newMesh(),
        map().cellMap(),
        map().faceMap(),
        addedPatches
    );


    const polyBoundaryMesh& newPatches = newMesh().boundaryMesh();
    newPatches.checkParallelSync(true);

    // Delete empty patches
    // ~~~~~~~~~~~~~~~~~~~~

    // Create reordering list to move patches-to-be-deleted to end
    labelList oldToNew(newPatches.size(), -1);
    DynamicList<label> sharedPatches(newPatches.size());
    label newI = 0;

    Info<< "Deleting empty patches" << endl;

    // Assumes all non-proc boundaries are on all processors!
    forAll(newPatches, patchi)
    {
        const polyPatch& pp = newPatches[patchi];

        if (!isA<processorPolyPatch>(pp))
        {
            if (returnReduce(pp.size(), sumOp<label>()) > 0)
            {
                oldToNew[patchi] = newI;
                if (!addedPatches.found(patchi))
                {
                    sharedPatches.append(newI);
                }
                newI++;
            }
        }
    }

    // Same for processor patches (but need no reduction)
    forAll(newPatches, patchi)
    {
        const polyPatch& pp = newPatches[patchi];

        if (isA<processorPolyPatch>(pp) && pp.size())
        {
            oldToNew[patchi] = newI++;
        }
    }

    const label nNewPatches = newI;

    // Move all deleteable patches to the end
    forAll(oldToNew, patchi)
    {
        if (oldToNew[patchi] == -1)
        {
            oldToNew[patchi] = newI++;
        }
    }

    //reorderPatches(newMesh(), oldToNew, nNewPatches);
    fvMeshTools::reorderPatches(newMesh(), oldToNew, nNewPatches, true);

    // Rename shared patches with region name
    if (prefixRegion)
    {
        Info<< "Prefixing patches with region name" << endl;

        renamePatches(newMesh(), regionNames[regionI], sharedPatches);
    }


    Info<< "Writing new mesh" << endl;

    newMesh().setInstance(newMeshInstance);
    newMesh().write();

    // Write addressing files like decomposePar
    Info<< "Writing addressing to base mesh" << endl;

    labelIOList pointProcAddressing
    (
        IOobject
        (
            "pointRegionAddressing",
            newMesh().facesInstance(),
            newMesh().meshSubDir,
            newMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        map().pointMap()
    );
    Info<< "Writing map " << pointProcAddressing.name()
        << " from region" << regionI
        << " points back to base mesh." << endl;
    pointProcAddressing.write();

    labelIOList faceProcAddressing
    (
        IOobject
        (
            "faceRegionAddressing",
            newMesh().facesInstance(),
            newMesh().meshSubDir,
            newMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        newMesh().nFaces()
    );
    forAll(faceProcAddressing, facei)
    {
        // face + turning index. (see decomposePar)
        // Is the face pointing in the same direction?
        label oldFacei = map().faceMap()[facei];

        if
        (
            map().cellMap()[newMesh().faceOwner()[facei]]
         == mesh.faceOwner()[oldFacei]
        )
        {
            faceProcAddressing[facei] = oldFacei+1;
        }
        else
        {
            faceProcAddressing[facei] = -(oldFacei+1);
        }
    }
    Info<< "Writing map " << faceProcAddressing.name()
        << " from region" << regionI
        << " faces back to base mesh." << endl;
    faceProcAddressing.write();

    labelIOList cellProcAddressing
    (
        IOobject
        (
            "cellRegionAddressing",
            newMesh().facesInstance(),
            newMesh().meshSubDir,
            newMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        map().cellMap()
    );
    Info<< "Writing map " <<cellProcAddressing.name()
        << " from region" << regionI
        << " cells back to base mesh." << endl;
    cellProcAddressing.write();

    labelIOList boundaryProcAddressing
    (
        IOobject
        (
            "boundaryRegionAddressing",
            newMesh().facesInstance(),
            newMesh().meshSubDir,
            newMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        labelList(nNewPatches, -1)
    );
    forAll(oldToNew, i)
    {
        if (!addedPatches.found(i))
        {
            label newI = oldToNew[i];
            if (newI >= 0 && newI < nNewPatches)
            {
                boundaryProcAddressing[oldToNew[i]] = i;
            }
        }
    }
    Info<< "Writing map " << boundaryProcAddressing.name()
        << " from region" << regionI
        << " boundary back to base mesh." << endl;
    boundaryProcAddressing.write();
}


//splitMeshRegions
Foam::labelList MeshManipulationFoam::addRegionPatches
(
    Foam::fvMesh& mesh,
    const Foam::wordList& regionNames,
    const Foam::edgeList& interfaces,
    const Foam::List<Foam::Pair<Foam::word>>& interfaceNames
)
{
  using namespace Foam;
    Info<< nl << "Adding patches" << nl << endl;

    labelList interfacePatches(interfaces.size());

    forAll(interfaces, interI)
    {
        const edge& e = interfaces[interI];
        const Pair<word>& names = interfaceNames[interI];

        //Info<< "For interface " << interI
        //    << " between regions " << e
        //    << " trying to add patches " << names << endl;


        mappedWallPolyPatch patch1
        (
            names[0],
            0,                  // overridden
            0,                  // overridden
            0,                  // overridden
            regionNames[e[1]],  // sampleRegion
            mappedPatchBase::NEARESTPATCHFACE,
            names[1],           // samplePatch
            point::zero,        // offset
            mesh.boundaryMesh()
        );

        interfacePatches[interI] = fvMeshTools::addPatch
        (
            mesh,
            patch1,
            dictionary(),   //optional per field value
            calculatedFvPatchField<scalar>::typeName,
            true            //validBoundary
        );

        mappedWallPolyPatch patch2
        (
            names[1],
            0,
            0,
            0,
            regionNames[e[0]],  // sampleRegion
            mappedPatchBase::NEARESTPATCHFACE,
            names[0],
            point::zero,        // offset
            mesh.boundaryMesh()
        );
        fvMeshTools::addPatch
        (
            mesh,
            patch2,
            dictionary(),   //optional per field value
            calculatedFvPatchField<scalar>::typeName,
            true            //validBoundary
        );

        Info<< "For interface between region " << regionNames[e[0]]
            << " and " << regionNames[e[1]] << " added patches" << endl
            << "    " << interfacePatches[interI]
            << "\t" << mesh.boundaryMesh()[interfacePatches[interI]].name()
            << endl
            << "    " << interfacePatches[interI]+1
            << "\t" << mesh.boundaryMesh()[interfacePatches[interI]+1].name()
            << endl;
    }
    return interfacePatches;
}


//splitMeshRegions
Foam::label MeshManipulationFoam::findCorrespondingRegion
(
    const Foam::labelList& existingZoneID,    // per cell the (unique) zoneID
    const Foam::labelList& cellRegion,
    const Foam::label nCellRegions,
    const Foam::label zoneI,
    const Foam::label minOverlapSize
)
{
  using namespace Foam;
    // Per region the number of cells in zoneI
    labelList cellsInZone(nCellRegions, 0);

    forAll(cellRegion, celli)
    {
        if (existingZoneID[celli] == zoneI)
        {
            cellsInZone[cellRegion[celli]]++;
        }
    }

    Pstream::listCombineGather(cellsInZone, plusEqOp<label>());
    Pstream::listCombineScatter(cellsInZone);

    // Pick region with largest overlap of zoneI
    label regionI = findMax(cellsInZone);


    if (cellsInZone[regionI] < minOverlapSize)
    {
        // Region covers too little of zone. Not good enough.
        regionI = -1;
    }
    else
    {
        // Check that region contains no cells that aren't in cellZone.
        forAll(cellRegion, celli)
        {
            if (cellRegion[celli] == regionI && existingZoneID[celli] != zoneI)
            {
                // celli in regionI but not in zoneI
                regionI = -1;
                break;
            }
        }
        // If one in error, all should be in error. Note that branch gets taken
        // on all procs.
        reduce(regionI, minOp<label>());
    }

    return regionI;
}


//splitMeshRegions
void MeshManipulationFoam::getZoneID
(
    const Foam::polyMesh& mesh,
    const Foam::cellZoneMesh& cellZones,
    Foam::labelList& zoneID,
    Foam::labelList& neiZoneID
)
{
  using namespace Foam;
    // Existing zoneID
    zoneID.setSize(mesh.nCells());
    zoneID = -1;

    forAll(cellZones, zoneI)
    {
        const cellZone& cz = cellZones[zoneI];

        forAll(cz, i)
        {
            label celli = cz[i];
            if (zoneID[celli] == -1)
            {
                zoneID[celli] = zoneI;
            }
            else
            {
                FatalErrorInFunction
                    << "Cell " << celli << " with cell centre "
                    << mesh.cellCentres()[celli]
                    << " is multiple zones. This is not allowed." << endl
                    << "It is in zone " << cellZones[zoneID[celli]].name()
                    << " and in zone " << cellZones[zoneI].name()
                    << exit(FatalError);
            }
        }
    }

    // Neighbour zoneID.
    neiZoneID.setSize(mesh.nFaces()-mesh.nInternalFaces());

    forAll(neiZoneID, i)
    {
        neiZoneID[i] = zoneID[mesh.faceOwner()[i+mesh.nInternalFaces()]];
    }
    syncTools::swapBoundaryFaceList(mesh, neiZoneID);
}


//splitMeshRegions
//void MeshManipulationFoam::matchRegions
int MeshManipulationFoam::matchRegions
(
    const bool sloppyCellZones,
    const Foam::polyMesh& mesh,

    const Foam::label nCellRegions,
    const Foam::labelList& cellRegion,

    Foam::labelList& regionToZone,
    Foam::wordList& regionNames,
    Foam::labelList& zoneToRegion
)
{
  using namespace Foam;
    const cellZoneMesh& cellZones = mesh.cellZones();

    regionToZone.setSize(nCellRegions, -1);
    regionNames.setSize(nCellRegions);
    zoneToRegion.setSize(cellZones.size(), -1);

    // Get current per cell zoneID
    labelList zoneID(mesh.nCells(), -1);
    labelList neiZoneID(mesh.nFaces()-mesh.nInternalFaces());
    getZoneID(mesh, cellZones, zoneID, neiZoneID);

    // Sizes per cellzone
    labelList zoneSizes(cellZones.size(), 0);
    {
        List<wordList> zoneNames(Pstream::nProcs());
        zoneNames[Pstream::myProcNo()] = cellZones.names();
        Pstream::gatherList(zoneNames);
        Pstream::scatterList(zoneNames);

        forAll(zoneNames, proci)
        {
            if (zoneNames[proci] != zoneNames[0])
            {
                FatalErrorInFunction
                    << "cellZones not synchronised across processors." << endl
                    << "Master has cellZones " << zoneNames[0] << endl
                    << "Processor " << proci
                    << " has cellZones " << zoneNames[proci]
                    << exit(FatalError);
            }
        }

        forAll(cellZones, zoneI)
        {
            zoneSizes[zoneI] = returnReduce
            (
                cellZones[zoneI].size(),
                sumOp<label>()
            );
        }
    }


    if (sloppyCellZones)
    {
        Info<< "Trying to match regions to existing cell zones;"
            << " region can be subset of cell zone." << nl << endl;

        forAll(cellZones, zoneI)
        {
            label regionI = findCorrespondingRegion
            (
                zoneID,
                cellRegion,
                nCellRegions,
                zoneI,
                label(0.5*zoneSizes[zoneI]) // minimum overlap
            );

            if (regionI != -1)
            {
                Info<< "Sloppily matched region " << regionI
                    //<< " size " << regionSizes[regionI]
                    << " to zone " << zoneI << " size " << zoneSizes[zoneI]
                    << endl;
                zoneToRegion[zoneI] = regionI;
                regionToZone[regionI] = zoneI;
                regionNames[regionI] = cellZones[zoneI].name();
            }
        }
    }
    else
    {
        Info<< "Trying to match regions to existing cell zones." << nl << endl;

        forAll(cellZones, zoneI)
        {
            label regionI = findCorrespondingRegion
            (
                zoneID,
                cellRegion,
                nCellRegions,
                zoneI,
                1               // minimum overlap
            );

            if (regionI != -1)
            {
                zoneToRegion[zoneI] = regionI;
                regionToZone[regionI] = zoneI;
                regionNames[regionI] = cellZones[zoneI].name();
            }
        }
    }
  
  int caughtU;
  int countU = 0;
    // Allocate region names for unmatched regions.
    forAll(regionToZone, regionI)
    {
        if (regionToZone[regionI] == -1)
        {
            regionNames[regionI] = "domain" + Foam::name(regionI);
        }
    else{
      caughtU = countU;
    }
    countU++;
    }

return caughtU;
}



//splitMeshRegions
void MeshManipulationFoam::writeCellToRegion
(
  const Foam::fvMesh& mesh, const Foam::labelList& cellRegion
)
{
  using namespace Foam;
    // Write to manual decomposition option
    {
        labelIOList cellToRegion
        (
            IOobject
            (
                "cellToRegion",
                mesh.facesInstance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            cellRegion
        );
        cellToRegion.write();

        Info<< "Writing region per cell file (for manual decomposition) to "
            << cellToRegion.objectPath() << nl << endl;
    }
    // Write for postprocessing
    {
        volScalarField cellToRegion
        (
            IOobject
            (
                "cellToRegion",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0),
            zeroGradientFvPatchScalarField::typeName
        );
        forAll(cellRegion, celli)
        {
            cellToRegion[celli] = cellRegion[celli];
        }
        cellToRegion.write();

        Info<< "Writing region per cell as volScalarField to "
            << cellToRegion.objectPath() << nl << endl;
    }
}


//mergeMeshes
void MeshManipulationFoam::getRootCase
(
  Foam::fileName& casePath
)
{
  using namespace Foam;
    casePath.clean();

    if (casePath.empty() || casePath == ".")
    {
        // handle degenerate form and '.'
        casePath = cwd();
    }
    else if (casePath[0] != '/' && casePath.name() == "..")
    {
        // avoid relative cases ending in '..' - makes for very ugly names
        casePath = cwd()/casePath;
        casePath.clean();
    }
}



//createPatch
void MeshManipulationFoam::changePatchID
(
  const Foam::polyMesh& mesh,
  const Foam::label faceID, 
  const Foam::label patchID,
  Foam::polyTopoChange& meshMod
)
{
  using namespace Foam;
  const label zoneID = mesh.faceZones().whichZone(faceID);

    bool zoneFlip = false;

    if (zoneID >= 0)
    {
        const faceZone& fZone = mesh.faceZones()[zoneID];

        zoneFlip = fZone.flipMap()[fZone.whichFace(faceID)];
    }

    meshMod.setAction
    (
        polyModifyFace
        (
            mesh.faces()[faceID],               // face
            faceID,                             // face ID
            mesh.faceOwner()[faceID],           // owner
            -1,                                 // neighbour
            false,                              // flip flux
            patchID,                            // patch ID
            false,                              // remove from zone
            zoneID,                             // zone ID
            zoneFlip                            // zone flip
        )
    );
}


void MeshManipulationFoam::filterPatches
(
  Foam::polyMesh& mesh,
  const Foam::HashSet<Foam::word>& addedPatchNames
)
{
  using namespace Foam;
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Patches to keep
    DynamicList<polyPatch*> allPatches(patches.size());

    label nOldPatches = returnReduce(patches.size(), sumOp<label>());

    // Copy old patches.
    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        // Note: reduce possible since non-proc patches guaranteed in same order
        if (!isA<processorPolyPatch>(pp))
        {

            // Add if
            // - non zero size
            // - or added from the createPatchDict
            // - or cyclic (since referred to by other cyclic half or
            //   proccyclic)

            if
            (
                addedPatchNames.found(pp.name())
             || returnReduce(pp.size(), sumOp<label>()) > 0
             || isA<coupledPolyPatch>(pp)
            )
            {
                allPatches.append
                (
                    pp.clone
                    (
                        patches,
                        allPatches.size(),
                        pp.size(),
                        pp.start()
                    ).ptr()
                );
            }
            else
            {
                Info<< "Removing zero-sized patch " << pp.name()
                    << " type " << pp.type()
                    << " at position " << patchi << endl;
            }
        }
    }
    // Copy non-empty processor patches
    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (isA<processorPolyPatch>(pp))
        {
            if (pp.size())
            {
                allPatches.append
                (
                    pp.clone
                    (
                        patches,
                        allPatches.size(),
                        pp.size(),
                        pp.start()
                    ).ptr()
                );
            }
            else
            {
                Info<< "Removing empty processor patch " << pp.name()
                    << " at position " << patchi << endl;
            }
        }
    }

    label nAllPatches = returnReduce(allPatches.size(), sumOp<label>());
    if (nAllPatches != nOldPatches)
    {
        Info<< "Removing patches." << endl;
        allPatches.shrink();
        mesh.removeBoundary();
        mesh.addPatches(allPatches);
    }
    else
    {
        Info<< "No patches removed." << endl;
        forAll(allPatches, i)
        {
            delete allPatches[i];
        }
    }
}


void MeshManipulationFoam::dumpCyclicMatch
(
  const Foam::fileName& prefix,
  const Foam::polyMesh& mesh
)
{
  using namespace Foam;
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchi)
    {
        if
        (
            isA<cyclicPolyPatch>(patches[patchi])
         && refCast<const cyclicPolyPatch>(patches[patchi]).owner()
        )
        {
            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(patches[patchi]);

            // Dump patches
            {
                OFstream str(prefix+cycPatch.name()+".obj");
                Pout<< "Dumping " << cycPatch.name()
                    << " faces to " << str.name() << endl;
                meshTools::writeOBJ
                (
                    str,
                    cycPatch,
                    cycPatch.points()
                );
            }

            const cyclicPolyPatch& nbrPatch = cycPatch.neighbPatch();
            {
                OFstream str(prefix+nbrPatch.name()+".obj");
                Pout<< "Dumping " << nbrPatch.name()
                    << " faces to " << str.name() << endl;
                meshTools::writeOBJ
                (
                    str,
                    nbrPatch,
                    nbrPatch.points()
                );
            }


            // Lines between corresponding face centres
            OFstream str(prefix+cycPatch.name()+nbrPatch.name()+"_match.obj");
            label vertI = 0;

            Pout<< "Dumping cyclic match as lines between face centres to "
                << str.name() << endl;

            forAll(cycPatch, facei)
            {
                const point& fc0 = mesh.faceCentres()[cycPatch.start()+facei];
                meshTools::writeOBJ(str, fc0);
                vertI++;
                const point& fc1 = mesh.faceCentres()[nbrPatch.start()+facei];
                meshTools::writeOBJ(str, fc1);
                vertI++;

                str<< "l " << vertI-1 << ' ' << vertI << nl;
            }
        }
    }
}


void MeshManipulationFoam::separateList
(
    const Foam::vectorField& separation,
    Foam::UList<Foam::vector>& field
)
{
  using namespace Foam;
    if (separation.size() == 1)
    {
        // Single value for all.

        forAll(field, i)
        {
            field[i] += separation[0];
        }
    }
    else if (separation.size() == field.size())
    {
        forAll(field, i)
        {
            field[i] += separation[i];
        }
    }
    else
    {
        FatalErrorInFunction
            << "Sizes of field and transformation not equal. field:"
            << field.size() << " transformation:" << separation.size()
            << abort(FatalError);
    }
}


template<class CombineOp>
void MeshManipulationFoam::syncPoints
(
    const Foam::polyMesh& mesh,
    Foam::pointField& points,
    const CombineOp& cop,
    const Foam::point& nullValue
)
{
  using namespace Foam;
    if (points.size() != mesh.nPoints())
    {
        FatalErrorInFunction
            << "Number of values " << points.size()
            << " is not equal to the number of points in the mesh "
            << mesh.nPoints() << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Is there any coupled patch with transformation?
    bool hasTransformation = false;

    if (Pstream::parRun())
    {
        // Send

        forAll(patches, patchi)
        {
            const polyPatch& pp = patches[patchi];

            if
            (
                isA<processorPolyPatch>(pp)
             && pp.nPoints() > 0
             && refCast<const processorPolyPatch>(pp).owner()
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(pp);

                // Get data per patchPoint in neighbouring point numbers.
                pointField patchInfo(procPatch.nPoints(), nullValue);

                const labelList& meshPts = procPatch.meshPoints();
                const labelList& nbrPts = procPatch.neighbPoints();

                forAll(nbrPts, pointi)
                {
                    label nbrPointi = nbrPts[pointi];
                    if (nbrPointi >= 0 && nbrPointi < patchInfo.size())
                    {
                        patchInfo[nbrPointi] = points[meshPts[pointi]];
                    }
                }

                OPstream toNbr
                (
                    Pstream::commsTypes::blocking,
                    procPatch.neighbProcNo()
                );
                toNbr << patchInfo;
            }
        }


        // Receive and set.

        forAll(patches, patchi)
        {
            const polyPatch& pp = patches[patchi];

            if
            (
                isA<processorPolyPatch>(pp)
             && pp.nPoints() > 0
             && !refCast<const processorPolyPatch>(pp).owner()
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(pp);

                pointField nbrPatchInfo(procPatch.nPoints());
                {
                    // We do not know the number of points on the other side
                    // so cannot use Pstream::read.
                    IPstream fromNbr
                    (
                        Pstream::commsTypes::blocking,
                        procPatch.neighbProcNo()
                    );
                    fromNbr >> nbrPatchInfo;
                }
                // Null any value which is not on neighbouring processor
                nbrPatchInfo.setSize(procPatch.nPoints(), nullValue);

                if (!procPatch.parallel())
                {
                    hasTransformation = true;
                    transformList(procPatch.forwardT(), nbrPatchInfo);
                }
                else if (procPatch.separated())
                {
                    hasTransformation = true;
                    separateList(-procPatch.separation(), nbrPatchInfo);
                }

                const labelList& meshPts = procPatch.meshPoints();

                forAll(meshPts, pointi)
                {
                    label meshPointi = meshPts[pointi];
                    points[meshPointi] = nbrPatchInfo[pointi];
                }
            }
        }
    }

    // Do the cyclics.
    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if
        (
            isA<cyclicPolyPatch>(pp)
         && refCast<const cyclicPolyPatch>(pp).owner()
        )
        {
            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(pp);

            const edgeList& coupledPoints = cycPatch.coupledPoints();
            const labelList& meshPts = cycPatch.meshPoints();
            const cyclicPolyPatch& nbrPatch = cycPatch.neighbPatch();
            const labelList& nbrMeshPts = nbrPatch.meshPoints();

            pointField half0Values(coupledPoints.size());

            forAll(coupledPoints, i)
            {
                const edge& e = coupledPoints[i];
                label point0 = meshPts[e[0]];
                half0Values[i] = points[point0];
            }

            if (!cycPatch.parallel())
            {
                hasTransformation = true;
                transformList(cycPatch.reverseT(), half0Values);
            }
            else if (cycPatch.separated())
            {
                hasTransformation = true;
                separateList(cycPatch.separation(), half0Values);
            }

            forAll(coupledPoints, i)
            {
                const edge& e = coupledPoints[i];
                label point1 = nbrMeshPts[e[1]];
                points[point1] = half0Values[i];
            }
        }
    }

    //- Note: hasTransformation is only used for warning messages so
    //  reduction not strictly nessecary.
    //reduce(hasTransformation, orOp<bool>());

    // Synchronize multiple shared points.
    const globalMeshData& pd = mesh.globalData();

    if (pd.nGlobalPoints() > 0)
    {
        if (hasTransformation)
        {
            WarningInFunction
                << "There are decomposed cyclics in this mesh with"
                << " transformations." << endl
                << "This is not supported. The result will be incorrect"
                << endl;
        }


        // Values on shared points.
        pointField sharedPts(pd.nGlobalPoints(), nullValue);

        forAll(pd.sharedPointLabels(), i)
        {
            label meshPointi = pd.sharedPointLabels()[i];
            // Fill my entries in the shared points
            sharedPts[pd.sharedPointAddr()[i]] = points[meshPointi];
        }

        // Combine on master.
        Pstream::listCombineGather(sharedPts, cop);
        Pstream::listCombineScatter(sharedPts);

        // Now we will all have the same information. Merge it back with
        // my local information.
        forAll(pd.sharedPointLabels(), i)
        {
            label meshPointi = pd.sharedPointLabels()[i];
            points[meshPointi] = sharedPts[pd.sharedPointAddr()[i]];
        }
    }
}
