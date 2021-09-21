#include <iostream>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "getDicts.H"
#include "splitMeshRegions.H"

using namespace Foam;

splitMeshRegions::splitMeshRegions() {}

splitMeshRegions::~splitMeshRegions() {}

std::pair<std::vector<int>, std::vector<std::string>> splitMeshRegions::execute(
    const bool useCellZones) {
  // Create controlDict
  bool writeDicts = false;
  std::unique_ptr<getDicts> initFoam;
  initFoam = std::unique_ptr<getDicts>(new getDicts());
  auto controlDict_ = initFoam->createControlDict(writeDicts);
  auto fvSchemes_ = initFoam->createFvSchemes(writeDicts);
  auto fvSolution_ = initFoam->createFvSolution(writeDicts);

  // Create time class without reading controlDict
  Foam::Time runTime(controlDict_.get(), ".", ".");
  Foam::argList::noParallel();

  Foam::word regionName = Foam::fvMesh::defaultRegion;

  Foam::fvMesh mesh(Foam::IOobject(regionName, runTime.timeName(), runTime,
                                   Foam::IOobject::MUST_READ));

  const word oldInstance = mesh.pointsInstance();
  word blockedFacesName;

  int impVar;

  const bool makeCellZones = false;
  const bool largestOnly = false;
  const bool insidePoint = false;
  const bool useCellZonesOnly = false;
  const bool useCellZonesFile = false;
  const bool overwrite = true;
  const bool detectOnly = false;
  const bool sloppyCellZones = false;
  const bool useFaceZones = false;
  const bool prefixRegion = false;

  if ((useCellZonesOnly || useCellZonesFile) &&
      (useCellZones || blockedFacesName.size())) {
    FatalErrorInFunction
        << "You cannot specify both -cellZonesOnly or -cellZonesFileOnly"
        << " (which specify complete split)"
        << " in combination with -blockedFaces or -cellZones"
        << " (which imply a split based on topology)" << exit(FatalError);
  }

  if (useFaceZones) {
    Info << "Using current faceZones to divide inter-region interfaces"
         << " into multiple patches." << nl << endl;
  } else {
    Info << "Creating single patch per inter-region interface." << nl << endl;
  }

  if (insidePoint && largestOnly) {
    FatalErrorInFunction << "You cannot specify both -largestOnly"
                         << " (keep region with most cells)"
                         << " and -insidePoint (keep region containing point)"
                         << exit(FatalError);
  }

  const cellZoneMesh& cellZones = mesh.cellZones();

  // Existing zoneID
  labelList zoneID(mesh.nCells(), -1);
  // Neighbour zoneID.
  labelList neiZoneID(mesh.nFaces() - mesh.nInternalFaces());
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
  if (useCellZonesOnly) {
    Info << "Using current cellZones to split mesh into regions."
         << " This requires all"
         << " cells to be in one and only one cellZone." << nl << endl;

    // label unzonedCelli = findIndex(zoneID, -1);
    label unzonedCelli = zoneID.find(-1);
    if (unzonedCelli != -1) {
      FatalErrorInFunction
          << "For the cellZonesOnly option all cells "
          << "have to be in a cellZone." << endl
          << "Cell " << unzonedCelli << " at"
          << mesh.cellCentres()[unzonedCelli]
          << " is not in a cellZone. There might be more unzoned cells."
          << exit(FatalError);
    }

    cellRegion = zoneID;
    nCellRegions = gMax(cellRegion) + 1;
    regionToZone.setSize(nCellRegions);
    regionNames.setSize(nCellRegions);
    zoneToRegion.setSize(cellZones.size(), -1);

    for (label regionI = 0; regionI < nCellRegions; regionI++) {
      regionToZone[regionI] = regionI;
      zoneToRegion[regionI] = regionI;
      regionNames[regionI] = cellZones[regionI].name();
    }
  }

  // Will be implemented in future
  else if (useCellZonesFile) {
    // const word zoneFile(args["cellZonesFileOnly"]);
    const word zoneFile = "";
    Info << "Reading split from cellZones file " << zoneFile << endl
         << "This requires all"
         << " cells to be in one and only one cellZone." << nl << endl;

    cellZoneMesh newCellZones(
        IOobject(zoneFile, mesh.facesInstance(), polyMesh::meshSubDir, mesh,
                 IOobject::MUST_READ, IOobject::NO_WRITE, false),
        mesh);

    labelList newZoneID(mesh.nCells(), -1);
    labelList newNeiZoneID(mesh.nFaces() - mesh.nInternalFaces());
    getZoneID(mesh, newCellZones, newZoneID, newNeiZoneID);

    // label unzonedCelli = findIndex(newZoneID, -1);
    label unzonedCelli = newZoneID.find(-1);
    if (unzonedCelli != -1) {
      FatalErrorInFunction
          << "For the cellZonesFileOnly option all cells "
          << "have to be in a cellZone." << endl
          << "Cell " << unzonedCelli << " at"
          << mesh.cellCentres()[unzonedCelli]
          << " is not in a cellZone. There might be more unzoned cells."
          << exit(FatalError);
    }
    cellRegion = newZoneID;
    nCellRegions = gMax(cellRegion) + 1;
    zoneToRegion.setSize(newCellZones.size(), -1);
    regionToZone.setSize(nCellRegions);
    regionNames.setSize(nCellRegions);

    for (label regionI = 0; regionI < nCellRegions; regionI++) {
      regionToZone[regionI] = regionI;
      zoneToRegion[regionI] = regionI;
      regionNames[regionI] = newCellZones[regionI].name();
    }
  } else {
    // Determine connected regions
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Mark additional faces that are blocked
    boolList blockedFace;

    // Read from faceSet
    if (blockedFacesName.size()) {
      faceSet blockedFaceSet(mesh, blockedFacesName);
      Info << "Read " << returnReduce(blockedFaceSet.size(), sumOp<label>())
           << " blocked faces from set " << blockedFacesName << nl << endl;

      blockedFace.setSize(mesh.nFaces(), false);

      forAllConstIter(faceSet, blockedFaceSet, iter) {
        blockedFace[iter.key()] = true;
      }
    }

    // Imply from differing cellZones
    if (useCellZones) {
      blockedFace.setSize(mesh.nFaces(), false);

      for (label facei = 0; facei < mesh.nInternalFaces(); facei++) {
        label own = mesh.faceOwner()[facei];
        label nei = mesh.faceNeighbour()[facei];

        if (zoneID[own] != zoneID[nei]) {
          blockedFace[facei] = true;
        }
      }

      // Different cellZones on either side of processor patch.
      forAll(neiZoneID, i) {
        label facei = i + mesh.nInternalFaces();

        if (zoneID[mesh.faceOwner()[facei]] != neiZoneID[i]) {
          blockedFace[facei] = true;
        }
      }
    }

    // Do a topological walk to determine regions
    regionSplit regions(mesh, blockedFace);
    nCellRegions = regions.nRegions();
    cellRegion.transfer(regions);

    // Make up region names. If possible match them to existing zones.
    impVar = matchRegions(sloppyCellZones, mesh, nCellRegions, cellRegion,

                          regionToZone, regionNames, zoneToRegion);

    // Override any default region names if single region selected
    if (largestOnly || insidePoint) {
      forAll(regionToZone, regionI) {
        if (regionToZone[regionI] == -1) {
          if (overwrite) {
            regionNames[regionI] = polyMesh::defaultRegion;
          } else if (insidePoint) {
            regionNames[regionI] = "insidePoint";
          } else if (largestOnly) {
            regionNames[regionI] = "largestOnly";
          }
        }
      }
    }
  }

  Info << endl << "Number of regions:" << nCellRegions << nl << endl;

  // Write decomposition to file
  writeCellToRegion(mesh, cellRegion);

  // Sizes per region
  // ~~~~~~~~~~~~~~~~

  labelList regionSizes(nCellRegions, 0);

  forAll(cellRegion, celli) { regionSizes[cellRegion[celli]]++; }
  forAll(regionSizes, regionI) { reduce(regionSizes[regionI], sumOp<label>()); }

  Info << "Region\tCells" << nl << "------\t-----" << endl;

  forAll(regionSizes, regionI) {
    Info << regionI << '\t' << regionSizes[regionI] << nl;
  }
  Info << endl;

  // Print region to zone
  Info << "Region\tZone\tName" << nl << "------\t----\t----" << endl;
  forAll(regionToZone, regionI) {
    Info << regionI << '\t' << regionToZone[regionI] << '\t'
         << regionNames[regionI] << nl;
  }
  Info << endl;

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

  getInterfaceSizes(mesh, useFaceZones, cellRegion, regionNames,

                    interfaces, interfaceNames, interfaceSizes,
                    faceToInterface);

  Info << "Sizes of interfaces between regions:" << nl << nl
       << "Interface\tRegion\tRegion\tFaces" << nl
       << "---------\t------\t------\t-----" << endl;

  forAll(interfaces, interI) {
    const edge& e = interfaces[interI];

    Info << interI << "\t\t" << e[0] << '\t' << e[1] << '\t'
         << interfaceSizes[interI] << nl;
  }
  Info << endl;

  if (detectOnly) {
    exit(1);  // return 0;
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

  Info << endl;

  // Remove any demand-driven fields ('S', 'V' etc)
  mesh.clearOut();

  if (nCellRegions == 1) {
    Info << "Only one region. Doing nothing." << endl;
  } else if (makeCellZones) {
    Info << "Putting cells into cellZones instead of splitting mesh." << endl;

    // Check if region overlaps with existing zone. If so keep.

    for (label regionI = 0; regionI < nCellRegions; regionI++) {
      label zoneI = regionToZone[regionI];

      if (zoneI != -1) {
        Info << "    Region " << regionI << " : corresponds to existing"
             << " cellZone " << zoneI << ' ' << cellZones[zoneI].name() << endl;
      } else {
        // Create new cellZone.
        labelList regionCells = findIndices(cellRegion, regionI);

        word zoneName = "region" + Foam::name(regionI);

        zoneI = cellZones.findZoneID(zoneName);

        if (zoneI == -1) {
          zoneI = cellZones.size();
          mesh.cellZones().setSize(zoneI + 1);
          mesh.cellZones().set(zoneI,
                               new cellZone(zoneName,     // name
                                            regionCells,  // addressing
                                            zoneI,        // index
                                            cellZones     // cellZoneMesh
                                            ));
        } else {
          mesh.cellZones()[zoneI].clearAddressing();
          mesh.cellZones()[zoneI] = regionCells;
        }
        Info << "    Region " << regionI << " : created new cellZone " << zoneI
             << ' ' << cellZones[zoneI].name() << endl;
      }
    }

    mesh.cellZones().writeOpt() = IOobject::AUTO_WRITE;

    if (!overwrite) {
      runTime++;
      mesh.setInstance(runTime.timeName());
    } else {
      mesh.setInstance(oldInstance);
    }

    Info << "Writing cellZones as new mesh to time " << runTime.timeName() << nl
         << endl;

    mesh.write();

    // Write cellSets for convenience
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Info << "Writing cellSets corresponding to cellZones." << nl << endl;

    forAll(cellZones, zoneI) {
      const cellZone& cz = cellZones[zoneI];

      cellSet(mesh, cz.name(), cz).write();
    }
  } else {
    // Add patches for interfaces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Add all possible patches. Empty ones get filtered later on.
    Info << nl << "Adding patches" << nl << endl;

    labelList interfacePatches(
        addRegionPatches(mesh, regionNames, interfaces, interfaceNames));

    if (!overwrite) {
      runTime++;
    }

    // Create regions
    // ~~~~~~~~~~~~~~

    if (insidePoint) {
      // Will be implemented in future
      // const point insidePoint = args.get<point>("insidePoint");
      const point insidePoint = Foam::Vector<double>();

      label regionI = -1;

      (void)mesh.tetBasePtIs();

      label celli = mesh.findCell(insidePoint);

      Info << nl << "Found point " << insidePoint << " in cell " << celli
           << endl;

      if (celli != -1) {
        regionI = cellRegion[celli];
      }

      reduce(regionI, maxOp<label>());

      Info << nl << "Subsetting region " << regionI << " containing point "
           << insidePoint << endl;

      if (regionI == -1) {
        FatalErrorInFunction
            << "Point " << insidePoint << " is not inside the mesh." << nl
            << "Bounding box of the mesh:" << mesh.bounds() << exit(FatalError);
      }

      createAndWriteRegion(mesh, cellRegion, regionNames, prefixRegion,
                           faceToInterface, interfacePatches, regionI,
                           (overwrite ? oldInstance : runTime.timeName()));
    } else if (largestOnly) {
      label regionI = findMax(regionSizes);

      Info << nl << "Subsetting region " << regionI << " of size "
           << regionSizes[regionI] << " as named region "
           << regionNames[regionI] << endl;

      createAndWriteRegion(mesh, cellRegion, regionNames, prefixRegion,
                           faceToInterface, interfacePatches, regionI,
                           (overwrite ? oldInstance : runTime.timeName()));
    } else {
      // Split all
      for (label regionI = 0; regionI < nCellRegions; regionI++) {
        Info << nl << "Region " << regionI << nl << "-------- " << endl;

        createAndWriteRegion(mesh, cellRegion, regionNames, prefixRegion,
                             faceToInterface, interfacePatches, regionI,
                             (overwrite ? oldInstance : runTime.timeName()));
      }
    }
  }

  Info << "End\n" << endl;

  label largestReg = findMax(regionSizes);
  std::string largestRegMsh = regionNames[largestReg];
  std::vector<std::string> allRegs;
  allRegs.push_back(largestRegMsh);

  // Iterate through wordlist regionNames. First element is the largest
  // region in mesh
  for (int i = 0; i < regionNames.size(); i++) {
    if (regionNames[i] != largestRegMsh) {
      std::string tmpstr = regionNames[i];
      allRegs.push_back(tmpstr);
    }
  }

  std::vector<int> rtrnVec;
  rtrnVec.push_back(impVar);
  rtrnVec.push_back(nCellRegions);
  return std::make_pair(rtrnVec, allRegs);
}

void splitMeshRegions::renamePatches(Foam::fvMesh& mesh,
                                     const Foam::word& prefix,
                                     const Foam::labelList& patchesToRename) {
  polyBoundaryMesh& polyPatches =
      const_cast<polyBoundaryMesh&>(mesh.boundaryMesh());
  forAll(patchesToRename, i) {
    label patchi = patchesToRename[i];
    polyPatch& pp = polyPatches[patchi];

    if (isA<coupledPolyPatch>(pp)) {
      WarningInFunction << "Encountered coupled patch " << pp.name()
                        << ". Will only rename the patch itself,"
                        << " not any referred patches."
                        << " This might have to be done by hand." << endl;
    }

    pp.name() = prefix + '_' + pp.name();
  }
  // Recalculate any demand driven data (e.g. group to name lookup)
  polyPatches.updateMesh();
}

// splitMeshByRegions
template <class GeoField>
void splitMeshRegions::subsetVolFields(const Foam::fvMesh& mesh,
                                       const Foam::fvMesh& subMesh,
                                       const Foam::labelList& cellMap,
                                       const Foam::labelList& faceMap,
                                       const Foam::labelHashSet& addedPatches) {
  const labelList patchMap(identity(mesh.boundaryMesh().size()));

  HashTable<const GeoField*> fields(
      mesh.objectRegistry::lookupClass<GeoField>());
  forAllConstIter(typename HashTable<const GeoField*>, fields, iter) {
    const GeoField& fld = *iter();

    Info << "Mapping field " << fld.name() << endl;

    tmp<GeoField> tSubFld(
        fvMeshSubset::interpolate(fld, subMesh, patchMap, cellMap, faceMap));

    // Hack: set value to 0 for introduced patches (since don't
    //       get initialised.
    forAll(tSubFld().boundaryField(), patchi) {
      if (addedPatches.found(patchi)) {
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

// splitMeshByRegions
template <class GeoField>
void splitMeshRegions::subsetSurfaceFields(
    const Foam::fvMesh& mesh, const Foam::fvMesh& subMesh,
    const Foam::labelList& cellMap, const Foam::labelList& faceMap,
    const Foam::labelHashSet& addedPatches) {
  const labelList patchMap(identity(mesh.boundaryMesh().size()));

  HashTable<const GeoField*> fields(
      mesh.objectRegistry::lookupClass<GeoField>());
  forAllConstIter(typename HashTable<const GeoField*>, fields, iter) {
    const GeoField& fld = *iter();

    Info << "Mapping field " << fld.name() << endl;

    tmp<GeoField> tSubFld(
        fvMeshSubset::interpolate(fld, subMesh, patchMap, cellMap, faceMap));

    // Hack: set value to 0 for introduced patches (since don't
    //       get initialised.
    forAll(tSubFld().boundaryField(), patchi) {
      if (addedPatches.found(patchi)) {
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

// splitMeshByRegions
Foam::labelList splitMeshRegions::getNonRegionCells(
    const Foam::labelList& cellRegion, const Foam::label regionI) {
  DynamicList<label> nonRegionCells(cellRegion.size());
  forAll(cellRegion, celli) {
    if (cellRegion[celli] != regionI) {
      nonRegionCells.append(celli);
    }
  }
  return nonRegionCells.shrink();
}

// splitMeshByRegions
void splitMeshRegions::addToInterface(
    const Foam::polyMesh& mesh, const Foam::label zoneID,
    const Foam::label ownRegion, const Foam::label neiRegion,
    Foam::EdgeMap<Foam::Map<Foam::label>>& regionsToSize) {
  edge interface(min(ownRegion, neiRegion), max(ownRegion, neiRegion));

  EdgeMap<Map<label>>::iterator iter = regionsToSize.find(interface);

  if (iter != regionsToSize.end()) {
    // Check if zone present
    Map<label>::iterator zoneFnd = iter().find(zoneID);
    if (zoneFnd != iter().end()) {
      // Found zone. Increment count.
      zoneFnd()++;
    } else {
      // New or no zone. Insert with count 1.
      iter().insert(zoneID, 1);
    }
  } else {
    // Create new interface of size 1.
    Map<label> zoneToSize;
    zoneToSize.insert(zoneID, 1);
    regionsToSize.insert(interface, zoneToSize);
  }
}

// splitMeshByRegions
void splitMeshRegions::getInterfaceSizes(
    const Foam::polyMesh& mesh, const bool useFaceZones,
    const Foam::labelList& cellRegion, const Foam::wordList& regionNames,

    Foam::edgeList& interfaces,
    Foam::List<Foam::Pair<Foam::word>>& interfaceNames,
    Foam::labelList& interfaceSizes, Foam::labelList& faceToInterface) {
  // From region-region to faceZone (or -1) to number of faces.

  EdgeMap<Map<label>> regionsToSize;

  // Internal faces
  // ~~~~~~~~~~~~~~

  forAll(mesh.faceNeighbour(), facei) {
    label ownRegion = cellRegion[mesh.faceOwner()[facei]];
    label neiRegion = cellRegion[mesh.faceNeighbour()[facei]];

    if (ownRegion != neiRegion) {
      addToInterface(mesh,
                     (useFaceZones ? mesh.faceZones().whichZone(facei) : -1),
                     ownRegion, neiRegion, regionsToSize);
    }
  }

  // Boundary faces
  // ~~~~~~~~~~~~~~

  // Neighbour cellRegion.
  labelList coupledRegion(mesh.nFaces() - mesh.nInternalFaces());

  forAll(coupledRegion, i) {
    label celli = mesh.faceOwner()[i + mesh.nInternalFaces()];
    coupledRegion[i] = cellRegion[celli];
  }
  syncTools::swapBoundaryFaceList(mesh, coupledRegion);

  forAll(coupledRegion, i) {
    label facei = i + mesh.nInternalFaces();
    label ownRegion = cellRegion[mesh.faceOwner()[facei]];
    label neiRegion = coupledRegion[i];

    if (ownRegion != neiRegion) {
      addToInterface(mesh,
                     (useFaceZones ? mesh.faceZones().whichZone(facei) : -1),
                     ownRegion, neiRegion, regionsToSize);
    }
  }

  if (Pstream::parRun()) {
    if (Pstream::master()) {
      // Receive and add to my sizes
      for (int slave = Pstream::firstSlave(); slave <= Pstream::lastSlave();
           slave++) {
        IPstream fromSlave(Pstream::commsTypes::blocking, slave);

        EdgeMap<Map<label>> slaveSizes(fromSlave);

        forAllConstIter(EdgeMap<Map<label>>, slaveSizes, slaveIter) {
          EdgeMap<Map<label>>::iterator masterIter =
              regionsToSize.find(slaveIter.key());

          if (masterIter != regionsToSize.end()) {
            // Same inter-region
            const Map<label>& slaveInfo = slaveIter();
            Map<label>& masterInfo = masterIter();

            forAllConstIter(Map<label>, slaveInfo, iter) {
              label zoneID = iter.key();
              label slaveSize = iter();

              Map<label>::iterator zoneFnd = masterInfo.find(zoneID);
              if (zoneFnd != masterInfo.end()) {
                zoneFnd() += slaveSize;
              } else {
                masterInfo.insert(zoneID, slaveSize);
              }
            }
          } else {
            regionsToSize.insert(slaveIter.key(), slaveIter());
          }
        }
      }
    } else {
      // Send to master
      {
        OPstream toMaster(Pstream::commsTypes::blocking, Pstream::masterNo());
        toMaster << regionsToSize;
      }
    }
  }

  // Rework

  Pstream::scatter(regionsToSize);

  // Now we have the global sizes of all inter-regions.
  // Invert this on master and distribute.
  label nInterfaces = 0;
  forAllConstIter(EdgeMap<Map<label>>, regionsToSize, iter) {
    const Map<label>& info = iter();
    nInterfaces += info.size();
  }

  interfaces.setSize(nInterfaces);
  interfaceNames.setSize(nInterfaces);
  interfaceSizes.setSize(nInterfaces);
  EdgeMap<Map<label>> regionsToInterface(nInterfaces);

  nInterfaces = 0;
  forAllConstIter(EdgeMap<Map<label>>, regionsToSize, iter) {
    const edge& e = iter.key();
    const word& name0 = regionNames[e[0]];
    const word& name1 = regionNames[e[1]];

    const Map<label>& info = iter();
    forAllConstIter(Map<label>, info, infoIter) {
      interfaces[nInterfaces] = iter.key();
      label zoneID = infoIter.key();
      if (zoneID == -1) {
        interfaceNames[nInterfaces] =
            Pair<word>(name0 + "_to_" + name1, name1 + "_to_" + name0);
      } else {
        const word& zoneName = mesh.faceZones()[zoneID].name();
        interfaceNames[nInterfaces] =
            Pair<word>(zoneName + "_" + name0 + "_to_" + name1,
                       zoneName + "_" + name1 + "_to_" + name0);
      }
      interfaceSizes[nInterfaces] = infoIter();

      if (regionsToInterface.found(e)) {
        regionsToInterface[e].insert(zoneID, nInterfaces);
      } else {
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

  forAll(mesh.faceNeighbour(), facei) {
    label ownRegion = cellRegion[mesh.faceOwner()[facei]];
    label neiRegion = cellRegion[mesh.faceNeighbour()[facei]];

    if (ownRegion != neiRegion) {
      label zoneID = -1;
      if (useFaceZones) {
        zoneID = mesh.faceZones().whichZone(facei);
      }

      edge interface(min(ownRegion, neiRegion), max(ownRegion, neiRegion));

      faceToInterface[facei] = regionsToInterface[interface][zoneID];
    }
  }
  forAll(coupledRegion, i) {
    label facei = i + mesh.nInternalFaces();
    label ownRegion = cellRegion[mesh.faceOwner()[facei]];
    label neiRegion = coupledRegion[i];

    if (ownRegion != neiRegion) {
      label zoneID = -1;
      if (useFaceZones) {
        zoneID = mesh.faceZones().whichZone(facei);
      }

      edge interface(min(ownRegion, neiRegion), max(ownRegion, neiRegion));

      faceToInterface[facei] = regionsToInterface[interface][zoneID];
    }
  }
}

// splitMeshByRegions
Foam::autoPtr<Foam::mapPolyMesh> splitMeshRegions::createRegionMesh(
    const Foam::fvMesh& mesh,
    // Region info
    const Foam::labelList& cellRegion, const Foam::label regionI,
    const Foam::word& regionName,
    // Interface info
    const Foam::labelList& interfacePatches,
    const Foam::labelList& faceToInterface,

    Foam::autoPtr<Foam::fvMesh>& newMesh) {
  // Create dummy system/fv*
  {
    IOobject io("fvSchemes", mesh.time().system(), regionName, mesh,
                IOobject::NO_READ, IOobject::NO_WRITE, false);

    Info << "Testing:" << io.objectPath() << endl;

    if (!io.typeHeaderOk<IOdictionary>(true))
    // if (!exists(io.objectPath()))
    {
      Info << "Writing dummy " << regionName / io.name() << endl;
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
    IOobject io("fvSolution", mesh.time().system(), regionName, mesh,
                IOobject::NO_READ, IOobject::NO_WRITE, false);

    if (!io.typeHeaderOk<IOdictionary>(true))
    // if (!exists(io.objectPath()))
    {
      Info << "Writing dummy " << regionName / io.name() << endl;
      dictionary dummyDict;
      IOdictionary(io, dummyDict).regIOobject::write();
    }
  }

  // Neighbour cellRegion.
  labelList coupledRegion(mesh.nFaces() - mesh.nInternalFaces());

  forAll(coupledRegion, i) {
    label celli = mesh.faceOwner()[i + mesh.nInternalFaces()];
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
  forAll(exposedFaces, i) {
    label facei = exposedFaces[i];
    label interfacei = faceToInterface[facei];

    label ownRegion = cellRegion[mesh.faceOwner()[facei]];
    label neiRegion = -1;

    if (mesh.isInternalFace(facei)) {
      neiRegion = cellRegion[mesh.faceNeighbour()[facei]];
    } else {
      neiRegion = coupledRegion[facei - mesh.nInternalFaces()];
    }

    // Check which side is being kept - determines which of the two
    // patches will be used.

    label otherRegion = -1;

    if (ownRegion == regionI && neiRegion != regionI) {
      otherRegion = neiRegion;
    } else if (ownRegion != regionI && neiRegion == regionI) {
      otherRegion = ownRegion;
    } else {
      FatalErrorInFunction << "Exposed face:" << facei
                           << " fc:" << mesh.faceCentres()[facei]
                           << " has owner region " << ownRegion
                           << " and neighbour region " << neiRegion
                           << " when handling region:" << regionI
                           << exit(FatalError);
    }

    // Find the patch.
    if (regionI < otherRegion) {
      exposedPatchIDs[i] = interfacePatches[interfacei];
    } else {
      exposedPatchIDs[i] = interfacePatches[interfacei] + 1;
    }
  }

  // Remove faces
  cellRemover.setRefinement(cellsToRemove, exposedFaces, exposedPatchIDs,
                            meshMod);

  autoPtr<mapPolyMesh> map =
      meshMod.makeMesh(newMesh,
                       IOobject(regionName, mesh.time().timeName(), mesh.time(),
                                IOobject::NO_READ, IOobject::AUTO_WRITE),
                       mesh);

  return map;
}

// splitMeshRegions
// void splitMeshRegions::createAndWriteRegion
void splitMeshRegions::createAndWriteRegion(
    const Foam::fvMesh& mesh, const Foam::labelList& cellRegion,
    const Foam::wordList& regionNames, const bool prefixRegion,
    const Foam::labelList& faceToInterface,
    const Foam::labelList& interfacePatches, const Foam::label regionI,
    const Foam::word& newMeshInstance) {
  Info << "Creating mesh for region " << regionI << ' ' << regionNames[regionI]
       << endl;

  autoPtr<fvMesh> newMesh;
  autoPtr<mapPolyMesh> map =
      createRegionMesh(mesh, cellRegion, regionI, regionNames[regionI],
                       interfacePatches, faceToInterface, newMesh);

  // Make map of all added patches
  labelHashSet addedPatches(2 * interfacePatches.size());
  forAll(interfacePatches, interfacei) {
    addedPatches.insert(interfacePatches[interfacei]);
    addedPatches.insert(interfacePatches[interfacei] + 1);
  }

  Info << "Mapping fields" << endl;

  // Map existing fields
  newMesh().updateMesh(map());

  // Add subsetted fields
  subsetVolFields<volScalarField>(mesh, newMesh(), map().cellMap(),
                                  map().faceMap(), addedPatches);
  subsetVolFields<volVectorField>(mesh, newMesh(), map().cellMap(),
                                  map().faceMap(), addedPatches);
  subsetVolFields<volSphericalTensorField>(mesh, newMesh(), map().cellMap(),
                                           map().faceMap(), addedPatches);
  subsetVolFields<volSymmTensorField>(mesh, newMesh(), map().cellMap(),
                                      map().faceMap(), addedPatches);
  subsetVolFields<volTensorField>(mesh, newMesh(), map().cellMap(),
                                  map().faceMap(), addedPatches);

  subsetSurfaceFields<surfaceScalarField>(mesh, newMesh(), map().cellMap(),
                                          map().faceMap(), addedPatches);
  subsetSurfaceFields<surfaceVectorField>(mesh, newMesh(), map().cellMap(),
                                          map().faceMap(), addedPatches);
  subsetSurfaceFields<surfaceSphericalTensorField>(
      mesh, newMesh(), map().cellMap(), map().faceMap(), addedPatches);
  subsetSurfaceFields<surfaceSymmTensorField>(mesh, newMesh(), map().cellMap(),
                                              map().faceMap(), addedPatches);
  subsetSurfaceFields<surfaceTensorField>(mesh, newMesh(), map().cellMap(),
                                          map().faceMap(), addedPatches);

  const polyBoundaryMesh& newPatches = newMesh().boundaryMesh();
  newPatches.checkParallelSync(true);

  // Delete empty patches
  // ~~~~~~~~~~~~~~~~~~~~

  // Create reordering list to move patches-to-be-deleted to end
  labelList oldToNew(newPatches.size(), -1);
  DynamicList<label> sharedPatches(newPatches.size());
  label newI = 0;

  Info << "Deleting empty patches" << endl;

  // Assumes all non-proc boundaries are on all processors!
  forAll(newPatches, patchi) {
    const polyPatch& pp = newPatches[patchi];

    if (!isA<processorPolyPatch>(pp)) {
      if (returnReduce(pp.size(), sumOp<label>()) > 0) {
        oldToNew[patchi] = newI;
        if (!addedPatches.found(patchi)) {
          sharedPatches.append(newI);
        }
        newI++;
      }
    }
  }

  // Same for processor patches (but need no reduction)
  forAll(newPatches, patchi) {
    const polyPatch& pp = newPatches[patchi];

    if (isA<processorPolyPatch>(pp) && pp.size()) {
      oldToNew[patchi] = newI++;
    }
  }

  const label nNewPatches = newI;

  // Move all deleteable patches to the end
  forAll(oldToNew, patchi) {
    if (oldToNew[patchi] == -1) {
      oldToNew[patchi] = newI++;
    }
  }

  // reorderPatches(newMesh(), oldToNew, nNewPatches);
  fvMeshTools::reorderPatches(newMesh(), oldToNew, nNewPatches, true);

  // Rename shared patches with region name
  if (prefixRegion) {
    Info << "Prefixing patches with region name" << endl;

    renamePatches(newMesh(), regionNames[regionI], sharedPatches);
  }

  Info << "Writing new mesh" << endl;

  newMesh().setInstance(newMeshInstance);
  newMesh().write();

  // Write addressing files like decomposePar
  Info << "Writing addressing to base mesh" << endl;

  labelIOList pointProcAddressing(
      IOobject("pointRegionAddressing", newMesh().facesInstance(),
               newMesh().meshSubDir, newMesh(), IOobject::NO_READ,
               IOobject::NO_WRITE, false),
      map().pointMap());
  Info << "Writing map " << pointProcAddressing.name() << " from region"
       << regionI << " points back to base mesh." << endl;
  pointProcAddressing.write();

  labelIOList faceProcAddressing(
      IOobject("faceRegionAddressing", newMesh().facesInstance(),
               newMesh().meshSubDir, newMesh(), IOobject::NO_READ,
               IOobject::NO_WRITE, false),
      newMesh().nFaces());
  forAll(faceProcAddressing, facei) {
    // face + turning index. (see decomposePar)
    // Is the face pointing in the same direction?
    label oldFacei = map().faceMap()[facei];

    if (map().cellMap()[newMesh().faceOwner()[facei]] ==
        mesh.faceOwner()[oldFacei]) {
      faceProcAddressing[facei] = oldFacei + 1;
    } else {
      faceProcAddressing[facei] = -(oldFacei + 1);
    }
  }
  Info << "Writing map " << faceProcAddressing.name() << " from region"
       << regionI << " faces back to base mesh." << endl;
  faceProcAddressing.write();

  labelIOList cellProcAddressing(
      IOobject("cellRegionAddressing", newMesh().facesInstance(),
               newMesh().meshSubDir, newMesh(), IOobject::NO_READ,
               IOobject::NO_WRITE, false),
      map().cellMap());
  Info << "Writing map " << cellProcAddressing.name() << " from region"
       << regionI << " cells back to base mesh." << endl;
  cellProcAddressing.write();

  labelIOList boundaryProcAddressing(
      IOobject("boundaryRegionAddressing", newMesh().facesInstance(),
               newMesh().meshSubDir, newMesh(), IOobject::NO_READ,
               IOobject::NO_WRITE, false),
      labelList(nNewPatches, -1));
  forAll(oldToNew, i) {
    if (!addedPatches.found(i)) {
      label newI = oldToNew[i];
      if (newI >= 0 && newI < nNewPatches) {
        boundaryProcAddressing[oldToNew[i]] = i;
      }
    }
  }
  Info << "Writing map " << boundaryProcAddressing.name() << " from region"
       << regionI << " boundary back to base mesh." << endl;
  boundaryProcAddressing.write();
}

// splitMeshRegions
Foam::labelList splitMeshRegions::addRegionPatches(
    Foam::fvMesh& mesh, const Foam::wordList& regionNames,
    const Foam::edgeList& interfaces,
    const Foam::List<Foam::Pair<Foam::word>>& interfaceNames) {
  Info << nl << "Adding patches" << nl << endl;

  labelList interfacePatches(interfaces.size());

  forAll(interfaces, interI) {
    const edge& e = interfaces[interI];
    const Pair<word>& names = interfaceNames[interI];

    // Info<< "For interface " << interI
    //    << " between regions " << e
    //    << " trying to add patches " << names << endl;

    mappedWallPolyPatch patch1(names[0],
                               0,                  // overridden
                               0,                  // overridden
                               0,                  // overridden
                               regionNames[e[1]],  // sampleRegion
                               mappedPatchBase::NEARESTPATCHFACE,
                               names[1],     // samplePatch
                               point::zero,  // offset
                               mesh.boundaryMesh());

    interfacePatches[interI] =
        fvMeshTools::addPatch(mesh, patch1,
                              dictionary(),  // optional per field value
                              calculatedFvPatchField<scalar>::typeName,
                              true  // validBoundary
        );

    mappedWallPolyPatch patch2(names[1], 0, 0, 0,
                               regionNames[e[0]],  // sampleRegion
                               mappedPatchBase::NEARESTPATCHFACE, names[0],
                               point::zero,  // offset
                               mesh.boundaryMesh());
    fvMeshTools::addPatch(mesh, patch2,
                          dictionary(),  // optional per field value
                          calculatedFvPatchField<scalar>::typeName,
                          true  // validBoundary
    );

    Info << "For interface between region " << regionNames[e[0]] << " and "
         << regionNames[e[1]] << " added patches" << endl
         << "    " << interfacePatches[interI] << "\t"
         << mesh.boundaryMesh()[interfacePatches[interI]].name() << endl
         << "    " << interfacePatches[interI] + 1 << "\t"
         << mesh.boundaryMesh()[interfacePatches[interI] + 1].name() << endl;
  }
  return interfacePatches;
}

// splitMeshRegions
Foam::label splitMeshRegions::findCorrespondingRegion(
    const Foam::labelList& existingZoneID,  // per cell the (unique) zoneID
    const Foam::labelList& cellRegion, const Foam::label nCellRegions,
    const Foam::label zoneI, const Foam::label minOverlapSize) {
  // Per region the number of cells in zoneI
  labelList cellsInZone(nCellRegions, 0);

  forAll(cellRegion, celli) {
    if (existingZoneID[celli] == zoneI) {
      cellsInZone[cellRegion[celli]]++;
    }
  }

  Pstream::listCombineGather(cellsInZone, plusEqOp<label>());
  Pstream::listCombineScatter(cellsInZone);

  // Pick region with largest overlap of zoneI
  label regionI = findMax(cellsInZone);

  if (cellsInZone[regionI] < minOverlapSize) {
    // Region covers too little of zone. Not good enough.
    regionI = -1;
  } else {
    // Check that region contains no cells that aren't in cellZone.
    forAll(cellRegion, celli) {
      if (cellRegion[celli] == regionI && existingZoneID[celli] != zoneI) {
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

// splitMeshRegions
void splitMeshRegions::getZoneID(const Foam::polyMesh& mesh,
                                 const Foam::cellZoneMesh& cellZones,
                                 Foam::labelList& zoneID,
                                 Foam::labelList& neiZoneID) {
  // Existing zoneID
  zoneID.setSize(mesh.nCells());
  zoneID = -1;

  forAll(cellZones, zoneI) {
    const cellZone& cz = cellZones[zoneI];

    forAll(cz, i) {
      label celli = cz[i];
      if (zoneID[celli] == -1) {
        zoneID[celli] = zoneI;
      } else {
        FatalErrorInFunction
            << "Cell " << celli << " with cell centre "
            << mesh.cellCentres()[celli]
            << " is multiple zones. This is not allowed." << endl
            << "It is in zone " << cellZones[zoneID[celli]].name()
            << " and in zone " << cellZones[zoneI].name() << exit(FatalError);
      }
    }
  }

  // Neighbour zoneID.
  neiZoneID.setSize(mesh.nFaces() - mesh.nInternalFaces());

  forAll(neiZoneID, i) {
    neiZoneID[i] = zoneID[mesh.faceOwner()[i + mesh.nInternalFaces()]];
  }
  syncTools::swapBoundaryFaceList(mesh, neiZoneID);
}

// splitMeshRegions
// void splitMeshRegions::matchRegions
int splitMeshRegions::matchRegions(const bool sloppyCellZones,
                                   const Foam::polyMesh& mesh,

                                   const Foam::label nCellRegions,
                                   const Foam::labelList& cellRegion,

                                   Foam::labelList& regionToZone,
                                   Foam::wordList& regionNames,
                                   Foam::labelList& zoneToRegion) {
  const cellZoneMesh& cellZones = mesh.cellZones();

  regionToZone.setSize(nCellRegions, -1);
  regionNames.setSize(nCellRegions);
  zoneToRegion.setSize(cellZones.size(), -1);

  // Get current per cell zoneID
  labelList zoneID(mesh.nCells(), -1);
  labelList neiZoneID(mesh.nFaces() - mesh.nInternalFaces());
  getZoneID(mesh, cellZones, zoneID, neiZoneID);

  // Sizes per cellzone
  labelList zoneSizes(cellZones.size(), 0);
  {
    List<wordList> zoneNames(Pstream::nProcs());
    zoneNames[Pstream::myProcNo()] = cellZones.names();
    Pstream::gatherList(zoneNames);
    Pstream::scatterList(zoneNames);

    forAll(zoneNames, proci) {
      if (zoneNames[proci] != zoneNames[0]) {
        FatalErrorInFunction
            << "cellZones not synchronised across processors." << endl
            << "Master has cellZones " << zoneNames[0] << endl
            << "Processor " << proci << " has cellZones " << zoneNames[proci]
            << exit(FatalError);
      }
    }

    forAll(cellZones, zoneI) {
      zoneSizes[zoneI] = returnReduce(cellZones[zoneI].size(), sumOp<label>());
    }
  }

  if (sloppyCellZones) {
    Info << "Trying to match regions to existing cell zones;"
         << " region can be subset of cell zone." << nl << endl;

    forAll(cellZones, zoneI) {
      label regionI = findCorrespondingRegion(
          zoneID, cellRegion, nCellRegions, zoneI,
          label(0.5 * zoneSizes[zoneI])  // minimum overlap
      );

      if (regionI != -1) {
        Info << "Sloppily matched region "
             << regionI
             //<< " size " << regionSizes[regionI]
             << " to zone " << zoneI << " size " << zoneSizes[zoneI] << endl;
        zoneToRegion[zoneI] = regionI;
        regionToZone[regionI] = zoneI;
        regionNames[regionI] = cellZones[zoneI].name();
      }
    }
  } else {
    Info << "Trying to match regions to existing cell zones." << nl << endl;

    forAll(cellZones, zoneI) {
      label regionI =
          findCorrespondingRegion(zoneID, cellRegion, nCellRegions, zoneI,
                                  1  // minimum overlap
          );

      if (regionI != -1) {
        zoneToRegion[zoneI] = regionI;
        regionToZone[regionI] = zoneI;
        regionNames[regionI] = cellZones[zoneI].name();
      }
    }
  }

  int caughtU;
  int countU = 0;
  // Allocate region names for unmatched regions.
  forAll(regionToZone, regionI) {
    if (regionToZone[regionI] == -1) {
      regionNames[regionI] = "domain" + Foam::name(regionI);
    } else {
      caughtU = countU;
    }
    countU++;
  }

  return caughtU;
}

// splitMeshRegions
void splitMeshRegions::writeCellToRegion(const Foam::fvMesh& mesh,
                                         const Foam::labelList& cellRegion) {
  // Write to manual decomposition option
  {
    labelIOList cellToRegion(
        IOobject("cellToRegion", mesh.facesInstance(), mesh, IOobject::NO_READ,
                 IOobject::NO_WRITE, false),
        cellRegion);
    cellToRegion.write();

    Info << "Writing region per cell file (for manual decomposition) to "
         << cellToRegion.objectPath() << nl << endl;
  }
  // Write for postprocessing
  {
    volScalarField cellToRegion(
        IOobject("cellToRegion", mesh.time().timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE, false),
        mesh, dimensionedScalar("zero", dimless, 0),
        zeroGradientFvPatchScalarField::typeName);
    forAll(cellRegion, celli) { cellToRegion[celli] = cellRegion[celli]; }
    cellToRegion.write();

    Info << "Writing region per cell as volScalarField to "
         << cellToRegion.objectPath() << nl << endl;
  }
}