#if defined(_MSC_VER) && !defined(_USE_MATH_DEFINES)
#  define _USE_MATH_DEFINES
#endif

#include "Mesh/foamGeoMesh.H"

#include <iostream>
#include <set>
#include <utility>
#include <string>
#include <boost/filesystem.hpp>

#include <fvCFD.H>
#include <fvMesh.H>
#include <fileName.H>
#include <cellModeller.H>

#include <gmsh.h>
#include <vtkPoints.h>
#include <vtkIdList.h>
#include <vtkGenericCell.h>
#include <vtkDoubleArray.h>
#include <vtkDataArray.h>
#include <vtkStringArray.h>
#include <vtkFieldData.h>

#include <foamVTKTopo.H>
#include <ZoneMesh.H>
#include <cellZone.H>
#include <topoSetSource.H>
#include <globalMeshData.H>
#include <timeSelector.H>
#include <IOobjectList.H>
#include <cellZoneSet.H>
#include <faceZoneSet.H>
#include <pointZoneSet.H>

#include "AuxiliaryFunctions.H"

namespace NEM {
namespace MSH {

vtkStandardNewMacro(foamGeoMesh)

foamGeoMesh *foamGeoMesh::Read(const std::string &fileName,
                               const std::string &geoArrayName) {
  auto foamGM = new foamGeoMesh();
  Foam::word regionName;
  if (fileName.empty()) {
    regionName = Foam::fvMesh::defaultRegion;
  } else {
    regionName = fileName;
  }

  foamGM->controlDict_ = std::unique_ptr<Foam::dictionary>(
                                           new Foam::dictionary("controlDict"));
  Foam::dictionary fmfle("FoamFile");
  fmfle.add("version", "2.0");
  fmfle.add("format", "ascii");
  fmfle.add("class", "dictionary");
  fmfle.add("location", "\"system\"");
  fmfle.add("object", "controlDict");
  foamGM->controlDict_->add("FoamFile",fmfle);
  foamGM->controlDict_->add("deltaT",1);
  foamGM->controlDict_->add("startTime",0);
  foamGM->controlDict_->add("writeInterval",1);
  
  foamGM->runTime_ = std::unique_ptr<Foam::Time>(
                          new Foam::Time(foamGM->controlDict_.get(), ".", "."));

  foamGM->fmesh_.reset(new Foam::fvMesh(
      Foam::IOobject(regionName, foamGM->runTime_->timeName(),
                     *foamGM->runTime_, Foam::IOobject::MUST_READ)));

  foamGM->setGeoMesh(foam2GM(foamGM->fmesh_.get()));

  return foamGM;
}

foamGeoMesh::foamGeoMesh()
    :foamGeoMesh(nullptr) {
  InitializeFoam();
}

foamGeoMesh::foamGeoMesh(Foam::fvMesh *foamMesh,
                         const std::string &phyGrpArrayName)
    :geoMeshBase(foam2GM(foamMesh,phyGrpArrayName)),fmesh_(foamMesh) {
  std::cout << "foamGeoMesh Constructed" << std::endl;
}

foamGeoMesh::~foamGeoMesh() {
  std::cout << "foamGeoMesh destructed" << std::endl;
}

geoMeshBase::GeoMesh foamGeoMesh::foam2GM(Foam::fvMesh *foamMesh,
                                          const std::string &phyGrpArrayName) {
  // create vtk database
  vtkSmartPointer<vtkUnstructuredGrid> vtkdataSet
      = vtkSmartPointer<vtkUnstructuredGrid>::New();

  // Check for cellsets and cellzones
  if (!foamMesh)
  {
    return {vtkdataSet, "", "", {}};
  } else {
    const cellZoneMesh& cellZones = foamMesh->cellZones();   // get cellZoneMesh
    label zoneI = cellZones.size();  
    
    // creating equivalent vtk topology from fvMesh
    // by default polyhedral cells will be decomposed to 
    // tets and pyramids. Additional points will be added
    // to underlying fvMesh.
    std::cout << "Performing topological decomposition.\n";
    foamVTKTopo topo(*foamMesh);

    // point data
    Foam::pointField pf = foamMesh->points();
    vtkSmartPointer<vtkPoints> points
        = vtkSmartPointer<vtkPoints>::New();
    for (int ipt=0; ipt<foamMesh->nPoints(); ipt++)
        points->InsertNextPoint(
                pf[ipt].x(),
                pf[ipt].y(),
                pf[ipt].z()
                );
    vtkdataSet->SetPoints(points);

    // cell data
    std::vector<int> pntIds;
    int nCelPnts = 0;
    for (int icl=0; icl<topo.vertLabels().size(); icl++) {
        nCelPnts = topo.vertLabels()[icl].size();
        pntIds.resize(nCelPnts, -1);
        for (int ip=0; ip< nCelPnts; ip++)
            pntIds[ip] = topo.vertLabels()[icl][ip];
        vtkSmartPointer<vtkIdList> vtkCellIds = vtkSmartPointer<vtkIdList>::New();
        vtkCellIds->SetNumberOfIds(pntIds.size());
        for (auto pit = pntIds.begin(); pit!= pntIds.end(); pit++)
            vtkCellIds->SetId(pit-pntIds.begin(), *pit);
        vtkdataSet->InsertNextCell(topo.cellTypes()[icl],vtkCellIds);
    }

    // Check for cellZones and extract physical groups
    List<scalar> physGrps(foamMesh->nCells(),0.0);
    if (zoneI > 0) {
      for (int i=0; i<zoneI; i++) {
        const cellZone& cz = cellZones[i];
        forAll(cz, j)
          physGrps[cz[j]] = i;
      }
      vtkSmartPointer<vtkDoubleArray> pg = vtkSmartPointer<vtkDoubleArray>::New();
      pg->SetName(phyGrpArrayName.c_str());
      pg->SetNumberOfComponents(1);
      for (int i = 0; i < foamMesh->nCells(); i++)
        pg->InsertNextTuple1(physGrps[i]);
      vtkdataSet->GetCellData()->AddArray(pg);
    }

    // Create boundary patches arrays and add to VTK
    const polyBoundaryMesh& patches = foamMesh->boundaryMesh();
    wordList patchNames = patches.names();
    wordList patchTypes = patches.types();
    const labelList &patchIds = patches.patchID();
    label patchSize = patches.size();
    std::vector<int> face_patch_map; // Which patch the current face is in
    label nFaces = foamMesh->nFaces();
    for (int i=0; i<nFaces; i++)
      face_patch_map.push_back(int(patches.whichPatch(i)));

    // Create sideSets
    vtkSmartPointer<vtkPolyData> sideSet = vtkSmartPointer<vtkPolyData>::New();
    faceList allFaces = foamMesh->faces();
    labelList faceOwners = foamMesh->faceOwner();
    pointField allPoints = foamMesh->points();

    // Filter internal faces out
    int nonInternalFaces = 0;
    for (int i=0; i<patchSize; i++) nonInternalFaces += patches[i].size();
    faceList sideSetOnlyFaces(nonInternalFaces);
    labelList sideSetOnlyFaceLabels(nonInternalFaces);
    labelList sideSetOnlyOwners(nonInternalFaces);
    int indx = 0;
    for (int i=0; i<face_patch_map.size(); i++) {
      if (face_patch_map[i] != -1) {
        sideSetOnlyFaces[indx] = allFaces[i];
        sideSetOnlyOwners[indx] = faceOwners[i];
        sideSetOnlyFaceLabels[indx] = i;
        indx++;
      }
    }

    sideSet->SetPoints(vtkdataSet->GetPoints());
    sideSet->Allocate();

    // Start adding cells (2D faces)
    std::vector<int> facePntIds;
    for (int i=0; i<sideSetOnlyFaces.size(); i++) {
      if (sideSetOnlyFaces[i].size() == 3) {
        // add triangle face (#5)
        vtkSmartPointer<vtkIdList> vtkCellIds = vtkSmartPointer<vtkIdList>::New();
        vtkCellIds->SetNumberOfIds(sideSetOnlyFaces[i].size());
        for (int j=0; j<sideSetOnlyFaces[i].size(); j++)
          vtkCellIds->SetId(j,sideSetOnlyFaces[i][j]);
        sideSet->InsertNextCell(5,vtkCellIds);
      } else if (sideSetOnlyFaces[i].size() == 4) {
        // add quad face (#9)
        vtkSmartPointer<vtkIdList> vtkCellIds = vtkSmartPointer<vtkIdList>::New();
        vtkCellIds->SetNumberOfIds(sideSetOnlyFaces[i].size());
        for (int j=0; j<sideSetOnlyFaces[i].size(); j++)
          vtkCellIds->SetId(j,sideSetOnlyFaces[i][j]);
        sideSet->InsertNextCell(9,vtkCellIds);
      } else if (sideSetOnlyFaces[i].size() > 4) { 
        //  add polygon face (#7)
        vtkSmartPointer<vtkIdList> vtkCellIds = vtkSmartPointer<vtkIdList>::New();
        vtkCellIds->SetNumberOfIds(sideSetOnlyFaces[i].size());
        for (int j=0; j<sideSetOnlyFaces[i].size(); j++)
          vtkCellIds->SetId(j,sideSetOnlyFaces[i][j]);
        sideSet->InsertNextCell(7,vtkCellIds);
      } else {
        Info << "Face type not supported yet!" << nl;
      }
    }

    auto sideSetEntities = vtkSmartPointer<vtkIntArray>::New();
    sideSetEntities->SetName("GeoEnt");

    auto sideSetOrigCellId = vtkSmartPointer<vtkIdTypeArray>::New();
    sideSetOrigCellId->SetName("OrigCellIds");

    auto sideSetCellFaceId = vtkSmartPointer<vtkIntArray>::New();
    sideSetCellFaceId->SetName("CellFaceIds");

    auto sideSetPatchId = vtkSmartPointer<vtkIntArray>::New();
    sideSetPatchId->SetName("PatchIds");

    auto sideSetPatchName = vtkSmartPointer<vtkStringArray>::New();
    sideSetPatchName->SetName("PatchNames");

    for (int i=0; i<sideSetOnlyFaces.size(); i++) {
      sideSetEntities->InsertNextValue(
          static_cast<int>(physGrps[sideSetOnlyOwners[i]]));
      sideSetOrigCellId->InsertNextValue(sideSetOnlyOwners[i]);
    }

    const cellShapeList& cellShapes = foamMesh->cellShapes();
    for (int i=0; i<sideSetOnlyOwners.size(); i++) {
      faceList fcs = cellShapes[sideSetOnlyOwners[i]].faces();
      for (int j=0; j<fcs.size(); j++) {
        if (fcs[j] == sideSetOnlyFaces[i])
          sideSetCellFaceId->InsertNextValue(j);
      }
    }

    for (int i=0; i<sideSetOnlyFaces.size(); i++) {
      sideSetPatchId->InsertNextValue(face_patch_map[sideSetOnlyFaceLabels[i]]);
      sideSetPatchName->InsertNextValue(patchNames[face_patch_map[sideSetOnlyFaceLabels[i]]]);
    }

    // Add arrays to sideSet
    sideSet->GetCellData()->AddArray(sideSetPatchId);
    sideSet->GetFieldData()->AddArray(sideSetPatchName);

    // Create sideSet struct
    auto sideSetStruct = SideSet(sideSet,sideSetEntities,sideSetOrigCellId,
                                 sideSetCellFaceId);

    std::string gmshMesh = "foamGeoMesh_" + nemAux::getRandomString(6);
    GmshInterface::Initialize();
    gmsh::model::add(gmshMesh);
    gmsh::model::setCurrent(gmshMesh);


    {  // Add geometric entities and physical groups
      std::set<std::pair<int, int>> dim_phyGrp;
      int idx = 0;
      auto phyGrpArray = vtkArrayDownCast<vtkDataArray>(
          vtkdataSet->GetCellData()->GetArray(phyGrpArrayName.c_str(),idx));
      vtkSmartPointer<vtkGenericCell> vtkGC =
          vtkSmartPointer<vtkGenericCell>::New();

      if (idx != -1) {
        // First sort
        for (vtkIdType i = 0; i < vtkdataSet->GetNumberOfCells(); ++i) {
          vtkGC->SetCellType(vtkdataSet->GetCellType(i));
          int dim = vtkGC->GetCellDimension();
          int phyGrp = static_cast<int>(phyGrpArray->GetComponent(i, 0));

          dim_phyGrp.insert({dim, phyGrp});
        }

        // then add. Each phyGrp gets its own geoEnt
        for (const auto &dp : dim_phyGrp) {
          gmsh::model::addDiscreteEntity(dp.first, dp.second);
          gmsh::model::addPhysicalGroup(dp.first, {dp.second}, dp.second);
        }
      }
    }

    return {vtkdataSet, gmshMesh, phyGrpArrayName, sideSetStruct};
  }
}

std::unique_ptr<Foam::fvMesh> foamGeoMesh::GM2foam(
      const GeoMesh &geoMesh,
      Foam::Time *runTime_,
      const std::string &phyGrpArrayName) {
  int numPoints = geoMesh.mesh->GetNumberOfPoints();
  int numCells = geoMesh.mesh->GetNumberOfCells();

  Foam::pointField pointData(numPoints);
  Foam::pointField pointData2(numPoints);
  Foam::cellShapeList cellShapeData(numCells);
  Foam::cellShapeList cellShapeData2(numCells);

  if (numPoints > 0) {
    // Fetch all points and coordinates
    std::vector<std::vector<double>> verts;
    verts.resize(numPoints);
    for (int ipt = 0; ipt < numPoints; ipt++) {
      std::vector<double> getPt = std::vector<double>(3);
      geoMesh.mesh->GetPoint(ipt, &getPt[0]);
      verts[ipt].resize(3);
      verts[ipt][0] = getPt[0];
      verts[ipt][1] = getPt[1];
      verts[ipt][2] = getPt[2];
    }

    // Gets Ids for cells
    std::vector<std::vector<int>> cellIds;
    cellIds.resize(numCells);

    // Gets celltypes for all cells in mesh
    std::vector<int> typeCell;
    typeCell.resize(numCells);

    for (int i=0; i<numPoints; i++) {
      pointData[i] = Foam::vector(verts[i][0],verts[i][1],verts[i][2]);
      pointData2[i] = Foam::vector(verts[i][0],verts[i][1],verts[i][2]);
    }

    for (int i=0; i<numCells; i++) {
      vtkIdList *ptIds = vtkIdList::New();
      geoMesh.mesh->GetCellPoints(i, ptIds);
      int numIds = ptIds->GetNumberOfIds();
      cellIds[i].resize(numIds);
      for (int j=0; j<numIds; j++)
        cellIds[i][j] = ptIds->GetId(j);
      Foam::labelList meshPoints(numIds);
      for (int k=0; k<numIds; k++)
        meshPoints[k] = cellIds[i][k];
      typeCell[i] = geoMesh.mesh->GetCellType(i);
      if (typeCell[i] == 12) {
        cellShapeData[i] = cellShape("hex", meshPoints, true);
        cellShapeData2[i] = cellShape("hex", meshPoints, true);
      }
      else if (typeCell[i] == 14) {
        cellShapeData[i] = cellShape("pyr", meshPoints, true);
        cellShapeData2[i] = cellShape("pyr", meshPoints, true);
      }
      else if (typeCell[i] == 10) {
        cellShapeData[i] = cellShape("tet", meshPoints, true);
        cellShapeData2[i] = cellShape("tet", meshPoints, true);
      }
      else {
        std::cerr << "Only Hexahedral, Tetrahedral," 
                  << " and Pyramid cells are supported in foamGeoMesh!"
                  << std::endl;
        throw;
      }
    }
  }

  Foam::word regName = Foam::fvMesh::defaultRegion;

  // Get patch names from mesh
  int idx1 = 0;
  int idx2 = 0;
  std::map<int,std::string> ptchNameIDMap;
  auto ptchIds = vtkIntArray::FastDownCast(
      geoMesh.sideSet.sides->GetCellData()->GetArray("PatchIds", idx1));
  auto ptchNms = vtkStringArray::SafeDownCast(
      geoMesh.sideSet.sides->GetFieldData()->GetAbstractArray("PatchNames",
                                                              idx2));
  Foam::faceListList bndryFaces(0);
  Foam::wordList BndryPatchNames(0);
  Foam::wordList BndryPatchTypes(0);
  Foam::wordList boundaryPatchPhysicalTypes(0);
  int totalPatches = 0;
  std::vector<int> startIds;
  std::vector<int> nFacesInPatch;
  Foam::PtrList<Foam::dictionary> dictList(0);
  Foam::PtrList<Foam::polyPatch> patchPtrList(0);
  Foam::cellList cellLst2(0);
  Foam::faceList faceLst2(0);

  if (idx1 != -1 && idx2 != -1) {
    // Patches exist
    for (int i = 0; i < ptchIds->GetNumberOfTuples(); ++i) {
      // Get patch Id
      auto cc = ptchIds->GetTypedComponent(i, 0);
      int id_p = cc;
      
      // Get patch name
      std::string nm = ptchNms->GetValue(i);
      
      // Create maps
      ptchNameIDMap[id_p] = nm;
    }

    totalPatches = ptchNameIDMap.size();
    std::cout << "Total Patches --> " << totalPatches << std::endl;
    std::vector<int> faceCounts(totalPatches,0);   // Counts faces for each patch
    std::vector<std::vector<int>> facePatchIds;
    facePatchIds.resize(totalPatches);
    for (int i = 0; i < ptchIds->GetNumberOfTuples(); ++i) {
      // Get patch Id
      auto cc = ptchIds->GetTypedComponent(i, 0);
      int id_p = cc;
      faceCounts[id_p] += 1;
      facePatchIds[id_p].push_back(i);
    }

    BndryPatchNames.setSize(totalPatches);
    BndryPatchTypes.setSize(totalPatches);

    auto iter = ptchNameIDMap.begin();
    int r_indx = 0;
    while (iter != ptchNameIDMap.end()) {
      BndryPatchNames[r_indx] = iter->second;
      BndryPatchTypes[r_indx] = "patch";
      iter++;
      r_indx++;
    }

    // Start pulling out cell (2D faces) point order from sideSets and create
    // faceList. Append those into faceListList at the end of every patch.
    bndryFaces.clear();
    for (int i=0; i<totalPatches; i++) {
      faceList thisPatch(0);
      thisPatch.clear();
      for (int j = 0; j < facePatchIds[i].size(); j++) {
        vtkIdList *ptIds = vtkIdList::New();
        geoMesh.sideSet.sides->GetCellPoints(facePatchIds[i][j], ptIds);
        int numIds = ptIds->GetNumberOfIds();
        labelList lstLbl(0);
        lstLbl.clear();
        for (int k = 0; k < numIds; k++) {
          lstLbl.append(ptIds->GetId(k));
        }
        if (!lstLbl.empty()) thisPatch.append(face(lstLbl));
      }
      bndryFaces.append(thisPatch);
      nFacesInPatch.push_back(thisPatch.size());
    }
    boundaryPatchPhysicalTypes.resize(totalPatches,polyPatch::typeName);
  }

  {
    Foam::polyMesh tmpfm(
      IOobject
      (
        regName,
        runTime_->timeName(),
        *runTime_,
        Foam::IOobject::NO_READ,
        Foam::IOobject::AUTO_WRITE
      ),
      std::move(pointData),        // Vertices
      cellShapeData,               // Cell shape and points
      bndryFaces,                  // Boundary faces
      BndryPatchNames,             // Boundary Patch Names
      BndryPatchTypes,             // Boundary Patch Types
      "defaultPatch",              // Default Patch Name
      Foam::polyPatch::typeName,   // Default Patch Type
      Foam::wordList()
    );

    const polyBoundaryMesh& patches2 = tmpfm.boundaryMesh();

    for (int i=0; i<patches2.size(); i++) {
      Foam::polyPatch *tmpPtch = new Foam::polyPatch(patches2[i]);
      patchPtrList.append(tmpPtch);
    }

    // Get pointField (points), faceList (faces), and cellList (cells).
    const Foam::cellList cellLst = tmpfm.cells();
    cellLst2 = cellLst;
    const Foam::faceList faceLst = tmpfm.faces();
    faceLst2 = faceLst;
  }

  auto fm = std::unique_ptr<Foam::fvMesh>(new Foam::fvMesh(
    IOobject
    (
      regName,
      runTime_->timeName(),
      *runTime_,
      Foam::IOobject::NO_READ,
      Foam::IOobject::AUTO_WRITE
    ),
    std::move(pointData2),
    std::move(faceLst2),
    std::move(cellLst2)
  ));

  // Add boundary patches to fvMesh using patchPtrList
  fm->addFvPatches(patchPtrList,true);

  // Get Physical Groups Vector from GeoMesh
  std::vector<double> allGroups;
  int idx;
  vtkSmartPointer<vtkDataArray> cd =
      geoMesh.mesh->GetCellData()->GetArray(phyGrpArrayName.c_str(), idx);
  if (idx != -1) {
    //allGroups.resize(cd->GetNumberOfTuples());
    double x[1];
    for (nemId_t i = 0; i < cd->GetNumberOfTuples(); ++i) {
      cd->GetTuple(i, x);
     allGroups.push_back(x[0]);
    }

    // Figure out number of physical groups
    std::vector<double> scratchVec = allGroups;
    std::sort(scratchVec.begin(), scratchVec.end());
    scratchVec.erase(unique(scratchVec.begin(), scratchVec.end()), scratchVec.end());
    int totalGroups = scratchVec.size();

    // Create cellZones
    label zoneI = 0;
    for (int i=0; i<totalGroups; i++) {
      Foam::labelList currentGroupList;
      for (int j=0; j < allGroups.size(); j++)
        if (int(allGroups[j]) == i)
          currentGroupList.append(j);

      const Foam::cellZoneMesh& cellZones = fm->cellZones();
      zoneI = cellZones.size();
      fm->cellZones().setSize(zoneI+1);
      Foam::word nameSet = "pg_"+std::to_string(i);
      auto clzn = new Foam::cellZone(nameSet,currentGroupList,zoneI,cellZones);
      fm->cellZones().set(zoneI,clzn);
    }
    fm->cellZones().writeOpt() = IOobject::AUTO_WRITE;
  }
  return fm;
}

void foamGeoMesh::write(const std::string &fileName) {
  // Create system directory
  Foam::word w = "system";
  if (!Foam::isDir(w))
    Foam::mkDir(w);

  // Write ControlDict
  Foam::fileName fcontrolDict_ = "system/controlDict";
  if (!Foam::exists(fcontrolDict_)) {
    Foam::OFstream outcontrolDict_(fcontrolDict_);
    IOobject::writeBanner(outcontrolDict_);
    controlDict_->write(outcontrolDict_,false);
  }

  // Write fvSchemes
  Foam::fileName ffvSchemes_ = "system/fvSchemes";
  if (!Foam::exists(ffvSchemes_)) {
    Foam::OFstream outfvSchemes_(ffvSchemes_);
    IOobject::writeBanner(outfvSchemes_);
    fvSchemes_->write(outfvSchemes_,false);
  }

  // Write fvSolution
  Foam::fileName ffvSolution_ = "system/fvSolution";
  if (!Foam::exists(ffvSolution_)) {
    Foam::OFstream outfvSolution_(ffvSolution_);
    IOobject::writeBanner(outfvSolution_);
    fvSolution_->write(outfvSolution_,false);
  }

  // Region to write in
  Foam::fileName newReg = "";
  if (fileName.empty())
    newReg = Foam::polyMesh::defaultRegion;
  else
    newReg = fileName;

  // Write dictionaries to the specific region too
  if(!Foam::isDir("system/"+newReg)) Foam::mkDir("system/"+newReg);

  Foam::fileName fcontrolDict_2 = "system/"+newReg+"/controlDict";
  if (!Foam::exists(fcontrolDict_2)) {
    Foam::cp(fcontrolDict_, fcontrolDict_2);
  }

  // Write fvSchemes
  Foam::fileName ffvSchemes_2 = "system/"+newReg+"/fvSchemes";
  if (!Foam::exists(ffvSchemes_2)) {
    Foam::cp(ffvSchemes_, ffvSchemes_2);
  }

  // Write fvSolution
  Foam::fileName ffvSolution_2 = "system/"+newReg+"/fvSolution";
  if (!Foam::exists(ffvSolution_2)) {
    Foam::cp(ffvSolution_, ffvSolution_2);
  }

  // Write mesh in fileName
  Foam::mkDir("0");

  // Write cellZones and cellSets
  const cellZoneMesh& czMesh = fmesh_->cellZones();
  label numZones = czMesh.size();

  fmesh_->setInstance("0");
  fmesh_->write();

  if (numZones > 0) {
    for (int i=0; i<czMesh.size(); i++) {
      const Foam::cellZone& cz = czMesh[i];
      Foam::cellSet(*fmesh_, cz.name(), cz).write();
    }
  }
  
  Foam::fileName currentDir = "0";
  Foam::fileName currentReg = fmesh_->dbDir();
  Foam::fileName newDir = "constant";
  Foam::fileName finalMesh = "";

  if (currentReg.empty()) {
    if (!Foam::isDir(currentDir+"/"+newReg))
      Foam::mkDir(currentDir+"/"+newReg);
    Foam::mv(currentDir + "/polyMesh", currentDir + "/" + newReg + "/polyMesh");
  } else {
    Foam::mv(currentDir+"/"+currentReg,currentDir+"/"+newReg);
  }
  finalMesh = currentDir + "/" + newReg;
  if (!Foam::isDir(newDir+"/"+newReg))
    Foam::mkDir(newDir+"/"+newReg);
  Foam::mv(finalMesh,newDir+"/"+newReg);
}

void foamGeoMesh::report(std::ostream &out) const { geoMeshBase::report(out); }

void foamGeoMesh::setFoamMesh(std::unique_ptr<Foam::fvMesh> foamMesh) {
  // Setting class parameters
  fmesh_ = std::unique_ptr<Foam::fvMesh>(foamMesh.get());

  // Parent
  setGeoMesh(foam2GM(foamMesh.get()));
}

void foamGeoMesh::resetNative() {
  // Clear
  controlDict_.reset();
  fvSchemes_.reset();
  fvSolution_.reset();
  runTime_.reset();
  fmesh_.reset();

  // Initialize required variables
  InitializeFoam();

  // Reset fvMesh
  fmesh_ = GM2foam(getGeoMesh(),runTime_.get());
}

void foamGeoMesh::InitializeFoam() {
  // Create controlDict
  controlDict_ = std::unique_ptr<Foam::dictionary>(
                                           new Foam::dictionary("controlDict"));
  Foam::dictionary fmfle("FoamFile");
  fmfle.add("version", "2.0");
  fmfle.add("format", "ascii");
  fmfle.add("class", "dictionary");
  fmfle.add("location", "\"system\"");
  fmfle.add("object", "controlDict");
  controlDict_->add("FoamFile",fmfle);
  controlDict_->add("deltaT",1);
  controlDict_->add("startTime",0);
  controlDict_->add("writeInterval",1);

  // Create fvSchemes
  fvSchemes_ = std::unique_ptr<Foam::dictionary>(
                                            new Foam::dictionary("fvSchemes"));
  Foam::dictionary gradSchemes("gradSchemes");
  Foam::dictionary divSchemes("divSchemes");
  Foam::dictionary laplacianSchemes("laplacianSchemes");
  gradSchemes.add("default","Gauss linear");
  gradSchemes.add("grad(p)","Gauss linear");
  divSchemes.add("default","none");
  divSchemes.add("div(phi,U)","Gauss linear");
  laplacianSchemes.add("default","none");
  laplacianSchemes.add("laplacian(nu,U)","Gauss linear corrected");
  laplacianSchemes.add("laplacian((1|A(U)),p)","Gauss linear corrected");
  Foam::dictionary fmfleFvscheme("FoamFile");
  fmfleFvscheme.add("version", "2.0");
  fmfleFvscheme.add("format", "ascii");
  fmfleFvscheme.add("class", "dictionary");
  fmfleFvscheme.add("location", "\"system\"");
  fmfleFvscheme.add("object", "fvSchemes");
  fvSchemes_->add("FoamFile",fmfleFvscheme);
  fvSchemes_->add("gradSchemes",gradSchemes);
  fvSchemes_->add("divSchemes",divSchemes);
  fvSchemes_->add("laplacianSchemes",laplacianSchemes);

  // Create fvSolution
  fvSolution_ = std::unique_ptr<Foam::dictionary>(
                                            new Foam::dictionary("fvSolution"));
  Foam::dictionary fmfleFvsol("FoamFile");
  fmfleFvsol.add("version", "2.0");
  fmfleFvsol.add("format", "ascii");
  fmfleFvsol.add("class", "dictionary");
  fmfleFvsol.add("location", "\"system\"");
  fmfleFvsol.add("object", "fvSolution");
  fvSolution_->add("FoamFile",fmfleFvsol);

  // Create time class without reading controlDict
  runTime_ = std::unique_ptr<Foam::Time>(new Foam::Time(controlDict_.get(), ".", "."));
  Foam::argList::noParallel();
}

} // namespace MSH
} // namespace NEM
