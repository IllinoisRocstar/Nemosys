#if defined(_MSC_VER) && !defined(_USE_MATH_DEFINES)
#  define _USE_MATH_DEFINES
#endif

#include "geoMeshBase.H"

#include <algorithm>
#include <iostream>
#include <map>
#include <numeric>
#include <set>

#include <gmsh.h>
#include <vtkCellIterator.h>
#include <vtkConnectivityFilter.h>
#include <vtkDataSetTriangleFilter.h>
#include <vtkIdTypeArray.h>
#include <vtkIntArray.h>
#include <vtkPolyDataNormals.h>
#include <vtkStringArray.h>
#include <vtkThreshold.h>
#include <vtkTriangleFilter.h>
#include <vtkType.h>  // for vtkIdType

#include "AuxiliaryFunctions.H"
#include "dataSetRegionSurfaceFilter.H"
#include "gmshTypes.H"

namespace NEM {
namespace MSH {

int geoMeshBase::getGmshTypeFromVTKType(int vtkType) {
  switch (vtkType) {
    case VTK_VERTEX: return MSH_PNT;

    case VTK_LINE: return MSH_LIN_2;
    case VTK_TRIANGLE: return MSH_TRI_3;
    case VTK_QUAD: return MSH_QUA_4;
    case VTK_TETRA: return MSH_TET_4;
    case VTK_HEXAHEDRON: return MSH_HEX_8;
    case VTK_WEDGE: return MSH_PRI_6;
    case VTK_PYRAMID: return MSH_PYR_5;

    case VTK_QUADRATIC_EDGE: return MSH_LIN_3;
    case VTK_QUADRATIC_TRIANGLE: return MSH_TRI_6;
    case VTK_QUADRATIC_QUAD: return MSH_QUA_8;
    case VTK_BIQUADRATIC_QUAD: return MSH_QUA_9;
    case VTK_QUADRATIC_TETRA: return MSH_TET_10;
    case VTK_QUADRATIC_HEXAHEDRON: return MSH_HEX_20;
    case VTK_TRIQUADRATIC_HEXAHEDRON: return MSH_HEX_27;

    default: return -1;
  }
}

// Helper functions that implement geoMeshBase::reconstructGeo
namespace {

// Convert vtk dataset into gmsh and add it to the current gmodel. Assumes
// gmsh::model::setCurrent has been called. Returns gmsh tag of highest
// dimensional discrete entity added. Adds all elements of one dimension on a
// single discrete entity. Supports only one 0D, 1D, and 2D element set at the
// time. Adding 3D elements will error. If gmshTag2cellData provided, will
// assume dataSet has "OrigCellIds", "CellFaceIds", and "TwinIds" CellData and
// add entries to the map for each element added.
int vtkSurf2gmsh(vtkUnstructuredGrid *dataSet, int elemOffset,
                 std::map<std::size_t, std::tuple<vtkIdType, int, vtkIdType>>
                     *gmshTag2cellData,
                 int tag1 = -1, int tag2 = -1) {
  // This is a dummy entity to add nodes onto.
  int tag0 = gmsh::model::addDiscreteEntity(0);
  // tag1 and tag2 are to hold the tag of discrete entities of dimension 1 and
  // 2, respectively. Delay adding discrete entities of dimension 1 and 2 until
  // we need to.
  int nodeOffset;
  {  // Get number of nodes already used.
    double numNodes;
    gmsh::option::getNumber("Mesh.NbNodes", numNodes);
    nodeOffset = static_cast<int>(numNodes) + 1;
  }

  {  // Add points
    std::vector<std::size_t> nodeTags;
    std::vector<double> coord;

    vtkSmartPointer<vtkPoints> dataSetPoints = dataSet->GetPoints();
    coord.resize(3 * dataSetPoints->GetNumberOfPoints());
    nodeTags.resize(dataSetPoints->GetNumberOfPoints());
    for (vtkIdType ptId = 0; ptId < dataSetPoints->GetNumberOfPoints();
         ++ptId) {
      nodeTags[ptId] = ptId + nodeOffset;
      dataSetPoints->GetPoint(ptId, &coord[3 * ptId]);
    }

    gmsh::model::mesh::addNodes(0, tag0, nodeTags, coord);
  }

  {  // Add elements
    auto origCellId = vtkIdTypeArray::FastDownCast(
        dataSet->GetCellData()->GetArray("GlobalIds"));
    auto origCellFaceId = vtkIntArray::FastDownCast(
        dataSet->GetCellData()->GetArray("CellFaceIds"));
    auto origTwinId = vtkIdTypeArray::FastDownCast(
        dataSet->GetCellData()->GetArray("TwinIds"));

    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> nodeTags;

    vtkSmartPointer<vtkCellIterator> it = dataSet->NewCellIterator();
    vtkIdType cell_idx = 0;
    for (it->InitTraversal(); !it->IsDoneWithTraversal(); it->GoToNextCell()) {
      auto etit =
          std::find(elementTypes.begin(), elementTypes.end(),
                    geoMeshBase::getGmshTypeFromVTKType(it->GetCellType()));
      std::size_t etnum = etit - elementTypes.begin();
      if (etit == elementTypes.end()) {
        elementTypes.emplace_back(
            geoMeshBase::getGmshTypeFromVTKType(it->GetCellType()));
        elementTags.emplace_back();
        nodeTags.emplace_back();
      }

      //      elementTags[etnum].emplace_back(it->GetCellId() + 1);

      vtkSmartPointer<vtkIdList> ptIds = it->GetPointIds();
      for (vtkIdType pt = 0; pt < it->GetNumberOfPoints(); ++pt) {
        nodeTags[etnum].emplace_back(ptIds->GetId(pt) + nodeOffset);
      }
      if (elemOffset > 0) {
        std::size_t elemTag = cell_idx + elemOffset;
        elementTags[etnum].emplace_back(elemTag);
        if (gmshTag2cellData) {
          (*gmshTag2cellData)[elemTag] =
              std::make_tuple(origCellId->GetValue(cell_idx),
                              origCellFaceId->GetValue(cell_idx),
                              origTwinId ? origTwinId->GetValue(cell_idx) : -1);
        }
      }  // else leave nodeTags empty and gmsh will auto assign
      ++cell_idx;
    }
    it->Delete();

    for (std::size_t i = 0; i < elementTypes.size(); ++i) {
      std::string elementName;
      int eDim, eOrder, numNodes, numPrimaryNodes;
      std::vector<double> nodeCoord;
      gmsh::model::mesh::getElementProperties(elementTypes[i], elementName,
                                              eDim, eOrder, numNodes, nodeCoord,
                                              numPrimaryNodes);

      if (eDim == 0) {
        gmsh::model::mesh::addElementsByType(tag0, elementTypes[i],
                                             elementTags[i], nodeTags[i]);
      } else if (eDim == 1) {
        if (tag1 < 0) {
          tag1 = gmsh::model::addDiscreteEntity(1);
        }
        gmsh::model::mesh::addElementsByType(tag1, elementTypes[i],
                                             elementTags[i], nodeTags[i]);
      } else if (eDim == 2) {
        if (tag2 < 0) {
          tag2 = gmsh::model::addDiscreteEntity(2);
        }
        gmsh::model::mesh::addElementsByType(tag2, elementTypes[i],
                                             elementTags[i], nodeTags[i]);
      } else {
        std::cerr << "ERROR: 3D Elements received after skinning." << std::endl;
        throw;
      }
    }
  }
  return tag2 > 0 ? tag2 : tag1;
}

std::set<int> cdgfmReadGeo(vtkUnstructuredGrid *dataSet,
                           const std::string &geoName) {
  vtkSmartPointer<vtkDataArray> geoArray = vtkArrayDownCast<vtkDataArray>(
      dataSet->GetCellData()->GetAbstractArray(geoName.c_str()));

  // Count unique entities.
  std::set<int> geoEntities;
  if (geoArray) {
    for (vtkIdType i = 0; i < geoArray->GetNumberOfValues(); ++i) {
      // TODO: Should use vtk data array dispatchers.
      //  https://vtk.org/Wiki/VTK/Tutorials/DataArrays#Best_Practices_for_vtkDataArray_Post-7.1
      geoEntities.insert(static_cast<int>(geoArray->GetComponent(i, 0)));
    }
  } else {
    std::cout << "No geometry array found." << std::endl;
    geoEntities.insert(1);

    vtkSmartPointer<vtkIntArray> newGeoArray =
        vtkSmartPointer<vtkIntArray>::New();
    newGeoArray->SetNumberOfValues(dataSet->GetNumberOfCells());
    newGeoArray->SetName(geoName.c_str());
    for (vtkIdType i = 0; i < dataSet->GetNumberOfCells(); ++i) {
      newGeoArray->SetTypedComponent(i, 0, 1);
    }
    dataSet->GetCellData()->AddArray(newGeoArray);
  }
  return geoEntities;
}

std::map<int, std::vector<int>> cdgfmAssignConnectedRegionId(
    vtkUnstructuredGrid *dataSet, const std::string &geoName,
    const std::set<int> &geoEntities, vtkIntArray *materials,
    vtkStringArray *originals) {
  std::map<int, int> globalId2regionId;

  int regionOffset = 1;
  std::map<int, std::vector<int>> geoId2regionIds;

  // Define Global Ids explicitly since vtkConnectivityFilter does not preserve
  // the built-in GlobalIds attribute.
  vtkSmartPointer<vtkIdTypeArray> cellIds =
      vtkSmartPointer<vtkIdTypeArray>::New();
  cellIds->SetNumberOfValues(dataSet->GetNumberOfCells());
  for (vtkIdType i = 0; i < dataSet->GetNumberOfCells(); ++i)
    cellIds->SetValue(i, i);
  cellIds->SetName("GlobalIds");
  dataSet->GetCellData()->AddArray(cellIds);
  dataSet->GetCellData()->RemoveArray("RegionId");
  // dataSet->GetCellData()->SetGlobalIds(cellIds);

  // Set-up threshold
  vtkSmartPointer<vtkThreshold> tf = vtkSmartPointer<vtkThreshold>::New();

  tf->SetInputData(dataSet);
  tf->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS,
                             geoName.c_str());

  // Extract connected regions from the mesh to pass to the geometry filter
  vtkSmartPointer<vtkConnectivityFilter> cf =
      vtkSmartPointer<vtkConnectivityFilter>::New();

  cf->SetInputConnection(tf->GetOutputPort());
  cf->SetExtractionModeToAllRegions();
  cf->ColorRegionsOn();

  // For each phys group, run threshold, and color connected regions,
  for (const auto &geoEntity : geoEntities) {
    std::cout << "Processing " << geoEntity << std::endl;

    tf->ThresholdBetween(geoEntity, geoEntity);

    // Run all the filters.
    cf->Update();

    // bypass API based on VTK version
    vtkSmartPointer<vtkUnstructuredGrid> ug =
        vtkUnstructuredGrid::SafeDownCast(cf->GetOutput());

    // Get the coloring array from the connectivity filter, (1) convert it to
    // vtkIntArray for interface filter and (2) add an offset to keep region
    // ids unique after merging back into single dataset.
    vtkSmartPointer<vtkIdTypeArray> regionId = vtkIdTypeArray::SafeDownCast(
        ug->GetCellData()->GetAbstractArray("RegionId"));
    vtkSmartPointer<vtkIdTypeArray> globalId = vtkIdTypeArray::SafeDownCast(
        ug->GetCellData()->GetAbstractArray("GlobalIds"));
    // ug->GetCellData()->GetGlobalIds());

    for (vtkIdType i = 0; i < ug->GetNumberOfCells(); ++i)
      globalId2regionId[globalId->GetTypedComponent(i, 0)] =
          regionId->GetTypedComponent(i, 0) + regionOffset;

    // Create a map from phys group ids to vector of region ids it was split
    // into by the connectivity filter. This will be used to define physical
    // groups later.
    std::vector<int> regionIds(cf->GetNumberOfExtractedRegions());
    std::iota(regionIds.begin(), regionIds.end(), regionOffset);
    geoId2regionIds[geoEntity] = regionIds;

    for (const auto &region : regionIds) {
      materials->InsertNextValue(region);
      originals->InsertNextValue(std::to_string(region));
    }

    // Advance the offset counter and merge the dataset
    regionOffset += cf->GetNumberOfExtractedRegions();
  }

  vtkSmartPointer<vtkIntArray> regionId = vtkSmartPointer<vtkIntArray>::New();
  regionId->SetName("RegionId");
  regionId->SetNumberOfValues(dataSet->GetNumberOfCells());
  for (vtkIdType i = 0; i < dataSet->GetNumberOfCells(); ++i)
    regionId->SetTypedComponent(i, 0, globalId2regionId[i]);
  dataSet->GetCellData()->AddArray(regionId);
  return geoId2regionIds;
}

void cdgfmAddBoundaryGmsh(
    vtkPolyDataAlgorithm *inputAlg, char *RegionArrayName,
    vtkIntArray *materials, std::map<int, std::vector<int>> &new2gmshTags,
    std::map<std::size_t, std::tuple<vtkIdType, int, vtkIdType>>
        &gmshElem2cellId) {
  // Loop over all surfaces and add to gmsh.
  int elemOffset = 1;  // Because mesh is empty, assign first element to tag 1

  for (vtkIdType i = 0; i < materials->GetNumberOfValues(); ++i) {
    auto material = materials->GetValue(i);
    vtkSmartPointer<vtkThreshold> t = vtkSmartPointer<vtkThreshold>::New();

    t->SetInputConnection(inputAlg->GetOutputPort());
    t->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS,
                              RegionArrayName);
    t->ThresholdBetween(material, material);

    t->Update();
    int gmshTag = vtkSurf2gmsh(t->GetOutput(), elemOffset, &gmshElem2cellId);
    elemOffset += t->GetOutput()->GetNumberOfCells();

    // Create map from original geo zones to vector of enclosing surfaces
    new2gmshTags[material].emplace_back(gmshTag);
  }
}

void cdgfmCreateSideSet(
    vtkUnstructuredGrid *dataSet, vtkPolyData *boundary,
    std::map<std::size_t, std::tuple<vtkIdType, int, vtkIdType>>
        &gmshElem2cellData,
    int max_dim, vtkPolyData *sideSet) {
  gmsh::vectorpair dimTags;
  gmsh::model::getEntities(dimTags, max_dim - 1);

  sideSet->SetPoints(dataSet->GetPoints());

  auto numSideSetCells = boundary->GetNumberOfCells();
  // Entity data needs to be added for lower dimensional cells.
  auto sideSetEntities = vtkSmartPointer<vtkIntArray>::New();
  sideSetEntities->SetName("GeoEnt");
  sideSetEntities->Allocate(numSideSetCells);
  auto sideSetOrigCellId = vtkSmartPointer<vtkIdTypeArray>::New();
  sideSetOrigCellId->SetName("OrigCellIds");
  sideSetOrigCellId->Allocate(numSideSetCells);
  auto sideSetFaceId = vtkSmartPointer<vtkIntArray>::New();
  sideSetFaceId->SetName("CellFaceIds");
  sideSetFaceId->Allocate(numSideSetCells);
  auto sideSetTwinId = vtkSmartPointer<vtkIdTypeArray>::New();
  sideSetTwinId->SetName("TwinIds");
  sideSetTwinId->Allocate(numSideSetCells);
  auto sideSetCellArr = vtkSmartPointer<vtkCellArray>::New();
  sideSetCellArr->Allocate(numSideSetCells);

  auto origCellIdArr = vtkIdTypeArray::FastDownCast(
      boundary->GetCellData()->GetArray("GlobalIds"));
  auto origCellFaceIdArr = vtkIntArray::FastDownCast(
      boundary->GetCellData()->GetArray("CellFaceIds"));

  // Get the index of the output of boundary filter from cell and cell face.
  // With this and gmshElem2cellId, can recover gmsh entity to vtk cell
  // relationship.
  std::map<std::pair<vtkIdType, int>, vtkIdType> cellFace2idx;
  // So we only add a cell once (in case it is represented by multiple tri
  // in gmsh mesh).
  std::set<std::pair<vtkIdType, int>> addedCellFace2ds;
  vtkIdType numCellsAdded = 0;
  // Parse by entity
  for (const auto &dimTag : dimTags) {
    // Get GMSH elements
    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> elemNodeTags;
    gmsh::model::mesh::getElements(elementTypes, elementTags, elemNodeTags,
                                   dimTag.first, dimTag.second);
    // Add VTK cell by GMSH element
    for (std::size_t iety = 0; iety < elementTypes.size(); ++iety) {
      // Get element type information from GMSH
      for (auto elem_tag : elementTags[iety]) {
        auto origCellId = std::get<0>(gmshElem2cellData[elem_tag]);
        auto cellFaceId = std::get<1>(gmshElem2cellData[elem_tag]);
        auto twin = std::get<2>(gmshElem2cellData[elem_tag]);
        if (addedCellFace2ds.find({origCellId, cellFaceId}) ==
            addedCellFace2ds.end()) {
          if (max_dim == 2) {
            auto idList = vtkSmartPointer<vtkIdList>::New();
            idList->Allocate(2);
            auto edge = dataSet->GetCell(origCellId)->GetEdge(cellFaceId);
            idList->InsertNextId(edge->GetPointId(0));
            idList->InsertNextId(edge->GetPointId(1));
            sideSetCellArr->InsertNextCell(idList);
          } else {  // max_dim == 3
            // TODO: Check higher order cells.
            sideSetCellArr->InsertNextCell(dataSet->GetCell(origCellId)
                                               ->GetFace(cellFaceId)
                                               ->GetPointIds());
          }
          // Add face to dataSet
          // geoEnt: Add to geometric entity array
          sideSetEntities->InsertNextValue(dimTag.second);
          sideSetOrigCellId->InsertNextValue(origCellId);
          sideSetFaceId->InsertNextValue(cellFaceId);
          ++numCellsAdded;
          sideSetTwinId->InsertNextValue(twin < 0 ? -1 : numCellsAdded);
          addedCellFace2ds.insert({origCellId, cellFaceId});
          if (twin >= 0) {
            auto twinCellId = origCellIdArr->GetValue(twin);
            auto twinFaceId = origCellFaceIdArr->GetValue(twin);
            if (addedCellFace2ds.find({twinCellId, twinFaceId}) ==
                addedCellFace2ds.end()) {
              if (max_dim == 2) {
                sideSetCellArr->InsertNextCell(dataSet->GetCell(twinCellId)
                                                   ->GetEdge(twinFaceId)
                                                   ->GetPointIds());
              } else {  // max_dim == 3
                sideSetCellArr->InsertNextCell(dataSet->GetCell(twinCellId)
                                                   ->GetFace(twinFaceId)
                                                   ->GetPointIds());
              }
              sideSetEntities->InsertNextValue(dimTag.second);
              sideSetOrigCellId->InsertNextValue(twinCellId);
              sideSetFaceId->InsertNextValue(twinFaceId);
              sideSetTwinId->InsertNextValue(numCellsAdded - 1);
              addedCellFace2ds.insert({twinCellId, twinFaceId});
              ++numCellsAdded;
            }
          }
        }
      }
    }
  }
  if (max_dim == 2) {
    sideSet->SetLines(sideSetCellArr);
  } else {  // max_dim == 3
    sideSet->SetPolys(sideSetCellArr);
  }
  sideSet->GetCellData()->AddArray(sideSetEntities);
  sideSet->GetCellData()->AddArray(sideSetOrigCellId);
  sideSet->GetCellData()->AddArray(sideSetFaceId);
  sideSet->GetCellData()->AddArray(sideSetTwinId);
}

// Compute discrete geometry from dataSet and name of geometry cell array
// geoName, with angleThreshold dihedral angle (in radians) to detect surfaces.
// Also sets sideSet with lower dimensional entities on boundaries (incl between
// different materials)
void computeDiscreteGeoFromMsh(vtkUnstructuredGrid *dataSet,
                               const std::string &geoName,
                               double angleThreshold, vtkPolyData *sideSet) {
  //****************************************************************************
  // 1. Access the geo array and count number of phys groups
  //****************************************************************************
  auto geoEntities = cdgfmReadGeo(dataSet, geoName);

  //****************************************************************************
  // 2. For each group, split into disjoint volumes.
  //****************************************************************************

  // Add each region as a material and a material property. These are
  // preserved by the interface filter and will be used to identify the
  // original phys group ids later on.
  vtkSmartPointer<vtkIntArray> materials = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkStringArray> originals =
      vtkSmartPointer<vtkStringArray>::New();
  materials->Initialize();
  originals->Initialize();

  auto geoId2regionIds = cdgfmAssignConnectedRegionId(
      dataSet, geoName, geoEntities, materials, originals);

  //****************************************************************************
  // 3. Run the interface filter (dataSetRegionSurfaceFilter) to extract the
  //    surfaces as PolyData. Point ordering is not checked to guarantee
  //    consistent normals so pass through another filter (vtkPolyDataNormals)
  //    to fix normals.
  //****************************************************************************

  // The DSRSF filter should only be applied to one dimension at a time.
  vtkSmartPointer<vtkIntArray> cellDim = vtkSmartPointer<vtkIntArray>::New();
  cellDim->SetName("CellDim");
  cellDim->SetNumberOfValues(dataSet->GetNumberOfCells());
  vtkSmartPointer<vtkGenericCell> gc = vtkSmartPointer<vtkGenericCell>::New();
  int max_dim = -1;
  for (vtkIdType i = 0; i < dataSet->GetNumberOfCells(); ++i) {
    dataSet->GetCell(i, gc);
    cellDim->SetTypedComponent(i, 0, gc->GetCellDimension());
    max_dim = std::max(gc->GetCellDimension(), max_dim);
  }
  dataSet->GetCellData()->AddArray(cellDim);
  auto tf = vtkSmartPointer<vtkThreshold>::New();
  tf->SetInputData(dataSet);
  tf->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS,
                             cellDim->GetName());
  tf->ThresholdBetween(max_dim, max_dim);
  tf->Update();

  vtkSmartPointer<dataSetRegionSurfaceFilter> dsrsf =
      vtkSmartPointer<dataSetRegionSurfaceFilter>::New();
  // Set-up interface filter
  materials->SetName(dsrsf->GetMaterialIDsName());
  originals->SetName(dsrsf->GetMaterialPropertiesName());
  dataSet->GetFieldData()->AddArray(materials);
  dataSet->GetFieldData()->AddArray(originals);
  dsrsf->SetRegionArrayName("RegionId");
  dsrsf->SetSingleSided(false);
  dsrsf->SetDimension(max_dim);
  dsrsf->SetInputConnection(tf->GetOutputPort());

  // Triangulated because gmsh can only classify surfaces of triangular meshes.
  // Redundant if max_dim == 2.
  vtkSmartPointer<vtkPolyDataAlgorithm> boundaryTriF;
  if (max_dim == 2) {
    dsrsf->Update();
    boundaryTriF = dsrsf;
  } else {  // max_dim == 3
    // The DSRSF filter doesn't guarantee consistent normals
    vtkSmartPointer<vtkPolyDataNormals> pdn =
        vtkSmartPointer<vtkPolyDataNormals>::New();

    // Set-up normals filter
    pdn->SetInputConnection(dsrsf->GetOutputPort());
    pdn->ConsistencyOn();
    pdn->SplittingOff();

    // Gmsh geometry algorithm requires triangles.
    auto trif = vtkSmartPointer<vtkTriangleFilter>::New();

    // Set-up triangulation filter.
    trif->SetInputConnection(pdn->GetOutputPort());

    // Run all the filters
    trif->Update();
    boundaryTriF = trif;
  }

  //****************************************************************************
  // 4. Loop through each interface created by DSRSF and add it as a gmsh
  //    surface/curve.
  //****************************************************************************

  // Loop over all surfaces and add to gmsh.
  std::map<int, std::vector<int>> new2gmshTags;
  // Map from gmsh element tag to (OrigCellIds, CellFaceIds, TwinIds) from
  // dataSet CellData; used to add the gmsh elements by entity later
  std::map<std::size_t, std::tuple<vtkIdType, int, vtkIdType>>
      gmshElem2cellData;
  cdgfmAddBoundaryGmsh(boundaryTriF, dsrsf->GetRegionArrayName(), materials,
                       new2gmshTags, gmshElem2cellData);

  //****************************************************************************
  // 5. Add the volumes and physical groups to gmsh. The maps created in
  //    previous steps are used here.
  //****************************************************************************

  std::map<int, int> dummyEntity2geoEntity;
  if (max_dim == 2) {
    // If max_dim == 2, we need to add surface elements to gmsh to run
    // classifySurfaces, which in turn requires only triangular elements.
    vtkSmartPointer<vtkDataSetTriangleFilter> trif =
        vtkSmartPointer<vtkDataSetTriangleFilter>::New();

    trif->SetInputConnection(tf->GetOutputPort());

    trif->Update();

    std::map<int, std::map<int, int>> region2surface;
    std::map<int, std::vector<int>> surfIds;
    for (const auto &geoEntity : geoEntities) {
      std::vector<int> surfIdVec;
      for (const auto &region : geoId2regionIds[geoEntity]) {
        int curveLoop = gmsh::model::geo::addCurveLoop(new2gmshTags[region]);
        int surface = gmsh::model::geo::addPlaneSurface({curveLoop});
        region2surface[geoEntity][region] = surface;
        surfIdVec.emplace_back(surface);
      }
      surfIds[geoEntity] = surfIdVec;
    }
    gmsh::model::geo::synchronize();
    for (const auto &geoEntity : geoEntities) {
      for (const auto &region : geoId2regionIds[geoEntity]) {
        auto surface = region2surface[geoEntity][region];
        auto regionThreshold = vtkSmartPointer<vtkThreshold>::New();
        regionThreshold->SetInputConnection(trif->GetOutputPort());
        regionThreshold->SetInputArrayToProcess(
            0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "RegionId");
        regionThreshold->ThresholdBetween(region, region);
        regionThreshold->Update();
        vtkSurf2gmsh(regionThreshold->GetOutput(), -1, nullptr, -1, surface);
      }
      // Because classifySurfaces changes the tags of 2d entities, adding
      // physical groups now is pointless: the new tags won't be assigned to the
      // correct physical groups. The workaround: create fake entities of
      // dimension 3 with the boundary being the surfaces that are supposed to
      // belong to a physical group. After classifySurfaces, read back the
      // boundary of these fake entities and assign the new 2d entities to the
      // appropriate physical group.
      auto dummyEntity =
          gmsh::model::addDiscreteEntity(3, geoEntity, surfIds[geoEntity]);
      dummyEntity2geoEntity[dummyEntity] = geoEntity;
    }
  } else {  // max_dim == 3
    // Loop over phys groups first
    for (const auto &geoEntity : geoEntities) {
      std::vector<int> volIds;

      // For each phys group, loop over volumes contained.
      for (const auto &region : geoId2regionIds[geoEntity]) {
        int surfLoop = gmsh::model::geo::addSurfaceLoop(new2gmshTags[region]);
        volIds.emplace_back(gmsh::model::geo::addVolume({surfLoop}));
      }

      gmsh::model::setPhysicalName(3, gmsh::model::addPhysicalGroup(3, volIds),
                                   std::to_string(geoEntity));
    }

    // Synchronize the added geo volumes
    gmsh::model::geo::synchronize();
  }

  // Remove duplicate nodes that are generated between nearby interfaces.
  gmsh::model::mesh::removeDuplicateNodes();

  //****************************************************************************
  // 6. Compute the discrete geometry from the mesh and volume data.
  //****************************************************************************

  constexpr bool includeBoundary = false;
  // constexpr bool forReparameterization = true;

  // Color the surfaces of the mesh. Acts on triangles and quadrangles.
  gmsh::model::mesh::classifySurfaces(angleThreshold, includeBoundary, false,
                                      angleThreshold);
  gmsh::model::mesh::classifySurfaces(angleThreshold, includeBoundary, true,
                                      angleThreshold);
  // gmsh::model::mesh::classifySurfaces(angleThreshold, includeBoundary,
  //                                     forReparameterization);

  //****************************************************************************
  // 7. Add lower dimensional elements to the sideSet. sideSet should have 4
  //    cell data: (1) "GeoEnt", which hold which geo entity the cell is on,
  //    (2) "OrigCellIds" and (3) "CellFaceIds", which hold which cell in
  //    dataSet the face (or in 2d case, edge) is from, and (4) "TwinIds",
  //    which, if positive, holds the index of the twin across an interface.
  //****************************************************************************

  cdgfmCreateSideSet(dataSet, dsrsf->GetOutput(), gmshElem2cellData, max_dim,
                     sideSet);

  // For each color, compute parametrizations of the discrete surfaces.
  // Note: Acts ONLY on triangles!
  gmsh::model::mesh::createGeometry();

  // In the case max_dim == 2, assign physical groups corresponding to
  // geoEntities based on the new tags.
  std::map<int, std::vector<int>> dummyEntity2boundary;
  if (max_dim == 2) {
    gmsh::vectorpair dimTags;
    gmsh::model::getEntities(dimTags, 3);
    for (auto dimTag : dimTags) {
      gmsh::vectorpair outDimTags;
      gmsh::model::getBoundary({dimTag}, outDimTags, false, false, false);
      std::vector<int> boundary;
      boundary.reserve(outDimTags.size());
      for (auto &outDimTag : outDimTags) {
        boundary.emplace_back(outDimTag.second);
      }
      dummyEntity2boundary[dimTag.second] = boundary;
    }
    for (const auto &boundary : dummyEntity2boundary) {
      gmsh::model::setPhysicalName(
          2, gmsh::model::addPhysicalGroup(2, boundary.second),
          std::to_string(dummyEntity2geoEntity[boundary.first]));
    }
    gmsh::model::removeEntities(dimTags);
  }
  // Clean-up
  gmsh::model::mesh::clear();
  dataSet->GetCellData()->RemoveArray(cellDim->GetName());
  dataSet->GetCellData()->RemoveArray("GlobalIds");
  dataSet->GetFieldData()->RemoveArray(dsrsf->GetMaterialIDsName());
  dataSet->GetFieldData()->RemoveArray(dsrsf->GetMaterialPropertiesName());
}

}  // namespace

void geoMeshBase::reconstructGeo() {
  if (_geoMesh.geo.empty()) {
    _geoMesh.geo = "geoMesh_" + nemAux::getRandomString(6);
  } else {
    gmsh::model::setCurrent(_geoMesh.geo);
    gmsh::model::remove();
  }
  if (_geoMesh.link.empty()) {
    _geoMesh.link = GEO_ENT_DEFAULT_NAME;
  }
  gmsh::model::add(_geoMesh.geo);
  gmsh::model::setCurrent(_geoMesh.geo);
  _geoMesh.sideSet = vtkSmartPointer<vtkPolyData>::New();
  computeDiscreteGeoFromMsh(_geoMesh.mesh, _geoMesh.link, _angleThreshold,
                            _geoMesh.sideSet);
}

}  // namespace MSH
}  // namespace NEM
