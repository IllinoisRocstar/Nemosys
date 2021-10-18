#if defined(_MSC_VER) && !defined(_USE_MATH_DEFINES)
#  define _USE_MATH_DEFINES
#endif

#include "Mesh/geoMeshBase.H"

#include <algorithm>
#include <iostream>
#include <map>
#include <numeric>
#include <set>

#include <vtkConnectivityFilter.h>
#include <vtkIdTypeArray.h>
#include <vtkIntArray.h>
#include <vtkMath.h>
#include <vtkPolyDataNormals.h>
#include <vtkThreshold.h>
#include <vtkDataSetTriangleFilter.h>
#include <vtkType.h>  // for vtkIdType

#include "AuxiliaryFunctions.H"
#include "MeshOperation/dataSetRegionBoundaryFilter.H"
#include "Mesh/gmshTypes.H"

#ifdef HAVE_GMSH
#include <gmsh.h>
#endif

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

// Setup for vtkConnectivityFilter
void addGlobalIds(vtkDataSet *dataSet) {
  vtkNew<vtkIdTypeArray> cellIds{};
  cellIds->SetNumberOfTuples(dataSet->GetNumberOfCells());
  for (vtkIdType i = 0; i < dataSet->GetNumberOfCells(); ++i) {
    cellIds->SetValue(i, i);
  }
  cellIds->SetName("GlobalIds");
  dataSet->GetCellData()->AddArray(cellIds);
  // These lines don't seem to affect vtkConnectivityFilter:
  // dataSet->GetCellData()->SetGlobalIds(cellIds);
  // dataSet->GetCellData()->SetCopyGlobalIds(true);
}

int local2GlobalRegions(vtkIdTypeArray *globalIds,
                        vtkIdTypeArray *localRegionId,
                        vtkIntArray *globalRegionId, int oldRegion,
                        int regionOffset, int numNewRegions,
                        std::map<int, std::vector<int>> &old2NewRegion) {
  for (vtkIdType i = 0; i < globalIds->GetNumberOfTuples(); ++i) {
    globalRegionId->SetTypedComponent(
        globalIds->GetTypedComponent(i, 0), 0,
        localRegionId->GetTypedComponent(i, 0) + regionOffset);
  }
  auto newRegionIter = old2NewRegion.emplace(
      std::piecewise_construct, std::forward_as_tuple(oldRegion),
      std::forward_as_tuple(numNewRegions));
  auto &newRegionVec = newRegionIter.first->second;
  std::iota(newRegionVec.begin(), newRegionVec.end(), regionOffset);
  return regionOffset + numNewRegions;
}

std::map<int, std::vector<int>> getConnectedEnts(vtkUnstructuredGrid *dataSet,
                                                 const std::string &geoName) {
  vtkSmartPointer<vtkDataArray> geoArray = vtkArrayDownCast<vtkDataArray>(
      dataSet->GetCellData()->GetAbstractArray(geoName.c_str()));

  // Count unique entities.
  std::set<int> geoEntities;
  if (geoArray) {
    for (vtkIdType i = 0; i < geoArray->GetNumberOfValues(); ++i) {
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

  vtkNew<vtkIntArray> connectedRegionId{};
  connectedRegionId->SetNumberOfValues(dataSet->GetNumberOfCells());

  int regionOffset = 1;
  std::map<int, std::vector<int>> geoId2regionIds;

  // Define Global Ids explicitly since vtkConnectivityFilter does not preserve
  // the built-in GlobalIds attribute.
  addGlobalIds(dataSet);

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
    auto componentId = vtkArrayDownCast<vtkIdTypeArray>(
        ug->GetCellData()->GetAbstractArray("RegionId"));
    auto globalId = vtkArrayDownCast<vtkIdTypeArray>(
        ug->GetCellData()->GetAbstractArray("GlobalIds"));

    // Create a map from phys group ids to vector of region ids it was split
    // into by the connectivity filter. This will be used to define physical
    // groups later.
    regionOffset = local2GlobalRegions(
        globalId, componentId, connectedRegionId, geoEntity, regionOffset,
        cf->GetNumberOfExtractedRegions(), geoId2regionIds);
  }

  connectedRegionId->SetName(geoName.c_str());
  dataSet->GetCellData()->AddArray(connectedRegionId);
  return geoId2regionIds;
}

// Mimic the workflow of (vtkPolyDataNormals with splitting on,
// vtkConnectivityFilter) for vtkDataSet with only lines
int colorLinesByAngle(vtkDataSet *lines, double angleThreshold) {
  vtkNew<vtkIdTypeArray> colorArr;
  colorArr->SetName("RegionId");
  colorArr->SetNumberOfTuples(lines->GetNumberOfCells());
  colorArr->FillTypedComponent(0, -1);
  vtkIdType nextRegion = 0;
  std::vector<vtkIdType> cellsToVisit;
  vtkNew<vtkIdList> neighborCells;
  vtkNew<vtkGenericCell> tempCell;
  for (vtkIdType i = 0; i < lines->GetNumberOfCells(); ++i) {
    if (colorArr->GetTypedComponent(i, 0) == -1) {
      cellsToVisit.clear();
      cellsToVisit.emplace_back(i);
      colorArr->SetTypedComponent(i, 0, nextRegion);
      while (!cellsToVisit.empty()) {
        auto cellIdx = cellsToVisit.back();
        cellsToVisit.pop_back();
        lines->GetCell(cellIdx, tempCell);
        auto points = tempCell->GetPointIds();
        assert(points->GetNumberOfIds() == 2);
        std::array<double, 3> vec{};
        {
          lines->GetPoint(points->GetId(1), vec.data());
          auto u = lines->GetPoint(points->GetId(0));
          std::transform(vec.begin(), vec.end(), u, vec.begin(),
                         std::minus<double>());
        }
        for (int j = 0; j < 2; ++j) {
          if (j == 1) {
            std::transform(vec.begin(), vec.end(), vec.begin(),
                           [](double x) { return -1 * x; });
          }
          auto sharedPointId = points->GetId(j);
          lines->GetPointCells(sharedPointId, neighborCells);
          for (vtkIdType k = 0; k < neighborCells->GetNumberOfIds(); ++k) {
            auto otherCellIdx = neighborCells->GetId(k);
            if (colorArr->GetTypedComponent(otherCellIdx, 0) == -1) {
              auto otherCellPoints =
                  lines->GetCell(otherCellIdx)->GetPointIds();
              auto otherPointId = otherCellPoints->GetId(0) == sharedPointId
                                      ? otherCellPoints->GetId(1)
                                      : otherCellPoints->GetId(0);
              std::array<double, 3> otherVec{};
              lines->GetPoint(otherPointId, otherVec.data());
              auto u = lines->GetPoint(sharedPointId);
              std::transform(otherVec.begin(), otherVec.end(), u,
                             otherVec.begin(), std::minus<double>());
              auto angle =
                  vtkMath::AngleBetweenVectors(vec.data(), otherVec.data());
              if (std::abs(angle - vtkMath::Pi()) < angleThreshold) {
                cellsToVisit.emplace_back(otherCellIdx);
                colorArr->SetTypedComponent(otherCellIdx, 0, nextRegion);
              }
            }
          }
        }
      }
      ++nextRegion;
    }
  }
  lines->GetCellData()->AddArray(colorArr);
  return nextRegion;
}

// Iterate through each of the surfaces (indexed from 1..numBoundEnts in the
// array given by surfaceArrayName), splitting by angleThreshold. Return map
// from old to new surfaces.
std::map<int, std::vector<int>> splitSurfaces(
    vtkPolyData *surfaces, const std::string &surfaceArrayName, int max_dim,
    double angleThreshold, int numBoundEnts, vtkIntArray *sideGeoEntArr) {
  std::map<int, std::vector<int>> old2newSurfs{};
  vtkNew<vtkThreshold> selectBoundary{};
  selectBoundary->SetInputData(surfaces);
  selectBoundary->SetInputArrayToProcess(0, 0, 0,
                                         vtkDataObject::FIELD_ASSOCIATION_CELLS,
                                         surfaceArrayName.c_str());

  vtkSmartPointer<vtkPolyDataNormals> pdn;
  vtkSmartPointer<vtkPolyData> pdnInput;
  vtkSmartPointer<vtkConnectivityFilter> pdConn{};
  if (max_dim == 3) {
    pdn = vtkSmartPointer<vtkPolyDataNormals>::New();
    // The custom boundary filter should leave normals consistent
    pdn->SetConsistency(false);
    pdn->SetSplitting(true);
    pdn->SetFeatureAngle(vtkMath::DegreesFromRadians(angleThreshold));
    pdnInput = vtkSmartPointer<vtkPolyData>::New();
    pdConn = vtkSmartPointer<vtkConnectivityFilter>::New();
    pdConn->SetExtractionModeToAllRegions();
    pdConn->ColorRegionsOn();
  }
  // Start from 1 because gmsh want strictly positive entity tags
  int regionOffset = 1;
  for (int i = 0; i < numBoundEnts; ++i) {
    vtkSmartPointer<vtkDataSet> coloredDS;
    int numRegions;
    selectBoundary->ThresholdBetween(i, i);
    if (max_dim == 2) {
      selectBoundary->Update();
      coloredDS = selectBoundary->GetOutput();
      numRegions = colorLinesByAngle(coloredDS, angleThreshold);
    } else {  // max_dim == 3
      selectBoundary->Update();
      auto thresholdOut = selectBoundary->GetOutput();
      pdnInput->SetPoints(thresholdOut->GetPoints());
      pdnInput->SetPolys(thresholdOut->GetCells());
      pdnInput->GetCellData()->ShallowCopy(thresholdOut->GetCellData());

      pdn->SetInputData(pdnInput);
      pdConn->SetInputConnection(pdn->GetOutputPort());
      pdConn->Update();
      coloredDS = pdConn->GetOutput();
      numRegions = pdConn->GetNumberOfExtractedRegions();
    }
    auto globalId = vtkArrayDownCast<vtkIdTypeArray>(
        coloredDS->GetCellData()->GetArray("GlobalIds"));
    auto splitSurfId = vtkArrayDownCast<vtkIdTypeArray>(
        coloredDS->GetCellData()->GetArray("RegionId"));
    regionOffset = local2GlobalRegions(globalId, splitSurfId, sideGeoEntArr, i,
                                       regionOffset, numRegions, old2newSurfs);
  }
  surfaces->GetCellData()->RemoveArray(surfaceArrayName.c_str());
  surfaces->GetCellData()->AddArray(sideGeoEntArr);
  return old2newSurfs;
}

#ifdef HAVE_GMSH
void constructGeoGmsh(vtkDataSet *dataSet, const int dim,
                      vtkIntArray *sideGeoEntArr,
                      const std::map<int, std::vector<int>> &vol2surf,
                      const std::map<int, std::vector<int>> &material2GeoEnt) {
  std::vector<double> coords(3 * dataSet->GetNumberOfPoints());
  for (vtkIdType i = 0; i < dataSet->GetNumberOfPoints(); ++i) {
    dataSet->GetPoint(i, coords.data() + 3 * i);
  }
  std::map<int, std::vector<std::size_t>> splitSurf2NodeTags;
  const int numPoints = dim;  // Assumes lines/triangles!
  for (vtkIdType i = 0; i < dataSet->GetNumberOfCells(); ++i) {
    auto &nodeTags = splitSurf2NodeTags[sideGeoEntArr->GetTypedComponent(i, 0)];
    auto pointsBegin = dataSet->GetCell(i)->GetPointIds()->GetPointer(0);
    std::transform(pointsBegin, pointsBegin + numPoints,
                   std::back_inserter(nodeTags),
                   [](vtkIdType x) { return x + 1; });
  }
  bool firstSurf = true;
  for (const auto &surf : splitSurf2NodeTags) {
    auto gmshTag = gmsh::model::addDiscreteEntity(dim - 1, surf.first);
    if (firstSurf) {
      // Assume (gmsh node tag) == (1 + vtk point id)
      gmsh::model::mesh::addNodes(dim - 1, gmshTag, {}, coords);
      firstSurf = false;
    }
    gmsh::model::mesh::addElementsByType(
        gmshTag, dim == 2 ? MSH_LIN_2 : MSH_TRI_3, {}, surf.second);
  }
  gmsh::model::mesh::reclassifyNodes();
  // Create curves and points to prepare for createGeometry
  gmsh::model::mesh::createTopology();
  gmsh::model::mesh::createGeometry();
  for (const auto &vol : vol2surf) {
    std::vector<int> postClassifySurfs;
    if (dim == 2) {
      auto curveLoop = gmsh::model::geo::addCurveLoop(vol.second);
      gmsh::model::geo::addPlaneSurface({curveLoop}, vol.first);
    } else {  // dim == 3
      auto surfLoop = gmsh::model::geo::addSurfaceLoop(vol.second);
      gmsh::model::geo::addVolume({surfLoop}, vol.first);
    }
  }
  gmsh::model::geo::synchronize();
  for (const auto &phyGrp : material2GeoEnt) {
    gmsh::model::addPhysicalGroup(dim, phyGrp.second, phyGrp.first);
  }
}
#endif

// Compute discrete geometry from dataSetCopy and name of geometry cell array
// geoName, with angleThreshold dihedral angle (in radians) to detect surfaces.
// Also sets sideSet with lower dimensional entities on boundaries (incl between
// different materials)
vtkSmartPointer<vtkPolyData> computeDiscreteGeoFromMsh(
    vtkUnstructuredGrid *dataSetOriginal, const std::string &geoName,
    double angleThreshold) {
  //****************************************************************************
  // 1. Access the geo array and count number of phys groups
  //****************************************************************************
  auto material2GeoEnt = getConnectedEnts(dataSetOriginal, geoName);

  // Prevent modifying the mesh any more
  vtkNew<vtkUnstructuredGrid> dataSet;
  dataSet->ShallowCopy(dataSetOriginal);

  //****************************************************************************
  // 2. Run the interface filter (dataSetRegionBoundaryFilter) to extract the
  //    surfaces as PolyData. For each resulting surface, split surfaces by
  //    angle.
  //****************************************************************************

  // The boundary filter should only be applied to one dimension.
  int max_dim = -1;
  for (vtkIdType i = 0; i < dataSet->GetNumberOfCells(); ++i) {
    max_dim = std::max(dataSet->GetCell(i)->GetCellDimension(), max_dim);
  }

  vtkNew<dataSetRegionBoundaryFilter> dsrbf;
  // Set-up interface filter
  dsrbf->SetMaterialArrayName(geoName);
  dsrbf->SetDimension(max_dim);
  dsrbf->SetInputData(dataSet);
  dsrbf->Update();

  vtkSmartPointer<vtkPolyData> sidePD = dsrbf->GetOutput();
  auto boundary2MeshRegion =
      vtkArrayDownCast<vtkIntArray>(sidePD->GetFieldData()->GetAbstractArray(
          dsrbf->GetRegionToMaterialArrayName().c_str()));
  addGlobalIds(sidePD);
  vtkNew<vtkIntArray> sideGeoEntArr;
  sideGeoEntArr->SetName("GeoEnt");
  sideGeoEntArr->SetNumberOfTuples(sidePD->GetNumberOfCells());
  int numBoundEnts = boundary2MeshRegion->GetNumberOfTuples();
  auto surf2ClassifiedSurfs =
      splitSurfaces(sidePD, dsrbf->GetBoundaryRegionArrayName(), max_dim,
                    angleThreshold, numBoundEnts, sideGeoEntArr);
  std::map<int, std::vector<int>> vol2surf;
  for (vtkIdType i = 0; i < boundary2MeshRegion->GetNumberOfTuples(); ++i) {
    auto &splitSurfs = surf2ClassifiedSurfs[i];
    for (vtkIdType j = 0; j < 2; ++j) {
      auto vol = boundary2MeshRegion->GetTypedComponent(i, j);
      if (vol != -1) {
        auto &boundary = vol2surf[vol];
        boundary.insert(boundary.end(), splitSurfs.begin(), splitSurfs.end());
      }
    }
  }
  sidePD->GetFieldData()->RemoveArray(boundary2MeshRegion->GetName());

  //****************************************************************************
  // 3. Pass to gmsh to reconstruct geometry.
  //****************************************************************************

#ifdef HAVE_GMSH
  if (max_dim == 2) {
    constructGeoGmsh(sidePD, max_dim, sideGeoEntArr, vol2surf, material2GeoEnt);
  } else {
    vtkNew<vtkDataSetTriangleFilter> trif;
    trif->SetInputData(sidePD);
    trif->Update();
    auto triangleSurf = trif->GetOutput();
    constructGeoGmsh(triangleSurf, max_dim,
                     vtkArrayDownCast<vtkIntArray>(
                         triangleSurf->GetCellData()->GetAbstractArray(
                             sideGeoEntArr->GetName())),
                     vol2surf, material2GeoEnt);
  }
#endif
  return sidePD;
}

}  // namespace

void geoMeshBase::reconstructGeo() {
#ifdef HAVE_GMSH
  if (_geoMesh.geo.empty()) {
    _geoMesh.geo = "geoMesh_" + nemAux::getRandomString(6);
  } else {
    gmsh::model::setCurrent(_geoMesh.geo);
    gmsh::model::remove();
  }
  gmsh::model::add(_geoMesh.geo);
  gmsh::model::setCurrent(_geoMesh.geo);
#endif
  if (_geoMesh.link.empty()) {
    _geoMesh.link = GEO_ENT_DEFAULT_NAME;
  }
  _geoMesh.sideSet = SideSet{
      computeDiscreteGeoFromMsh(_geoMesh.mesh, _geoMesh.link, _angleThreshold)};
  assert(_geoMesh.mesh->GetPoints() == _geoMesh.sideSet.sides->GetPoints());
}

}  // namespace MSH
}  // namespace NEM
