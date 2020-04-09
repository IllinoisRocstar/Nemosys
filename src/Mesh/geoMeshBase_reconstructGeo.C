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
#include <vtkDataSetRegionSurfaceFilter.h>
#include <vtkIdTypeArray.h>
#include <vtkIntArray.h>
#include <vtkPolyDataNormals.h>
#include <vtkStringArray.h>
#include <vtkThreshold.h>
#include <vtkTriangleFilter.h>
#include <vtkType.h>  // for vtkIdType

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

int geoMeshBase::vtkSurf2gmsh(vtkUnstructuredGrid *dataSet) {
  int tag0 = gmsh::model::addDiscreteEntity(0);
  int tag1 = gmsh::model::addDiscreteEntity(1);
  int tag2 = gmsh::model::addDiscreteEntity(2);
  int nodeOffset;
  {  // Get number of nodes already used.
    double nodeOffsetD;
    gmsh::option::getNumber("Mesh.NbNodes", nodeOffsetD);
    nodeOffset = static_cast<int>(nodeOffsetD);
  }

  {  // Add points
    std::vector<std::size_t> nodeTags;
    std::vector<double> coord;

    vtkSmartPointer<vtkPoints> dataSetPoints = dataSet->GetPoints();
    coord.resize(3 * dataSetPoints->GetNumberOfPoints());
    nodeTags.resize(dataSetPoints->GetNumberOfPoints());
    for (vtkIdType ptId = 0; ptId < dataSetPoints->GetNumberOfPoints();
         ++ptId) {
      nodeTags[ptId] = ptId + 1 + nodeOffset;
      dataSetPoints->GetPoint(ptId, &coord[3 * ptId]);
    }

    gmsh::model::mesh::addNodes(0, tag0, nodeTags, coord);
  }

  {  // Add elements
    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> nodeTags;

    vtkSmartPointer<vtkCellIterator> it = dataSet->NewCellIterator();
    for (it->InitTraversal(); !it->IsDoneWithTraversal(); it->GoToNextCell()) {
      auto etit = std::find(elementTypes.begin(), elementTypes.end(),
                            getGmshTypeFromVTKType(it->GetCellType()));
      std::size_t etnum = etit - elementTypes.begin();
      if (etit == elementTypes.end()) {
        elementTypes.emplace_back(getGmshTypeFromVTKType(it->GetCellType()));
        elementTags.emplace_back();
        nodeTags.emplace_back();
      }

      //      elementTags[etnum].emplace_back(it->GetCellId() + 1);

      vtkSmartPointer<vtkIdList> ptIds = it->GetPointIds();
      for (vtkIdType pt = 0; pt < it->GetNumberOfPoints(); ++pt) {
        nodeTags[etnum].emplace_back(ptIds->GetId(pt) + 1 + nodeOffset);
      }
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
        gmsh::model::mesh::addElementsByType(tag1, elementTypes[i],
                                             elementTags[i], nodeTags[i]);
      } else if (eDim == 2) {
        gmsh::model::mesh::addElementsByType(tag2, elementTypes[i],
                                             elementTags[i], nodeTags[i]);
      } else {
        std::cerr << "ERROR: 3D Elements received after skinning." << std::endl;
        throw;
      }
    }
  }

  return tag2;
}

void geoMeshBase::computeDiscreteGeoFromMsh(vtkUnstructuredGrid *dataSet,
                                            const std::string &geoName,
                                            double angleThreshold) {
  //****************************************************************************
  // 1. Access the geo array and count number of phys groups
  //****************************************************************************

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

  //****************************************************************************
  // 2. For each group, split into disjoint volumes.
  //****************************************************************************

  std::map<int, int> globalId2regionId;

  // Add each region as a material and a material property. These are
  // preserved by the interface filter and will be used to identify the
  // original phys group ids later on.
  vtkSmartPointer<vtkIntArray> materials = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkStringArray> originals =
      vtkSmartPointer<vtkStringArray>::New();
  materials->Initialize();
  originals->Initialize();

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
    std::cout << "Processing " << geoEntity << " of " << *geoEntities.end()
              << std::endl;

    tf->ThresholdBetween(geoEntity, geoEntity);

    // Run all the filters.
    cf->Update();
    vtkSmartPointer<vtkUnstructuredGrid> ug = cf->GetOutput();

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

  //****************************************************************************
  // 3. Run the interface filter (vtkDataSetRegionSurfaceFilter) to extract the
  //    surfaces as PolyData. Point ordering is not checked to guarantee
  //    consistent normals so pass through another filter (vtkPolyDataNormals)
  //    to fix normals.
  //****************************************************************************

  // The DSRSF filter acts only on 3D elements. Providing a mixed mesh treats
  // each dimension individually, causing duplicate MaterialIDs that do not
  // match across dims.
  vtkSmartPointer<vtkIntArray> cellDim = vtkSmartPointer<vtkIntArray>::New();
  cellDim->SetName("CellDim");
  cellDim->SetNumberOfValues(dataSet->GetNumberOfCells());
  vtkSmartPointer<vtkGenericCell> gc = vtkSmartPointer<vtkGenericCell>::New();
  for (vtkIdType i = 0; i < dataSet->GetNumberOfCells(); ++i) {
    dataSet->GetCell(i, gc);
    cellDim->SetTypedComponent(i, 0, gc->GetCellDimension());
  }
  dataSet->GetCellData()->AddArray(cellDim);

  tf->SetInputData(dataSet);
  tf->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS,
                             cellDim->GetName());
  tf->ThresholdBetween(3, 3);

  // Initialize interface filter
  vtkSmartPointer<vtkDataSetRegionSurfaceFilter> dsrsf =
      vtkSmartPointer<vtkDataSetRegionSurfaceFilter>::New();

  // Set-up interface filter
  materials->SetName(dsrsf->GetMaterialIDsName());
  originals->SetName(dsrsf->GetMaterialPropertiesName());
  dataSet->GetFieldData()->AddArray(materials);
  dataSet->GetFieldData()->AddArray(originals);
  dsrsf->SetRegionArrayName("RegionId");
  dsrsf->SetInputConnection(tf->GetOutputPort());

  /*
  // DEBUG
  vtkSmartPointer<vtkXMLDataSetWriter> writer =
      vtkSmartPointer<vtkXMLDataSetWriter>::New();
  writer->SetFileName("nucmesh_sample/nucmesh_1.vtu");
  dataSet->GetPointData()->SetActiveGlobalIds("GlobalIds");
  writer->SetInputData(dataSet);
  writer->Write();
  */

  // The DSRSF filter doesn't guarantee consistent normals
  vtkSmartPointer<vtkPolyDataNormals> pdn =
      vtkSmartPointer<vtkPolyDataNormals>::New();

  // Set-up normals filter
  pdn->SetInputConnection(dsrsf->GetOutputPort());
  pdn->ConsistencyOn();
  pdn->SplittingOff();

  // Gmsh geometry algorithm requires triangles.
  vtkSmartPointer<vtkTriangleFilter> trif =
      vtkSmartPointer<vtkTriangleFilter>::New();

  // Set-up triangulation filter.
  trif->SetInputConnection(pdn->GetOutputPort());

  // Run all the filters
  trif->Update();

  //****************************************************************************
  // 4. Loop through each interface created by DSRSF and add it as a gmsh
  //    surface.
  //****************************************************************************

  // Access field data produced by DSRSF to track what phys groups were on
  // each side of the interface.
  vtkSmartPointer<vtkIntArray> ancestors;

  materials = vtkIntArray::FastDownCast(
      trif->GetOutput()->GetFieldData()->GetAbstractArray(
          dsrsf->GetMaterialIDsName()));
  ancestors = vtkIntArray::FastDownCast(
      trif->GetOutput()->GetFieldData()->GetAbstractArray(
          dsrsf->GetMaterialPIDsName()));
  originals = vtkStringArray::SafeDownCast(
      trif->GetOutput()->GetFieldData()->GetAbstractArray(
          dsrsf->GetMaterialPropertiesName()));

  // Loop over all surfaces and add to gmsh.
  std::map<int, int> orig2new;
  std::map<int, std::vector<int>> new2gmshTags;

  for (vtkIdType i = 0; i < materials->GetNumberOfValues(); ++i) {
    int material = materials->GetTypedComponent(i, 0);
    int ancestor0 = ancestors->GetTypedComponent(i, 0);
    int ancestor1 = ancestors->GetTypedComponent(i, 1);

    vtkSmartPointer<vtkThreshold> t = vtkSmartPointer<vtkThreshold>::New();

    t->SetInputConnection(trif->GetOutputPort());
    t->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS,
                              dsrsf->GetRegionArrayName());
    t->ThresholdBetween(material, material);

    t->Update();
    int gmshTag = vtkSurf2gmsh(t->GetOutput());

    // Create map from original geo zones to vector of enclosing surfaces
    if (ancestor1 == -1) {
      int original = std::stoi(originals->GetValue(i));
      orig2new[original] = material;
      new2gmshTags[ancestor0].emplace_back(gmshTag);
    } else {
      new2gmshTags[ancestor0].emplace_back(gmshTag);
      new2gmshTags[ancestor1].emplace_back(gmshTag);
    }
  }

  // Remove duplicate nodes that are generated between nearby interfaces.
  gmsh::model::mesh::removeDuplicateNodes();

  //****************************************************************************
  // 5. Add the volumes and physical groups to gmsh. The maps created in
  //    previous steps are used here.
  //****************************************************************************

  // Loop over phys groups first
  for (const auto &geoEntity : geoEntities) {
    std::vector<int> volIds;

    // For each phys group, loop over volumes contained.
    for (const auto &region : geoId2regionIds[geoEntity]) {
      int surfLoop =
          gmsh::model::geo::addSurfaceLoop(new2gmshTags[orig2new[region]]);
      volIds.emplace_back(gmsh::model::geo::addVolume({surfLoop}));
    }

    gmsh::model::setPhysicalName(3, gmsh::model::addPhysicalGroup(3, volIds),
                                 std::to_string(geoEntity));
  }

  // Synchronize the added geo volumes
  gmsh::model::geo::synchronize();

  //****************************************************************************
  // 6. Compute the discrete geometry from the mesh and volume data.
  //****************************************************************************

  constexpr bool includeBoundary = false;
  // constexpr bool forReparameterization = true;

  // Color the surfaces of the mesh. Acts on triangles and quadrangles.
  gmsh::model::mesh::classifySurfaces(angleThreshold, includeBoundary, false);
  gmsh::model::mesh::classifySurfaces(angleThreshold, includeBoundary, true);
  // gmsh::model::mesh::classifySurfaces(angleThreshold, includeBoundary,
  //                                     forReparameterization);

  // For each color, compute parametrizations of the discrete surfaces.
  // Note: Acts ONLY on triangles!
  gmsh::model::mesh::createGeometry();
}

void geoMeshBase::reconstructGeo() {
  gmsh::model::setCurrent(_geoMesh.geo);
  gmsh::model::remove();
  gmsh::model::add(_geoMesh.geo);
  computeDiscreteGeoFromMsh(_geoMesh.mesh, _geoMesh.link, _angleThreshold);
}

}  // namespace MSH
}  // namespace NEM
