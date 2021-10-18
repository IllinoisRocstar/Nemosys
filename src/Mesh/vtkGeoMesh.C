#if defined(_MSC_VER) && !defined(_USE_MATH_DEFINES)
#  define _USE_MATH_DEFINES
#endif

#include "Mesh/vtkGeoMesh.H"

#include <iostream>
#include <set>
#include <utility>

#include <vtkGenericDataObjectReader.h>
#include <vtkGenericDataObjectWriter.h>
#include <vtkSTLReader.h>
#include <vtkXMLDataSetWriter.h>
#include <vtkXMLGenericDataObjectReader.h>

#include "AuxiliaryFunctions.H"

#ifdef HAVE_GMSH
#include <gmsh.h>
#endif

namespace NEM {
namespace MSH {

vtkStandardNewMacro(vtkGeoMesh)

vtkGeoMesh *vtkGeoMesh::Read(const std::string &fileName,
                             const std::string &geoArrayName) {
  vtkSmartPointer<vtkUnstructuredGrid> inUnstructuredGrid;

  std::string fileExt = nemAux::find_ext(fileName);
  if (fileExt == ".vtk") {  // Legacy type
    vtkSmartPointer<vtkGenericDataObjectReader> reader =
        vtkSmartPointer<vtkGenericDataObjectReader>::New();
    reader->SetFileName(fileName.c_str());
    reader->Update();

    inUnstructuredGrid = vtkUnstructuredGrid::SafeDownCast(reader->GetOutput());
  } else if (fileExt.substr(0, 3) == ".vt") {  // XML type
    vtkSmartPointer<vtkXMLGenericDataObjectReader> reader =
        vtkSmartPointer<vtkXMLGenericDataObjectReader>::New();
    reader->SetFileName(fileName.c_str());
    reader->Update();

    inUnstructuredGrid = vtkUnstructuredGrid::SafeDownCast(reader->GetOutput());
    /*
  } else if (fileExt == ".stl") {  // Allow STLs
    vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName(fileName.c_str());
    reader->Update();

    // TODO: This cast fails. STL reader outputs vtkPolyData!
    inUnstructuredGrid = vtkUnstructuredGrid::SafeDownCast(reader->GetOutput());
    */
  } else {  // Not a VTK type
    std::cerr << "ERROR: File extension " << fileExt
              << " is not supported by vtkGeoMesh." << std::endl;
    exit(1);
  }

  if (!inUnstructuredGrid) {  // Not a mesh
    std::cerr << "ERROR: Only VTK Unstructured Grid (vtu) is supported."
              << std::endl;
    exit(1);
  }

  return new vtkGeoMesh(inUnstructuredGrid, geoArrayName);
}

vtkGeoMesh::vtkGeoMesh()
    : vtkGeoMesh(vtkSmartPointer<vtkUnstructuredGrid>::New()) {}

vtkGeoMesh::vtkGeoMesh(vtkUnstructuredGrid *inUnstructuredGrid,
                       const std::string &geoArrayName)
    : geoMeshBase(vtk2GM(inUnstructuredGrid, geoArrayName)) {
  std::cout << "vtkGeoMesh constructed" << std::endl;
}

// vtkGeoMesh::vtkGeoMesh(geoMeshBase *inGmb) : geoMeshBase(inGmb),
//                                            _vtkMesh() {
//}

vtkGeoMesh::~vtkGeoMesh() { std::cout << "vtkGeoMesh destructed" << std::endl; }

void vtkGeoMesh::write(const std::string &fileName) {
  auto mesh = getGeoMesh().mesh;

  std::string fileExt = nemAux::find_ext(fileName);
  if (fileExt == ".vtk") {  // Legacy type
    vtkSmartPointer<vtkGenericDataObjectWriter> writer =
        vtkSmartPointer<vtkGenericDataObjectWriter>::New();
    writer->SetFileName(fileName.c_str());
    writer->SetInputData(mesh);
    writer->Write();
  } else {  // XML type
    vtkSmartPointer<vtkXMLDataSetWriter> writer =
        vtkSmartPointer<vtkXMLDataSetWriter>::New();
    writer->SetFileName(fileName.c_str());
    writer->SetInputData(mesh);
    writer->Write();
  }
}

void vtkGeoMesh::report(std::ostream &out) const { geoMeshBase::report(out); }

void vtkGeoMesh::getVtkMesh(vtkUnstructuredGrid *dest) {
  dest->DeepCopy(getGeoMesh().mesh);
}

void vtkGeoMesh::setVtkMesh(vtkUnstructuredGrid *vtkMesh) {
  setGeoMesh(vtk2GM(vtkMesh));
}

geoMeshBase::GeoMesh vtkGeoMesh::vtk2GM(vtkUnstructuredGrid *vtkMesh,
                                        const std::string &phyGrpArrayName) {
  std::string gmshMesh = "vtkGeoMesh_" + nemAux::getRandomString(6);
#ifdef HAVE_GMSH
  GmshInterface::Initialize();
  gmsh::model::add(gmshMesh);
  gmsh::model::setCurrent(gmshMesh);
#endif

  if (!vtkMesh ||                                                 // No mesh
      vtkMesh->GetNumberOfPoints() == 0 ||                        // No points
      !vtkMesh->GetCellData()->HasArray(phyGrpArrayName.c_str())  // No geometry
  )
    return {vtkMesh, "", "", {}};

  {  // Add geometric entities and physical groups
    std::set<std::pair<int, int>> dim_phyGrp;
    auto phyGrpArray = vtkArrayDownCast<vtkDataArray>(
        vtkMesh->GetCellData()->GetArray(phyGrpArrayName.c_str()));
    vtkSmartPointer<vtkGenericCell> vtkGC =
        vtkSmartPointer<vtkGenericCell>::New();

    // First sort
    for (vtkIdType i = 0; i < vtkMesh->GetNumberOfCells(); ++i) {
      vtkGC->SetCellType(vtkMesh->GetCellType(i));
      int dim = vtkGC->GetCellDimension();
      int phyGrp = static_cast<int>(phyGrpArray->GetComponent(i, 0));

      dim_phyGrp.insert({dim, phyGrp});
    }

#ifdef HAVE_GMSH
    // then add. Each phyGrp gets its own geoEnt
    for (const auto &dp : dim_phyGrp) {
      gmsh::model::addDiscreteEntity(dp.first, dp.second);
      gmsh::model::addPhysicalGroup(dp.first, {dp.second}, dp.second);
    }
#endif
  }

  /*
  {  // DEBUG
    gmsh::vectorpair dimTags;
    gmsh::model::getEntities(dimTags);
    for (const auto &dimTag : dimTags) {
      std::cout << "geoEnt   dim: " << dimTag.first
                << "  tag: " << dimTag.second << std::endl;
    }

    gmsh::model::getPhysicalGroups(dimTags);
    for (const auto &dimTag : dimTags) {
      std::cout << "phyGrp   dim: " << dimTag.first
                << "  tag: " << dimTag.second << std::endl;
    }
  }
  */

  return {vtkMesh, gmshMesh, phyGrpArrayName, {}};
}

void vtkGeoMesh::resetNative() {}

}  // namespace MSH
}  // namespace NEM
