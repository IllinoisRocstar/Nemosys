#include "oshGeoMesh.H"

#include <array>

#include <vtkCellTypes.h>
#include <vtkCharArray.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>
#include <vtkIntArray.h>
#include <vtkLongArray.h>
#include <vtkLongLongArray.h>
#include <vtkShortArray.h>
#include <vtkSignedCharArray.h>
#include <vtkTypeFloat64Array.h>
#include <vtkTypeInt32Array.h>
#include <vtkTypeInt64Array.h>
#include <vtkTypeInt8Array.h>
#include <Omega_h_build.hpp>
#include <Omega_h_element.hpp>
#include <Omega_h_file.hpp>

namespace NEM {
namespace MSH {

oshGeoMesh *oshGeoMesh::Read(const std::string &fileName) {
  auto lib = std::make_shared<Omega_h::Library>();

  return oshGeoMesh::Read(fileName, lib.get());
}

oshGeoMesh *oshGeoMesh::Read(const std::string &fileName,
                             Omega_h::Library *lib) {
  Omega_h::Mesh oshMesh = Omega_h::read_mesh_file(fileName, lib->world());

  return new oshGeoMesh(&oshMesh, lib);
}

oshGeoMesh::oshGeoMesh()
    : oshGeoMesh(std::make_shared<Omega_h::Mesh>().get()) {}

oshGeoMesh::oshGeoMesh(Omega_h::Mesh *oshMesh)
    : oshGeoMesh(oshMesh, (oshMesh->library()
                               ? oshMesh->library()
                               : std::make_shared<Omega_h::Library>().get())) {}

oshGeoMesh::oshGeoMesh(Omega_h::Mesh *oshMesh, Omega_h::Library *lib)
    : geoMeshBase({osh2vtk(oshMesh), "", ""}), _oshLibrary(*lib), _oshMesh() {
  _oshMesh = *oshMesh;
  _oshMesh.set_library(&_oshLibrary);
  std::cout << "oshGeoMesh constructed" << std::endl;
}

oshGeoMesh::~oshGeoMesh() { std::cout << "oshGeoMesh destructed" << std::endl; }

void oshGeoMesh::write(const std::string &fileName) {}
void oshGeoMesh::report(std::ostream &out) const {}

VTKCellType getVTKTypeFromOmega_hFamilyDim(Omega_h_Family family,
                                           Omega_h::Int dim) {
  return (
      family == OMEGA_H_SIMPLEX
          ? (dim == 3 ? VTK_TETRA
                      : (dim == 2 ? VTK_TRIANGLE
                                  : (dim == 1 ? VTK_LINE
                                              : (dim == 0 ? VTK_VERTEX
                                                          : VTK_EMPTY_CELL))))
          : (dim == 3 ? VTK_HEXAHEDRON
                      : (dim == 2 ? VTK_QUAD
                                  : (dim == 1 ? VTK_LINE
                                              : (dim == 0 ? VTK_VERTEX
                                                          : VTK_EMPTY_CELL)))));
}

Omega_h_Family getOmega_hFamilyFromVTKType(VTKCellType vtkCellType) {
  switch (vtkCellType) {
    case VTK_EMPTY_CELL:
    case VTK_VERTEX:
    case VTK_LINE:
    case VTK_TRIANGLE:
    case VTK_TETRA: return OMEGA_H_SIMPLEX;

    case VTK_QUAD:
    case VTK_HEXAHEDRON: return OMEGA_H_HYPERCUBE;

    default: {
      std::cerr << "Omega_h supports only linear simplices and hypercubes. "
                   "VTKCellType "
                << vtkCellType << " is not supported." << std::endl;
      exit(1);
    }
  }
}

Omega_h::Int getOmega_hDimFromVTKType(VTKCellType vtkCellType) {
  switch (vtkCellType) {
    case VTK_TETRA:
    case VTK_HEXAHEDRON: return 3;

    case VTK_TRIANGLE:
    case VTK_QUAD: return 2;

    case VTK_LINE: return 1;

    case VTK_VERTEX:
    case VTK_EMPTY_CELL:
    default: {
      std::cerr << "Omega_h supports only linear simplices and hypercubes. "
                   "VTKCellType "
                << vtkCellType << " is not supported." << std::endl;
      exit(1);
    }
  }
}

template <typename OT, typename VT>
void Otag2Varray(const Omega_h::TagBase *tagBase, VT *vtkArray) {
  const Omega_h::Tag<OT> *tag = Omega_h::as<OT>(tagBase);
  Omega_h::Int ncomps = tagBase->ncomps();
  Omega_h::HostRead<OT> tagArray(tag->array());

  vtkArray->SetName(tagBase->name().c_str());
  vtkArray->SetNumberOfComponents(ncomps);
  vtkArray->SetNumberOfValues(tagArray.size());

  for (Omega_h::LO i = 0; i < tagArray.size(); ++i)
    vtkArray->SetValue(i, tagArray.get(i));
}

void getVtkDataArrayFromOmega_hTag(const Omega_h::TagBase *tagBase,
                                   vtkSmartPointer<vtkDataArray> &vtkArray) {
  switch (tagBase->type()) {
    case OMEGA_H_I8: {
      vtkSmartPointer<vtkTypeInt8Array> vtkArray_tmp =
          vtkSmartPointer<vtkTypeInt8Array>::New();
      Otag2Varray<Omega_h::I8, vtkTypeInt8Array>(tagBase, vtkArray_tmp);
      vtkArray = vtkArray_tmp;
      break;
    }
    case OMEGA_H_I32: {
      vtkSmartPointer<vtkTypeInt32Array> vtkArray_tmp =
          vtkSmartPointer<vtkTypeInt32Array>::New();
      Otag2Varray<Omega_h::I32, vtkTypeInt32Array>(tagBase, vtkArray_tmp);
      vtkArray = vtkArray_tmp;
      break;
    }
    case OMEGA_H_I64: {
      vtkSmartPointer<vtkTypeInt64Array> vtkArray_tmp =
          vtkSmartPointer<vtkTypeInt64Array>::New();
      Otag2Varray<Omega_h::I64, vtkTypeInt64Array>(tagBase, vtkArray_tmp);
      vtkArray = vtkArray_tmp;
      break;
    }
    case OMEGA_H_F64: {
      vtkSmartPointer<vtkTypeFloat64Array> vtkArray_tmp =
          vtkSmartPointer<vtkTypeFloat64Array>::New();
      Otag2Varray<Omega_h::Real, vtkTypeFloat64Array>(tagBase, vtkArray_tmp);
      vtkArray = vtkArray_tmp;
      break;
    }
  }
}

vtkUnstructuredGrid *oshGeoMesh::osh2vtk(Omega_h::Mesh *oshMesh) {
  vtkUnstructuredGrid *vtkMesh = vtkUnstructuredGrid::New();

  if (!oshMesh->is_valid()) return vtkMesh;

  // Omega_h dimension is important to read its data.
  //
  // It stores 2D mesh without a z-coordinate while VTK always requires 3D. A
  // "2D mesh" in VTK will have z-coordinates set to zero.
  //
  // Also, an Omega_h mesh contains the main element and all its lower-dimension
  // constituents. E.g., a 3D mesh containing tetrahedrons will also contain the
  // four triangles, six edges, and four vertices. When converting to VTK, we
  // will discard the lower-dimensional elements.
  Omega_h::Int dim = oshMesh->dim();

  // An Omega_h mesh cannot have mixed elements. The family determines if the
  // mesh is composed of simplices (tetrahedrons and triangles) or hypercubes
  // (hexahedrons and quadrangles).
  Omega_h_Family family = oshMesh->family();
  VTKCellType vtk_type = getVTKTypeFromOmega_hFamilyDim(family, dim);

  {  // Add points
    vtkSmartPointer<vtkPoints> points =
        vtkSmartPointer<vtkPoints>::Take(vtkPoints::New(VTK_DOUBLE));
    points->Allocate(oshMesh->nverts());

    Omega_h::HostRead<Omega_h::Real> h_coord(oshMesh->coords());

    for (Omega_h::LO i = 0; i < oshMesh->nverts(); ++i) {
      points->InsertNextPoint(h_coord[i * dim + 0], h_coord[i * dim + 1],
                              dim == 2 ? 0 : h_coord[i * dim + 2]);
    }

    vtkMesh->SetPoints(points);
  }

  {  // Add cells
    Omega_h::Int deg = Omega_h::element_degree(family, dim, OMEGA_H_VERT);
    Omega_h::HostRead<Omega_h::LO> h_ev2v(oshMesh->ask_elem_verts());

    vtkMesh->Allocate(oshMesh->nelems());

    for (Omega_h::LO nelem = 0; nelem < oshMesh->nelems(); ++nelem) {
      vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
      ptIds->Allocate(deg);

      for (Omega_h::Int nvert = 0; nvert < deg; ++nvert)
        ptIds->InsertNextId(h_ev2v[nelem * deg + nvert]);

      vtkMesh->InsertNextCell(vtk_type, ptIds);
    }
  }

  {  // Add point data
    for (Omega_h::Int itag = 0; itag < oshMesh->ntags(Omega_h::VERT); ++itag) {
      vtkSmartPointer<vtkDataArray> vtkArray;

      getVtkDataArrayFromOmega_hTag(oshMesh->get_tag(Omega_h::VERT, itag),
                                    vtkArray);

      vtkMesh->GetPointData()->AddArray(vtkArray);
    }
  }

  {  // Add cell data
    for (Omega_h::Int itag = 0; itag < oshMesh->ntags(dim); ++itag) {
      vtkSmartPointer<vtkDataArray> vtkArray;

      getVtkDataArrayFromOmega_hTag(oshMesh->get_tag(dim, itag), vtkArray);

      vtkMesh->GetCellData()->AddArray(vtkArray);
    }
  }

  return vtkMesh;
}

template <typename VT, typename OT>
void Varray2Otag(Omega_h::Mesh *oshMesh, VT *vtkArray, Omega_h::Int oshDim) {
  Omega_h::HostWrite<OT> h_oshArray(vtkArray->GetNumberOfValues(),
                                    vtkArray->GetClassName());

  for (vtkIdType i = 0; i < vtkArray->GetNumberOfValues(); ++i)
    h_oshArray.set(i, vtkArray->GetValue(i));

  Omega_h::Read<OT> oshArray(h_oshArray);

  oshMesh->add_tag(oshDim, vtkArray->GetName(),
                   vtkArray->GetNumberOfComponents(), oshArray);
}

template <typename VT>
void Varray2Otag2(Omega_h::Mesh *oshMesh, VT *vtkArray, Omega_h::Int oshDim) {
  Varray2Otag<
      VT, typename std::conditional<
              sizeof(typename VT::ValueType) <= 1, Omega_h::I8,
              typename std::conditional<
                  sizeof(typename VT::ValueType) <= 4, Omega_h::I32,
                  typename std::enable_if<sizeof(typename VT::ValueType) <= 8,
                                          Omega_h::I64>::type>::type>::type>(
      oshMesh, vtkArray, oshDim);
}

void getOmega_hArrayFromVtkDataArray(Omega_h::Mesh *oshMesh,
                                     vtkSmartPointer<vtkDataArray> &vtkArray,
                                     Omega_h::Int oshDim) {
  switch (vtkArray->GetDataType()) {
    case VTK_FLOAT: {
      vtkSmartPointer<vtkFloatArray> vtkArray_tmp =
          vtkFloatArray::FastDownCast(vtkArray);
      Varray2Otag<vtkFloatArray, Omega_h::Real>(oshMesh, vtkArray_tmp, oshDim);
      vtkArray = vtkArray_tmp;
      break;
    }
    case VTK_DOUBLE: {
      vtkSmartPointer<vtkDoubleArray> vtkArray_tmp =
          vtkDoubleArray::FastDownCast(vtkArray);
      Varray2Otag<vtkDoubleArray, Omega_h::Real>(oshMesh, vtkArray_tmp, oshDim);
      vtkArray = vtkArray_tmp;
      break;
    }
    case VTK_CHAR: {
      vtkSmartPointer<vtkCharArray> vtkArray_tmp =
          vtkCharArray::FastDownCast(vtkArray);
      Varray2Otag2<vtkCharArray>(oshMesh, vtkArray_tmp, oshDim);
      vtkArray = vtkArray_tmp;
      break;
    }
    case VTK_SIGNED_CHAR: {
      vtkSmartPointer<vtkSignedCharArray> vtkArray_tmp =
          vtkSignedCharArray::FastDownCast(vtkArray);
      Varray2Otag2<vtkSignedCharArray>(oshMesh, vtkArray_tmp, oshDim);
      vtkArray = vtkArray_tmp;
      break;
    }
    case VTK_SHORT: {
      vtkSmartPointer<vtkShortArray> vtkArray_tmp =
          vtkShortArray::FastDownCast(vtkArray);
      Varray2Otag2<vtkShortArray>(oshMesh, vtkArray_tmp, oshDim);
      vtkArray = vtkArray_tmp;
      break;
    }
    case VTK_INT: {
      vtkSmartPointer<vtkIntArray> vtkArray_tmp =
          vtkIntArray::FastDownCast(vtkArray);
      Varray2Otag2<vtkIntArray>(oshMesh, vtkArray_tmp, oshDim);
      vtkArray = vtkArray_tmp;
      break;
    }
    case VTK_LONG: {
      vtkSmartPointer<vtkLongArray> vtkArray_tmp =
          vtkLongArray::FastDownCast(vtkArray);
      Varray2Otag2<vtkLongArray>(oshMesh, vtkArray_tmp, oshDim);
      vtkArray = vtkArray_tmp;
      break;
    }
    case VTK_ID_TYPE: {
      vtkSmartPointer<vtkIdTypeArray> vtkArray_tmp =
          vtkIdTypeArray::FastDownCast(vtkArray);
      Varray2Otag2<vtkIdTypeArray>(oshMesh, vtkArray_tmp, oshDim);
      vtkArray = vtkArray_tmp;
      break;
    }
    case VTK_LONG_LONG: {
      vtkSmartPointer<vtkLongLongArray> vtkArray_tmp =
          vtkLongLongArray::FastDownCast(vtkArray);
      Varray2Otag2<vtkLongLongArray>(oshMesh, vtkArray_tmp, oshDim);
      vtkArray = vtkArray_tmp;
      break;
    }
    /*
    case VTK_VOID:
    case VTK_BIT:
    case VTK_UNSIGNED_CHAR:
    case VTK_UNSIGNED_SHORT:
    case VTK_UNSIGNED_INT:
    case VTK_UNSIGNED_LONG:
    case VTK_STRING:
    case VTK_OPAQUE:
    case VTK_UNSIGNED_LONG_LONG:
    case VTK___INT64:
    case VTK_UNSIGNED___INT64:
    case VTK_VARIANT:
    case VTK_OBJECT:
    case VTK_UNICODE_STRING:
    */
    default: {
      std::cerr << "VTK type " << vtkArray->GetDataTypeAsString()
                << " is not supported as an Omega_h type." << std::endl;
      exit(1);
    }
  }
}

Omega_h::Mesh *oshGeoMesh::vtk2osh(vtkUnstructuredGrid *vtkMesh,
                                   Omega_h::Library *lib) {
  auto *oshMesh = new Omega_h::Mesh(lib);

  vtkSmartPointer<vtkCellTypes> ct = vtkSmartPointer<vtkCellTypes>::New();
  vtkMesh->GetCellTypes(ct);

  if (ct->GetNumberOfTypes() != 1) {
    std::cerr << "Error: Omega_h mesh does not support mixed meshes. Received "
                 "mesh with "
              << ct->GetNumberOfTypes() << " different cell types."
              << std::endl;
    exit(1);
  }

  Omega_h_Family oshFamily =
      getOmega_hFamilyFromVTKType(static_cast<VTKCellType>(ct->GetCellType(0)));
  Omega_h::Int oshDim =
      getOmega_hDimFromVTKType(static_cast<VTKCellType>(ct->GetCellType(0)));
  Omega_h::Int oshDeg =
      Omega_h::element_degree(oshFamily, oshDim, Omega_h::VERT);

  {  // Add points and cells (needs both for build_from_elems_and_coords)
    // Add points
    Omega_h::HostWrite<Omega_h::Real> h_oshCoords(oshDim *
                                                  vtkMesh->GetNumberOfPoints());
    vtkSmartPointer<vtkPoints> vtkPoints = vtkMesh->GetPoints();

    std::array<double, 3> point{};
    for (vtkIdType i = 0; i < vtkPoints->GetNumberOfPoints(); ++i) {
      vtkPoints->GetPoint(i, point.data());

      h_oshCoords[i * oshDim] = point[0];
      h_oshCoords[i * oshDim + 1] = point[1];
      if (oshDim >= 3) h_oshCoords[i * oshDim + 2] = point[2];
    }

    Omega_h::Reals oshCoords(h_oshCoords);

    // Add cells
    Omega_h::HostWrite<Omega_h::LO> h_oshEv2v(vtkMesh->GetNumberOfCells() *
                                              oshDeg);

    vtkSmartPointer<vtkIdList> ids = vtkSmartPointer<vtkIdList>::New();
    for (vtkIdType i = 0; i < vtkMesh->GetNumberOfCells(); ++i) {
      vtkMesh->GetCellPoints(i, ids);

      for (Omega_h::Int d = 0; d < oshDeg; ++d)
        h_oshEv2v.set(i * oshDeg + d, ids->GetId(d));
    }

    Omega_h::LOs oshEv2v(h_oshEv2v);

    Omega_h::build_from_elems_and_coords(oshMesh, oshFamily, oshDim, oshEv2v,
                                         oshCoords);
  }

  {  // Add point data
    for (int a = 0; a < vtkMesh->GetPointData()->GetNumberOfArrays(); ++a) {
      vtkSmartPointer<vtkDataArray> da = vtkMesh->GetPointData()->GetArray(a);

      Omega_h::HostWrite<Omega_h::Real> h_oshArray(da->GetNumberOfValues());

      getOmega_hArrayFromVtkDataArray(oshMesh, da, Omega_h::VERT);
    }
  }

  {  // Add cell data
    for (int a = 0; a < vtkMesh->GetCellData()->GetNumberOfArrays(); ++a) {
      vtkSmartPointer<vtkDataArray> da = vtkMesh->GetCellData()->GetArray(a);

      Omega_h::HostWrite<Omega_h::Real> h_oshArray(da->GetNumberOfValues());

      getOmega_hArrayFromVtkDataArray(oshMesh, da, oshDim);
    }
  }

  return oshMesh;
}

}  // namespace MSH
}  // namespace NEM
