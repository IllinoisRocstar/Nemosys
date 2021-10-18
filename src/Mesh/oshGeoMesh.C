#if defined(_MSC_VER) && !defined(_USE_MATH_DEFINES)
#  define _USE_MATH_DEFINES
#endif

#include "Mesh/oshGeoMesh.H"

#include <array>

#include <vtkCellIterator.h>
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
#include <Omega_h_align.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_class.hpp>
#include <Omega_h_element.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_mesh.hpp>

#include "AuxiliaryFunctions.H"

namespace NEM {
namespace MSH {

std::shared_ptr<Omega_h::Library> OmegaHInterface::GetLibrary() {
  return instance->library;
}

OmegaHInterface::OmegaHInterface()
    : library(std::make_shared<Omega_h::Library>()) {}

std::shared_ptr<OmegaHInterface> OmegaHInterface::instance{new OmegaHInterface};

vtkStandardNewMacro(oshGeoMesh)

oshGeoMesh *oshGeoMesh::Read(const std::string &fileName,
                             Omega_h::Library *lib) {
  if (!lib) {
    lib = OmegaHInterface::GetLibrary().get();
  }
  Omega_h::Mesh oshMesh = Omega_h::read_mesh_file(fileName, lib->world());

  return new oshGeoMesh(&oshMesh, lib);
}

oshGeoMesh::oshGeoMesh() : oshGeoMesh(nullptr) {}

oshGeoMesh::oshGeoMesh(Omega_h::Mesh *oshMesh, Omega_h::Library *lib)
    : geoMeshBase(osh2GM(oshMesh)), _oshMesh() {
  if (!lib) lib = OmegaHInterface::GetLibrary().get();
  std::cout << "oshGeoMesh constructed" << std::endl;
  _oshMesh = std::unique_ptr<Omega_h::Mesh>(
      oshMesh ? new Omega_h::Mesh(*oshMesh) : new Omega_h::Mesh());
  _oshMesh->set_library(lib);
}

oshGeoMesh::~oshGeoMesh() { std::cout << "oshGeoMesh destructed" << std::endl; }

void oshGeoMesh::write(const std::string &fileName) {
  std::string fileExt = nemAux::find_ext(fileName);

  if (fileExt == ".osh") {
    Omega_h::binary::write(fileName, _oshMesh.get());
  } else if (fileExt == ".msh") {
    Omega_h::gmsh::write(fileName, _oshMesh.get());
  } else if (fileExt == ".pvtu") {
    Omega_h::vtk::write_parallel(fileName, _oshMesh.get());
  } else if (fileExt == ".vtu") {
    Omega_h::vtk::write_vtu(fileName, _oshMesh.get());
  } else {
    std::cerr << "Omega_h library does not support writing " << fileExt
              << " format." << std::endl;
  }
}

void oshGeoMesh::report(std::ostream &out) const {
  geoMeshBase::report(out);
  out << "Family:\t"
      << (_oshMesh->family() == OMEGA_H_SIMPLEX ? "SIMPLEX" : "HYPERCUBE")
      << '\n';
  out << "Dim:\t" << _oshMesh->dim() << '\n';
}

void oshGeoMesh::takeGeoMesh(geoMeshBase *otherGeoMesh) {
  auto otherOshGM = dynamic_cast<oshGeoMesh *>(otherGeoMesh);
  if (otherOshGM) {
    setGeoMesh(otherOshGM->getGeoMesh());
    otherOshGM->setGeoMesh(
        {vtkSmartPointer<vtkUnstructuredGrid>::New(), {}, {}, {}});
    _oshMesh = std::move(otherOshGM->_oshMesh);
    otherOshGM->resetNative();
  } else {
    geoMeshBase::takeGeoMesh(otherGeoMesh);
  }
}

void oshGeoMesh::reconstructGeo() {
  geoMeshBase::reconstructGeo();
  auto mesh = getGeoMesh().mesh;
  auto link = getGeoMesh().link;
  auto sideSet = getGeoMesh().sideSet;
  auto oshFamily = _oshMesh->family();
  auto oshDim = _oshMesh->dim();
  if (sideSet.sides && sideSet.sides->GetNumberOfCells() > 0) {
    for (int i = 0; i <= oshDim; ++i) {
      if (_oshMesh->has_tag(i, "class_dim")) {
        _oshMesh->remove_tag(i, "class_dim");
      }
      if (_oshMesh->has_tag(i, "class_id")) {
        _oshMesh->remove_tag(i, "class_id");
      }
    }
    {  // Set class_dim and class_id for dimension oshDim
      Omega_h::HostWrite<Omega_h::ClassId> h_class_id(mesh->GetNumberOfCells());
      auto linkArr = mesh->GetCellData()->GetArray(link.c_str());
      if (linkArr) {
        for (Omega_h::LO i = 0; i < h_class_id.size(); ++i) {
          h_class_id[i] =
              static_cast<Omega_h::ClassId>(linkArr->GetComponent(i, 0));
        }
        // Note we pass dummy because classify_equal_order just assumes that
        // we pass in h_class_id in same order that we create cells.
        // This is fine because reconstructGeo doesn't change the mesh.
        auto dummyEqv2v = Omega_h::LOs(
            mesh->GetNumberOfCells() *
                Omega_h::element_degree(oshFamily, oshDim, Omega_h::VERT),
            0);
        Omega_h::classify_equal_order(_oshMesh.get(), oshDim, dummyEqv2v,
                                      h_class_id.write());
      }
    }
    {  // Set class_dim and class_id for dimension oshDim - 1
      auto it = vtkSmartPointer<vtkCellIterator>::Take(
          sideSet.sides->NewCellIterator());
      auto numNodes =
          Omega_h::element_degree(oshFamily, oshDim - 1, Omega_h::VERT);
      Omega_h::HostWrite<Omega_h::LO> ev2v(numNodes *
                                           sideSet.sides->GetNumberOfCells());
      for (int i = (it->InitTraversal(), 0); !it->IsDoneWithTraversal();
           it->GoToNextCell(), ++i) {
        auto ids = it->GetPointIds();
        for (Omega_h::Int d = 0; d < numNodes; ++d) {
          ev2v[numNodes * i + d] = ids->GetId(d);
        }
      }
      Omega_h::LOs eqv2v(ev2v);
      // If we have a non-empty sideSet, we should have a non-empty geo ent
      // array
      Omega_h::HostWrite<Omega_h::LO> h_class_id(
          sideSet.sides->GetNumberOfCells());
      auto geoEntArr = sideSet.getGeoEntArr();
      for (Omega_h::LO i = 0; i < h_class_id.size(); ++i) {
        h_class_id[i] = geoEntArr->GetValue(i);
      }
      Omega_h::classify_equal_order(_oshMesh.get(), oshDim - 1, eqv2v,
                                    h_class_id.write());
    }
  }
  Omega_h::finalize_classification(_oshMesh.get());
}

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

    case VTK_VERTEX: return 0;
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

int oshFace2vtkFace(int oshFace, Omega_h_Family oshFamily,
                    Omega_h::Int oshDim) {
  switch (oshDim) {
    case 2: return oshFace;
    case 3:
      switch (oshFamily) {
        case OMEGA_H_SIMPLEX: return (oshFace + 3) % 4;
        case OMEGA_H_HYPERCUBE:
          switch (oshFace) {
            case 0: return 4;
            case 1: return 2;
            case 2: return 1;
            case 3: return 3;
            case 4: return 0;
            case 5: return 5;
            default: return -1;
          }
      }
    default: return -1;
  }
}

geoMeshBase::GeoMesh oshGeoMesh::osh2GM(Omega_h::Mesh *oshMesh,
                                        const std::string &geo,
                                        const std::string &link) {
  auto vtkMesh = vtkSmartPointer<vtkUnstructuredGrid>::New();

  if (!oshMesh || !oshMesh->is_valid()) return {vtkMesh, "", "", {}};
  // Omega_h dimension is important to read its data.
  //
  // It stores 2D mesh without a z-coordinate while VTK always requires 3D. A
  // "2D mesh" in VTK will have z-coordinates set to zero.
  //
  // Also, an Omega_h mesh contains the main element and all its lower-dimension
  // constituents. E.g., a 3D mesh containing tetrahedrons will also contain the
  // four triangles, six edges, and four vertices. When converting to VTK, we
  // will only add the lower-dimensional elements if the classification of
  // elements appears valid.
  Omega_h::Int dim = oshMesh->dim();

  // An Omega_h mesh cannot have mixed elements. The family determines if the
  // mesh is composed of simplices (tetrahedrons and triangles) or hypercubes
  // (hexahedrons and quadrangles).
  Omega_h_Family family = oshMesh->family();
  VTKCellType vtk_type = getVTKTypeFromOmega_hFamilyDim(family, dim);

  auto linkName = link;
  if (link.empty()) {
    linkName = GEO_ENT_DEFAULT_NAME;
  }
  bool validSideSet = true;

  {  // Add points
    vtkSmartPointer<vtkPoints> points =
        vtkSmartPointer<vtkPoints>::Take(vtkPoints::New(VTK_DOUBLE));
    points->Allocate(oshMesh->nverts());
    std::vector<double> nodes;

    Omega_h::HostRead<Omega_h::Real> h_coord(oshMesh->coords());

    for (Omega_h::LO i = 0; i < oshMesh->nverts(); ++i) {
      points->InsertNextPoint(h_coord[i * dim + 0], h_coord[i * dim + 1],
                              dim == 2 ? 0 : h_coord[i * dim + 2]);
    }
    // Because we add them in order, the indices should be the same.
    vtkMesh->SetPoints(points);
    std::vector<size_t> nodeTags;
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
      auto oshTagName = oshMesh->get_tag(Omega_h::VERT, itag)->name();
      if (oshTagName != "class_dim" && oshTagName != "class_id" &&
          oshTagName != "coordinates" && oshTagName != "global" &&
          oshTagName != "local") {
        getVtkDataArrayFromOmega_hTag(oshMesh->get_tag(Omega_h::VERT, itag),
                                      vtkArray);
      }

      vtkMesh->GetPointData()->AddArray(vtkArray);
    }
  }

  {  // Add cell data
    for (Omega_h::Int itag = 0; itag < oshMesh->ntags(dim); ++itag) {
      vtkSmartPointer<vtkDataArray> vtkArray;
      auto oshTagName = oshMesh->get_tag(dim, itag)->name();
      if (oshTagName != "class_dim" && oshTagName != "class_id" &&
          oshTagName != "global" && oshTagName != "local") {
        getVtkDataArrayFromOmega_hTag(oshMesh->get_tag(dim, itag), vtkArray);
        vtkMesh->GetCellData()->AddArray(vtkArray);
      }
    }
  }

  // Add boundary elements (of one lower dimension), if they exist, so we don't
  // have to reconstruct geometry.
  auto sideSetPD = vtkSmartPointer<vtkPolyData>::New();
  sideSetPD->SetPoints(vtkMesh->GetPoints());
  sideSetPD->Allocate();
  auto sideSetEntities = vtkSmartPointer<vtkIntArray>::New();
  auto sideSetOrigCellId = vtkSmartPointer<vtkIdTypeArray>::New();
  sideSetOrigCellId->SetNumberOfComponents(2);
  auto sideSetCellFaceId = vtkSmartPointer<vtkIntArray>::New();
  sideSetCellFaceId->SetNumberOfComponents(2);
  // Map from highest dimensional entities to bounding entities of one lower
  // dimension (by "class_id").
  std::map<Omega_h::ClassId, std::set<Omega_h::ClassId>> entities;
  // List of bounding entities of one lower dimension (by "class_id").
  std::vector<Omega_h::ClassId> entities_boundary;
  Omega_h::HostRead<Omega_h::ClassId> arr_id;
  {  // Add all entities of max dimension
    if (oshMesh->has_tag(dim, "class_id")) {
      arr_id = oshMesh->get_array<Omega_h::ClassId>(dim, "class_id");
    } else {
      validSideSet = false;
    }
  }
  // Find boundary elements of dimension dim - 1. Also recover the bounding
  // entities (of dimension dim - 1) for the entities of highest dimension
  if (validSideSet && oshMesh->has_tag(dim - 1, "class_id") &&
      oshMesh->has_tag(dim - 1, "class_dim")) {
    Omega_h::HostRead<Omega_h::ClassId> arr_id_boundary =
        oshMesh->get_array<Omega_h::ClassId>(dim - 1, "class_id");
    Omega_h::HostRead<Omega_h::I8> arr_dim =
        oshMesh->get_array<Omega_h::I8>(dim - 1, "class_dim");
    auto parent = oshMesh->ask_up(dim - 1, dim);
    Omega_h::HostRead<Omega_h::LO> parent_a2ab = parent.a2ab;
    Omega_h::HostRead<Omega_h::LO> parent_ab2b = parent.ab2b;
    Omega_h::HostRead<Omega_h::I8> parent_codes = parent.codes;
    for (Omega_h::LO i = 0; i < arr_dim.size(); ++i) {
      if (arr_dim[i] == dim - 1) {
        bool newEntityBoundary = false;
        auto oshEntity = arr_id_boundary[i];
        auto entityIter = std::find(entities_boundary.begin(),
                                entities_boundary.end(), oshEntity);
        if (entityIter == entities_boundary.end()) {
          entities_boundary.emplace_back(oshEntity);
          newEntityBoundary = true;
        }
        std::array<vtkIdType, 2> origCell{-1, -1};
        std::array<int, 2> cellFace{-1, -1};
        Omega_h::ClassId firstParentEnt = -1;
        for (int j = 0; j < parent_a2ab[i + 1] - parent_a2ab[i]; ++j) {
          auto ab = parent_ab2b[i + j];
          auto b = parent_ab2b[ab];
          if (j == 0) {
            firstParentEnt = arr_id[b];
          }
          if (newEntityBoundary) {
            entities[arr_id[b]].emplace(oshEntity);
          }
          auto code = parent_codes[ab];
          auto oshFace = Omega_h::code_which_down(code);
          if (j == 0 || (j == 1 && arr_id[b] < firstParentEnt)) {
            origCell[1] = origCell[0];
            cellFace[1] = cellFace[0];
            origCell[0] = b;
            cellFace[0] = oshFace2vtkFace(oshFace, family, dim);
          } else {
            origCell[1] = b;
            cellFace[1] = oshFace2vtkFace(oshFace, family, dim);
          }
        }
        auto side = dim == 2
                        ? vtkMesh->GetCell(origCell[0])->GetEdge(cellFace[0])
                        : vtkMesh->GetCell(origCell[0])->GetFace(cellFace[0]);
        sideSetPD->InsertNextCell(side->GetCellType(), side->GetPointIds());
        sideSetOrigCellId->InsertNextTypedTuple(origCell.data());
        sideSetCellFaceId->InsertNextTypedTuple(cellFace.data());
        sideSetEntities->InsertNextValue(oshEntity);
      }
    }
    Omega_h::HostRead<Omega_h::LO> point_id = oshMesh->ask_verts_of(dim - 1);
    // It might be more efficient to get the underlying array pointer
    // and use std::copy.
    auto vtkEntities = vtkSmartPointer<vtkIntArray>::New();
    vtkEntities->SetName(linkName.c_str());
    vtkEntities->SetNumberOfValues(vtkMesh->GetNumberOfCells());
    for (vtkIdType i = 0; i < vtkMesh->GetNumberOfCells(); ++i) {
      // The physical group name matches the class_id
      vtkEntities->SetValue(i, arr_id[i]);
    }
    vtkMesh->GetCellData()->AddArray(vtkEntities);
    return {vtkMesh,
            {},
            linkName,
            {sideSetPD, sideSetEntities, sideSetOrigCellId, sideSetCellFaceId}};
  } else {
    return {vtkMesh, {}, "", {}};
  }
}

template <typename VT, typename OT>
void Varray2Otag(Omega_h::Mesh *oshMesh, VT *vtkArray, Omega_h::Int oshDim) {
  Omega_h::HostWrite<OT> h_oshArray;
  h_oshArray = Omega_h::HostWrite<OT>(vtkArray->GetNumberOfValues(),
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

Omega_h::Mesh *oshGeoMesh::GM2osh(const GeoMesh &geoMesh,
                                  Omega_h::Library *lib) {
  if (!lib) {
    lib = OmegaHInterface::GetLibrary().get();
  }
  auto *oshMesh = new Omega_h::Mesh(lib);

  auto vtkMesh = geoMesh.mesh;
  if (vtkMesh->GetNumberOfCells() == 0) {
    return oshMesh;
  }
  auto link = geoMesh.link;
  auto sideSet = geoMesh.sideSet;
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

    {  // Add cells of dimension oshDim
      auto it = vtkMesh->NewCellIterator();
      auto numNodes = Omega_h::element_degree(oshFamily, oshDim, Omega_h::VERT);
      Omega_h::HostWrite<Omega_h::LO> ev2v(numNodes *
                                           vtkMesh->GetNumberOfCells());
      int i = 0;
      for (it->InitTraversal(); !it->IsDoneWithTraversal();
           it->GoToNextCell()) {
        auto ids = it->GetPointIds();
        for (Omega_h::Int d = 0; d < numNodes; ++d) {
          ev2v[numNodes * i + d] = ids->GetId(d);
        }
        ++i;
      }
      it->Delete();
      Omega_h::LOs eqv2v(ev2v);
      Omega_h::build_from_elems_and_coords(oshMesh, oshFamily, oshDim, eqv2v,
                                           Omega_h::Reals(h_oshCoords));
      if (!link.empty()) {
        Omega_h::HostWrite<Omega_h::ClassId> h_class_id(
            vtkMesh->GetNumberOfCells());
        auto linkArr = vtkMesh->GetCellData()->GetArray(link.c_str());
        if (linkArr) {
          for (i = 0; i < h_class_id.size(); ++i) {
            h_class_id[i] =
                static_cast<Omega_h::ClassId>(linkArr->GetComponent(i, 0));
          }
          Omega_h::classify_equal_order(oshMesh, oshDim, eqv2v,
                                        h_class_id.write());
        }
      }
    }
    // Add cells of dimension oshDim - 1 from the sideSet
    if (sideSet.sides && sideSet.sides->GetNumberOfCells() > 0) {
      auto it = sideSet.sides->NewCellIterator();
      auto numNodes =
          Omega_h::element_degree(oshFamily, oshDim - 1, Omega_h::VERT);
      Omega_h::HostWrite<Omega_h::LO> ev2v(numNodes *
                                           sideSet.sides->GetNumberOfCells());
      int i = 0;
      for (it->InitTraversal(); !it->IsDoneWithTraversal();
           it->GoToNextCell()) {
        auto ids = it->GetPointIds();
        for (Omega_h::Int d = 0; d < numNodes; ++d) {
          ev2v[numNodes * i + d] = ids->GetId(d);
        }
        ++i;
      }
      it->Delete();
      Omega_h::LOs eqv2v(ev2v);
      // If we have a non-empty sideSet, we should have a non-empty geo ent
      // array
      Omega_h::HostWrite<Omega_h::LO> h_class_id(
          sideSet.sides->GetNumberOfCells());
      auto geoEntArr = geoMesh.sideSet.getGeoEntArr();
      for (i = 0; i < h_class_id.size(); ++i) {
        h_class_id[i] = geoEntArr->GetValue(i);
      }
      Omega_h::classify_equal_order(oshMesh, oshDim - 1, eqv2v,
                                    h_class_id.write());
    }
  }

  {  // Add point data
    for (int a = 0; a < vtkMesh->GetPointData()->GetNumberOfArrays(); ++a) {
      vtkSmartPointer<vtkDataArray> da = vtkMesh->GetPointData()->GetArray(a);
      if (da) {
        getOmega_hArrayFromVtkDataArray(oshMesh, da, Omega_h::VERT);
      }
    }
  }

  {  // Add cell data
    for (int a = 0; a < vtkMesh->GetCellData()->GetNumberOfArrays(); ++a) {
      vtkSmartPointer<vtkDataArray> da = vtkMesh->GetCellData()->GetArray(a);
      if (da && strcmp(da->GetName(), link.c_str()) != 0) {
        getOmega_hArrayFromVtkDataArray(oshMesh, da, oshDim);
      }
    }
  }
  Omega_h::finalize_classification(oshMesh);
  return oshMesh;
}

void oshGeoMesh::setOshMesh(Omega_h::Mesh *oshMesh) {
  this->setGeoMesh(osh2GM(oshMesh, getGeoMesh().geo, getGeoMesh().link));
  _oshMesh = std::unique_ptr<Omega_h::Mesh>(new Omega_h::Mesh(*oshMesh));
}

void oshGeoMesh::resetNative() {
  _oshMesh = std::unique_ptr<Omega_h::Mesh>(GM2osh(getGeoMesh()));
}

}  // namespace MSH
}  // namespace NEM
