/*******************************************************************************
* Promesh                                                                      *
* Copyright (C) 2022, IllinoisRocstar LLC. All rights reserved.                *
*                                                                              *
* Promesh is the property of IllinoisRocstar LLC.                              *
*                                                                              *
* IllinoisRocstar LLC                                                          *
* Champaign, IL                                                                *
* www.illinoisrocstar.com                                                      *
* promesh@illinoisrocstar.com                                                  *
*******************************************************************************/
/*******************************************************************************
* This file is part of Promesh                                                 *
*                                                                              *
* This version of Promesh is free software: you can redistribute it and/or     *
* modify it under the terms of the GNU Lesser General Public License as        *
* published by the Free Software Foundation, either version 3 of the License,  *
* or (at your option) any later version.                                       *
*                                                                              *
* Promesh is distributed in the hope that it will be useful, but WITHOUT ANY   *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS    *
* FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more *
* details.                                                                     *
*                                                                              *
* You should have received a copy of the GNU Lesser General Public License     *
* along with this program. If not, see <https://www.gnu.org/licenses/>.        *
*                                                                              *
*******************************************************************************/
#include "Mesh/geoMeshFactory.H"

#include "AuxiliaryFunctions.H"

#include "Mesh/exoGeoMesh.H"
#include "Mesh/inpGeoMesh.H"
#include "Mesh/oshGeoMesh.H"
#include "Mesh/vtkGeoMesh.H"

#ifdef HAVE_CFMSH
#  include "Mesh/foamGeoMesh.H"
#endif

#ifdef HAVE_GMSH
#  include "Mesh/gmshGeoMesh.H"
#endif

#ifdef HAVE_OCC
#  include "Mesh/smeshGeoMesh.H"
#endif

namespace NEM {
namespace MSH {

namespace {

MeshType MeshTypeFromFilename(const std::string &fileName) {
  std::string fileExt = nemAux::find_ext(fileName);

  if (fileExt == ".vtu" || fileExt == ".vtk") {
    return MeshType::VTK_GEO_MESH;
  } else if (fileExt == ".msh" || fileExt == ".geo") {
    return MeshType::GMSH_GEO_MESH;
  } else if (fileExt == ".osh") {
    return MeshType::OSH_GEO_MESH;
  } else if (fileExt == ".exo" || fileExt == ".e" || fileExt == ".gen" ||
             fileExt == ".g") {
    return MeshType::EXO_GEO_MESH;
  } else if (fileExt == ".inp") {
    return MeshType::INP_GEO_MESH;
  } else if (fileExt == ".foam") {
    return MeshType::FOAM_GEO_MESH;
  } else {
    std::cerr << "File extension " << fileExt << " is not supported."
              << std::endl;
    exit(1);
  }
}

MeshType MeshTypeFromFormatName(const std::string &formatName) {
  // convert to lowercase to limit the number of options
  std::string formatname = formatName;
  nemAux::toLower(formatname);

  if (formatname == "vtk") {
    return MeshType::VTK_GEO_MESH;
  } else if (formatname == "gmsh") {
    return MeshType::GMSH_GEO_MESH;
  } else if (formatname == "omega_h" || formatname == "omegah" ||
             formatname == "omega h") {
    return MeshType::OSH_GEO_MESH;
  } else if (formatname == "exodus" || formatname == "exodus ii") {
    return MeshType::EXO_GEO_MESH;
  } else if (formatname == "calculix" || formatname == "abaqus") {
    return MeshType::INP_GEO_MESH;
  } else if (formatname == "open foam" || formatname == "openfoam" ||
             formatname == "foam") {
    return MeshType::FOAM_GEO_MESH;
  } else {
    std::cerr << "Format name " << formatName << " is not supported."
              << std::endl;
    exit(1);
  }
}

}  // namespace

geoMeshBase *Read(const std::string &fileName) {
  return Read(fileName, MeshTypeFromFilename(fileName));
}

geoMeshBase *Read(const std::string &fileName, const std::string &formatName) {
  if (formatName.empty()) {
    return Read(fileName);
  } else {
    return Read(fileName, MeshTypeFromFormatName(formatName));
  }
}

geoMeshBase *Read(const std::string &fileName, MeshType meshType) {
  switch (meshType) {
    case MeshType::VTK_GEO_MESH: return vtkGeoMesh::Read(fileName);
    case MeshType::GMSH_GEO_MESH: {
#ifdef HAVE_GMSH
      return gmshGeoMesh::Read(fileName);
#else
      std::cerr << "Please build NEMoSys with ENABLE_GMSH=ON to use this"
                << " feature!" << std::endl;
      throw;
#endif
    }
    case MeshType::OSH_GEO_MESH: return oshGeoMesh::Read(fileName);
    case MeshType::EXO_GEO_MESH: return exoGeoMesh::Read(fileName);
    case MeshType::INP_GEO_MESH: return inpGeoMesh::Read(fileName);
    case MeshType::FOAM_GEO_MESH: {
#ifdef HAVE_CFMSH
      std::string fname = fileName;
      size_t dot_loc = fname.find_last_of('.');
      if (dot_loc != std::string::npos) {
        fname.erase(dot_loc);
      }
      return foamGeoMesh::Read(fname);
#else
      std::cerr << "Please build NEMoSys with ENABLE_CFMSH=ON to use this"
                << " feature!" << std::endl;
      throw;
#endif
    }
    case MeshType::SMESH_GEO_MESH: {
      std::cerr << "Unsupported.\n";
      throw;
    }
  }
  // Don't put inside switch statement so static analyzers can detect unhandled
  // cases
  return nullptr;
}

geoMeshBase *New(MeshType meshType) {
  switch (meshType) {
    case MeshType::VTK_GEO_MESH: return vtkGeoMesh::New();
    case MeshType::GMSH_GEO_MESH: {
#ifdef HAVE_GMSH
      return gmshGeoMesh::New();
#else
      std::cerr << "Please build NEMoSys with ENABLE_GMSH=ON to use this"
                << " feature!" << std::endl;
      throw;
#endif
    }
    case MeshType::OSH_GEO_MESH: return oshGeoMesh::New();
    case MeshType::EXO_GEO_MESH: return exoGeoMesh::New();
    case MeshType::INP_GEO_MESH: return inpGeoMesh::New();
    case MeshType::FOAM_GEO_MESH: {
#ifdef HAVE_CFMSH
      return foamGeoMesh::New();
#else
      std::cerr << "Please build NEMoSys with ENABLE_CFMSH=ON to use this"
                << " feature!" << std::endl;
      throw;
#endif
    }
    case MeshType::SMESH_GEO_MESH: {
#ifdef HAVE_OCC
      return smeshGeoMesh::New();
#else
      std::cerr << "Requires ENABLE_OPENCASCADE=ON\n";
      throw;
#endif
    }
  }
  return nullptr;
}

geoMeshBase *New(const std::string &fileName) {
  return New(MeshTypeFromFilename(fileName));
}

geoMeshBase *New(const std::string &fileName, const std::string &formatName) {
  if (formatName.empty()) {
    return New(fileName);
  } else {
    return New(MeshTypeFromFormatName(formatName));
  }
}

}  // namespace MSH
}  // namespace NEM
