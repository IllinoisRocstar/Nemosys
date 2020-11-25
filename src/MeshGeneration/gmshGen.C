#if defined(_MSC_VER) && !defined(_USE_MATH_DEFINES)
#  define _USE_MATH_DEFINES
#endif

#include "AuxiliaryFunctions.H"
#include "gmshGen.H"
#include "gmshParams.H"
#include "gmshTypes.H"
#include "meshGen.H"

#include <gmsh.h>

#include <iostream>
#include <utility>
#include <vector>

#include <math.h>

namespace NEM {

namespace GEN {

// Default Constructor
gmshGen::gmshGen() {
  // Default meshing parameters
  meshParams = new gmshParams();
  meshParams->minSize = 0.01;
  meshParams->maxSize = 50.0;
  meshParams->algo2D = "Frontal";
  meshParams->algo3D = "HXT";
  meshParams->extSizeFromBoundary = true;
  meshParams->sizeFromCurvature = false;
  meshParams->minElePer2Pi = 6;
  meshParams->optimize = false;
  meshParams->optimizeThreshold = 0.3;
}

gmshGen::gmshGen(gmshParams *params) : meshParams(params) {}

gmshGen::~gmshGen() {}

int gmshGen::createMeshFromSTL(const char *fname) {
  int ret;
  ret = createMeshFromSTEP(fname);
  if (ret == 0)
    return 0;
  else
    return 1;
}

int gmshGen::createMeshFromSTEP(const char *fname) {
  std::cout << "Mimimum mesh size " << meshParams->minSize << "\n"
            << "Maximum mesh size " << meshParams->maxSize << "\n"
            << "Surface algorithm " << meshParams->algo2D << "\n"
            << "Volume algorithm " << meshParams->algo3D << std::endl;
  // Initialize Gmsh
  gmsh::initialize();
  // Set global options from JSON
  globalOptions();

  std::string file = fname;
  std::vector<std::pair<int, int>> imported;
  std::ifstream f;
  f.open(file);
  if (!f.good()) {
    std::cerr << "Error: Input Geometry File " << file << " not found."
              << std::endl;
    exit(-1);
  }
  f.close();
  gmsh::model::occ::importShapes(file, imported, false);
  gmsh::model::occ::synchronize();

  // getGeomNames();

  // Apply mesh size fields
  if (meshParams->mSizeField) meshSizeFields();

  gmsh::model::mesh::generate(3);
  gmsh::model::mesh::removeDuplicateNodes();
  // By default, ouput MSH file
  std::string newname = nemAux::trim_fname(fname, ".msh");
  gmsh::write(newname);
  // Also output another format
  std::string ext = nemAux::find_ext(meshParams->ofname);
  if (ext != ".vtu") {
    auto pos = meshParams->meshExtensions.find(ext);
    if (pos == meshParams->meshExtensions.end()) {
      std::cerr << "Error: Output file extension " << ext
                << " not supported in this mesh engine." << std::endl;
      std::cout << "Supported formats are ";
      for (auto it = meshParams->meshExtensions.begin();
           it != meshParams->meshExtensions.end(); ++it) {
        std::cout << *it << " ";
      }
      std::cout << "\n" << std::endl;
    } else {
      // IF save format is INP or UNV, save physical groups to blocks/nodesets
      if (ext == ".inp" || ext == ".unv") {
        gmsh::option::setNumber("Mesh.SaveGroupsOfElements", 1);
        gmsh::option::setNumber("Mesh.SaveGroupsOfNodes", 1);
      }
      gmsh::write(meshParams->ofname);
    }
  }
  // gmsh::fltk::run();
  gmsh::finalize();
  return 0;
}

//********************************************************************//
// Sets the global geometry and mesh options
//********************************************************************//
void gmshGen::globalOptions() {
  //------------------------ General Options ----------------------//
  //---------------------------------------------------------------//
  gmsh::option::setNumber("General.Terminal", 1);

  //----------------------- Geometry Options ----------------------//
  //---------------------------------------------------------------//
  gmsh::option::setNumber("Geometry.OCCImportLabels", 1);

  //------------------------- Mesh Options ------------------------//
  //---------------------------------------------------------------//
  gmsh::option::setNumber("Mesh.MshFileVersion", 2.2);
  gmsh::option::setNumber("Mesh.CharacteristicLengthMin", meshParams->minSize);
  gmsh::option::setNumber("Mesh.CharacteristicLengthMax", meshParams->maxSize);
  if (meshParams->extSizeFromBoundary)
    gmsh::option::setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 1);
  else
    gmsh::option::setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 0);
  if (meshParams->sizeFromCurvature)
    gmsh::option::setNumber("Mesh.CharacteristicLengthFromCurvature", 1);
  else
    gmsh::option::setNumber("Mesh.CharacteristicLengthFromCurvature", 0);
  gmsh::option::setNumber("Mesh.MinimumElementsPerTwoPi",
                          meshParams->minElePer2Pi);
  if (meshParams->optimize)
    gmsh::option::setNumber("Mesh.Optimize", 1);
  else
    gmsh::option::setNumber("Mesh.Optimize", 0);
  gmsh::option::setNumber("Mesh.OptimizeThreshold",
                          meshParams->optimizeThreshold);
  gmsh::option::setNumber("Mesh.SaveAll", 1);
  gmsh::option::setNumber("Mesh.ElementOrder", 1);
  gmsh::option::setNumber("Mesh.RecombineAll", 0);
  gmsh::option::setNumber("Mesh.Smoothing", 1);
  gmsh::option::setNumber("Mesh.CharacteristicLengthFactor", 1);
  gmsh::option::setNumber("Mesh.CharacteristicLengthFromPoints", 1);
  gmsh::option::setNumber("Mesh.LcIntegrationPrecision", 1e-09);
  gmsh::option::setNumber("Mesh.MaxIterDelaunay3D", 0);
  gmsh::option::setNumber("Mesh.RandomFactor", 1e-09);
  gmsh::option::setNumber("Mesh.RandomFactor3D", 1e-12);

  std::string algo = meshParams->algo2D;
  int a = ALGO_SURF_FRONTAL;  // default: "Frontal"
  if (algo == "Frontal")
    a = ALGO_SURF_FRONTAL;
  else if (algo == "MeshAdapt")
    a = ALGO_SURF_MESHADAPT;
  else if (algo == "Automatic")
    a = ALGO_SURF_AUTO;
  else if (algo == "Delaunay")
    a = ALGO_SURF_DELAUNAY;
  else if (algo == "Frontal Quads" || algo == "Frontal Quad")
    a = ALGO_SURF_FRONTALQUAD;
  else if (algo == "Packing of Parallelograms")
    a = ALGO_SURF_POP;
  else {
    std::cerr << "Warning: Surface algorithm " << algo << " not supported. \n"
              << "Using default Frontal algorithm." << std::endl;
  }
  gmsh::option::setNumber("Mesh.Algorithm", a);

  algo = meshParams->algo3D;
  a = ALGO_VOL_DELAUNAY;
  if (algo == "Delaunay")
    a = ALGO_VOL_DELAUNAY;
  else if (algo == "Frontal")
    a = ALGO_VOL_FRONTAL;
  else if (algo == "HXT")
    a = ALGO_VOL_HXT;
  else {
    std::cerr << "Warning: Volume algorithm " << algo << " not supported. \n"
              << "Using default Delaunay algorithm." << std::endl;
  }
  gmsh::option::setNumber("Mesh.Algorithm3D", a);
}

//********************************************************************//
// Gets the surface colors from geometry
//********************************************************************//
void gmshGen::getGeomNames() {
  // Get all surfaces
  std::vector<std::pair<int, int>> surfaces;
  gmsh::model::getEntities(surfaces, 2);
  std::string name = "";
  for (const std::pair<int, int> &s : surfaces) {
    int dim = s.first;
    int tag = s.second;
    gmsh::model::getEntityName(dim, tag, name);
    std::cout << "Surface " << tag << " name: " << name << std::endl;
  }
  gmsh::model::getEntityName(3, 1, name);
  std::cout << "Volume " << 1 << " name: " << name << std::endl;
}

//********************************************************************//
// Gets the surface colors from geometry
//********************************************************************//
void gmshGen::getSurfaceColors() {
  // Get all surfaces
  std::vector<std::pair<int, int>> surfaces;
  gmsh::model::getEntities(surfaces, 2);
  for (const std::pair<int, int> &s : surfaces) {
    int dim = s.first;
    int tag = s.second;
    int r, g, b, a;
    gmsh::model::getColor(dim, tag, r, g, b, a);
    std::cout << "Surface " << tag << " color: " << r << "," << g << "," << b
              << std::endl;
  }
}

//********************************************************************//
// Applies mesh size fields
//********************************************************************//
void gmshGen::meshSizeFields() {
  bool bgfAssigned = false;
  int highest_id = -1;
  int id = -1;

  // Loop through fields in search of Frustum.IField option
  for (auto &sf : meshParams->sizeFields) {
    std::string type = sf.type;
    if (sf.type == "Frustum") {
      double thick = 0.0;
      for (auto &prm : sf.params) {
        if (prm.first == "Thickness") {
          thick = prm.second;
          break;
        }
      }
      for (auto &prm : sf.params) {
        std::string key = prm.first;
        double val = prm.second;
        if (key == "IField") {
          for (auto &sf2 : meshParams->sizeFields) {
            if (sf2.type == "Cylinder" && sf2.id == val) {
              double radius = 0;
              double thick = 0;
              double vin = 10.0;
              double vout = 10.0;
              double xc = 0.0;
              double yc = 0.0;
              double zc = 0.0;
              double xa = 0.0;
              double ya = 0.0;
              double za = 0.0;
              // Get Cylinder field parameters
              for (auto &prm2 : sf2.params) {
                if (prm2.first == "Radius") radius = prm2.second;
                if (prm2.first == "Thickness") thick = prm2.second;
                if (prm2.first == "VIn") vin = prm2.second;
                if (prm2.first == "VOut") vout = prm2.second;
                if (prm2.first == "XCenter") xc = prm2.second;
                if (prm2.first == "YCenter") yc = prm2.second;
                if (prm2.first == "ZCenter") zc = prm2.second;
                if (prm2.first == "XAxis") xa = prm2.second;
                if (prm2.first == "YAxis") ya = prm2.second;
                if (prm2.first == "ZAxis") za = prm2.second;
              }

              // Clear the Frustum parameters
              sf.params.clear();

              // Calculate Frustum parameters from Cylinder params and insert
              std::pair<std::string, double> p;
              p = {"R1_inner", radius};
              sf.params.push_back(p);
              p = {"R1_outer", radius + thick};
              sf.params.push_back(p);
              p = {"R2_inner", radius};
              sf.params.push_back(p);
              p = {"R2_outer", radius + thick};
              sf.params.push_back(p);

              p = {"V1_inner", vin};
              sf.params.push_back(p);
              p = {"V2_inner", vin};
              sf.params.push_back(p);
              p = {"V1_outer", vout};
              sf.params.push_back(p);
              p = {"V2_outer", vout};
              sf.params.push_back(p);

              double mag = sqrt(xa * xa + ya * ya + za * za);

              p = {"X1", xc - xa / 2 - xa * (thick / mag)};
              sf.params.push_back(p);
              p = {"Y1", yc - ya / 2 - ya * (thick / mag)};
              sf.params.push_back(p);
              p = {"Z1", zc - za / 2 - za * (thick / mag)};
              sf.params.push_back(p);

              p = {"X2", xc + xa / 2 + xa * (thick / mag)};
              sf.params.push_back(p);
              p = {"Y2", yc + ya / 2 + ya * (thick / mag)};
              sf.params.push_back(p);
              p = {"Z2", zc + za / 2 + za * (thick / mag)};
              sf.params.push_back(p);
            }
          }
          break;
        }
      }
    }
  }

  // Loop through all size fields
  for (auto sf = (meshParams->sizeFields).begin();
       sf != (meshParams->sizeFields).end(); sf++) {
    if (sf->id == meshParams->bgField) {
      bgfAssigned = true;
      id = sf->id;
    }

    gmsh::model::mesh::field::add(sf->type, sf->id);
    for (auto prm = (sf->params).begin(); prm != (sf->params).end(); prm++) {
      if (prm->first != "")
        gmsh::model::mesh::field::setNumber(sf->id, prm->first, prm->second);
    }

    // Keep this code for later when surface names are used
    // for (auto prm = (sf->strg_list_params).begin();
    // 					prm != (sf->strg_list_params).end();
    // prm++)
    // {
    // 	gmsh::model::mesh::field::setNumbers(sf->id, prm->first, prm->second);
    // 	//std::cout << sf->id << " " << prm->first << " " << prm->second <<
    // std::endl;
    // }

    for (auto prm = (sf->num_list_params).begin();
         prm != (sf->num_list_params).end(); prm++) {
      gmsh::model::mesh::field::setNumbers(sf->id, prm->first, prm->second);
    }

    highest_id = sf->id;
  }

  // Assign the Background Field (size field)
  if (bgfAssigned) {
    gmsh::model::mesh::field::setAsBackgroundMesh(meshParams->bgField);
    std::cout << "Using size field ID = " << id << std::endl;
  } else {
    gmsh::model::mesh::field::setAsBackgroundMesh(highest_id);
    std::cout << "Warning: BackgroundField ID not found. "
              << "Using size field with highest ID." << std::endl;
  }
  gmsh::model::occ::synchronize();
}

}  // namespace GEN

}  // namespace NEM
