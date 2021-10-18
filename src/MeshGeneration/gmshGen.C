#if defined(_MSC_VER) && !defined(_USE_MATH_DEFINES)
#define _USE_MATH_DEFINES
#endif

#include "AuxiliaryFunctions.H"
#include "MeshGeneration/gmshGen.H"
#include "MeshGeneration/gmshParams.H"
#include "Mesh/gmshTypes.H"

#include <gmsh.h>

#include <iostream>
#include <utility>
#include <vector>

#include <cmath>

namespace NEM {

namespace GEN {

// Default Constructor
gmshGen::gmshGen() {
  // Default meshing parameters
  meshParams = new gmshParams();
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
    gmsh::finalize();
    exit(-1);
  }
  f.close();
  gmsh::model::occ::importShapes(file, imported, false);
  gmsh::model::occ::synchronize();

  // getGeomNames();

  if(meshParams->fragmentAll) {
    // get all volumes
    gmsh::vectorpair volumes;
    gmsh::model::getEntities(volumes, 3);
    // fragment
    gmsh::vectorpair outTags;
    std::vector<gmsh::vectorpair> outDimTagsMap;
    gmsh::model::occ::fragment(volumes, {}, outTags, outDimTagsMap);
    gmsh::model::occ::synchronize();
  }

  // Apply mesh size fields
  if (!meshParams->sizeFields.empty()) meshSizeFields();

  // Apply mesh color naming
  if (!meshParams->color2groupMap.empty()) applyColorNames();

  // Apply transfinite volumes
  if (!meshParams->transfiniteBlocks.empty()) applyTransfiniteVolumes();

  gmsh::model::mesh::generate(3);
  gmsh::model::mesh::removeDuplicateNodes();
  // By default, ouput MSH file
  mshFname = nemAux::trim_fname(fname, ".msh");
  gmsh::write(mshFname);

  std::string ext = nemAux::find_ext(meshParams->ofname);

  if (ext != ".msh" && ext != ".vtu") {
    auto &meshExtensions = gmshParams::getMeshExtensions();
    auto pos = std::find(meshExtensions.begin(), meshExtensions.end(), ext);
    if (pos == meshExtensions.end()) {
      std::cerr << "Error: Output file extension " << ext
                << " not supported in this mesh engine." << std::endl;
      std::cout << "Supported formats are ";
      for (auto meshExtension : meshExtensions) {
        std::cout << meshExtension << " ";
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
  gmsh::option::setNumber("Mesh.SaveAll", meshParams->saveAll);
  gmsh::option::setNumber("Mesh.RecombineAll", 0);
  gmsh::option::setNumber("Mesh.Smoothing", 1);
  gmsh::option::setNumber("Mesh.CharacteristicLengthFactor", 1);
  gmsh::option::setNumber("Mesh.CharacteristicLengthFromPoints", 1);
  gmsh::option::setNumber("Mesh.LcIntegrationPrecision", 1e-09);
  gmsh::option::setNumber("Mesh.MaxIterDelaunay3D", 0);
  gmsh::option::setNumber("Mesh.RandomFactor", 1e-09);
  gmsh::option::setNumber("Mesh.RandomFactor3D", 1e-12);
  gmsh::option::setNumber("Mesh.ElementOrder",
                          meshParams->elementOrder);
  gmsh::option::setNumber("Mesh.SubdivisionAlgorithm",
                          meshParams->subdivisionAlg);

  std::string algo = meshParams->algo2D;
  int a = ALGO_2D_FRONTAL;  // default: "Frontal"
  if (algo == "Frontal")
    a = ALGO_2D_FRONTAL;
  else if (algo == "MeshAdapt")
    a = ALGO_2D_MESHADAPT;
  else if (algo == "Automatic")
    a = ALGO_2D_AUTO;
  else if (algo == "Delaunay")
    a = ALGO_2D_DELAUNAY;
  else if (algo == "Frontal Quads" || algo == "Frontal Quad")
    a = ALGO_2D_FRONTAL_QUAD;
  else if (algo == "Packing of Parallelograms")
    a = ALGO_2D_PACK_PRLGRMS;
  else {
    std::cerr << "Warning: Surface algorithm " << algo << " not supported. \n"
              << "Using default Frontal algorithm." << std::endl;
  }
  gmsh::option::setNumber("Mesh.Algorithm", a);

  algo = meshParams->algo3D;
  a = ALGO_3D_DELAUNAY;
  if (algo == "Delaunay")
    a = ALGO_3D_DELAUNAY;
  else if (algo == "Frontal")
    a = ALGO_3D_FRONTAL;
  else if (algo == "HXT")
    a = ALGO_3D_HXT;
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
          std::cout << "thickness = " << thick << std::endl;
          break;
        }
      }
      for (auto &prm : sf.params) {
        std::string key = prm.first;
        int val = static_cast<int>(prm.second);
        if (key == "IField") {
          for (auto &sf2 : meshParams->sizeFields) {
            if (sf2.type == "Cylinder" && sf2.id == val) {
              double radius = 0;
              //double thick = 0;
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
                //if (prm2.first == "Thickness") thick = prm2.second;
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

void gmshGen::applyColorNames() {
  // Get all volumes
  std::vector<std::pair<int, int>> volumes;
  gmsh::model::getEntities(volumes, 3);

  std::vector<int> volumeTags;
  for(const std::pair<int,int> &v : volumes) {
    // int dim = v.first;
    int tag = v.second;
    volumeTags.push_back(tag);
  }
  gmsh::model::addPhysicalGroup(3, volumeTags);

  // Get all surfaces
  std::vector<std::pair<int, int>> surfaces;
  gmsh::model::getEntities(surfaces, 2);

  /*
   * 1. Find all surfaces with given color.
   * 2. Name each surface with unique name based on group name.
   *    e.g. if name is <group_name>, surfaces will be
   *    named <group_name>_0, <group_name>_1, ...
   * 3. Give physical group the assigned name <group_name>
   */

  // loop over color, group names and find constituent surfaces
  for(auto iter = meshParams->color2groupMap.begin();
      iter != meshParams->color2groupMap.end(); ++iter) {
    const auto &groupColor = iter->first;
    std::string groupColorStr = std::to_string(groupColor.at(0)) + "," +
                                std::to_string(groupColor.at(1)) + "," +
                                std::to_string(groupColor.at(2));
    std::string groupName = iter->second;
    std::vector<int> surfTags;

    int id = 0;
    for(const std::pair<int,int> &s : surfaces) {
      int dim = s.first;
      int tag = s.second;
      int r, g, b, a;
      gmsh::model::getColor(dim, tag, r, g, b, a);
      if (groupColor != std::array<int, 3>{r, g, b}) continue;
      std::string surfName = groupName + "_" + std::to_string(id);
      gmsh::model::setEntityName(dim, tag, surfName);
      surfTags.push_back(tag);
      ++id;
    }
    if(surfTags.size() == 0) {
      std::cout << "NO SURFACES WITH COLOR " << groupColorStr << " FOUND."
                << std::endl;
    } else {
      std::cout << "Found "
                << surfTags.size()
                << ((surfTags.size() == 1) ? " surface" : " surfaces")
                << " in group "
                << groupName
                << " mapped by color "
                << groupColorStr
                << std::endl;
      std::cout << "Adding physical group " << groupName << std::endl;
      int groupTag = gmsh::model::addPhysicalGroup(2, surfTags);
      gmsh::model::setPhysicalName(2, groupTag, groupName.c_str());
    }
  }
  gmsh::model::occ::synchronize();
}

void gmshGen::applyTransfiniteVolumes() {
  gmsh::vectorpair volumes;
  gmsh::model::getEntities(volumes, 3);

  for(const std::pair<int,int> &v : volumes) {
    // int volumeDim = v.first;
    int volumeTag = v.second;

    // Use lower_bound because TransfiniteBlock::operator< is based on id
    auto blockIter = std::lower_bound(
        meshParams->transfiniteBlocks.begin(),
        meshParams->transfiniteBlocks.end(), volumeTag,
        [](const TransfiniteBlock &block, int tag) { return block.id < tag; });
    // if not a transfinite volume, continue
    if (blockIter == meshParams->transfiniteBlocks.end() ||
        blockIter->id != volumeTag) {
      continue;
    }

    const auto &block = *blockIter;

    gmsh::vectorpair volume = { v };
    gmsh::vectorpair surfaces;
    gmsh::vectorpair curves;
    // get all entities (combined - false, oriented - false, recursive - false)
    gmsh::model::getBoundary(volume, surfaces, false, false, false);
    gmsh::model::getBoundary(surfaces, curves, false, false, false);

    // loop over curves c, ignore duplicates
    gmsh::vectorpair uniqueCurves;
    for(const std::pair<int,int> &c : curves) {
      auto iter = std::find(uniqueCurves.begin(), uniqueCurves.end(), c);
      if (iter == uniqueCurves.end()) {
        uniqueCurves.push_back(c);
      }
    }

    if(uniqueCurves.size() != 12) {
      std::cout << "Found " << uniqueCurves.size() << " unique curves in volume " << volumeTag <<
                ". There should be 12." << std::endl;
      std::cout << "Error : only hexahedral transfinite volumes are supported." << std::endl;
      throw;
    }

    std::cout << "Found transfinite block : " << volumeTag << std::endl;

    struct CurveData {
      int axis;
      std::pair<int,int> curve;
      bool reverse;
    };

    std::vector<CurveData> orientedCurves;

    for(const std::pair<int,int> &c : uniqueCurves) {
      // get points in curve, getBoundary requires vector of entities
      gmsh::vectorpair curves = { c };
      gmsh::vectorpair points;
      gmsh::model::getBoundary(curves, points, false, true, true);

      std::pair<int,int> a = points[0];
      std::pair<int,int> b = points[1];

      std::vector<double> a_vert, b_vert;
      gmsh::model::getValue(a.first, a.second, {}, a_vert);
      gmsh::model::getValue(b.first, b.second, {}, b_vert);

      // normalize
      double ab[] = { b_vert[0]-a_vert[0],
                      b_vert[1]-a_vert[1],
                      b_vert[2]-a_vert[2] };
      double ab_den = std::sqrt(ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2]);
      ab[0] = ab[0]/ab_den;
      ab[1] = ab[1]/ab_den;
      ab[2] = ab[2]/ab_den;

      // get corresponding direction
      // i=0 -> x-axis
      // i=1 -> y-axis
      // i=2 -> z-axis
      for(int axis = 0; axis < 3; ++axis) {
        double dot = ab[0]*block.axis[axis][0] +
                     ab[1]*block.axis[axis][1] +
                     ab[2]*block.axis[axis][2];
        double theta = std::acos(dot);
        if(theta < M_PI/4 || theta > 3*M_PI/4) {
          bool reverse = theta >= M_PI / 4;
          CurveData data;
          data.axis = axis;
          data.curve = c;
          data.reverse = reverse;
          orientedCurves.push_back(data);
        }
      }
    }

    for(auto &oc : orientedCurves) {
      // curve data
      std::pair<int,int> curve = oc.curve;
      bool reverse = oc.reverse;

      // transfinite data
      int vert = block.vert[oc.axis];
      double coef = block.coef[oc.axis];
      std::string type = block.type[oc.axis];

      if(type == "Progression") {
        // take reciprocal of coefficient if curve is "reversed" w.r.t axis
        coef = reverse ? 1/coef : coef;
      }

      gmsh::model::mesh::setTransfiniteCurve(curve.second, vert, type, coef);
    }

    std::cout << "Setting transfinite surfaces..." << std::endl;
    for(const std::pair<int,int> &s : surfaces) {
      gmsh::model::mesh::setTransfiniteSurface(s.second);
      gmsh::model::mesh::setRecombine(s.first, s.second);
    }

    std::cout << "Setting transfinite volume..." << std::endl;
    gmsh::model::mesh::setTransfiniteVolume(v.second);
  }
}

}  // namespace GEN

}  // namespace NEM
