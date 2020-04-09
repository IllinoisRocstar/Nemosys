#define _USE_MATH_DEFINES
#include "NucMeshDriver.H"

#include <gmsh.h>
#include <iomanip>

#include "AuxiliaryFunctions.H"
#include "circles.H"
#include "polygon.H"

nemAux::Timer Tgeom, Tboolean, Textrude, Tmesh, Tpause, Tconserve;

NucMeshDriver::~NucMeshDriver() {
  std::cout << "NucMeshDriver destroyed" << std::endl;
}

NucMeshDriver *NucMeshDriver::readJSON(const jsoncons::json &inputjson) {
  std::cout << "Reading Input JSON File." << std::endl;

  if (inputjson.contains("Geometry and Mesh")) {
    auto nucmeshobj = new NucMeshDriver(inputjson);
    return nucmeshobj;
  } else {
    std::cerr << "Error: 'Geometry and Mesh' keyword not found, expected after "
                 "'Output File Name'"
              << std::endl;
    exit(-1);
  }
}

NucMeshDriver::NucMeshDriver(jsoncons::json inputjson) {
  id = 1;
  physTag = 1;
  gui = false;

  std::cout << "NucMeshDriver created" << std::endl;

  if (inputjson.contains("Output File Name")) {
    ofname = inputjson["Output File Name"].as_string();
  } else {
    std::cerr << "Error: 'Output File Name' keyword not found, expected after "
                 "'Program Type'"
              << std::endl;
    exit(-1);
  }

  jsoncons::json shapes = inputjson["Geometry and Mesh"];

  Tgeom.start();
  gmsh::initialize();

  gmsh::model::add("NucMesh");
  gmsh::option::setNumber("General.NumThreads", 10);
  gmsh::option::setNumber("Geometry.OCCBooleanPreserveNumbering", 1);
  gmsh::option::setNumber("Geometry.OCCParallel", 1);
  gmsh::option::setNumber("Mesh.FlexibleTransfinite", 0);
  gmsh::option::setNumber("Mesh.MshFileVersion", 2.2);
  gmsh::option::setNumber("Mesh.CharacteristicLengthFromCurvature", 0);

  //--------------------------------------------------------------------//
  // Main JSON file geometry parsing
  //--------------------------------------------------------------------//
  parseGeomAndMesh(shapes);
  Tgeom.stop();
  std::cout << gui << std::endl;
  if (gui)
    gmsh::fltk::initialize();

  //--------------------------------------------------------------------//
  // Perform Boolean Fragments
  //--------------------------------------------------------------------//
  Tboolean.start();
  std::vector<std::pair<int, int>> surf;
  std::vector<std::pair<int, int>> outSurf;
  std::vector<std::vector<std::pair<int, int>>> outSurfMap;
  // Get all surfaces
  gmsh::model::getEntities(surf, 2);
  // Boolean Fragments
  gmsh::model::occ::fragment(surf, surf, outSurf, outSurfMap, -1, true, true);
  gmsh::model::occ::synchronize();

  int oldID, newID;
  std::pair<int, int> oldNew;                  // holds the old/new id pair
  std::vector<std::pair<int, int>> oldNew_vec; // vector for id pair

  //--------------------------------------------------------------------//
  // Gather surface ids that have changed for Boolean Fragments
  //--------------------------------------------------------------------//
  for (auto i = 0; i < outSurfMap.size() / 2; ++i) {
    int size = outSurfMap[i].size();
    if (size > 1) {
      oldID = surf[i].second;
      newID = outSurfMap[i][0].second;
      if (newID != oldID) {
        std::cout << oldID << " became " << newID << std::endl;
        oldNew.first = oldID;
        oldNew.second = newID;
        oldNew_vec.push_back(oldNew);
      }
    }
  }

  //--------------------------------------------------------------------//
  // Update all surfaces with new ID's
  //--------------------------------------------------------------------//
  for (auto itr = shape_map.begin(); itr != shape_map.end(); ++itr)
    itr->second->updateSurfaces(oldNew_vec);

  std::cout << "Updated surface IDs\n";

  //--------------------------------------------------------------------//
  // Apply meshType and get physical surfaces
  //--------------------------------------------------------------------//
  // iterate through shape_map to apply meshType and collect physical Surfs
  std::map<int, int> physSurf_map;
  for (auto &&itr : shape_map) {
    itr.second->applyMeshType();
    physSurf_map = itr.second->getPhysSurf(phystag_map, physSurf_map);
  }

  std::cout << "Mesh type applied\n";

  //--------------------------------------------------------------------//
  // Make physical groups for surfaces
  //--------------------------------------------------------------------//
  for (const auto &itr : phystag_map) {
    std::vector<int> s; // vector of surfaces
    std::string name;   // regions name
    int tag = -1;       // physical group tag
    for (const auto &itr2 : physSurf_map) {
      if (itr2.second == itr.second) {
        tag = itr2.second;
        name = itr.first;
        s.push_back(itr2.first);
      }
    }
    if (tag == -1) {
      std::cerr << "Error: Physical group cannot be applied." << std::endl;
      exit(-1);
    }
    if (!s.empty()) {
      gmsh::model::addPhysicalGroup(2, s, tag);
      gmsh::model::setPhysicalName(2, tag, name);
    }
  }
  std::cout << "2D Physical groups applied\n";
  gmsh::model::occ::synchronize();
  Tboolean.stop();

  //--------------------------------------------------------------------//
  // Extrude Surfaces if '3D' option is set to true
  //--------------------------------------------------------------------//
  // Extrude surfaces based on Height
  if (extrude) {
    Textrude.start();
    std::map<int, std::vector<int>> volPhysTag_map;
    std::vector<std::pair<int, int>> entities, outEntities, volumes,
        volSurfaces;
    std::vector<int> pTag;
    // Collect all surfaces
    gmsh::model::getEntities(entities, 2);
    // Normalize the Heights vector for gmsh
    double last = heights.back();
    if (last == 0.0) {
      std::cerr << "Error: The last value in 'Heights' vector cannot be 0. "
                   "Cannot perform extrusion."
                << std::endl;
      exit(-1);
    } else {
      for (int i = 0; i < heights.size(); i++) {
        heights[i] = heights[i] / last;
      }
      // Extrude all surfaces
      gmsh::model::occ::extrude(entities, 0, 0, last, outEntities, layers,
                                heights, true);
      gmsh::model::occ::synchronize();
      std::cout << "Extruded to 3D\n";
      Textrude.stop();
    }

    //--------------------------------------------------------------------//
    // For each new volume, go through it's surfaces and obtain physical
    // tag. Put physical tag and volume ID into 'volPhysTag_map'.
    //--------------------------------------------------------------------//
    // Gather the volumes
    gmsh::model::getEntities(volumes, 3);
    // Loop through volumes and get surfaces
    for (int i = 0; i < (int)volumes.size(); ++i) {
      volSurfaces.clear();
      gmsh::model::getBoundary({volumes[i]}, volSurfaces, false, false, false);
      // Loop through surfaces and find physical tags
      for (int j = 0; j < (int)volSurfaces.size(); ++j) {
        gmsh::model::getPhysicalGroupsForEntity(volSurfaces[j].first,
                                                volSurfaces[j].second, pTag);
        // If phys tag exists, populate vol phys tag map with tags and volumes
        if (pTag.size() > 0) {
          auto search = volPhysTag_map.find(pTag[0]);
          if (search == volPhysTag_map.end()) {
            volPhysTag_map.insert(
                std::pair<int, std::vector<int>>(pTag[0], {volumes[i].second}));
            continue;
          } else {
            search->second.push_back(volumes[i].second);
            continue;
          }
        }
      }
    }
    volumes.clear();
    entities.clear();
    outEntities.clear();
    volSurfaces.clear();

    //--------------------------------------------------------------------//
    // Apply physical tags/names from original surfaces to corresponding
    // extruded volume. Then remove 2D/surface physical tags.
    //--------------------------------------------------------------------//
    std::string name;
    // Iterate through vol phys tag map and apply volume physical groups
    for (const auto &itr : volPhysTag_map) {
      // Add physical group
      int tag = 0;
      tag = gmsh::model::addPhysicalGroup(3, itr.second, -1);
      // Get physical name for physical tag of surface
      gmsh::model::getPhysicalName(2, itr.first, name);
      // Set physical name of Volume
      gmsh::model::setPhysicalName(3, tag, name);
    }
    std::cout << "3D Physical groups applied\n";

    // Remove the surface physical tag
    for (const auto &itr : volPhysTag_map) {
      gmsh::model::removePhysicalGroups({{2, itr.first}});
    }
    std::cout << "2D Physical groups removed\n";
    gmsh::model::occ::synchronize();
  } else {
    //--------------------------------------------------------------------//
    // For 2D/Surface case we need all the mesh normals to be pointing
    // the same direction in Proteus. Thus, we can use the
    //'setOutwardOrientation'method in Gmsh. However, this method is only
    // used on volumes. So we extrude the surfaces and set the mesh
    // orientation. When the mesh is saved, only phyiscal entities will be
    // saved, that is, only the 2D/surfaces.
    //--------------------------------------------------------------------//
    std::vector<std::pair<int, int>> entities, outEntities, volumes;
    // Collect all surfaces
    gmsh::model::getEntities(entities, 2);
    // Extrude all surfaces
    gmsh::model::occ::extrude(entities, 0, 0, -0.1, outEntities, {1}, {1},
                              false);
    gmsh::model::occ::synchronize();
    // Collect all volumes
    gmsh::model::getEntities(volumes, 3);
    int volTag;
    // Set the mesh orientation for those volumes
    for (int i = 0; i < volumes.size(); i++) {
      volTag = volumes[i].second;
      gmsh::model::mesh::setOutwardOrientation(volTag);
    }
    std::cout << "Surface mesh orientation set\n";
  }

  // Output region ids and names
  std::cout << "Region   Name" << std::endl;
  for (const auto &itr : phystag_map) {
    std::cout << "  " << itr.second << "       " << itr.first << "\n"
              << std::endl;
  }

  // Apply size field to certain physical group
  //	int tag = phystag_map["C"];
  //	std::vector<int> surf_tags_int;
  //	std::vector<double> surf_tags;
  //	std::vector<std::pair<int,int>> surfs;
  //	std::vector<std::pair<int,int>> lines;
  //	gmsh::model::getEntitiesForPhysicalGroup(2,tag,surf_tags_int);
  //
  //	for(int i = 0; i < surf_tags_int.size(); i++){
  //		surfs.push_back({2,surf_tags_int[i]});
  //		surf_tags.push_back(surf_tags_int[i]*1.0);
  //	}
  //
  //	gmsh::model::getBoundary(surfs, lines, false, false, false);
  //
  //	std::vector<double> lin;
  //	for (int i = 0; i < lines.size(); i++){
  //		lin.push_back(lines[i].second*1.0);
  //	}
  //
  //	gmsh::model::mesh::field::add("Distance", 1);
  //	gmsh::model::mesh::field::setNumber(1, "NNodesByEdge", 100);
  //	gmsh::model::mesh::field::setNumbers(1, "EdgesList", lin);
  //
  //	gmsh::model::mesh::field::add("Threshold", 2);
  //  gmsh::model::mesh::field::setNumber(2, "IField", 1);
  //  gmsh::model::mesh::field::setNumber(2, "LcMin", 0.1);
  //  gmsh::model::mesh::field::setNumber(2, "LcMax", 0.4);
  //  gmsh::model::mesh::field::setNumber(2, "DistMin", 1);
  //  gmsh::model::mesh::field::setNumber(2, "DistMax", 3);
  //
  //  gmsh::model::mesh::field::add("Restrict", 3);
  //  gmsh::model::mesh::field::setNumber(3, "IField", 2);
  //  gmsh::model::mesh::field::setNumbers(3, "FacesList", surf_tags);
  //  gmsh::model::mesh::field::setNumbers(3, "EdgesList", lin);
  //
  //  gmsh::model::mesh::field::setAsBackgroundMesh(3);
  //  gmsh::model::occ::synchronize();

  if (gui)
    openGUI();

  std::cout << "Meshing...\n" << std::endl;
  Tmesh.start();
  // Mesh 1D to get wireframe, save as vtk
  gmsh::model::mesh::generate(1);
  gmsh::option::setNumber("Mesh.SaveAll", 1);
  std::cout << "Writing line mesh to " << ofname << "_lines.vtk" << std::endl;
  gmsh::write(ofname + ".vtk");

  // Mesh 2D, save as msh and vtk
  if (!extrude) {
    gmsh::model::mesh::generate(2);
    gmsh::model::mesh::removeDuplicateNodes();
    gmsh::option::setNumber("Mesh.SaveAll", 0);
    std::cout << "Writing surface mesh to " << ofname << ".vtk" << std::endl;
    gmsh::write(ofname + ".msh");
    gmsh::write(ofname + ".vtk");
  }

  // Mesh 3D, save as msh and vtk
  if (extrude) {
    gmsh::model::mesh::generate(3);
    gmsh::option::setNumber("Mesh.SaveAll", 0);
    gmsh::write(ofname + ".msh");
    gmsh::write(ofname + ".vtk");
  }
  std::cout << "Meshing...Complete\n" << std::endl;
  Tmesh.stop();

  double s = 1000.0;
  std::cout << std::setprecision(3) << std::fixed
            << "Mesh conservation time = " << (Tconserve.elapsed() / s)
            << " seconds\n";
  std::cout << "Geometry creation time = "
            << (Tgeom.elapsed() - Tconserve.elapsed()) / s << " seconds\n";
  std::cout << "Geometry coherence time = " << Tboolean.elapsed() / s
            << " seconds\n";
  if (extrude)
    std::cout << "Extrusion time = " << Textrude.elapsed() / s << " seconds\n";
  std::cout << "Meshing time = " << Tmesh.elapsed() / s << " seconds\n";
  std::cout << "Total run time = "
            << (Tgeom.elapsed() + Tboolean.elapsed() + Textrude.elapsed() +
                Tmesh.elapsed()) /
                   s
            << " seconds\n"
            << std::endl;

  // gmsh::fltk::run();
  if (gui)
    openGUI();
  gmsh::finalize();
}

//********************************************************************//
// Parses the Geometry and Mesh main section
//********************************************************************//
void NucMeshDriver::parseGeomAndMesh(jsoncons::json shapes) {
  // Loop through Geometry and Mesh to get shapes
  std::cout << "Geometry and Mesh" << std::endl;
  for (const auto &obj : shapes.array_range()) {
    if (obj.contains("Global Options"))
      parseOptions(obj["Global Options"]);

    if (obj.contains("Saved Objects"))
      parseSavedObjects(obj["Saved Objects"]);

    if (obj.contains("Visible"))
      if (obj["Visible"] == false)
        continue;

    if (obj.contains("Show"))
      if (obj["Show"] == false)
        continue;

    if (obj.contains("Name"))
      std::cout << "Creating '" << obj["Name"].as_string() << "'" << std::endl;
    // if (obj.contains("Include")) parseGeomAndMesh(obj["Include"]);

    if (obj.contains("Circles"))
      makeCircles(obj["Circles"]);

    if (obj.contains("Circle"))
      makeCircles(obj["Circle"]);

    if (obj.contains("Hexagon") || obj.contains("Polygon")) {
      if (obj.contains("Hexagon"))
        makePolygons(obj["Hexagon"], 6);
      else
        makePolygons(obj["Polygon"]);
    }
    if (obj.contains("BREAK"))
      if (obj["BREAK"] == true)
        openGUI();

    if (obj.contains("Array"))
      makeArray(obj);
  }
}
//********************************************************************//
// Parses the global options data of json
//********************************************************************//
void NucMeshDriver::parseOptions(jsoncons::json opts) {
  // Iterate through options
  for (const auto &it : opts.array_range()) {
    if (it.contains("Open GUI"))
      gui = it["Open GUI"].as_bool();
    if (it.contains("3D")) {
      extrude = it["3D"].as_bool();
      if (extrude == true) {
        int numLayersSize;
        if (it.contains("Number of Layers")) {
          numLayersSize = it["Number of Layers"].size();
          if (numLayersSize > 0) {
            for (int i = 0; i < numLayersSize; ++i) {
              int layerValue = it["Number of Layers"][i].as<int>();
              if (layerValue == 0) {
                std::cerr << "Error: Value of 0 (zero) found in 'Number of "
                             "Layers' array."
                          << std::endl;
                exit(-1);
              } else
                layers.push_back(layerValue);
            }
          } else {
            std::cerr << "Error: 'Number of Layers' array is empty"
                      << std::endl;
            exit(-1);
          }
        } else {
          std::cerr << "Error: 'Number of Layers' keyword not found"
                    << std::endl;
          exit(-1);
        }
        if (it.contains("Heights")) {
          if (it["Heights"].size() > numLayersSize) {
            std::cerr
                << "Error: 'Heights' vector larger than 'Number of Layers'";
            exit(-1);
          }
          for (int i = 0; i < it["Heights"].size(); ++i) {
            heights.push_back(it["Heights"][i].as_double());
          }
        }
        gmsh::option::setNumber("Mesh.SurfaceEdges", 0);
        gmsh::option::setNumber("Mesh.SurfaceFaces", 0);
        gmsh::option::setNumber("Mesh.VolumeEdges", 1);
        gmsh::option::setNumber("Mesh.VolumeFaces", 1);
      }

    } else {
      std::cerr << "'3D' keyword not found" << std::endl;
      exit(-1);
    }

    if (it.contains("Min Mesh Size"))
      gmsh::option::setNumber("Mesh.CharacteristicLengthMin",
                              it["Min Mesh Size"].as<double>());

    if (it.contains("Max Mesh Size"))
      gmsh::option::setNumber("Mesh.CharacteristicLengthMax",
                              it["Max Mesh Size"].as<double>());

    if (it.contains("Meshing Algorithm")) {
      std::string algo = it["Meshing Algorithm"].as<std::string>();
      int a = 6; // default: "Frontal"
      if (algo == "Frontal")
        a = 6;
      else if (algo == "MeshAdapt")
        a = 1;
      else if (algo == "Automatic")
        a = 2;
      else if (algo == "Delaunay")
        a = 5;
      else if (algo == "Frontal Quads" || algo == "Frontal Quad")
        a = 8;
      else if (algo == "Packing of Parallelograms")
        a = 9;
      gmsh::option::setNumber("Mesh.Algorithm", a);
    }
    if (it.contains("Recombine Algorithm")) {
      std::string combine = it["Recombine Algorithm"].as<std::string>();
      int a = 1; // default: "Blossom"
      if (combine == "Simple")
        a = 0;
      else if (combine == "Blossom")
        a = 1;
      else if (combine == "Simple Full-Quad")
        a = 2;
      else if (combine == "Blossom Full-Quad")
        a = 3;
      gmsh::option::setNumber("Mesh.RecombinationAlgorithm", a);
    }
    if (it.contains("Mesh Smoothing Iters"))
      gmsh::option::setNumber("Mesh.Smoothing",
                              it["Mesh Smoothing Iters"].as<int>());
    if (it.contains("Extend from Boundary")) {
      bool extend_from_bnd = it["Extend from Boundary"].as_bool();
      int e = 1;
      if (!extend_from_bnd)
        e = 0;
      gmsh::option::setNumber("Mesh.CharacteristicLengthExtendFromBoundary", e);
    }
  }
}

//********************************************************************//
// Parses the saved objects data of json
//********************************************************************//
void NucMeshDriver::parseSavedObjects(jsoncons::json savedObj) {
  std::cout << "Parsing Saved Objects..." << std::endl;

  // Iterate through Saved Objects key
  for (const auto &it : savedObj.array_range()) {
    std::string alias;
    std::string shapetype;

    if (it.contains("Alias")) {
      alias = it["Alias"].as<std::string>();
      std::cout << "Alias '" << alias << "' found\n";
    } else {
      std::cerr << "Error: 'Alias' keyword not found in 'Saved Objects'"
                << std::endl;
      exit(-1);
    }

    if (it.contains("Circles"))
      shapetype = "Circles";
    if (it.contains("Circle"))
      shapetype = "Circle";
    if (it.contains("Polygon"))
      shapetype = "Polygon";

    bool visible = true;
    // If meshed area conservation is required and visible, change the radii
    if (it.contains("Conserve Area") && it["Conserve Area"] == true) {
      for (const auto &c : it["Circles"].array_range()) {
        if (c.contains("Visible")) visible = !(c["Visible"] == false);
      }
      jsoncons::json j;
      if (visible) {
        double tolerance = 1e-8;
        if (it.contains("Tolerance"))
          tolerance = it["Tolerance"].as_double();
        j = correctMeshArea(it[shapetype], shapetype, tolerance);
      } else
        j = it[shapetype];

      savedobj_map.insert(
          std::pair<std::string, std::pair<std::string, jsoncons::json>>(
              alias, {shapetype, j}));
    } else {
      savedobj_map.insert(
          std::pair<std::string, std::pair<std::string, jsoncons::json>>(
              alias, {shapetype, it[shapetype]}));
    }
    std::cout << "Saved Objects Parsed\n" << std::endl;
    ;
  }
}

//********************************************************************//
// Calculates the area of meshed shape and adjusts the radii
//********************************************************************//
jsoncons::json NucMeshDriver::correctMeshArea(jsoncons::json obj,
                                              std::string shapetype,
                                              double tolerance) {
  std::cout << "\nConserving shape meshed area " << std::endl;
  Tconserve.start();

  // start a new model to not disrupt main model
  gmsh::model::add("Preserve");

  std::vector<double> radii;
  double lastRadius = 0.0, area;
  int nRad;

  if (shapetype == "Circles" || shapetype == "Circle") {
    if (obj[0].contains("Radii")) {
      nRad = obj[0]["Radii"].size();
      if (nRad <= 0) {
        std::cerr << "Error: No 'Radii' values found for aliased object."
                  << std::endl;
        exit(-1);
      }
      for (int i = 0; i < nRad; ++i) {
        radii.push_back(obj[0]["Radii"][i].as_double());
        lastRadius = obj[0]["Radii"][i].as_double();
      }
    } else {
      std::cerr << "Error: No 'Radii' found for aliased object." << std::endl;
      exit(-1);
    }

    double eps = tolerance;
    double residual = 10;
    std::vector<double> r_a = radii;
    std::vector<double> r_b;
    std::vector<double> r_c(nRad);
    for (int i = 0; i < nRad; i++) {
      r_b.push_back(radii[i] * 1.5);
    }
    area = M_PI * lastRadius * lastRadius;
    double /*A_b, */A_c;
    int iter = 0;

    while (residual > eps) {
      iter++;
      for (int i = 0; i < nRad; i++) {
        r_c[i] = ((r_a[i] + r_b[i]) / 2);
      }

      obj[0].insert_or_assign("Radii", r_c);
      makeCircles(obj, true);
      gmsh::model::mesh::generate(2);

      gmsh::plugin::setNumber("MeshVolume", "Dimension", 2);
      gmsh::plugin::run("MeshVolume");

      std::vector<int> viewTags;
      gmsh::view::getTags(viewTags);
      if (viewTags.size() == 1) {
        std::vector<std::vector<double>> data;
        std::vector<std::string> dataTypes;
        std::vector<int> numElements;
        gmsh::view::getListData(viewTags[0], dataTypes, numElements, data);
        A_c = data[0][3];
      } else {
        std::cerr << "Error: View data not found." << std::endl;
        exit(-1);
      }
      viewTags.clear();
      gmsh::clear();

      // residual = std::abs(r_b[nRad-1]-r_a[nRad-1]);
      residual = std::abs(A_c - area);
      std::cout << std::scientific << std::setprecision(4)
                << "Residual = " << residual << std::endl;

      if (A_c - area > 0)
        r_b = r_c;
      else
        r_a = r_c;
    }
    std::cout << std::setprecision(8) << std::fixed
              << "New outer radius = " << r_b[nRad - 1] << std::endl;
    std::cout << "Mesh area converged in " << iter << " iterations"
              << std::endl;
  } else {
    std::cout << "Warning: Mesh area conservation only enabled for circles."
              << std::endl;
  }
  Tconserve.stop();
  return obj;
}

//********************************************************************//
// Parses json for polygon data then creates polygon object
//********************************************************************//
void NucMeshDriver::makePolygons(jsoncons::json poly, int ns) {
  // Iterate through Polygon/Hexagon key
  for (const auto &it : poly.array_range()) {
    bool visible = true;
    if (it.contains("Visible")) {
      if (it["Visible"] == true)
        visible = true;
      else
        visible = false;
    }
    if (visible) {
      if (it.contains("Name"))
        std::cout << "Creating '" << it["Name"].as_string() << "'" << std::endl;
      // If there is an alias check for new data under the alias
      if (it.contains("Alias")) {
        // Get the alias name
        std::string alias = it["Alias"].as_string();

        std::vector<double> orig_center, orig_radii;
        std::vector<std::string> orig_meshType, orig_names;
        std::vector<std::pair<int, int>> orig_elems;
        double orig_rot;
        // Get the original alias data
        auto find = savedobj_map.find(alias);
        if (find != savedobj_map.end()) {
          // the alias was found in map
          if (find->second.first == "Polygon" ||
              find->second.first == "Hexagon") {
            jsoncons::json j = find->second.second[0];
            if (j.contains("Center")) {
              orig_center.push_back(j["Center"][0].as<double>());
              orig_center.push_back(j["Center"][1].as<double>());
              orig_center.push_back(j["Center"][2].as<double>());
              if (j["Center"].size() == 5) {
                orig_center.push_back(j["Center"][3].as<double>());
                orig_center.push_back(j["Center"][4].as<double>());
              }
            }
            if (j.contains("Rotation")) {
              orig_rot = j["Rotation"].as<double>();
            }
            if (j.contains("Radii")) {
              int nRad = j["Radii"].size();
              for (int i = 0; i < nRad; ++i)
                orig_radii.push_back(j["Radii"][i].as<double>());
            }
            if (j.contains("Mesh Type")) {
              int nType = j["Mesh Type"].size();
              for (int i = 0; i < nType; ++i)
                orig_meshType.push_back(j["Mesh Type"][i].as<std::string>());
            }
            if (j.contains("Number of Elems")) {
              int nE = j["Number of Elems"].size();
              for (int i = 0; i < nE; ++i) {
                std::pair<int, int> p;
                p.first = j["Number of Elems"][i][0].as<int>();
                p.second = j["Number of Elems"][i][1].as<int>();
                orig_elems.push_back(p);
              }
            }
            if (j.contains("Region Names")) {
              int nNames = j["Region Names"].size();
              for (int i = 0; i < nNames; ++i) {
                orig_names.push_back(j["Region Names"][i].as<std::string>());
              }
            }
          }
        }

        // Check if any new parameters have been defined for alias
        std::vector<double> cen;
        if (it.contains("Center")) {
          // Get the center coords
          cen.push_back(it["Center"][0].as<double>());
          cen.push_back(it["Center"][1].as<double>());
          cen.push_back(it["Center"][2].as<double>());
          if (it["Center"].size() == 5) {
            cen.push_back(it["Center"][3].as<double>());
            cen.push_back(it["Center"][4].as<double>());
          }
        }

        double rot;
        if (it.contains("Rotation")) {
          // Get the rotation
          rot = it["Rotation"].as<double>();
        }

        std::vector<double> rad;
        if (it.contains("Radii")) {
          // Number of radii
          int nRad = it["Radii"].size();

          // Loop through radii
          for (int i = 0; i < nRad; ++i)
            rad.push_back(it["Radii"][i].as<double>());
        }

        std::vector<std::string> type;
        if (it.contains("Mesh Type")) {
          // Number of mesh types
          int nType = it["Mesh Type"].size();

          // Loop through mesh types
          for (int i = 0; i < nType; ++i)
            type.push_back(it["Mesh Type"][i].as<std::string>());
        }

        std::vector<std::pair<int, int>> elems;
        if (it.contains("Number of Elems")) {
          // Number of element pairs
          int nE = it["Number of Elems"].size();

          // Loop through element pairs
          for (int i = 0; i < nE; ++i) {
            std::pair<int, int> p;
            p.first = it["Number of Elems"][i][0].as<int>();
            p.second = it["Number of Elems"][i][1].as<int>();

            elems.push_back(p);
          }
        }

        std::vector<std::string> name;
        if (it.contains("Region Names")) {
          // Number of region names
          int nNames = it["Region Names"].size();
          // Loop through regions names
          for (int i = 0; i < nNames; ++i) {
            name.push_back(it["Region Names"][i].as<std::string>());

            // check if name is in physical name map
            auto iter = phystag_map.find(name[i]);
            // if name is not in map, add it
            if (iter == phystag_map.end()) {
              phystag_map.insert({name[i], physTag});
              physTag++;
            }
          }
        }

        auto search = savedobj_map.find(alias);
        if (search != savedobj_map.end()) {
          // the alias was found in map
          if (search->second.first == "Polygon" ||
              search->second.first == "Hexagon") {
            if (it.contains("Center"))
              search->second.second[0].insert_or_assign("Center", cen);
            if (it.contains("Rotation"))
              search->second.second[0].insert_or_assign("Rotation", rot);
            if (it.contains("Radii"))
              search->second.second[0].insert_or_assign("Radii", rad);
            if (it.contains("Mesh Type"))
              search->second.second[0].insert_or_assign("Mesh Type", type);
            if (it.contains("Number of Elems"))
              search->second.second[0].insert_or_assign("Number of Elems",
                                                        elems);
            if (it.contains("Region Names"))
              search->second.second[0].insert_or_assign("Region Names", name);
            if (it.contains("BREAK"))
              if (it["BREAK"] == true)
                search->second.second[0].insert_or_assign("BREAK", true);
            makePolygons(search->second.second);

            // Reset the alias map back to original data
            search->second.second[0].insert_or_assign("Center", orig_center);
            search->second.second[0].insert_or_assign("Rotation", orig_rot);
            search->second.second[0].insert_or_assign("Radii", orig_radii);
            search->second.second[0].insert_or_assign("Mesh Type",
                                                      orig_meshType);
            search->second.second[0].insert_or_assign("Number of Elems",
                                                      orig_elems);
            search->second.second[0].insert_or_assign("Region Names",
                                                      orig_names);
            search->second.second[0].insert_or_assign("BREAK", false);
          }
        } else {
          std::cerr << alias << " not found in map" << std::endl;
        }
      } else {
        int nsides = 0;
        if (ns == 6)
          nsides = ns;
        else
          // Get the number of sides
          nsides = it["Number of Sides"].as<int>();

        // Get the center coords
        std::vector<double> cen;
        cen.push_back(it["Center"][0].as<double>());
        cen.push_back(it["Center"][1].as<double>());
        cen.push_back(it["Center"][2].as<double>());
        if (it["Center"].size() == 5) {
          cen.push_back(it["Center"][3].as<double>());
          cen.push_back(it["Center"][4].as<double>());
        }

        // Number of radii
        int nRad = it["Radii"].size();

        // Loop through radii
        std::vector<double> rad;
        for (int i = 0; i < nRad; ++i)
          rad.push_back(it["Radii"][i].as<double>());

        // Number of mesh types
        int nType = it["Mesh Type"].size();

        // Loop through mesh types
        std::vector<std::string> type;
        for (int i = 0; i < nType; ++i)
          type.push_back(it["Mesh Type"][i].as<std::string>());

        // Number of element pairs
        int nE = it["Number of Elems"].size();

        // Loop through element pairs
        std::vector<std::pair<int, int>> elems;
        for (int i = 0; i < nE; ++i) {
          std::pair<int, int> p;
          p.first = it["Number of Elems"][i][0].as<int>();
          p.second = it["Number of Elems"][i][1].as<int>();

          elems.push_back(p);
        }

        // Number of region names
        int nNames = it["Region Names"].size();

        // Loop through regions names
        std::vector<std::string> name;
        for (int i = 0; i < nNames; ++i) {
          name.push_back(it["Region Names"][i].as<std::string>());

          // check if name is in physical name map
          auto iter = phystag_map.find(name[i]);
          if (iter == phystag_map.end()) {
            // std::cout << "checking phys tags in map" << std::endl;
            phystag_map.insert({name[i], physTag});
            physTag++;
          }
        }

        // Get rotation angle
        auto rot = it["Rotation"].as<double>();

        NEM::GEO::shape *p =
            new NEM::GEO::polygon(nsides, cen, rad, type, elems, name, rot);
        p->draw();

        shape_map.insert(std::pair<int, NEM::GEO::shape *>(id, p));
        id++;

        nsides = 0;
        type.clear();
        elems.clear();
      }
      if (it.contains("BREAK")) {
        if (it["BREAK"] == true)
          openGUI();
      }
    }
  }
}

//********************************************************************//
// Parses json for circle data then creates circle object
//********************************************************************//
void NucMeshDriver::makeCircles(jsoncons::json circ, bool conserving) {
  // Iterate through Circles key
  for (const auto &it : circ.array_range()) {
    bool visible = true;
    if (it.contains("Visible")) {
      if (it["Visible"] == true)
        visible = true;
      else
        visible = false;
    }

    if (visible) {
      if (it.contains("Name"))
        std::cout << "Creating '" << it["Name"].as_string() << "'" << std::endl;

      if (it.contains("Alias")) {
        std::string alias = it["Alias"].as_string();

        std::vector<double> orig_center, orig_radii;
        std::vector<std::string> orig_meshType, orig_names;
        std::vector<std::pair<int, int>> orig_elems;
        double orig_rot;
        // Get the original alias data
        auto find = savedobj_map.find(alias);
        if (find != savedobj_map.end()) {
          // the alias was found in map
          if (find->second.first == "Circle" ||
              find->second.first == "Circles") {
            jsoncons::json j = find->second.second[0];
            if (j.contains("Center")) {
              orig_center.push_back(j["Center"][0].as<double>());
              orig_center.push_back(j["Center"][1].as<double>());
              orig_center.push_back(j["Center"][2].as<double>());
              if (j["Center"].size() == 5) {
                orig_center.push_back(j["Center"][3].as<double>());
                orig_center.push_back(j["Center"][4].as<double>());
              }
            }
            if (j.contains("Rotation")) {
              orig_rot = j["Rotation"].as<double>();
            }
            if (j.contains("Radii")) {
              int nRad = j["Radii"].size();
              for (int i = 0; i < nRad; ++i)
                orig_radii.push_back(j["Radii"][i].as<double>());
            }
            if (j.contains("Mesh Type")) {
              int nType = j["Mesh Type"].size();
              for (int i = 0; i < nType; ++i)
                orig_meshType.push_back(j["Mesh Type"][i].as<std::string>());
            }
            if (j.contains("Number of Elems")) {
              int nE = j["Number of Elems"].size();
              for (int i = 0; i < nE; ++i) {
                std::pair<int, int> p;
                p.first = j["Number of Elems"][i][0].as<int>();
                p.second = j["Number of Elems"][i][1].as<int>();
                orig_elems.push_back(p);
              }
            }
            if (j.contains("Region Names")) {
              int nNames = j["Region Names"].size();
              for (int i = 0; i < nNames; ++i) {
                orig_names.push_back(j["Region Names"][i].as<std::string>());
              }
            }
          }
        }

        // Check if any new parameters have been defined for alias
        std::vector<double> cen;
        if (it.contains("Center")) {
          // Get the center coords
          cen.push_back(it["Center"][0].as<double>());
          cen.push_back(it["Center"][1].as<double>());
          cen.push_back(it["Center"][2].as<double>());
          if (it["Center"].size() == 5) {
            cen.push_back(it["Center"][3].as<double>());
            cen.push_back(it["Center"][4].as<double>());
          }
        }

        std::vector<double> rad;
        if (it.contains("Radii")) {
          // Number of radii
          int nRad = it["Radii"].size();

          // Loop through radii
          for (int i = 0; i < nRad; ++i)
            rad.push_back(it["Radii"][i].as<double>());
        }

        std::vector<std::string> type;
        if (it.contains("Mesh Type")) {
          // Number of mesh types
          int nType = it["Mesh Type"].size();

          // Loop through mesh types
          for (int i = 0; i < nType; ++i)
            type.push_back(it["Mesh Type"][i].as<std::string>());
        }

        std::vector<std::pair<int, int>> elems;
        if (it.contains("Number of Elems")) {
          // Number of element pairs
          int nE = it["Number of Elems"].size();

          // Loop through element pairs
          for (int i = 0; i < nE; ++i) {
            std::pair<int, int> p;
            p.first = it["Number of Elems"][i][0].as<int>();
            p.second = it["Number of Elems"][i][1].as<int>();

            elems.push_back(p);
          }
        }

        std::vector<std::string> name;
        if (it.contains("Region Names")) {
          // Number of region names
          int nNames = it["Region Names"].size();
          // Loop through regions names
          for (int i = 0; i < nNames; ++i) {
            name.push_back(it["Region Names"][i].as<std::string>());

            // check if name is in physical name map
            auto iter = phystag_map.find(name[i]);
            // if name is not in map, add it
            if (iter == phystag_map.end()) {
              phystag_map.insert({name[i], physTag});
              physTag++;
            }
          }
        }

        auto search = savedobj_map.find(alias);
        if (search != savedobj_map.end()) {
          // the alias was found in map
          if (search->second.first == "Circles" ||
              search->second.first == "Circle") {
            if (it.contains("Center"))
              search->second.second[0].insert_or_assign("Center", cen);
            if (it.contains("Radii"))
              search->second.second[0].insert_or_assign("Radii", rad);
            if (it.contains("Mesh Type"))
              search->second.second[0].insert_or_assign("Mesh Type", type);
            if (it.contains("Number of Elems"))
              search->second.second[0].insert_or_assign("Number of Elems",
                                                        elems);
            if (it.contains("Region Names"))
              search->second.second[0].insert_or_assign("Region Names", name);
            if (it.contains("BREAK"))
              if (it["BREAK"] == true)
                search->second.second[0].insert_or_assign("BREAK", true);

            makeCircles(search->second.second);

            // Reset the alias map back to original data
            search->second.second[0].insert_or_assign("Center", orig_center);
            search->second.second[0].insert_or_assign("Rotation", orig_rot);
            search->second.second[0].insert_or_assign("Radii", orig_radii);
            search->second.second[0].insert_or_assign("Mesh Type",
                                                      orig_meshType);
            search->second.second[0].insert_or_assign("Number of Elems",
                                                      orig_elems);
            search->second.second[0].insert_or_assign("Region Names",
                                                      orig_names);
            search->second.second[0].insert_or_assign("BREAK", false);
          }
        } else {
          std::cout << alias << " not found in map" << std::endl;
        }
      } else {
        // Get the center coords
        std::vector<double> cen;
        cen.push_back(it["Center"][0].as<double>());
        cen.push_back(it["Center"][1].as<double>());
        cen.push_back(it["Center"][2].as<double>());
        if (it["Center"].size() == 5) {
          cen.push_back(it["Center"][3].as<double>());
          cen.push_back(it["Center"][4].as<double>());
        }

        // Number of radii
        int nRad = it["Radii"].size();

        // Loop through radii
        std::vector<double> rad;
        for (int i = 0; i < nRad; ++i)
          rad.push_back(it["Radii"][i].as<double>());

        // Number of mesh types
        int nType = it["Mesh Type"].size();

        // Loop through mesh types
        std::vector<std::string> type;
        for (int i = 0; i < nType; ++i)
          type.push_back(it["Mesh Type"][i].as<std::string>());

        // Number of element pairs
        int nE = it["Number of Elems"].size();

        // Loop through element pairs
        std::vector<std::pair<int, int>> elems;
        for (int i = 0; i < nE; ++i) {
          std::pair<int, int> p;
          p.first = it["Number of Elems"][i][0].as<int>();
          p.second = it["Number of Elems"][i][1].as<int>();

          elems.push_back(p);
        }

        // Number of region names
        int nNames = it["Region Names"].size();

        // Loop through regions names
        std::vector<std::string> name;
        for (int i = 0; i < nNames; ++i) {
          name.push_back(it["Region Names"][i].as<std::string>());

          // check if name is in physical name map
          auto iter = phystag_map.find(name[i]);
          if (iter == phystag_map.end()) {
            phystag_map.insert({name[i], physTag});
            physTag++;
          }
        }

        NEM::GEO::shape *c = new NEM::GEO::circles(cen, rad, type, elems, name);
        c->draw();

        if (!conserving) {
          shape_map.insert(std::pair<int, NEM::GEO::shape *>(id, c));
          id++;
        }
      }
      if (it.contains("BREAK")) {
        if (it["BREAK"] == true)
          openGUI();
      }
    }
  }
}

//********************************************************************//
// Parses json for array data then creates array of shape objects
//********************************************************************//
void NucMeshDriver::makeArray(jsoncons::json arr) {
  std::string arrType = arr["Array"].as_string();

  if (arrType == "Rectangular")
    std::cout << "Rectangular Array declared" << std::endl;
  if (arrType == "Polar")
    std::cout << "Polar Array declared" << std::endl;
  if (arrType == "Hexagonal")
    std::cout << "Hexagonal Array declared" << std::endl;

  // if (arr.contains("Name"))
  //    std::cout << "Creating '" << arr["Name"].as_string() << "'" <<
  //    std::endl;

  //--------------- Rectangular Array --------------//
  //------------------------------------------------//
  if (arrType == "Rectangular") {
    int nx = arr["NX"].as<int>();
    int ny = arr["NY"].as<int>();
    double dx = arr["DX"].as<double>();
    double dy = arr["DY"].as<double>();

    jsoncons::json s = arr["Shapes"]; // the array of shapes
    jsoncons::json circ, poly;

    std::vector<double> cent, c_cen;
    std::vector<std::vector<double>> oc;

    // Iterate through shapes in array
    for (const auto &shapes : s.array_range()) {
      if (shapes.contains("Name"))
        std::cout << "Creating '" << shapes["Name"].as_string() << "'"
                  << std::endl;
      if (shapes.contains("Circles")) {
        circ = shapes["Circles"];
        for (auto i = 0; i < circ.size(); ++i) {
          cent.push_back(circ[i]["Center"][0].as_double());
          cent.push_back(circ[i]["Center"][1].as_double());
          cent.push_back(circ[i]["Center"][2].as_double());
          oc.push_back(cent);
          cent.clear();
        }

        // update the center to make array of shape
        for (int j = 0; j < ny; ++j) {
          for (int i = 0; i < nx; ++i) {
            int index = 0;
            for (const auto &c : circ.array_range()) {
              if (c.contains("Center")) {
                // update center
                c_cen = {oc[index][0] + i * dx, oc[index][1] + j * dy,
                         oc[index][2]};

                // update json object
                circ[index].insert_or_assign("Center", c_cen);
              }
              index++;
            }
            // make the circle
            makeCircles(circ);
          }
        }
        oc.clear();
      }

      if (shapes.contains("Polygon")) {
        poly = shapes["Polygon"];
        for (auto i = 0; i < poly.size(); ++i) {
          cent.push_back(poly[i]["Center"][0].as_double());
          cent.push_back(poly[i]["Center"][1].as_double());
          cent.push_back(poly[i]["Center"][2].as_double());
          oc.push_back(cent);
          cent.clear();
        }
        // update the center to make array of shape
        for (int j = 0; j < ny; ++j) {
          for (int i = 0; i < nx; ++i) {
            int index = 0;
            for (const auto &p : poly.array_range()) {
              if (p.contains("Center")) {
                // update center
                c_cen = {oc[index][0] + i * dx, oc[index][1] + j * dy,
                         oc[index][2]};

                // update json object
                poly[index].insert_or_assign("Center", c_cen);
              }
              index++;
            }
            // make the polygon
            makePolygons(poly);
          }
        }
        oc.clear();
      }
    }
  }

  //------------------ Polar Array -----------------//
  //------------------------------------------------//
  if (arrType == "Polar") {
    // get polar array parameters
    std::vector<double> center; // center of array
    if (arr.contains("Center")) {
      center.push_back(arr["Center"][0].as_double());
      center.push_back(arr["Center"][1].as_double());
      center.push_back(arr["Center"][2].as_double());
      if (arr["Center"].size() == 5) {
        center.push_back(arr["Center"][3].as_double());
        center.push_back(arr["Center"][4].as_double());
      }
    }

    if (center.size() == 5) {
      double angle = center[4] * M_PI / 180.0;
      center[0] = center[0] + center[3] * std::cos(angle);
      center[1] = center[1] + center[3] * std::sin(angle);
    }

    double radius = arr["Radius"].as_double(); // array radius
    int n = 0;                                 // number of array elements
    if (arr.contains("N"))
      n = arr["N"].as<int>();

    double start = arr["Start Angle"].as_double(); // start angle
    if (start > 360.0) {
      std::cerr << "Error: 'Start Angle' in polar array is greater than 360.0 "
                   "degrees."
                << std::endl;
      exit(-1);
    }
    double arc = arr["Arc"].as_double(); // arc angle
    if (arc > 360.0) {
      std::cerr << "Error: 'Arc' in polar array is greater than 360.0 degrees. "
                   "Limiting 'Arc' to 360.0 degrees."
                << std::endl;
      arc = 360.0;
    }

    int nPat = 0; // number of pattern repeats
    if (arr.contains("N Patterns"))
      nPat = arr["N Patterns"].as<int>();

    std::vector<int> pattern_vec; // pattern vector
    if (arr.contains("Pattern"))
      for (int j = 0; j < nPat; ++j)
        for (int i = 0; i < arr["Pattern"].size(); ++i)
          if (arr["Pattern"][i].as<int>() > 0)
            pattern_vec.push_back(arr["Pattern"][i].as<int>());

    bool rotWithArray = false; // rotate with array bool
    if (arr.contains("Rotate with Array")) {
      rotWithArray = arr["Rotate with Array"].as<bool>();
      if (rotWithArray)
        std::cout << "'Rotate with Array' set to true" << std::endl;
    }

    // map containing array pattern items and corresponding angles
    std::map<int, std::vector<double>> item_map;

    double inc; // angle increment
    if (n != 0 && nPat == 0) {
      inc = (arc) / n;
      // create vector of angles
      std::vector<double> angles;
      for (int i = 0; i < n; ++i) {
        double a = inc * i + start;
        angles.push_back(a);
      }
      // insert angle vector into map with 1 key
      item_map.insert({1, angles});
    } else if (n == 0 && nPat != 0) {
      inc = arc / pattern_vec.size();

      // loop through pattern vector and insert into map
      int size = pattern_vec.size();
      for (int i = 0; i < size; ++i) {
        double a = inc * i + start;

        // search for item in map
        auto search = item_map.find(pattern_vec[i]);
        if (search == item_map.end())
          item_map.insert({pattern_vec[i], {a}});
        else
          search->second.push_back(a);
      }
    } else {
      std::cerr
          << "Error: 'N Patterns' and 'N' cannot both be greater than zero.\n";
      exit(1);
    }

    jsoncons::json s = arr["Shapes"]; // the array of shapes
    jsoncons::json circ, poly;

    std::vector<double> c_cen;
    std::vector<std::vector<double>> oc;

    // Iterate through shapes in array
    for (const auto &shapes : s.array_range()) {
      if (shapes.contains("Name"))
        std::cout << "Creating '" << shapes["Name"].as_string() << "'"
                  << std::endl;
      if (shapes.contains("Circles")) {
        circ = shapes["Circles"];
        int index = 0;
        for (const auto &c : circ.array_range()) {
          if (c.contains("Item")) {
            int item = c["Item"].as<int>();
            std::vector<double> angles;

            // search for item in map
            auto search = item_map.find(item);
            if (search != item_map.end()) {
              angles = search->second;
              for (const auto &angle : angles) {
                // set the center
                c_cen = {center[0], center[1], center[2], radius, angle};

                // update json object center
                circ[index].insert_or_assign("Center", c_cen);

                makeCircles(circ);
              }
            } else {
              std::cout << "item " << item << " not found in item map\n";
            }
          } else {
            std::vector<double> angles;
            // search for item in map
            auto search = item_map.find(1);
            if (search != item_map.end()) {
              angles = search->second;
              for (const auto &angle : angles) {
                // set the center
                c_cen = {center[0], center[1], center[2], radius, angle};

                // update json object center
                circ[index].insert_or_assign("Center", c_cen);

                makeCircles(circ);
              }
            }
          }
        }
        index++;
      }

      if (shapes.contains("Polygon")) {
        poly = shapes["Polygon"];
        int index = 0;
        for (const auto &c : poly.array_range()) {
          if (c.contains("Item")) {
            int item = c["Item"].as<int>();
            std::vector<double> angles;

            // search for item in map
            auto search = item_map.find(item);
            if (search != item_map.end()) {
              angles = search->second;
              for (const auto &angle : angles) {
                // set the center
                c_cen = {center[0], center[1], center[2], radius, angle};

                // update json object
                poly[index].insert_or_assign("Center", c_cen);

                // if 'Rotate with Array', apply rotation
                if (rotWithArray)
                  poly[index].insert_or_assign("Rotation", angle);

                makePolygons(poly);
              }
            } else {
              std::cout << "item " << item << " not found in item map\n";
            }
          } else {
            std::vector<double> angles;

            // search for item in map
            auto search = item_map.find(1);
            if (search != item_map.end()) {
              angles = search->second;
              for (const auto &angle : angles) {
                // set the center
                c_cen = {center[0], center[1], center[2], radius, angle};

                // update json object
                poly[index].insert_or_assign("Center", c_cen);

                // if 'Rotate with Array', apply rotation
                if (rotWithArray)
                  poly[index].insert_or_assign("Rotation", angle);

                makePolygons(poly);
              }
            }
          }
        }
        index++;
      }
    }
  }

  //---------------- Hexagonal Array ---------------//
  //------------------------------------------------//
  if (arrType == "Hexagonal") {
    std::vector<double> cent, c_cen;
    std::vector<std::vector<double>> orig_center;

    std::vector<double> arr_center; // center of array
    // get hexagonal array parameters
    if (arr.contains("Center")) {
      arr_center.push_back(arr["Center"][0].as_double());
      arr_center.push_back(arr["Center"][1].as_double());
      arr_center.push_back(arr["Center"][2].as_double());
      if (arr["Center"].size() == 5) {
        arr_center.push_back(arr["Center"][3].as_double());
        arr_center.push_back(arr["Center"][4].as_double());
      }
    }

    if (arr_center.size() == 5) {
      double angle = arr_center[4] * M_PI / 180.0;
      arr_center[0] = arr_center[0] + arr_center[3] * std::cos(angle);
      arr_center[1] = arr_center[1] + arr_center[3] * std::sin(angle);
    }

    double radius;
    if (arr.contains("Radius")) {
      radius = arr["Radius"].as_double();  // array radius
    } else {
      std::cerr << "Error: No 'Radius' found for array parameters."
                << std::endl;
      exit(-1);
    }

    int n = 0; // number of array elements
    if (arr.contains("N"))
      n = arr["N"].as<int>();

    int type;
    if (arr.contains("Type"))
      type = arr["Type"].as<int>();
    else {
      type = 1;
      std::cout << "Warning: Hexagonal array 'Type' not set. Using 'Type = 1'."
                << std::endl;
    }
    double padding = 0;
    if (arr.contains("Padding"))
      padding = arr["Padding"].as_double(); // padding
    else {
      std::string name;
      if (arr.contains("Name")) {
        name = arr["Name"].as_string();
        std::cerr << "\nWarning: 'Padding' keyword in hexagonal array"
                  << " '" << name << "' not found. Using Padding = 0.0"
                  << std::endl;
      } else {
        std::cerr << "\nWarning: 'Padding' keyword in unnamed hexagonal array "
                     "not found. Using Padding = 0.0"
                  << std::endl;
      }
    }

    double start = 0.0;
    if (arr.contains("Start Angle")) {
      start = arr["Start Angle"].as_double(); // start angle
      if (start > 360.0) {
        std::cerr << "\nError: 'Start Angle' in hexagonal array is greater "
                     "than 360.0 degrees.\n"
                  << std::endl;
        exit(-1);
      }
    }
    double arc = 360.0;
    if (arr.contains("Arc")) {
      arc = arr["Arc"].as_double(); // arc angle
      if (arc > 360.0) {
        std::cerr << "Error: 'Arc' in hexagonal array is greater than 360.0 "
                     "degrees. Limiting 'Arc' to 360.0 degrees."
                  << std::endl;
        arc = 360.0;
      }
    }

    jsoncons::json circ, poly;
    int half = n / 2;
    double ang = 60.0 * M_PI / 180.0;
    double xOff = 0.0, yOff = 0.0;
    // If Circles are used and visible, offset to edge of circle
    jsoncons::json s = arr["Shapes"]; // the array of shapes
    for (const auto &shapes : s.array_range()) {
      if (shapes.contains("BREAK"))
        if (shapes["BREAK"] == true)
          openGUI();
      if (shapes.contains("Circles")) {
        circ = shapes["Circles"];
        for (const auto &it : circ.array_range()) {
          if (it.contains("Visible")) {
            if (it["Visible"] == true) {
              xOff = std::cos(ang) * (2 * radius + padding);
              yOff = std::sin(ang) * (2 * radius + padding);
            }
          }
        }
      }
      // Else offset to edge of hexagon
      else {
        xOff = std::cos(ang) * (std::sqrt(3) * radius + padding);
        yOff = std::sin(ang) * (std::sqrt(3) * radius + padding);
      }
    }

    double a = 0.0;
    if (type == 2) {
      a = 30 * M_PI / 180.0; // rad
    }

    double /*r, */theta;
    for (int row = 0; row < n; row++) {
      int cols = n - std::abs(row - half);
      for (int col = 0; col < cols; col++) {
        double x = (xOff * (col * 2 + 1 - cols));
        double y = (yOff * (row - half));

        // convert x,y to polar coordinates
        // r = std::sqrt(x * x + y * y);
        if (std::abs(x) > 1e-8) {
          if (x < 0 && y < 0)
            theta = std::atan(y / x) + M_PI;
          if (x < 0 && y > 0)
            theta = std::atan(y / x) + M_PI;
          if (x > 0 && y < 0)
            theta = std::atan(y / x) + 2 * M_PI;
          if (x > 0 && y > 0)
            theta = std::atan(y / x);
          if (y == 0 && x < 0)
            theta = M_PI;
          if (y == 0 && x > 0)
            theta = 0.0;
        } else if (x == 0 && y > 0)
          theta = M_PI / 2;
        else if (x == 0 && y < 0)
          theta = 3 * M_PI / 2;
        else
          theta = 0.0;

        theta = theta * 180.0 / M_PI;
        if (theta > arc || theta < start)
          continue;
        else {
          // Iterate through shapes in array
          double new_x, new_y, new_z;
          for (const auto &shapes : s.array_range()) {
            // If a 'Name' exists, output
            if (shapes.contains("Name")) {
              std::cout << "Creating '" << shapes["Name"].as_string() << "'"
                        << std::endl;
            }
            if (shapes.contains("Circles")) {
              circ = shapes["Circles"];
              int index = 0;
              for (const auto &c : circ.array_range()) {
                if (c.contains("Center")) {
                  cent.push_back(c["Center"][0].as_double());
                  cent.push_back(c["Center"][1].as_double());
                  cent.push_back(c["Center"][2].as_double());
                  orig_center.push_back(cent);
                  cent.clear();
                }
                // TODO: change for alias object, maybe
                else
                  orig_center.push_back({0, 0, 0});
                // set the shape center
                new_x = (arr_center[0] + orig_center[index][0] +
                         x * std::cos(a) - y * std::sin(a));
                new_y = (arr_center[1] + orig_center[index][1] +
                         x * std::sin(a) + y * std::cos(a));
                new_z = arr_center[2] + orig_center[index][2];
                c_cen = {new_x, new_y, new_z};
                // update json object center
                circ[index].insert_or_assign("Center", c_cen);
                index++;
              }
              makeCircles(circ);
            }
            orig_center.clear();

            if (shapes.contains("Polygon")) {
              poly = shapes["Polygon"];
              int index = 0;
              for (const auto &c : poly.array_range()) {
                // Set original center to zero, if otherwise, change it
                if (c.contains("Center")) {
                  cent.push_back(c["Center"][0].as_double());
                  cent.push_back(c["Center"][1].as_double());
                  cent.push_back(c["Center"][2].as_double());
                  orig_center.push_back(cent);
                  cent.clear();
                }
                // TODO: change for alias object, maybe
                else
                  orig_center.push_back({0, 0, 0});
                // set the center
                new_x = (arr_center[0] + orig_center[index][0] +
                         x * std::cos(a) - y * std::sin(a));
                new_y = (arr_center[1] + orig_center[index][1] +
                         x * std::sin(a) + y * std::cos(a));
                new_z = arr_center[2] + orig_center[index][2];
                c_cen = {new_x, new_y, new_z};

                // update json object
                poly[index].insert_or_assign("Center", c_cen);

                // Change rotation if hex array type is 1
                if (type == 1) {
                  if (c.contains("Number of Sides")) {
                    if (c["Number of Sides"] == 6)
                      poly[index].insert_or_assign("Rotation", 30.0);
                  } else {
                    if (c.contains("Alias")) {
                      std::string alias = c["Alias"].as_string();
                      auto search = savedobj_map.find(alias);
                      if (search != savedobj_map.end()) {
                        // the alias was found in map
                        if (search->second.first == "Polygon" ||
                            search->second.first == "Hexagon") {
                          if (search->second.second[0]["Number of Sides"] == 6)
                            poly[index].insert_or_assign("Rotation", 30.0);
                        }
                      }
                    }
                  }
                }
                index++;
              }
              makePolygons(poly);
            }
            orig_center.clear();
          }
        }
      }
    }
  }
  if (arr.contains("BREAK"))
    if (arr["BREAK"] == true)
      openGUI();
}

//********************************************************************//
// Causes the Gmsh GUI to wait for command line input
//********************************************************************//
void NucMeshDriver::openGUI() {
  // Tpause.start();
  gmsh::graphics::draw();
  while (1) {
    gmsh::fltk::awake();
    gmsh::fltk::wait();
    gmsh::fltk::update();
    std::cout << "\nPress '1' to continue, '0' to cancel" << std::endl;
    bool in;
    std::cin >> in;
    if (in)
      break;
    else
      exit(0);
  }
}
