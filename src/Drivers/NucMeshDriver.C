#define _USE_MATH_DEFINES
#include "NucMeshDriver.H"

#include <gmsh.h>
#include <cmath>
#include <iomanip>

#include "AuxiliaryFunctions.H"
#include "circles.H"
#include "circlesInPolys.H"
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
    throw;
  }
}

NucMeshDriver::NucMeshDriver(jsoncons::json inputjson) {
  id = 1;
  physTag = 1;
  gui = false;
  skipAll = false;
  parsingCount = 0;

  std::cout << "NucMeshDriver created" << std::endl;

  if (inputjson.contains("Output File Name"))
    ofname = inputjson["Output File Name"].as_string();
  else {
    std::cerr << "Error: 'Output File Name' keyword not found, expected after "
                 "'Program Type'"
              << std::endl;
    throw;
  }

  jsoncons::json shapes = inputjson["Geometry and Mesh"];

  Tgeom.start();
  gmsh::initialize();

  gmsh::model::add("NucMesh");
  gmsh::option::setNumber("General.NumThreads", 10);
  gmsh::option::setNumber("Geometry.OCCBooleanPreserveNumbering", 1);
  gmsh::option::setNumber("Geometry.OCCParallel", 1);
  gmsh::option::setNumber("Geometry.Tolerance", 1e-8);
  gmsh::option::setNumber("Geometry.MatchMeshTolerance", 1e-6);
  gmsh::option::setNumber("Mesh.FlexibleTransfinite", 0);
  gmsh::option::setNumber("Mesh.MshFileVersion", 2.2);
  gmsh::option::setNumber("Mesh.CharacteristicLengthFromCurvature", 0);

  //--------------------------------------------------------------------//
  // Main JSON file geometry parsing
  //--------------------------------------------------------------------//
  parseGeomAndMesh(shapes);
  Tgeom.stop();
  if (gui) gmsh::fltk::initialize();

  //--------------------------------------------------------------------//
  // Perform Boolean Fragments
  //--------------------------------------------------------------------//
  Tboolean.start();
  // gmsh::fltk::run();
  std::vector<std::pair<int, int>> surf;
  std::vector<std::pair<int, int>> outSurf;
  std::vector<std::vector<std::pair<int, int>>> outSurfMap;
  // Get all surfaces
  gmsh::model::getEntities(surf, 2);
  // Boolean Fragments
  gmsh::model::occ::fragment(surf, surf, outSurf, outSurfMap, -1, true, true);
  gmsh::model::occ::synchronize();

  int oldID, newID;
  std::pair<int, int> oldNew;                   // holds the old/new id pair
  std::vector<std::pair<int, int>> oldNew_vec;  // vector for id pair

  // gmsh::fltk::run();

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
  // gmsh::fltk::run();

  //--------------------------------------------------------------------//
  // Apply meshType and get physical surfaces
  //--------------------------------------------------------------------//
  // iterate through shape_map to apply meshType and collect physical Surfs
  std::map<int, int> physSurf_map;
  for (auto &&itr : shape_map) {
    itr.second->applyMeshType();
    physSurf_map = itr.second->getPhysSurf(phystag_map, physSurf_map);
  }

  // gmsh::fltk::run();
  std::cout << "Mesh type applied\n";

  //--------------------------------------------------------------------//
  // Make physical groups for surfaces
  //--------------------------------------------------------------------//
  for (const auto &itr : phystag_map) {
    std::vector<int> s;  // vector of surfaces
    std::string name;    // regions name
    int tag = -1;        // physical group tag
    for (const auto &itr2 : physSurf_map) {
      if (itr2.second == itr.second) {
        tag = itr2.second;
        name = itr.first;
        s.push_back(itr2.first);
      }
    }
    if (tag == -1) {
      std::cerr << "Error: Physical group cannot be applied." << std::endl;
      throw;
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
      throw;
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
    int tag = 1;
    for (const auto &itr : volPhysTag_map) {
      // Add physical group
      gmsh::model::addPhysicalGroup(3, itr.second, tag);
      // Get physical name for physical tag of surface
      gmsh::model::getPhysicalName(2, itr.first, name);
      // Set physical name of Volume
      gmsh::model::setPhysicalName(3, tag, name);
      tag++;
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
  if (gui) openGUI();

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
  if (gui) openGUI();
  gmsh::finalize();
}

//********************************************************************//
// Parses the Geometry and Mesh main section
//********************************************************************//
void NucMeshDriver::parseGeomAndMesh(jsoncons::json shapes) {
  // Loop through Geometry and Mesh to get shapes
  std::cout << "Geometry and Mesh" << std::endl;
  for (const auto &obj : shapes.array_range()) {
    if (obj.contains("Global Options")) parseOptions(obj["Global Options"]);

    if (obj.contains("Saved Objects")) parseSavedObjects(obj["Saved Objects"]);

    if (obj.contains("Visible"))
      if (obj["Visible"] == false) continue;

    if (obj.contains("Show"))
      if (obj["Show"] == false) continue;

    if (obj.contains("Name"))
      std::cout << "Creating '" << obj["Name"].as_string() << "'" << std::endl;

    if (obj.contains("Include")) {
      std::cout << "Include called." << std::endl;

      std::string fname = obj["Include"].as_string();
      std::cout << fname << std::endl;
      std::ifstream inputStream(fname);
      if (!inputStream.good() || nemAux::find_ext(fname) != ".json") {
        std::cerr << "Error opening file " << fname << std::endl;
        throw;
      }

      jsoncons::json includejson;
      inputStream >> includejson;
      if (includejson.is_array()) {
        parsingCount++;
        if (parsingCount > 200) {
          std::cerr << "Effor: Parsing JSON seems to be stuck in infinte loop."
                    << " Check included JSON files." << std::endl;
          throw;
        }
        parseGeomAndMesh(includejson);
      }
    }

    if (obj.contains("Circles")) makeCircles(obj["Circles"]);

    if (obj.contains("Circle")) makeCircles(obj["Circle"]);

    if (obj.contains("CirclesInPolys"))
      makeCirclesInPolys(obj["CirclesInPolys"]);

    if (obj.contains("Hexagon") || obj.contains("Polygon")) {
      if (obj.contains("Hexagon"))
        makePolygons(obj["Hexagon"], 6);
      else
        makePolygons(obj["Polygon"]);
    }
    if (obj.contains("BREAK"))
      if (obj["BREAK"] == true) openGUI();

    if (obj.contains("Array")) makeArray(obj);
  }
}

//********************************************************************//
// Parses the global options data of json
//********************************************************************//
void NucMeshDriver::parseOptions(jsoncons::json opts) {
  // Iterate through options
  for (const auto &it : opts.array_range()) {
    if (it.contains("Open GUI")) gui = it["Open GUI"].as_bool();
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
                throw;
              } else
                layers.push_back(layerValue);
            }
          } else {
            std::cerr << "Error: 'Number of Layers' array is empty"
                      << std::endl;
            throw;
          }
        } else {
          std::cerr << "Error: 'Number of Layers' keyword not found"
                    << std::endl;
          throw;
        }
        if (it.contains("Heights")) {
          if (it["Heights"].size() > numLayersSize) {
            std::cerr
                << "Error: 'Heights' vector larger than 'Number of Layers'";
            throw;
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
      throw;
    }

    if (it.contains("Min Mesh Size"))
      gmsh::option::setNumber("Mesh.CharacteristicLengthMin",
                              it["Min Mesh Size"].as<double>());

    if (it.contains("Max Mesh Size"))
      gmsh::option::setNumber("Mesh.CharacteristicLengthMax",
                              it["Max Mesh Size"].as<double>());

    if (it.contains("Meshing Algorithm")) {
      std::string algo = it["Meshing Algorithm"].as<std::string>();
      int a = 6;  // default: "Frontal"
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
      int a = 1;  // default: "Blossom"
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
      if (!extend_from_bnd) e = 0;
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
      throw;
    }

    if (it.contains("Circles")) shapetype = "Circles";
    if (it.contains("Circle")) shapetype = "Circle";
    if (it.contains("CirclesInPolys")) shapetype = "CirclesInPolys";
    if (it.contains("Polygon")) shapetype = "Polygon";
    if (it.contains("Array")) shapetype = "Array";

    bool visible = true;
    // If meshed area conservation is required and visible, change the radii
    if (it.contains("Conserve Area") && it["Conserve Area"] == true) {
      if (it.contains("Circles") || it.contains("Circle")) {
        for (const auto &c : it["Circles"].array_range()) {
          if (c.contains("Visible")) visible = !(c["Visible"] == false);
        }
        jsoncons::json j;
        if (visible) {
          double tolerance = 1e-8;
          if (it.contains("Tolerance")) tolerance = it["Tolerance"].as_double();
          j = correctMeshArea(it[shapetype], shapetype, tolerance);
        } else
          j = it[shapetype];

        savedobj_map.insert(
            std::pair<std::string, std::pair<std::string, jsoncons::json>>(
                alias, {shapetype, j}));
      }

      if (it.contains("CirclesInPolys")) {
        for (const auto &c : it["CirclesInPolys"].array_range()) {
          if (c.contains("Visible")) visible = !(c["Visible"] == false);
        }
        jsoncons::json j;
        if (visible) {
          double tolerance = 1e-8;
          if (it.contains("Tolerance")) tolerance = it["Tolerance"].as_double();
          std::cout << shapetype << std::endl;
          j = correctMeshArea(it[shapetype], shapetype, tolerance);
        } else
          j = it[shapetype];

        savedobj_map.insert(
            std::pair<std::string, std::pair<std::string, jsoncons::json>>(
                alias, {shapetype, j}));
      }
    } else {
      savedobj_map.insert(
          std::pair<std::string, std::pair<std::string, jsoncons::json>>(
              alias, {shapetype, it[shapetype]}));
    }
    std::cout << "Saved Object Parsed\n" << std::endl;
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

  int nRad;
  bool circPoly = false;

  // For CirclesInPolys objects, we need to only change the
  // radii of the circles and not the polygon. Therefore,
  // we are going to take the 'Circle Radii' and put the
  // values into a 'Circle' object, get the adujsted radii,
  // and put the new radii back to the CirclesInPolys object.

  if (shapetype == "Circles" || shapetype == "Circle" ||
      shapetype == "CirclesInPolys") {
    if (shapetype == "CirclesInPolys") circPoly = true;

    if (shapetype == "Circles") {
      if (obj[0].contains("Radii")) {
        nRad = obj[0]["Radii"].size();
        if (nRad <= 0) {
          std::cerr << "Error: No 'Radii' values found for aliased object."
                    << std::endl;
          throw;
        }
        for (int i = 0; i < nRad; ++i)
          radii.push_back(obj[0]["Radii"][i].as_double());
      } else {
        std::cerr << "Error: 'Radii' keyword not found for aliased object."
                  << std::endl;
        throw;
      }
    }

    jsoncons::json circle;
    if (shapetype == "CirclesInPolys") {
      // Instantiate new circle json object
      circle = obj;
      if (obj[0].contains("Circle Radii")) {
        nRad = obj[0]["Circle Radii"].size();
        if (nRad <= 0) {
          std::cerr
              << "Error: No 'Circle Radii' values found for aliased object."
              << std::endl;
          throw;
        }
        for (int i = 0; i < nRad; ++i)
          radii.push_back(obj[0]["Circle Radii"][i].as_double());
      } else {
        std::cerr
            << "Error: 'Circle Radii' keyword not found for aliased object."
            << std::endl;
        throw;
      }
      // Insert 'Radii' parameter to circle object
      circle[0].insert_or_assign("Radii", radii);
    }

    double lastRadius = 0.0, area;
    double eps = tolerance;
    double residual;
    std::vector<double> r_a, r_b, r_c;
    double A_c;
    int iter;

    // Conserving mesh area for each ring of circles,
    // so loop through each radii and iterate using
    // bisection method
    for (int i = 0; i < nRad; ++i) {
      iter = 0;
      r_a.push_back(radii[i]);
      r_b.push_back(radii[i] * 1.5);
      r_c.push_back((r_a[i] + r_b[i]) / 2);
      lastRadius = radii[i];
      area = M_PI * lastRadius * lastRadius;
      residual = 10;
      while (residual > eps) {
        if (iter > 50) {
          std::cerr << "Error: Iterations reached maximum." << std::endl;
          throw;
        }
        iter++;

        r_c[i] = (r_a[i] + r_b[i]) / 2;

        if (circPoly) {
          circle[0].insert_or_assign("Radii", r_c);
          makeCircles(circle, true);
        } else {
          obj[0].insert_or_assign("Radii", r_c);
          makeCircles(obj, true);
        }

        gmsh::model::mesh::generate(2);
        // Using MeshVolume (Area) plugin to get meshed area
        gmsh::plugin::setNumber("MeshVolume", "Dimension", 2);
        gmsh::plugin::run("MeshVolume");

        std::vector<int> viewTags;
        gmsh::view::getTags(viewTags);
        if (viewTags.size() == 1) {
          std::vector<std::vector<double>> data;
          std::vector<std::string> dataTypes;
          std::vector<int> numElements;
          gmsh::view::getListData(viewTags[0], dataTypes, numElements, data);
          // This is the meshed area data
          A_c = data[0][3];
        } else {
          std::cerr << "Error: View data not found." << std::endl;
          throw;
        }
        viewTags.clear();
        gmsh::clear();

        residual = std::abs(A_c - area);
        std::cout << std::scientific << std::setprecision(4)
                  << "Residual = " << residual << std::endl;

        if (A_c - area > 0)
          r_b = r_c;
        else
          r_a = r_c;
      }
      r_c = r_b;
      std::cout << std::setprecision(8) << std::fixed
                << "New radius = " << r_b[i] << std::endl;
      std::cout << "Mesh area converged in " << iter << " iterations"
                << std::endl;
    }

    // If CirclesInPolys object, insert new circle radii
    if (circPoly) obj[0].insert_or_assign("Circle Radii", r_b);  // CHECK THIS
  } else {
    std::cout << "Warning: Mesh area conservation only enabled for Circles and "
                 "CirclesInPolys."
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
            if (j.contains("Rotation"))
              orig_rot = j["Rotation"].as<double>();

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
        if (it["BREAK"] == true) openGUI();
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
            if (j.contains("Rotation"))
              orig_rot = j["Rotation"].as<double>();

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
        if (it["BREAK"] == true) openGUI();
      }
    }
  }
}

//********************************************************************//
// Parses json for circlesInPolys data then creates circlesInPoly object
//********************************************************************//
void NucMeshDriver::makeCirclesInPolys(jsoncons::json circInPoly,
                                       bool conserving) {
  // Iterate through key
  for (const auto &it : circInPoly.array_range()) {
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

        std::vector<double> orig_center, orig_circle_radii, orig_poly_radii;
        std::vector<std::string> orig_meshType, orig_names;
        std::vector<std::pair<int, int>> orig_elems;
        double orig_rot;
        int orig_nsides;
        // Get the original alias data
        auto find = savedobj_map.find(alias);
        if (find != savedobj_map.end()) {
          // the alias was found in map
          if (find->second.first == "CirclesInPolys") {
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
            if (j.contains("Rotation"))
              orig_rot = j["Rotation"].as<double>();

            if (j.contains("Number of Sides"))
              orig_nsides = j["Number of Sides"].as<int>();

            if (j.contains("Circle Radii")) {
              int nCircleRad = j["Circle Radii"].size();
              for (int i = 0; i < nCircleRad; ++i)
                orig_circle_radii.push_back(j["Circle Radii"][i].as<double>());
            }

            if (j.contains("Polygon Radii")) {
              int nPolyRad = j["Polygon Radii"].size();
              for (int i = 0; i < nPolyRad; ++i)
                orig_poly_radii.push_back(j["Polygon Radii"][i].as<double>());
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
        } else {
          std::cerr << "Alias '" << alias << "' was not found in map."
                    << std::endl;
          throw;
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

        int nsides;
        if (it.contains("Number of Sides")) {
          // Get the rotation
          nsides = it["Number of Sides"].as<int>();
        }

        std::vector<double> circle_rad;
        if (it.contains("Circle Radii")) {
          // Number of circle radii
          int nCircleRad = it["Circle Radii"].size();

          // Loop through circle radii
          for (int i = 0; i < nCircleRad; ++i)
            circle_rad.push_back(it["Circle Radii"][i].as<double>());
        }

        std::vector<double> poly_rad;
        if (it.contains("Polygon Radii")) {
          // Number of polys radii
          int nPolyRad = it["Polygon Radii"].size();

          // Loop through poly radii
          for (int i = 0; i < nPolyRad; ++i) {
            poly_rad.push_back(it["Polygon Radii"][i].as<double>());
            std::cout << it["Polygon Radii"][i].as<double>() << std::endl;
          }
        }

        // Check for overlapping circle/polygon
        if (circle_rad.size() > 0 && poly_rad.size() > 0) {
          double incircle = checkInscribedCircle(circle_rad, poly_rad, nsides);
          if (incircle < 1e10) {
            std::cerr << "Error: Circle and Polygon will overlap in alias '"
                      << alias << "." << std::endl;
            std::cerr << "       The largest circle radius allowed is "
                      << incircle * 0.97 << std::endl;
            throw;
          }
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
          if (search->second.first == "CirclesInPolys") {
            if (it.contains("Center"))
              search->second.second[0].insert_or_assign("Center", cen);
            if (it.contains("Rotation"))
              search->second.second[0].insert_or_assign("Rotation", rot);
            if (it.contains("Number of Sides"))
              search->second.second[0].insert_or_assign("Number of Sides",
                                                        nsides);
            if (it.contains("Circle Radii"))
              search->second.second[0].insert_or_assign("Circle Radii",
                                                        circle_rad);
            if (it.contains("Polygon Radii"))
              search->second.second[0].insert_or_assign("Polygon Radii",
                                                        poly_rad);
            if (it.contains("Mesh Type"))
              search->second.second[0].insert_or_assign("Mesh Type", type);
            if (it.contains("Number of Elems"))
              search->second.second[0].insert_or_assign("Number of Elems",
                                                        elems);
            if (it.contains("Region Names"))
              search->second.second[0].insert_or_assign("Region Names", name);

            makeCirclesInPolys(search->second.second);

            // Reset the alias map back to original data
            search->second.second[0].insert_or_assign("Center", orig_center);
            search->second.second[0].insert_or_assign("Rotation", orig_rot);
            search->second.second[0].insert_or_assign("Number of Sides",
                                                      orig_nsides);
            search->second.second[0].insert_or_assign("Circle Radii",
                                                      orig_circle_radii);
            search->second.second[0].insert_or_assign("Polygon Radii",
                                                      orig_poly_radii);
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

        // Number of circle radii
        int nCircleRad = it["Circle Radii"].size();
        // Number of poly radii
        int nPolyRad = it["Polygon Radii"].size();

        // Loop through circle radii
        std::vector<double> circle_rad;
        for (int i = 0; i < nCircleRad; ++i)
          circle_rad.push_back(it["Circle Radii"][i].as<double>());

        // Loop through polygon radii
        std::vector<double> poly_rad;
        for (int i = 0; i < nPolyRad; ++i)
          poly_rad.push_back(it["Polygon Radii"][i].as<double>());

        // Check for overlapping circle/polygon
        double incircle = checkInscribedCircle(circle_rad, poly_rad, nsides);
        if (incircle < 1e10) {
          std::cerr << "Error: Circle and Polygon will overlap." << std::endl;
          std::cerr << "       The largest circle radius allowed is "
                    << incircle * 0.97 << std::endl;
          std::cerr << jsoncons::pretty_print(it) << std::endl;
          throw;
        }
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

        NEM::GEO::shape *cp = new NEM::GEO::circlesInPolys(
            nsides, cen, circle_rad, poly_rad, type, elems, name, rot);

        cp->draw();

        if (!conserving) {
          shape_map.insert(std::pair<int, NEM::GEO::shape *>(id, cp));
          id++;
        }
        nsides = 0;
        type.clear();
        elems.clear();
      }
      if (it.contains("BREAK")) {
        if (it["BREAK"] == true) openGUI();
      }
    }
  }
}

//********************************************************************//
// Parses json for array data then creates array of shape objects
//********************************************************************//
void NucMeshDriver::makeArray(jsoncons::json arr) {
  std::string arrType = arr["Array"].as_string();

  if (arrType == "Rectangular") {
    std::cout << "Rectangular Array declared" << std::endl;
    makeRectangularArray(arr);
  }
  if (arrType == "Cartesian") {
    std::cout << "Cartesian Array declared" << std::endl;

    makeCartesianArray(arr);
  }
  if (arrType == "Polar") {
    std::cout << "Polar Array declared" << std::endl;
    makePolarArray(arr);
  }
  if (arrType == "Hexagonal") {
    std::cout << "Hexagonal Array declared" << std::endl;
    makeHexagonalArray(arr);
  }

  // if (arr.contains("Name"))
  //    std::cout << "Creating '" << arr["Name"].as_string() << "'" <<
  //    std::endl;
  if (arr.contains("BREAK"))
    if (arr["BREAK"] == true) openGUI();
}

//********************************************************************//
// Parses json for rectangular array then creates array of shape objects
//********************************************************************//
void NucMeshDriver::makeRectangularArray(jsoncons::json arr) {
  int nx = 0;
  int ny = 0;
  double dx = 0;
  double dy = 0;
  if (arr.contains("NX"))
    nx = arr["NX"].as<int>();
  else {
    std::cerr << "Error: Keyword 'NX' not found for rectangular array."
              << std::endl;
    throw;
  }
  if (arr.contains("NY"))
    ny = arr["NY"].as<int>();
  else {
    std::cerr << "Error: Keyword 'NY' not found for rectangular array."
              << std::endl;
    throw;
  }
  if (arr.contains("DX"))
    dx = arr["DX"].as<double>();
  else {
    std::cerr << "Error: Keyword 'DX' not found for rectangular array."
              << std::endl;
    throw;
  }
  if (arr.contains("DY"))
    dy = arr["DY"].as<double>();
  else {
    std::cerr << "Error: Keyword 'DY' not found for rectangular array."
              << std::endl;
    throw;
  }

  // map containing array pattern items and corresponding location
  std::map<int, std::vector<std::vector<double>>> item_map;
  int k = 0;
  double transX = 0;
  double transY = 0;
  if (arr.contains("Pattern")) {
    int pattern_size = arr["Pattern"].size();
    if (pattern_size != nx * ny) {
      std::cerr << "Error: Rectangular array 'Pattern' size not equal to NX*NY."
                << std::endl;
      throw;
    }
    for (int j = ny - 1; j > -1; --j) {
      for (int i = nx - 1; i > -1; --i) {
        int item_id = arr["Pattern"][k].as<int>();
        if (item_id > 0) {
          transX = (nx - i - 1) * dx;
          transY = j * dy;
          // search for item in map
          auto search = item_map.find(item_id);
          if (search == item_map.end())
            item_map.insert({item_id, {{transX, transY}}});
          else
            search->second.push_back({transX, transY});
        }
        k++;
      }
    }
  }

  jsoncons::json s = arr["Shapes"];  // the array of shapes
  jsoncons::json circ, poly, circInPoly;

  std::vector<double> cent, updated_center;
  std::vector<std::vector<double>> orig_center;

  // Iterate through shapes in array
  for (const auto &shapes : s.array_range()) {
    // if (shapes.contains("Name"))
    // std::cout << "Creating '" << shapes["Name"].as_string() << "'"
    //          << std::endl;

    if (shapes.contains("CirclesInPolys")) {
      circInPoly = shapes["CirclesInPolys"];
      int circlesSize = circInPoly.size();
      for (auto i = 0; i < circlesSize; ++i) {
        // get original center of circle
        cent.push_back(circInPoly[i]["Center"][0].as_double());
        cent.push_back(circInPoly[i]["Center"][1].as_double());
        cent.push_back(circInPoly[i]["Center"][2].as_double());
        orig_center.push_back(cent);
        cent.clear();
        // get the item number
        int item = 0;
        if (circInPoly[i].contains("Item")) {
          item = circInPoly[i]["Item"].as<int>();
          std::vector<std::vector<double>> translations;
          // search for item in map
          auto search = item_map.find(item);
          if (search != item_map.end()) {
            translations = search->second;
            for (const auto &transPair : translations) {
              // set the center
              updated_center = {orig_center[i][0] + transPair[0],
                                orig_center[i][1] + transPair[1],
                                orig_center[i][2]};

              // update the json object
              circInPoly[i].insert_or_assign("Center", updated_center);

              // make the circle
              makeCirclesInPolys(circInPoly);
            }
          }
        } else {
          // update the center to make array of shape
          for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
              int index = 0;
              for (const auto &p : circInPoly.array_range()) {
                if (p.contains("Center")) {
                  // update center
                  updated_center = {orig_center[index][0] + i * dx,
                                    orig_center[index][1] + j * dy,
                                    orig_center[index][2]};

                  // update json object
                  circInPoly[index].insert_or_assign("Center", updated_center);
                }
                index++;
              }
              // make the circleInPoly
              makeCirclesInPolys(circInPoly);
            }
          }
          orig_center.clear();
        }
      }
    }

    if (shapes.contains("Circles")) {
      circ = shapes["Circles"];
      int circlesSize = circ.size();
      for (auto i = 0; i < circlesSize; ++i) {
        // get original center of circle
        cent.push_back(circ[i]["Center"][0].as_double());
        cent.push_back(circ[i]["Center"][1].as_double());
        cent.push_back(circ[i]["Center"][2].as_double());
        orig_center.push_back(cent);
        cent.clear();
        // get the item number
        int item = 0;
        if (circ[i].contains("Item")) {
          item = circ[i]["Item"].as<int>();
          std::vector<std::vector<double>> translations;
          // search for item in map
          auto search = item_map.find(item);
          if (search != item_map.end()) {
            translations = search->second;
            for (const auto &transPair : translations) {
              // set the center
              updated_center = {orig_center[i][0] + transPair[0],
                                orig_center[i][1] + transPair[1],
                                orig_center[i][2]};

              // update the json object
              circ[i].insert_or_assign("Center", updated_center);

              // make the circle
              makeCircles(circ);
            }
          }
        } else {
          // update the center to make array of shape
          for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
              int index = 0;
              for (const auto &c : circ.array_range()) {
                if (c.contains("Center")) {
                  // update center
                  updated_center = {orig_center[index][0] + i * dx,
                                    orig_center[index][1] + j * dy,
                                    orig_center[index][2]};

                  // update json object
                  circ[index].insert_or_assign("Center", updated_center);
                }
                index++;
              }
              // make the circle
              makeCircles(circ);
            }
          }
          orig_center.clear();
        }
      }
    }

    if (shapes.contains("Polygon")) {
      poly = shapes["Polygon"];
      int polySize = poly.size();
      for (auto i = 0; i < polySize; ++i) {
        // get original center of circle
        cent.push_back(poly[i]["Center"][0].as_double());
        cent.push_back(poly[i]["Center"][1].as_double());
        cent.push_back(poly[i]["Center"][2].as_double());
        orig_center.push_back(cent);
        cent.clear();
        // get the item number
        int item = 0;
        if (poly[i].contains("Item")) {
          item = poly[i]["Item"].as<int>();
          std::vector<std::vector<double>> translations;
          // search for item in map
          auto search = item_map.find(item);
          if (search != item_map.end()) {
            translations = search->second;
            for (const auto &transPair : translations) {
              // set the center
              updated_center = {orig_center[i][0] + transPair[0],
                                orig_center[i][1] + transPair[1],
                                orig_center[i][2]};

              // update the json object
              poly[i].insert_or_assign("Center", updated_center);

              // make the polygon
              makePolygons(poly);
            }
          }
        } else {
          // update the center to make array of shape
          for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
              int index = 0;
              for (const auto &p : poly.array_range()) {
                if (p.contains("Center")) {
                  // update center
                  updated_center = {orig_center[index][0] + i * dx,
                                    orig_center[index][1] + j * dy,
                                    orig_center[index][2]};

                  // update json object
                  poly[index].insert_or_assign("Center", updated_center);
                }
                index++;
              }
              // make the polygon
              makePolygons(poly);
            }
          }
          orig_center.clear();
        }
      }
    }
  }
}

//********************************************************************//
// Parses json for cartesian array then creates array of shape objects
//********************************************************************//
void NucMeshDriver::makeCartesianArray(jsoncons::json arr) {
  int nx = 0;
  int ny = 0;

  if (arr.contains("NX")) {
    nx = arr["NX"].as<int>();
  } else {
    std::cerr << "Error: Keyword 'NX' not found for rectangular array."
              << std::endl;
    throw;
  }
  if (arr.contains("NY"))
    ny = arr["NY"].as<int>();
  else {
    std::cerr << "Error: Keyword 'NY' not found for rectangular array."
              << std::endl;
    throw;
  }

  // Get the center of the array, that is the "size"
  // of the objects to be repeated
  std::vector<double> arr_center(5, 0);  // center of array
  if (arr.contains("Center")) {
    arr_center[0] = arr["Center"][0].as_double();
    arr_center[1] = arr["Center"][1].as_double();
    arr_center[2] = arr["Center"][2].as_double();
    if (arr["Center"].size() == 5) {
      arr_center[3] = arr["Center"][3].as_double();
      arr_center[4] = arr["Center"][4].as_double();
    }
  } else {
    std::cout << "Warning: Cartesian array 'Center' not defined, using [0,0,0]"
              << std::endl;
  }
  if (arr_center.size() == 5) {
    double angle = arr_center[4] * M_PI / 180.0;
    arr_center[0] = arr_center[0] + arr_center[3] * std::cos(angle);
    arr_center[1] = arr_center[1] + arr_center[3] * std::sin(angle);
  }

  // Get the array radius
  double arr_radius = 0;
  if (arr.contains("Radius")) {
    arr_radius = arr["Radius"].as_double();
  }

  // Compute the pitch of the grid
  double dx = 0;
  double dy = 0;
  dx = 2 * arr_radius / std::sqrt(2);
  dy = 2 * arr_radius / std::sqrt(2);

  // Compute the first shapes center location
  double startX = 0;
  double startY = 0;
  double startZ = arr_center[2];
  if (nx % 2 == 0)
    startX = arr_center[0] - (nx / 2 - 1) * dx - dx / 2;
  else {
    int s = (int)(nx / 2);
    startX = arr_center[0] - s * dx;
  }
  if (ny % 2 == 0)
    startY = arr_center[1] - (ny / 2 - 1) * dy - dy / 2;
  else {
    int s = (int)(ny / 2);
    startY = arr_center[1] - s * dy;
  }

  // map containing array pattern items and corresponding location
  std::map<int, std::vector<std::vector<double>>> item_map;
  int k = 0;
  double transX = 0;
  double transY = 0;
  if (arr.contains("Pattern")) {
    int pattern_size = arr["Pattern"].size();
    if (pattern_size != nx * ny) {
      std::cerr << "Error: Rectangular array 'Pattern' size not equal to NX*NY."
                << std::endl;
      throw;
    }
    for (int j = ny - 1; j > -1; --j) {
      for (int i = nx - 1; i > -1; --i) {
        int item_id = arr["Pattern"][k].as<int>();
        if (item_id > 0) {
          transX = (nx - i - 1) * dx;
          transY = j * dy;
          // search for item in map
          auto search = item_map.find(item_id);
          if (search == item_map.end())
            item_map.insert({item_id, {{transX, transY}}});
          else
            search->second.push_back({transX, transY});
        }
        k++;
      }
    }
  }

  jsoncons::json s = arr["Shapes"];  // the array of shapes
  jsoncons::json circ, poly, circInPoly;

  std::vector<double> cent, updated_center;

  // Iterate through shapes in array
  for (const auto &shapes : s.array_range()) {
    if (shapes.contains("Name"))
      std::cout << "Creating '" << shapes["Name"].as_string() << "'"
                << std::endl;
    if (shapes.contains("CirclesInPolys")) {
      circInPoly = shapes["CirclesInPolys"];
      int circlesSize = circInPoly.size();
      for (auto i = 0; i < circlesSize; ++i) {
        // get the item number
        int item = 0;
        if (circInPoly[i].contains("Item")) {
          item = circInPoly[i]["Item"].as<int>();
          std::vector<std::vector<double>> translations;
          // search for item in map
          auto search = item_map.find(item);
          if (search != item_map.end()) {
            translations = search->second;
            for (const auto &transPair : translations) {
              // set the center
              updated_center = {startX + transPair[0], startY + transPair[1],
                                startZ};

              // update the json object
              circInPoly[i].insert_or_assign("Center", updated_center);

              // make the circle
              makeCirclesInPolys(circInPoly);
            }
          }
        } else {
          // update the center to make array of shape
          for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
              int index = 0;
              for (const auto &p : circInPoly.array_range()) {
                if (p.contains("Center")) {
                  // update center
                  updated_center = {startX + i * dx, startY + j * dy, startZ};

                  // update json object
                  circInPoly[index].insert_or_assign("Center", updated_center);
                }
                index++;
              }
              // make the circleInPoly
              makeCirclesInPolys(circInPoly);
            }
          }
        }
      }
    }

    if (shapes.contains("Circles")) {
      std::cout << "contains circles" << std::endl;
      circ = shapes["Circles"];
      int circlesSize = circ.size();
      for (auto i = 0; i < circlesSize; ++i) {
        // get the item number
        int item = 0;
        if (circ[i].contains("Item")) {
          item = circ[i]["Item"].as<int>();
          std::cout << "item = " << item << std::endl;
          std::vector<std::vector<double>> translations;
          // search for item in map
          auto search = item_map.find(item);
          if (search != item_map.end()) {
            std::cout << "found item " << item << " in map" << std::endl;
            translations = search->second;
            for (const auto &transPair : translations) {
              // set the center
              std::cout << "(" << transPair[0] << ", " << transPair[1] << ")"
                        << std::endl;
              updated_center = {startX + transPair[0], startY + transPair[1],
                                startZ};

              // update the json object
              circ[i].insert_or_assign("Center", updated_center);

              // make the circle
              makeCircles(circ);
            }
          }
        } else {
          // update the center to make array of shape
          for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
              int index = 0;
              for (const auto &c : circ.array_range()) {
                if (c.contains("Center")) {
                  // update center
                  updated_center = {startX + i * dx, startY + j * dy, startZ};

                  // update json object
                  circ[index].insert_or_assign("Center", updated_center);
                }
                index++;
              }
              // make the circle
              makeCircles(circ);
            }
          }
        }
      }
    }

    if (shapes.contains("Polygon")) {
      poly = shapes["Polygon"];
      int polySize = poly.size();
      for (auto i = 0; i < polySize; ++i) {
        // get the item number
        int item = 0;
        if (poly[i].contains("Item")) {
          item = poly[i]["Item"].as<int>();
          // std::cout << "item = " << item << std::endl;
          std::vector<std::vector<double>> translations;
          // search for item in map
          auto search = item_map.find(item);
          if (search != item_map.end()) {
            // std::cout << "found item " << item << " in map" << std::endl;
            translations = search->second;
            for (const auto &transPair : translations) {
              // set the center
              updated_center = {startX + transPair[0], startY + transPair[1],
                                startZ};

              // update the json object
              poly[i].insert_or_assign("Center", updated_center);

              // make the polygon
              makePolygons(poly);
            }
          }
        } else {
          // update the center to make array of shape
          for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
              int index = 0;
              for (const auto &p : poly.array_range()) {
                if (p.contains("Center")) {
                  // update center
                  updated_center = {startX + i * dx, startY + j * dy, startZ};

                  // update json object
                  poly[index].insert_or_assign("Center", updated_center);
                }
                index++;
              }
              // make the polygon
              makePolygons(poly);
            }
          }
        }
      }
    }
  }
}

//********************************************************************//
// Parses json for polar array data then creates array of shape objects
//********************************************************************//
void NucMeshDriver::makePolarArray(jsoncons::json arr) {
  // get polar array parameters
  std::vector<double> center;  // center of array
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

  double radius = arr["Radius"].as_double();  // array radius
  int n = 0;                                  // number of array elements
  if (arr.contains("N")) n = arr["N"].as<int>();

  double start = arr["Start Angle"].as_double();  // start angle
  if (start > 360.0) {
    std::cerr << "Error: 'Start Angle' in polar array is greater than 360.0 "
                 "degrees."
              << std::endl;
    throw;
  }
  double arc = arr["Arc"].as_double();  // arc angle
  if (arc > 360.0) {
    std::cerr << "Error: 'Arc' in polar array is greater than 360.0 degrees. "
                 "Limiting 'Arc' to 360.0 degrees."
              << std::endl;
    arc = 360.0;
  }

  int nPat = 0;  // number of pattern repeats
  if (arr.contains("N Patterns")) nPat = arr["N Patterns"].as<int>();

  std::vector<int> pattern_vec;  // pattern vector
  if (arr.contains("Pattern"))
    for (int j = 0; j < nPat; ++j)
      for (int i = 0; i < arr["Pattern"].size(); ++i)
        if (arr["Pattern"][i].as<int>() > 0)
          pattern_vec.push_back(arr["Pattern"][i].as<int>());

  bool rotWithArray = false;  // rotate with array bool
  if (arr.contains("Rotate with Array")) {
    rotWithArray = arr["Rotate with Array"].as<bool>();
    if (rotWithArray)
      std::cout << "'Rotate with Array' set to true" << std::endl;
  }

  // map containing array pattern items and corresponding angles
  std::map<int, std::vector<double>> item_map;

  double inc;  // angle increment
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
    throw;
  }

  jsoncons::json s = arr["Shapes"];  // the array of shapes
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
              if (rotWithArray) poly[index].insert_or_assign("Rotation", angle);

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
              if (rotWithArray) poly[index].insert_or_assign("Rotation", angle);

              makePolygons(poly);
            }
          }
        }
      }
      index++;
    }
  }
}

//********************************************************************//
// Parses json for hex array data then creates array of shape objects
//********************************************************************//
void NucMeshDriver::makeHexagonalArray(jsoncons::json arr) {
  std::vector<double> cent, c_cen;
  std::vector<std::vector<double>> orig_center;

  std::vector<double> arr_center;  // center of array
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
    std::cerr << "Error: No 'Radius' found for array parameters." << std::endl;
    throw;
  }

  int n = 0;  // number of array elements
  if (arr.contains("N")) n = arr["N"].as<int>();

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
    padding = arr["Padding"].as_double();  // padding
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
    start = arr["Start Angle"].as_double();  // start angle
    if (start > 360.0) {
      std::cerr << "\nError: 'Start Angle' in hexagonal array is greater "
                   "than 360.0 degrees.\n"
                << std::endl;
      throw;
    }
  }
  double arc = 360.0;
  if (arr.contains("Arc")) {
    arc = arr["Arc"].as_double();  // arc angle
    if (arc > 360.0) {
      std::cerr << "Error: 'Arc' in hexagonal array is greater than 360.0 "
                   "degrees. Limiting 'Arc' to 360.0 degrees."
                << std::endl;
      arc = 360.0;
    }
  }

  // Patterning
  // Hex array objects are created with the
  // following ID's with n = 5.
  //      17   18    19
  //    13  14    15    16
  //  8    9   10   11   12
  //     4   5    6    7
  //       1    2    3
  //
  // We are putting these into a map for patterning
  int k = 0;
  int numOnSide = (n + 1) / 2;  // number of objects on a side of hex array
  int numObj =
      3 * numOnSide * (numOnSide - 1) + 1;  // total number of objects in array
  int id = 0;                               // ID of object in array
  int begin = numObj;                       // the starting ID
  int count = 0;
  std::map<int, int> item_map;

  bool pattern = false;
  if (arr.contains("Pattern")) {
    pattern = true;
    int pattern_size = arr["Pattern"].size();
    if (pattern_size != numObj) {
      std::cerr << "Error: Hexagonal array 'Pattern' size not equal to "
                << numObj << std::endl;
      throw;
    }
    // Loop from top of array to bottom
    for (int i = n; i > 0; --i) {
      int cols = 0;  // the number of objects in a given row (number of columns)
      if (i > numOnSide - 1) {
        cols = numOnSide + k;
        k++;
      } else {
        cols = n + numOnSide - 1 - k;
        k++;
      }
      begin = begin - cols;
      id = begin;
      int item_id = 0;
      for (int j = cols - 1; j > -1; --j) {
        id++;
        item_id = arr["Pattern"][count].as<int>();
        if (item_id > 0) {
          // search for item in map
          auto search = item_map.find(id);
          if (search == item_map.end())
            item_map.insert({id, item_id});
          else
            std::cout << "Error: Item already in map." << std::endl;
        }
        count++;
      }
    }
  }

  jsoncons::json circ, poly, circInPoly;
  int half = n / 2;
  double ang = 60.0 * M_PI / 180.0;
  double xOff = 0.0, yOff = 0.0;
  // If Circles are used and visible, offset to edge of circle
  jsoncons::json s = arr["Shapes"];  // the array of shapes
  for (const auto &shapes : s.array_range()) {
    if (shapes.contains("BREAK"))
      if (shapes["BREAK"] == true) openGUI();
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
    a = 30 * M_PI / 180.0;  // rad
  }

  int cnt = 1;
  double /*r, */ theta, x, y;
  for (int row = 0; row < n; row++) {
    int cols = n - std::abs(row - half);
    for (int col = 0; col < cols; col++) {
      x = (xOff * (col * 2 + 1 - cols));
      y = (yOff * (row - half));

      // convert x,y to polar coordinates
      if (y >= 0)
        theta = std::atan2(y, x);
      else
        theta = std::atan2(y, x) + 2 * M_PI;

      theta = theta * 180.0 / M_PI;
      if (theta > arc || theta < start) {
        cnt++;
        continue;
      } else {
        // Iterate through shapes in array
        double new_x, new_y, new_z;
        for (const auto &shapes : s.array_range()) {
          // If a 'Name' exists, output
          if (shapes.contains("Name")) {
            std::cout << "Creating '" << shapes["Name"].as_string() << "'"
                      << std::endl;
          }
          //********** Circles *********//
          if (shapes.contains("Circles")) {
            circ = shapes["Circles"];
            int index = 0;
            for (const auto &c : circ.array_range()) {
              int item_id = 0;
              if (c.contains("Item")) {
                item_id = c["Item"].as<int>();
                auto search = item_map.find(cnt);
                if (search != item_map.end()) {
                  if (item_id == search->second) {
                    if (c.contains("Center")) {
                      cent.push_back(c["Center"][0].as_double());
                      cent.push_back(c["Center"][1].as_double());
                      cent.push_back(c["Center"][2].as_double());
                      orig_center.push_back(cent);
                      cent.clear();
                    } else if (c.contains("Alias")) {
                      std::string alias = c["Alias"].as_string();
                      auto search = savedobj_map.find(alias);
                      if (search != savedobj_map.end()) {
                        // the alias was found in map
                        if (search->second.first == "Circle" ||
                            search->second.first == "Circles") {
                          cent.push_back(search->second.second[0]["Center"][0]
                                             .as_double());
                          cent.push_back(search->second.second[0]["Center"][1]
                                             .as_double());
                          cent.push_back(search->second.second[0]["Center"][2]
                                             .as_double());
                          orig_center.push_back(cent);
                          cent.clear();
                        }
                      }
                    } else
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
                    makeCircles(circ);
                  }
                }
              } else {
                if (c.contains("Center")) {
                  cent.push_back(c["Center"][0].as_double());
                  cent.push_back(c["Center"][1].as_double());
                  cent.push_back(c["Center"][2].as_double());
                  orig_center.push_back(cent);
                  cent.clear();
                } else if (c.contains("Alias")) {
                  std::string alias = c["Alias"].as_string();
                  auto search = savedobj_map.find(alias);
                  if (search != savedobj_map.end()) {
                    // the alias was found in map
                    if (search->second.first == "Circle" ||
                        search->second.first == "Circles") {
                      cent.push_back(
                          search->second.second[0]["Center"][0].as_double());
                      cent.push_back(
                          search->second.second[0]["Center"][1].as_double());
                      cent.push_back(
                          search->second.second[0]["Center"][2].as_double());
                      orig_center.push_back(cent);
                      cent.clear();
                    }
                  }
                } else
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
                if (pattern) makeCircles(circ);
                index++;
              }
            }
            if (!pattern) makeCircles(circ);
          }
          orig_center.clear();
          //********** Polygons *********//
          if (shapes.contains("Polygon")) {
            poly = shapes["Polygon"];
            int index = 0;
            for (const auto &p : poly.array_range()) {
              int item_id = 0;
              if (p.contains("Item")) {
                item_id = p["Item"].as<int>();
                auto search = item_map.find(cnt);
                if (search != item_map.end()) {
                  if (item_id == search->second) {
                    // Set original center to zero, if otherwise, change it
                    if (p.contains("Center")) {
                      cent.push_back(p["Center"][0].as_double());
                      cent.push_back(p["Center"][1].as_double());
                      cent.push_back(p["Center"][2].as_double());
                      orig_center.push_back(cent);
                      cent.clear();
                    } else
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
                      if (p.contains("Number of Sides")) {
                        if (p["Number of Sides"] == 6)
                          poly[index].insert_or_assign("Rotation", 30.0);
                      } else {
                        if (p.contains("Alias")) {
                          std::string alias = p["Alias"].as_string();
                          auto search = savedobj_map.find(alias);
                          if (search != savedobj_map.end()) {
                            // the alias was found in map
                            if (search->second.first == "Polygon" ||
                                search->second.first == "Hexagon") {
                              if (search->second.second[0]["Number of Sides"] ==
                                  6)
                                poly[index].insert_or_assign("Rotation", 30.0);
                            }
                          }
                        }
                      }
                    }
                    makePolygons(poly);
                  }
                }
              } else {
                // Set original center to zero, if otherwise, change it
                if (p.contains("Center")) {
                  cent.push_back(p["Center"][0].as_double());
                  cent.push_back(p["Center"][1].as_double());
                  cent.push_back(p["Center"][2].as_double());
                  orig_center.push_back(cent);
                  cent.clear();
                } else if (p.contains("Alias")) {
                  std::string alias = p["Alias"].as_string();
                  auto search = savedobj_map.find(alias);
                  if (search != savedobj_map.end()) {
                    // the alias was found in map
                    if (search->second.first == "Polygon" ||
                        search->second.first == "Hexagon") {
                      cent.push_back(
                          search->second.second[0]["Center"][0].as_double());
                      cent.push_back(
                          search->second.second[0]["Center"][1].as_double());
                      cent.push_back(
                          search->second.second[0]["Center"][2].as_double());
                      orig_center.push_back(cent);
                      cent.clear();
                    }
                  }
                } else
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
                  if (p.contains("Number of Sides")) {
                    if (p["Number of Sides"] == 6)
                      poly[index].insert_or_assign("Rotation", 30.0);
                  } else {
                    if (p.contains("Alias")) {
                      std::string alias = p["Alias"].as_string();
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
                if (pattern) makePolygons(poly);
                index++;
              }
            }
            if (!pattern) makePolygons(poly);
          }
          orig_center.clear();

          //********** Circles in Polygons *********//
          if (shapes.contains("CirclesInPolys")) {
            circInPoly = shapes["CirclesInPolys"];
            int index = 0;
            for (const auto &cp : circInPoly.array_range()) {
              int item_id = 0;
              if (cp.contains("Item")) {
                item_id = cp["Item"].as<int>();
                auto search = item_map.find(cnt);
                if (search != item_map.end()) {
                  if (item_id == search->second) {
                    // Set original center to zero, if otherwise, change it
                    if (cp.contains("Center")) {
                      cent.push_back(cp["Center"][0].as_double());
                      cent.push_back(cp["Center"][1].as_double());
                      cent.push_back(cp["Center"][2].as_double());
                      orig_center.push_back(cent);
                      cent.clear();
                    } else
                      orig_center.push_back({0, 0, 0});
                    // set the center
                    new_x = (arr_center[0] + orig_center[index][0] +
                             x * std::cos(a) - y * std::sin(a));
                    new_y = (arr_center[1] + orig_center[index][1] +
                             x * std::sin(a) + y * std::cos(a));
                    new_z = arr_center[2] + orig_center[index][2];
                    c_cen = {new_x, new_y, new_z};

                    // update json object
                    circInPoly[index].insert_or_assign("Center", c_cen);

                    // Change rotation if hex array type is 1
                    if (type == 1) {
                      if (cp.contains("Number of Sides")) {
                        if (cp["Number of Sides"] == 6)
                          circInPoly[index].insert_or_assign("Rotation", 30.0);
                      } else {
                        if (cp.contains("Alias")) {
                          std::string alias = cp["Alias"].as_string();
                          auto search = savedobj_map.find(alias);
                          if (search != savedobj_map.end()) {
                            // the alias was found in map
                            if (search->second.first == "CirclesInPolys") {
                              if (search->second.second[0]["Number of Sides"] ==
                                  6)
                                circInPoly[index].insert_or_assign("Rotation",
                                                                   30.0);
                            }
                          }
                        }
                      }
                    }
                    makeCirclesInPolys(circInPoly);
                  }
                }
              } else {
                // Set original center to zero, if otherwise, change it
                if (cp.contains("Center")) {
                  cent.push_back(cp["Center"][0].as_double());
                  cent.push_back(cp["Center"][1].as_double());
                  cent.push_back(cp["Center"][2].as_double());
                  orig_center.push_back(cent);
                  cent.clear();
                } else if (cp.contains("Alias")) {
                  std::string alias = cp["Alias"].as_string();
                  auto search = savedobj_map.find(alias);
                  if (search != savedobj_map.end()) {
                    // the alias was found in map
                    if (search->second.first == "CirclesInPolys") {
                      cent.push_back(
                          search->second.second[0]["Center"][0].as_double());
                      cent.push_back(
                          search->second.second[0]["Center"][1].as_double());
                      cent.push_back(
                          search->second.second[0]["Center"][2].as_double());
                      orig_center.push_back(cent);
                      cent.clear();
                    }
                  }
                } else
                  orig_center.push_back({0, 0, 0});
                // set the center
                new_x = (arr_center[0] + orig_center[index][0] +
                         x * std::cos(a) - y * std::sin(a));
                new_y = (arr_center[1] + orig_center[index][1] +
                         x * std::sin(a) + y * std::cos(a));
                new_z = arr_center[2] + orig_center[index][2];
                c_cen = {new_x, new_y, new_z};

                // update json object
                circInPoly[index].insert_or_assign("Center", c_cen);

                // Change rotation if hex array type is 1
                if (type == 1) {
                  if (cp.contains("Number of Sides")) {
                    if (cp["Number of Sides"] == 6)
                      circInPoly[index].insert_or_assign("Rotation", 30.0);
                  } else {
                    if (cp.contains("Alias")) {
                      std::string alias = cp["Alias"].as_string();
                      auto search = savedobj_map.find(alias);
                      if (search != savedobj_map.end()) {
                        // the alias was found in map
                        if (search->second.first == "CirclesInPolys") {
                          if (search->second.second[0]["Number of Sides"] == 6)
                            circInPoly[index].insert_or_assign("Rotation",
                                                               30.0);
                        }
                      }
                    }
                  }
                }
                if (pattern) makeCirclesInPolys(circInPoly);
                index++;
              }
            }
            if (!pattern) makeCirclesInPolys(circInPoly);
          }
          orig_center.clear();
        }
      }
      cnt++;
    }
  }
  std::cout << "Hexagonal Array created" << std::endl;
}

//********************************************************************//
// Causes the Gmsh GUI to wait for command line input
//********************************************************************//
void NucMeshDriver::openGUI() {
  // Tpause.start();
  gmsh::graphics::draw();
  if (!skipAll) {
    while (1) {
      gmsh::fltk::awake();
      gmsh::fltk::wait();
      gmsh::fltk::update();
      std::cout
          << "\nPress '1' to continue, '2' to skip all breaks, '0' to cancel"
          << std::endl;
      int in;
      std::cin >> in;
      if (in == 1)
        break;
      else if (in == 2) {
        skipAll = true;
        break;
      } else if (in == 0)
        exit(0); 
      else {
        std::cout << "\nUnknown command, try again." << std::endl;
      }
    }
  }
}

//********************************************************************//
// Checks that polygon radii are larger than last circle radius
//********************************************************************//
double NucMeshDriver::checkInscribedCircle(std::vector<double> circle_radii,
                                           std::vector<double> poly_radii,
                                           int nSides) {
  double firstPolyRadius = poly_radii[0];
  int last = circle_radii.size() - 1;
  double lastCircleRadius = circle_radii[last];

  double degrees = M_PI / nSides;
  double incircle = firstPolyRadius * std::cos(degrees);

  if (incircle < lastCircleRadius)
    return incircle;
  else
    return 1e20;
}
