#define _USE_MATH_DEFINES
#include "NucMeshDriver.H"

#include <gmsh.h>

#include "AuxiliaryFunctions.H"
#include "circles.H"
#include "polygon.H"

NucMeshDriver::~NucMeshDriver() {
  std::cout << "NucMeshDriver destroyed" << std::endl;
}

NucMeshDriver *NucMeshDriver::readJSON(const jsoncons::json &inputjson) {
  std::cout << "Reading Input JSON File." << std::endl;

  if (inputjson.contains("Geometry and Mesh")) {
    // json shapes = inputjson["Geometry and Mesh"];

    auto nucmeshobj = new NucMeshDriver(inputjson);
    return nucmeshobj;
  } else {
    std::cerr << "'Geometry and Mesh' keyword not found." << std::endl;
    exit(-1);
  }
}

NucMeshDriver::NucMeshDriver(jsoncons::json inputjson) {
  id = 1;
  physTag = 1;

  std::cout << "NucMeshDriver created" << std::endl;

  if (inputjson.contains("Output File Name")) {
    ofname = inputjson["Output File Name"].as_string();
  } else {
    std::cerr << "'Output File Name' keyword not found." << std::endl;
    exit(-1);
  }

  if (inputjson.contains("Extension")) {
    extension = inputjson["Extension"].as_string();
    // TODO: add conditions for acceptable file extensions
  } else {
    std::cerr << "'Extension' keyword not found." << std::endl;
    exit(-1);
  }

  jsoncons::json shapes = inputjson["Geometry and Mesh"];

  nemAux::Timer Tgeom, Tboolean, Textrude, Tmesh;
  Tgeom.start();

  gmsh::initialize();
  gmsh::model::add("NucMesh");
  gmsh::option::setNumber("General.NumThreads", 10);
  gmsh::option::setNumber("Geometry.OCCBooleanPreserveNumbering", 1);
  gmsh::option::setNumber("Geometry.OCCParallel", 1);
  gmsh::option::setNumber("Mesh.FlexibleTransfinite", 0);
  gmsh::option::setNumber("Mesh.MshFileVersion", 2.2);

  // Loop through Geometry and Mesh to get shapes
  std::cout << "Geometry and Mesh" << std::endl;
  for (const auto &obj : shapes.array_range()) {
    if (obj.contains("Global Options")) parseOptions(obj["Global Options"]);

    if (obj.contains("Saved Objects")) parseSavedObjects(obj["Saved Objects"]);

    if (obj.contains("Visible"))
      if (obj["Visible"] == false) continue;

    if (obj.contains("Show"))
      if (obj["Show"] == false) continue;

    if (obj.contains("Circles")) makeCircles(obj["Circles"]);

    if (obj.contains("Circle")) makeCircles(obj["Circle"]);

    if (obj.contains("Hexagon") || obj.contains("Polygon")) {
      if (obj.contains("Hexagon"))
        makePolygons(obj["Hexagon"], 6);
      else
        makePolygons(obj["Polygon"]);
    }

    if (obj.contains("Array")) makeArray(obj);
  }
  Tgeom.stop();

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
  std::pair<int, int> oldNew;                   // holds the old/new id pair
  std::vector<std::pair<int, int>> oldNew_vec;  // vector for id pair

  // gather surface ids that have changed for Boolean Fragments
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

  for (auto itr = shape_map.begin(); itr != shape_map.end(); ++itr)
    itr->second->updateSurfaces(oldNew_vec);

  std::cout << "Updated surface IDs\n";

  // iterate through shape_map to apply meshType and collect physical Surfs
  std::map<int, int> physSurf_map;
  for (auto &&itr : shape_map) {
    itr.second->applyMeshType();

    physSurf_map = itr.second->getPhysSurf(phystag_map, physSurf_map);
  }

  std::cout << "Mesh type applied\n";

  // Make physical groups
  for (const auto &itr : phystag_map) {
    std::vector<int> s;  // vector of surfaces
    std::string name;    // regions name
    int tag;             // physical group tag
    for (const auto &itr2 : physSurf_map) {
      if (itr2.second == itr.second) {
        tag = itr2.second;
        name = itr.first;
        s.push_back(itr2.first);
      }
    }
    gmsh::model::addPhysicalGroup(2, s, tag);
    gmsh::model::setPhysicalName(2, tag, name);
  }
  std::cout << "Physical groups applied\n";

  gmsh::model::occ::synchronize();
  Tboolean.stop();

  Textrude.start();
  // std::vector<std::pair<int, int>> entities, outEntities;
  // std::vector<int> nLayers = {30};
  // std::vector<double> layers = {1};
  // gmsh::model::getEntities(entities, 2);
  // gmsh::model::occ::extrude(entities, 0, 0, 30, outEntities, nLayers, layers,
  //                          true);
  // gmsh::model::occ::synchronize();
  Textrude.stop();

  Tmesh.start();
  // Mesh 1D to get wireframe, save as vtk
  gmsh::model::mesh::generate(1);
  gmsh::option::setNumber("Mesh.SaveAll", 1);
  gmsh::write(ofname + ".vtk");

  // Mesh 2D, save as msh
  gmsh::model::mesh::generate(2);
  gmsh::model::mesh::removeDuplicateNodes();
  gmsh::option::setNumber("Mesh.SaveAll", 0);
  gmsh::write(ofname + ".msh");
  Tmesh.stop();

  // Convert mesh
  // convertMsh(ofname);

  double s = 1000.0;
  std::cout << "Geometry creation time = " << Tgeom.elapsed() / s
            << " seconds\n";
  std::cout << "Geometry coherence time = " << Tboolean.elapsed() / s
            << " seconds\n";
  std::cout << "Extrusion time = " << Textrude.elapsed() / s << " seconds\n";
  std::cout << "Meshing time = " << Tmesh.elapsed() / s << " seconds\n";
  std::cout << "Total run time = "
            << Tgeom.elapsed() / s + Tboolean.elapsed() / s +
                   Textrude.elapsed() / s + Tmesh.elapsed() / s
            << " seconds\n";

  // gmsh::fltk::run();
  gmsh::finalize();
}

// Parses the global options data of json
void NucMeshDriver::parseOptions(jsoncons::json opts) {
  // Iterate through options
  for (const auto &it : opts.array_range()) {
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
      else if (algo == "Frontal Quads" || "Frontal Quad")
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
    if (it.contains("Mesh Smoothing"))
      gmsh::option::setNumber("Mesh.Smoothing", it["Mesh Smoothing"].as<int>());
  }
}

// Parses the saved objects data of json
void NucMeshDriver::parseSavedObjects(jsoncons::json savedObj) {
  std::cout << "Parsing Saved Objects..." << std::endl;

  // Iterate through Saved Objects key
  for (const auto &it : savedObj.array_range()) {
    std::string alias;
    std::string shapetype;

    if (it.contains("Alias")) {
      alias = it["Alias"].as<std::string>();
      std::cout << "Alias '" << alias << "' found\n";
    }

    for (auto i = 1; i < it.size(); ++i) {
      if (it.contains("Circles")) shapetype = "Circles";
      if (it.contains("Circle")) shapetype = "Circle";
      if (it.contains("Polygon")) shapetype = "Polygon";

      savedobj_map.insert(
          std::pair<std::string, std::pair<std::string, jsoncons::json>>(
              alias, {shapetype, it[shapetype]}));
    }
  }
  std::cout << "Saved Objects Parsed\n";
}

// Parses json for polygon data then creates polygon object
void NucMeshDriver::makePolygons(jsoncons::json poly, int ns) {
  // Iterate through Polygon/Hexagon key
  for (const auto &it : poly.array_range()) {
    // If there is an alias check for new data under the alias
    if (it.contains("Alias")) {
      // Get the alias name
      std::string alias = it["Alias"].as_string();

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
      if (it.contains("Rotation"))
        // Get the rotation
        rot = it["Rotation"].as<double>();

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
          // std::cout << alias << " found in map" << std::endl;
          if (it.contains("Center"))
            search->second.second[0].insert_or_assign("Center", cen);
          if (it.contains("Rotation"))
            search->second.second[0].insert_or_assign("Rotation", rot);
          if (it.contains("Radii"))
            search->second.second[0].insert_or_assign("Radii", rad);
          if (it.contains("Mesh Type"))
            search->second.second[0].insert_or_assign("Mesh Type", type);
          if (it.contains("Number of Elems"))
            search->second.second[0].insert_or_assign("Number of Elems", elems);
          if (it.contains("Region Names"))
            search->second.second[0].insert_or_assign("Region Names", name);
          makePolygons(search->second.second);
        }
      } else {
        std::cerr << alias << " not found in map" << std::endl;
      }
    } else {
      if (it["Visible"] == true) {
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
    }
  }
}

// Parses json for circle data then creates circle object
void NucMeshDriver::makeCircles(jsoncons::json circ) {
  // Iterate through Circles key
  for (const auto &it : circ.array_range()) {
    if (it.contains("Alias")) {
      std::string alias = it["Alias"].as_string();

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
            search->second.second[0].insert_or_assign("Number of Elems", elems);
          if (it.contains("Region Names"))
            search->second.second[0].insert_or_assign("Region Names", name);

          makeCircles(search->second.second);
        }
      } else {
        std::cout << alias << " not found in map" << std::endl;
      }
    } else {
      if (it["Visible"] == true) {
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

        shape_map.insert(std::pair<int, NEM::GEO::shape *>(id, c));
        id++;
      }
    }
  }
}

// Parses json for array data then creates array of shape objects
void NucMeshDriver::makeArray(jsoncons::json arr) {
  std::string arrType = arr["Array"].as_string();

  //--------------- Rectangular Array --------------//
  //------------------------------------------------//
  if (arrType == "Rectangular") {
    std::cout << "Rectangular Array declared" << std::endl;

    int nx = arr["NX"].as<int>();
    int ny = arr["NY"].as<int>();
    double dx = arr["DX"].as<double>();
    double dy = arr["DY"].as<double>();

    jsoncons::json s = arr["Shapes"];  // the array of shapes
    jsoncons::json circ, poly;

    std::vector<double> cent, c_cen;
    std::vector<std::vector<double>> oc;

    // Iterate through shapes in array
    for (const auto &shapes : s.array_range()) {
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
    std::cout << "Polar Array declared" << std::endl;

    // get polar array parameters
    std::vector<double> center;  // center of array
    center.push_back(arr["Center"][0].as_double());
    center.push_back(arr["Center"][1].as_double());
    center.push_back(arr["Center"][2].as_double());

    double radius = arr["Radius"].as_double();  // array radius
    int n = 0;                                  // number of array elements
    if (arr.contains("N")) n = arr["N"].as<int>();

    double start = arr["Start Angle"].as_double();  // start angle
    double arc = arr["Arc"].as_double();            // arc angle
    // TODO: limit arc to 360 degrees

    int nPat = 0;  // number of pattern repeats
    if (arr.contains("N Patterns")) nPat = arr["N Patterns"].as<int>();

    std::vector<int> pattern_vec;  // pattern vector
    if (arr.contains("Pattern"))
      for (int j = 0; j < nPat; ++j)
        for (int i = 0; i < arr["Pattern"].size(); ++i)
          if (arr["Pattern"][i].as<int>() > 0)
            pattern_vec.push_back(arr["Pattern"][i].as<int>());

    bool rotWithArray = false;  // rotate with array bool
    if (arr.contains("Rotate with Array"))
      rotWithArray = arr["Rotate with Array"].as<bool>();

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
      std::cout << "'N Patterns' and 'N' cannot both be greater than zero.\n";
      exit(1);
    }

    jsoncons::json s = arr["Shapes"];  // the array of shapes
    jsoncons::json circ, poly;

    std::vector<double> cent, c_cen;
    std::vector<std::vector<double>> oc;

    // Iterate through shapes in array
    for (const auto &shapes : s.array_range()) {
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
          }
        }
        index++;
      }
    }
  }
}
