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
#define _USE_MATH_DEFINES
#include "Drivers/TemplateMeshDriver.H"

//#include <cstdio>
#include <gmsh.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <streambuf>

#include "AuxiliaryFunctions.H"

namespace NEM {
namespace DRV {

TemplateMeshDriver::Opts::Opts(std::string templateName)
    : templateName(std::move(templateName)) {}

TemplateMeshDriver::TemplateMeshDriver(Files file, Opts opts)
    : file_(std::move(file)), opts_(std::move(opts)) {
  std::cout << "TemplateMeshDriver created" << std::endl;
}

TemplateMeshDriver::TemplateMeshDriver()
    : TemplateMeshDriver(Files{{}}, Opts{{}}) {}

TemplateMeshDriver::~TemplateMeshDriver() {
  std::cout << "TemplateMeshDriver destroyed" << std::endl;
}

const TemplateMeshDriver::Files &TemplateMeshDriver::getFiles() const {
  return file_;
}

void TemplateMeshDriver::setFiles(Files file) { this->file_ = std::move(file); }

const TemplateMeshDriver::Opts &TemplateMeshDriver::getOpts() const {
  return opts_;
}

void TemplateMeshDriver::setOpts(Opts opts) { this->opts_ = std::move(opts); }

jsoncons::string_view TemplateMeshDriver::getProgramType() const {
  return programType;
}

void TemplateMeshDriver::execute() const {
  static const std::map<std::string, std::string> tpl_map = {
      {"Spiral Tape Pipe", "spiral_tape_pipe.tpl"}};

  std::string params_name{};
  // Generate parameters file for template
  if (this->opts_.templateName == "Spiral Tape Pipe") {
    params_name = spiralTapePipe(this->opts_.templateParams);
  } else {
    std::cerr << "Error: Template name not found" << std::endl;
    exit(-1);
  }

  gmsh::initialize();
  gmsh::option::setNumber("General.Terminal", 1);

  std::cout << "Generating mesh from template" << std::endl;
  nemAux::Timer T;
  T.start();
  // Execute the template script
  executeTemplate(tpl_map.at(this->opts_.templateName), this->file_.outputFile,
                  params_name);
  T.stop();
  std::cout << "Meshing time = " << (T.elapsed() / 1000.0) << " seconds\n";
}

void TemplateMeshDriver::executeTemplate(std::string tplName,
                                         std::string outName,
                                         std::string params) {
  // Decrypt the tpl and return the file name
  std::string name = decrypt(tplName);
  // Insert the name of the parameters file into the decrypted tpl
  insertParams(name, params);

  // Open and mesh the tpl
  gmsh::open(name);
  std::remove(name.c_str());
  // gmsh::fltk::run();
  gmsh::model::mesh::generate(2);
  gmsh::model::mesh::generate(3);
  gmsh::model::mesh::removeDuplicateNodes();
  gmsh::write(outName + ".msh");
  gmsh::finalize();
}

std::string TemplateMeshDriver::spiralTapePipe(jsoncons::json inputjson) {
  std::cout << "Reading Template Parameters" << std::endl;

  double rx = 1;                // Ellipse radius x-dir
  double ry = 1;                // Ellipse radius y-dir
  double thickness = 0.05;      // Thickness of spiral tape in y-dir
  double extrude_len = 3;       // Extrusion length
  double n_turns = 0.5;         // Number of spiral turns
  double width_percent = 0.85;  // Percentage of x-dir diameter for tape width
  double mSize_min = 0.048;     // Minimum mesh size, size a walls
  double mSize_max = 0.1;       // Maximum mesh size, size in bulk
  double dist_min = 0.05;       // Distance from walls with mSize_min
  double dist_max = 0.2;        // Distance from walls with mSize_max
  double bl_wall_n = 0.0038;    // Boundary layer mesh size normal to wall
  double bl_far = 0.08;         // Mesh size away from wall
  double bl_thickness = 0.02;   // Boundary layer mesh thickness
  double ratio = 1.3;           // Mesh size ratio normal to wall
  double extrude_layers = 20;   // Number of extruded elements during extrusion

  int element_order = 1;  // Finite element order
  int fan_points = 3;     // Number of fan points in boundary layer at corners

  if (inputjson.contains("rx")) rx = inputjson["rx"].as_double();
  if (inputjson.contains("ry")) ry = inputjson["ry"].as_double();
  if (inputjson.contains("thickness"))
    thickness = inputjson["thickness"].as_double();
  if (inputjson.contains("extrude_len"))
    extrude_len = inputjson["extrude_len"].as_double();
  if (inputjson.contains("n_turns")) n_turns = inputjson["n_turns"].as_double();
  if (inputjson.contains("width_percent"))
    width_percent = inputjson["width_percent"].as_double();

  if (inputjson.contains("dist_min"))
    dist_min = inputjson["dist_min"].as_double();
  if (inputjson.contains("dist_max"))
    dist_max = inputjson["dist_max"].as_double();
  if (inputjson.contains("mSize_min"))
    mSize_min = inputjson["mSize_min"].as_double();
  if (inputjson.contains("mSize_max"))
    mSize_max = inputjson["mSize_max"].as_double();
  if (inputjson.contains("bl_wall_n"))
    bl_wall_n = inputjson["bl_wall_n"].as_double();
  if (inputjson.contains("bl_far")) bl_far = inputjson["bl_far"].as_double();
  if (inputjson.contains("bl_thickness"))
    bl_thickness = inputjson["bl_thickness"].as_double();
  if (inputjson.contains("ratio")) ratio = inputjson["ratio"].as_double();
  if (inputjson.contains("fan_points"))
    fan_points = inputjson["fan_points"].as<int>();
  if (inputjson.contains("extrude_layers"))
    extrude_layers = inputjson["extrude_layers"].as_double();
  if (inputjson.contains("element_order"))
    element_order = inputjson["element_order"].as<int>();

  std::ofstream out;

  time_t t = time(0);
  struct tm *now = localtime(&t);
  char paramsName[80];
  strftime(paramsName, 80, "params_%Y-%m-%d-%H:%M:%S.txt", now);

  out.open(paramsName);
  out << "rx = " << rx << ";" << std::endl;
  out << "ry = " << ry << ";" << std::endl;
  out << "thickness = " << thickness << ";" << std::endl;
  out << "extrude_len = " << extrude_len << ";" << std::endl;
  out << "n_turns = " << n_turns << ";" << std::endl;
  out << "width_percent = " << width_percent << ";" << std::endl;

  out << "dist_min = " << dist_min << ";" << std::endl;
  out << "dist_max = " << dist_max << ";" << std::endl;
  out << "mSize_min = " << mSize_min << ";" << std::endl;
  out << "mSize_max = " << mSize_max << ";" << std::endl;
  out << "bl_wall_n = " << bl_wall_n << ";" << std::endl;
  out << "bl_far = " << bl_far << ";" << std::endl;
  out << "bl_thickness = " << bl_thickness << ";" << std::endl;
  out << "ratio = " << ratio << ";" << std::endl;
  out << "fan_points = " << fan_points << ";" << std::endl;
  out << "extrude_layers = " << extrude_layers << ";" << std::endl;
  out << "element_order = " << element_order << ";" << std::endl;
  out.close();

  return paramsName;
}

void TemplateMeshDriver::insertParams(std::string file,
                                      std::string params_name) {
  std::ifstream in(file);
  std::string str;
  std::vector<std::string> text;
  while (!in.eof()) {
    std::getline(in, str, '\n');

    text.push_back(str);
    std::string a = "SetFactory(\"built-in\");";
    std::string b = "SetFactory(\"OpenCASCADE\");";
    if (a.compare(str) == 0 || b.compare(str) == 0) {
      text.push_back("Include \"" + params_name + "\";");
    }
  }
  in.close();

  std::ofstream out(file);
  for (int i = 0; i < text.size(); ++i) {
    out << text[i] << std::endl;
  }
  out.close();
}

// TODO: This method should have empty implementation in public_hub
std::string TemplateMeshDriver::decrypt(std::string filename) {
  std::ifstream in(filename, std::ios_base::binary);
  if (!in.is_open()) {
    std::cerr << "Error: Template file " << filename << " is not open!"
              << std::endl;
    exit(-1);
  } else {
    unsigned vsize;
    in.read(reinterpret_cast<char *>(&vsize), sizeof(unsigned));
    std::vector<double> stuff(vsize);
    in.read(reinterpret_cast<char *>(&stuff[0]), vsize * sizeof(double));

    std::string result;
    for (double &d : stuff) {
      int inv_d = (int)std::sqrt(d);
      result += inv_d;
    }

    size_t lastindex = filename.find_last_of(".");
    std::string rawname = filename.substr(0, lastindex);
    std::string name = rawname + ".geo";
    std::ofstream ascii(name);
    ascii << result << std::endl;
    ascii.close();

    return name;
  }
}

// TODO: This method should have empty implementation in public_hub
void TemplateMeshDriver::encrypt(std::string inFile, std::string outFile) {
  std::ifstream file(inFile);
  std::string str;
  file.seekg(0, std::ios::end);
  str.reserve(file.tellg());
  file.seekg(0, std::ios::beg);
  str.assign((std::istreambuf_iterator<char>(file)),
             std::istreambuf_iterator<char>());

  int a;
  double b;
  std::vector<double> line;
  for (char const &c : str) {
    a = (int)c;
    b = std::pow(a, 2);
    line.push_back(b);
  }
  unsigned linesize = line.size();
  std::ofstream o(outFile, std::ios_base::binary);
  o.write(reinterpret_cast<char *>(&linesize), sizeof(unsigned));
  o.write(reinterpret_cast<char *>(&line[0]), line.size() * sizeof(double));
  o.close();
  line.clear();
}

}  // namespace DRV
}  // namespace NEM
