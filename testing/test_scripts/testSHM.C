#include <gtest/gtest.h>
#include <algorithm>
#include <boost/algorithm/string/replace.hpp>
#include <boost/filesystem.hpp>
#include <fstream>
#include <iterator>
#include <string>
#include <Drivers/MeshGen/SnappyMeshMeshGenDriver.H>
#include <Drivers/NemDriver.H>
#include <MeshGeneration/snappymeshGen.H>
#include <MeshGeneration/snappymeshParams.H>
#include <Mesh/vtkMesh.H>
namespace fs = boost::filesystem;

const char *inp_json;
meshBase *mesh;
meshBase *ref;
jsoncons::json inputjson;

// Aux functions
bool compareFiles(const std::string &p1, const std::string &p2) {
  std::ifstream f1(p1, std::ifstream::binary | std::ifstream::ate);
  std::ifstream f2(p2, std::ifstream::binary | std::ifstream::ate);

  if (f1.fail() || f2.fail()) {
    return false;  // file problem
  }

  if (f1.tellg() != f2.tellg()) {
    return false;  // size mismatch
  }

  // seek back to beginning and use std::equal to compare contents
  f1.seekg(0, std::ifstream::beg);
  f2.seekg(0, std::ifstream::beg);
  return std::equal(std::istreambuf_iterator<char>(f1.rdbuf()),
                    std::istreambuf_iterator<char>(),
                    std::istreambuf_iterator<char>(f2.rdbuf()));
}

void copyDirectoryRecursively(const fs::path &sourceDir,
                              const fs::path &destinationDir) {
  if (!fs::exists(sourceDir) || !fs::is_directory(sourceDir)) {
    throw std::runtime_error("Source directory " + sourceDir.string() +
                             " does not exist or is not a directory");
  }
  if (fs::exists(destinationDir)) {
    throw std::runtime_error("Destination directory " +
                             destinationDir.string() + " already exists");
  }
  if (!fs::create_directory(destinationDir)) {
    throw std::runtime_error("Cannot create destination directory " +
                             destinationDir.string());
  }

  for (const auto &dirEnt : fs::recursive_directory_iterator{sourceDir}) {
    const auto &path = dirEnt.path();
    auto relativePathStr = path.string();
    boost::algorithm::replace_first(relativePathStr, sourceDir.string(), "");
    fs::copy(path, destinationDir / relativePathStr);
  }
}

void snappyLogistics() {
  const char dir_path[] = "./constant/polyMesh";
  const char cpydir_path[] = "./constant/polyMesh.Orig";
  boost::filesystem::path dir1(dir_path);
  boost::filesystem::path dir2(cpydir_path);
  boost::filesystem::remove_all(dir1);
  copyDirectoryRecursively(dir2, dir1);

  return;
}

// Test implementations
int generate(const char *jsonF) {
  std::string fname(jsonF);
  std::ifstream inputStream(fname);
  if (!inputStream.good()) {
    std::cerr << "Error opening file " << jsonF << std::endl;
    exit(1);
  }

  jsoncons::json inputjson_tmp;
  inputStream >> inputjson_tmp;
  inputjson = inputjson_tmp[0];

  auto driver = std::unique_ptr<NEM::DRV::SnappyMeshMeshGenDriver>(
      dynamic_cast<NEM::DRV::SnappyMeshMeshGenDriver *>(
          NEM::DRV::NemDriver::readJSON(inputjson).release()));
  EXPECT_NE(driver, nullptr);

  auto paramsCopy = driver->getParams();
  snappymeshGen generator{&paramsCopy};
  // Parameter not used
  generator.createMeshFromSTL(nullptr);
  const auto &inputGeoFile = driver->getFiles().inputGeoFile;
  std::string newname =
      inputGeoFile.substr(0, inputGeoFile.find_last_of('.')) + ".vtu";
  mesh = meshBase::Create(generator.getDataSet(), newname);
  mesh->setFileName(driver->getFiles().outputMeshFile);
  mesh->report();
  mesh->write();

  return 0;
}

// Write a new name for snappyHexMesh
// TEST macros
TEST(snappyHexMesh, Generation) { EXPECT_EQ(0, generate(inp_json)); }

TEST(snappyHexMesh, NumberOfNodes) {
  if (ref) delete ref;
  ref = meshBase::Create(inputjson["Reference File"].as<std::string>());
  EXPECT_EQ(mesh->getNumberOfPoints(), ref->getNumberOfPoints());
}

TEST(snappyHexMesh, NumberOfCells) {
  if (ref) delete ref;
  ref = meshBase::Create(inputjson["Reference File"].as<std::string>());
  EXPECT_EQ(mesh->getNumberOfCells(), ref->getNumberOfCells());
}

// test constructor
int main(int argc, char **argv) {
  // IO
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc >= 1);
  inp_json = argv[1];

  if (!inp_json) {
    std::cerr << "No input file defined" << std::endl;
  }

  int res = RUN_ALL_TESTS();

  // clean up
  if (mesh) delete mesh;

  snappyLogistics();  // clean up and replace mesh

  return res;
}
