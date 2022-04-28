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
#include <gtest/gtest.h>

#include <Drivers/NucMeshDriver.H>
#include <NucMesh/CirclesAndPolys.H>
#include <NucMesh/FuelElement.H>
#include <NucMesh/HexagonalArray.H>
#include <NucMesh/PolarArray.H>

#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <vtkUnstructuredGrid.h>

#include <Mesh/diffMesh.H>
#include <Mesh/exoGeoMesh.H>
#include <Mesh/geoMeshFactory.H>

namespace {
std::string hex_array_pattern_test_json;
std::string hex_array_reference;
std::string fuel_element_reference;
std::vector<std::pair<std::string, std::string>> json_and_reference;
std::string threeD_json;
std::string threeD_reference;

std::unique_ptr<NEM::DRV::NucMeshDriver> getDriver(
    const std::string &json_file) {
  std::ifstream inputStream(json_file);
  jsoncons::json inputjson;
  inputStream >> inputjson;
  auto driver = NEM::DRV::NemDriver::readJSON(inputjson);
  // This was assigned but never accessed
  // auto nucMeshDriver = dynamic_cast<NEM::DRV::NucMeshDriver *>(driver.get());
  return std::unique_ptr<NEM::DRV::NucMeshDriver>{
      dynamic_cast<NEM::DRV::NucMeshDriver *>(driver.release())};
}

}  // namespace

// Test that JSON and C++ API match
TEST(NucMesh, JSON_deserializer) {
  vtkSmartPointer<NEM::MSH::geoMeshBase> vug1, vug2;
  {  // Use C++ API
    NEM::DRV::NucMeshDriver::Opts opts{};
    using namespace NEM::NUCMESH;
    CirclesAndPolys hex{
        6,
        {{PolyRing::ShapeType::POLYGON, 0.2, RingMeshOption::ApplyTriMesh(),
          30., "B"},
         {PolyRing::ShapeType::POLYGON, 0.3,
          RingMeshOption::ApplyStructuredMesh({2, 3}), 30., "C"}}};
    CirclesAndPolys hex2{
        6,
        {{PolyRing::ShapeType::POLYGON, 0.05, RingMeshOption::ApplyTriMesh(),
          30., "B"},
         {PolyRing::ShapeType::POLYGON, 0.3,
          RingMeshOption::ApplyStructuredMesh({2, 3}), 30., "D"}}};
    auto &arr = opts.makeShape<HexagonalArray>(3, 0.3 * std::sqrt(3) + 0.1);
    arr.setCenter({-1.5, 0., 0.});
    arr.insertPatternShape(1, std::move(hex));
    arr.insertPatternShape(2, std::move(hex2));
    for (int i = 0; i < 3; ++i) {
      arr.setPatternRowCol(0, i, 1);
      arr.setPatternRowCol(4, i, 1);
    }
    for (int i = 0; i < 4; ++i) {
      arr.setPatternRowCol(1, i, 2);
      arr.setPatternRowCol(3, i, 2);
    }
    arr.setPatternCoordCenter(0, 0, 0);
    arr.setPatternCoordCenter(1, 0, 2);
    arr.setPatternCoordCenter(-1, 0, 2);
    arr.setPatternCoordCenter(2, 0, 1);
    arr.setPatternCoordCenter(-2, 0, 1);
    NEM::DRV::NucMeshDriver driver{NEM::DRV::NucMeshDriver::Files{{}},
                                   std::move(opts)};
    vug1 = driver.draw();
  }
  {  // Read from JSON
    auto driver2 = getDriver(hex_array_pattern_test_json);
    ASSERT_NE(driver2, nullptr);
    driver2->execute();
    vug2 = vtkSmartPointer<NEM::MSH::geoMeshBase>::Take(
        NEM::MSH::Read(driver2->getFiles().outputFile));
  }
  EXPECT_EQ(vug1->getNumberOfPoints(), vug2->getNumberOfPoints());
  EXPECT_EQ(vug1->getNumberOfCells(), vug2->getNumberOfCells());
  EXPECT_EQ(vug1->getNumberOfCellDataArrays(),
            vug2->getNumberOfCellDataArrays());
  EXPECT_EQ(vug1->getNumberOfPointDataArrays(),
            vug2->getNumberOfPointDataArrays());
  auto gold_mesh = vtkSmartPointer<NEM::MSH::geoMeshBase>::Take(
      NEM::MSH::Read(hex_array_reference));
  ASSERT_NE(gold_mesh, nullptr);
  EXPECT_EQ(NEM::MSH::diffMesh(gold_mesh, vug2, 1.e-9, 1.e-6, 0.03, 0.03), 0);
}

TEST(NucMesh, FuelElement) {
  vtkSmartPointer<NEM::MSH::geoMeshBase> generated;
  {
    NEM::DRV::NucMeshDriver::Opts opts{};
    using namespace NEM::NUCMESH;
    std::map<std::string, PlateMeshOption> meshType_map = {
        {"End", PlateMeshOption::ApplyStructuredMesh({3, 1})},
        {"Side", PlateMeshOption::ApplyStructuredMesh({3, 2})},
        {"Plate", PlateMeshOption::ApplyStructuredMesh({1, 12})},
        {"Fuel", PlateMeshOption::ApplyStructuredMesh({1, 10})}};
    std::map<std::string, std::string> mat_map = {
        {"End", "A"}, {"Side", "B"}, {"Plate", "C"}, {"Fuel", "D"}};
    std::vector<Plate> plates19{
        {0.2032, 0.0508, 0.47498, 0.0813, 0.4953, meshType_map, mat_map},
        {0.19812, 0.0, 0.47498, 0.0813, 0.2413, meshType_map, mat_map},
        {0.127, 0.0508, 0.47498, 0.0813, 0.2413, meshType_map, mat_map},
        {0.19812, 0.0, 0.47498, 0.0813, 0.2413, meshType_map, mat_map},
        {0.127, 0.0508, 0.47498, 0.0813, 0.2413, meshType_map, mat_map},
        {0.19812, 0.0, 0.47498, 0.0813, 0.2413, meshType_map, mat_map},
        {0.127, 0.0508, 0.47498, 0.0813, 0.2413, meshType_map, mat_map},
        {0.19812, 0.0, 0.47498, 0.0813, 0.2413, meshType_map, mat_map},
        {0.127, 0.0508, 0.47498, 0.0813, 0.2413, meshType_map, mat_map},
        {0.19812, 0.0, 0.47498, 0.0813, 0.2413, meshType_map, mat_map},
        {0.127, 0.0508, 0.47498, 0.0813, 0.2413, meshType_map, mat_map},
        {0.19812, 0.0, 0.47498, 0.0813, 0.2413, meshType_map, mat_map},
        {0.127, 0.0508, 0.47498, 0.0813, 0.2413, meshType_map, mat_map},
        {0.19812, 0.0, 0.47498, 0.0813, 0.2413, meshType_map, mat_map},
        {0.127, 0.0508, 0.47498, 0.0813, 0.2413, meshType_map, mat_map},
        {0.19812, 0.0, 0.47498, 0.0813, 0.2413, meshType_map, mat_map},
        {0.127, 0.0508, 0.47498, 0.0813, 0.4445, meshType_map, mat_map},
        {0.19812, 0.0, 0.47498, 0.0813, 0.2413, meshType_map, mat_map},
        {0.254, 0.0508, 0.47498, 0.0813, 0.4953, meshType_map, mat_map}};

    auto &fuelelement = opts.makeShape<FuelElement>(7.6581, 45, plates19);
    NEM::DRV::NucMeshDriver driver{NEM::DRV::NucMeshDriver::Files{{}},
                                   std::move(opts)};
    generated = driver.draw();

    auto gold_mesh = vtkSmartPointer<NEM::MSH::geoMeshBase>::Take(
        NEM::MSH::Read(fuel_element_reference));
    ASSERT_NE(gold_mesh, nullptr);
    EXPECT_EQ(
        NEM::MSH::diffMesh(gold_mesh, generated, 1.e-9, 1.e-6, 0.03, 0.03), 0);
  }
}

/*TEST(NucMesh, FuelAssembly) {
  vtkSmartPointer<NEM::MSH::geoMeshBase> vug1;
  {
    NEM::DRV::NucMeshDriver::Opts opts{};
    using namespace NEM::NUCMESH;
    std::map<std::string, PlateMeshOption> meshType_map = {
        {"End", PlateMeshOption::ApplyStructuredMesh({3,1})},
        {"Side", PlateMeshOption::ApplyStructuredMesh({3,2})},
        {"Plate", PlateMeshOption::ApplyStructuredMesh({1,12})},
        {"Fuel", PlateMeshOption::ApplyStructuredMesh({1,10})}};
    std::map<std::string, std::string> mat_map = {
        {"End", "A"}, {"Side", "B"}, {"Plate", "C"}, {"Fuel", "D"}};
    std::vector<Plate> plates1{
        {0.2032, 0.0508, 0.47498, 0.0813, 0.2413, meshType_map, mat_map}};
    std::vector<Plate> plates2{
        {0.2032, 0.0508, 0.47498, 0.0813, 0.2413, meshType_map, mat_map},
        {0.15, 0.0508, 0.47498, 0.0813, 0.2413, meshType_map, mat_map}};
    std::vector<Plate> plates3{
        {0.2032, 0.0508, 0.47498, 0.0813, 0.2413, meshType_map, mat_map},
        {0.2032, 0.0, 0.47498, 0.0813, 0.2413, meshType_map, mat_map},
        {0.2032, 0.0508, 0.47498, 0.0813, 0.1413, meshType_map, mat_map}};
    std::vector<Plate> plates4{
        {0.2032, 0.0508, 0.47498, 0.0813, 0.2413, meshType_map, mat_map},
        {0.2032, 0.0, 0.47498, 0.0813, 0.2413, meshType_map, mat_map},
        {0.2032, 0.0508, 0.47498, 0.0813, 0.2413, meshType_map, mat_map},
        {0.2032, 0.0, 0.47498, 0.0813, 0.1413, meshType_map, mat_map}};
    std::vector<Plate> plates18{
        {0.2032, 0.0508, 0.47498, 0.0813, 0.2413, meshType_map, mat_map},
        {0.2032, 0.0, 0.47498, 0.0813, 0.2413, meshType_map, mat_map},
        {0.2032, 0.0508, 0.47498, 0.0813, 0.2413, meshType_map, mat_map},
        {0.2032, 0.0, 0.47498, 0.0813, 0.1413, meshType_map, mat_map},
        {0.2032, 0.0508, 0.47498, 0.0813, 0.2413, meshType_map, mat_map},
        {0.2032, 0.0, 0.47498, 0.0813, 0.1413, meshType_map, mat_map},
        {0.2032, 0.0508, 0.47498, 0.0813, 0.2413, meshType_map, mat_map},
        {0.2032, 0.0, 0.47498, 0.0813, 0.1413, meshType_map, mat_map},
        {0.2032, 0.0508, 0.47498, 0.0813, 0.2413, meshType_map, mat_map},
        {0.2032, 0.0, 0.47498, 0.0813, 0.1413, meshType_map, mat_map},
        {0.2032, 0.0508, 0.47498, 0.0813, 0.2413, meshType_map, mat_map},
        {0.2032, 0.0, 0.47498, 0.0813, 0.1413, meshType_map, mat_map},
        {0.2032, 0.0508, 0.47498, 0.0813, 0.2413, meshType_map, mat_map},
        {0.2032, 0.0, 0.47498, 0.0813, 0.1413, meshType_map, mat_map},
        {0.2032, 0.0508, 0.47498, 0.0813, 0.2413, meshType_map, mat_map},
        {0.2032, 0.0, 0.47498, 0.0813, 0.1413, meshType_map, mat_map},
        {0.2032, 0.0508, 0.47498, 0.0813, 0.2413, meshType_map, mat_map},
        {0.2032, 0.0, 0.47498, 0.0813, 0.1413, meshType_map, mat_map}};
    //FuelElement fe{7.5, 45, plates18 };
    auto &fuelelement = opts.makeShape<FuelElement>(7.5, 45, plates2);

    *//*auto &arr = opts.makeShape<PolarArray>(4, 0, 180, 0, true);
    arr.setCenter({0,0,0});
    arr.insertPatternShape(0, fe);
    arr.insertPatternShape(1, fe);
    arr.insertPatternShape(2, fe);
    arr.insertPatternShape(3, fe);*//*

    NEM::DRV::NucMeshDriver driver{NEM::DRV::NucMeshDriver::Files{{"fuel_assembly.vtu"}},
                                   std::move(opts)};
    //vug1 = driver.draw();
    //driver.execute();
  }
}*/

// Regression tests
TEST(NucMesh, Match_Reference) {
  for (auto &test_case : json_and_reference) {
    std::cout << "Running " << test_case.first << '\n';
    auto nucMeshDriver = getDriver(test_case.first);
    EXPECT_NE(nucMeshDriver, nullptr);
    if (!nucMeshDriver) { continue; }
    nucMeshDriver->execute();
    auto generated_mesh = vtkSmartPointer<NEM::MSH::geoMeshBase>::Take(
        NEM::MSH::Read(nucMeshDriver->getFiles().outputFile));
    EXPECT_NE(generated_mesh, nullptr);
    auto gold_mesh = vtkSmartPointer<NEM::MSH::geoMeshBase>::Take(
        NEM::MSH::Read(test_case.second));
    EXPECT_NE(gold_mesh, nullptr);
    if (gold_mesh && generated_mesh) {
      // The tolerance is so large because MEFISTO2 C and Fortran version
      // differences (despite being translated using f2c)
      EXPECT_EQ(
          NEM::MSH::diffMesh(gold_mesh, generated_mesh, 1.e-9, 1.e-6, 0.1, 0.1),
          0);
    }
  }
}

// The test case used relies on netgen for quad-dominant meshing
#ifdef HAVE_NGEN
TEST(NucMesh, Extrusion) {
  auto driver = getDriver(threeD_json);
  ASSERT_NE(driver, nullptr);
  driver->execute();
  auto generated_mesh = vtkSmartPointer<NEM::MSH::exoGeoMesh>::Take(
      NEM::MSH::exoGeoMesh::Read(driver->getFiles().outputFile));
  auto gold_mesh = vtkSmartPointer<NEM::MSH::exoGeoMesh>::Take(
      NEM::MSH::exoGeoMesh::Read(threeD_reference));
  EXPECT_EQ(gold_mesh->getNumSideSets(), generated_mesh->getNumSideSets());
  EXPECT_EQ(gold_mesh->getNumElemBlocks(), generated_mesh->getNumElemBlocks());
}
#endif

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc >= 3 && argc % 2 == 1);
  hex_array_pattern_test_json = argv[1];
  hex_array_reference = argv[2];
  fuel_element_reference = argv[3];
  for (int i = 4; i < argc - 1; i += 2) {
#ifdef HAVE_NGEN
    if (i == argc - 2) {
      threeD_json = argv[i];
      threeD_reference = argv[i + 1];
      continue;
    }
#endif
    json_and_reference.emplace_back(argv[i], argv[i + 1]);
  }
  int res = RUN_ALL_TESTS();

  return res;
}
