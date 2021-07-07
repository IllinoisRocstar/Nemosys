#if defined(_MSC_VER) && !defined(_USE_MATH_DEFINES)
#  define _USE_MATH_DEFINES
#endif
#include <gtest/gtest.h>
#include <gmsh.h>

#include <Mesh/geoMeshBase.H>

class testGeoMeshBase : public NEM::MSH::geoMeshBase {
  using geoMeshBase::geoMeshBase;

  void write(const std::string &) override {}
  void report(std::ostream &) const override {}

  void resetNative() override {}
};

TEST(geoMeshBase, ConstructorDefault) {
//  auto *pgmb = new testGeoMeshBase{};
  testGeoMeshBase gmb{};
//  delete pgmb;
}

TEST(geoMeshBase, ConstructorDefaultGMSHValidity) {
  testGeoMeshBase gmb{};
  { testGeoMeshBase gmb2{}; }
  gmsh::model::add("test");
}

TEST(geoMeshBase, ExternalGMSHValidity) {
  NEM::MSH::GmshInterface::Initialize();

  gmsh::model::add("test1");
  {
    testGeoMeshBase gmb{};
    gmsh::model::add("test2");
  }
  gmsh::model::add("test3");

  NEM::MSH::GmshInterface::Finalize();
}

TEST(geoMeshBase, GetSetGeoEntArrayName) {
  testGeoMeshBase gmb{};
  std::string geoEnt = "geoEnt";

  gmb.setGeoEntArrayName(geoEnt);

  ASSERT_EQ(gmb.getGeoEntArrayName(), geoEnt);
}

/*
TEST(geoMeshBase, computeDiscreteGeoFromMsh) {
  NEM::MSH::GmshInterface::Initialize();
  gmsh::option::setNumber("General.Terminal", 1.0);
  gmsh::model::add("geoMeshBase");

  vtkSmartPointer<vtkUnstructuredGrid> ug =
      vtkSmartPointer<vtkUnstructuredGrid>::New();
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  // 0-7: for a hex
  points->InsertNextPoint(0, 0, 0);  // 0
  points->InsertNextPoint(1, 0, 0);  // 1
  points->InsertNextPoint(1, 1, 0);
  points->InsertNextPoint(0, 1, 0);  // 2
  points->InsertNextPoint(0, 0, 1);  // 3
  points->InsertNextPoint(1, 0, 1);
  points->InsertNextPoint(1, 1, 1);
  points->InsertNextPoint(0, 1, 1);

  // 8-11: for qquad
  points->InsertNextPoint(0.5, 0, -0.01);
  points->InsertNextPoint(1, 0.5, -0.01);
  points->InsertNextPoint(0.5, 1, -0.01);
  points->InsertNextPoint(0, 0.5, -0.01);

  // 12-15: for qtet
  points->InsertNextPoint(0.5, 0.5, -0.01);
  points->InsertNextPoint(0, -0.01, 0.5);
  points->InsertNextPoint(0.5, -0.01, 0.5);
  points->InsertNextPoint(-0.01, 0.5, 0.5);

  // 16: for second tet
  points->InsertNextPoint(0, 0, -1);

  ug->SetPoints(points);

  // vtkIdType quad1[4] = {0, 1, 2, 3};
  // vtkIdType qquad1[8] = {0, 1, 2, 3, 8, 9, 10, 11};
  // vtkIdType hex1[8] = {0, 1, 2, 3, 4, 5, 6, 7};
  // vtkIdType qtet[10] = {0, 1, 3, 4, 8, 12, 11, 13, 14, 15};
  vtkIdType tet1[4] = {0, 1, 3, 4};
  vtkIdType tet2[4] = {0, 3, 1, 16};
  // ug->InsertNextCell(VTK_QUADRATIC_QUAD, 8, qquad1);
  // ug->InsertNextCell(VTK_HEXAHEDRON, 8, hex1);
  // ug->InsertNextCell(VTK_QUADRATIC_TETRA, 10, qtet);
  ug->InsertNextCell(VTK_TETRA, 4, tet1);
  ug->InsertNextCell(VTK_TETRA, 4, tet2);

  vtkSmartPointer<vtkIntArray> geoIds = vtkSmartPointer<vtkIntArray>::New();
  geoIds->SetName("GeoIds");
  geoIds->InsertNextValue(0);
  geoIds->InsertNextValue(1);
  ug->GetCellData()->AddArray(geoIds);

  NEM::MSH::geoMeshBase::computeDiscreteGeoFromMsh(ug, geoIds->GetName());
//  NEM::MSH::geoMeshBase::computeDiscreteGeoFromMsh(ug);

  gmsh::model::mesh::generate(3);

  gmsh::write("computeDiscreteGeoFromMsh.msh");
  NEM::MSH::GmshInterface::Finalize();
}
*/

int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);

  int res = RUN_ALL_TESTS();

  return res;
}
