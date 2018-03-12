#include <NemDriver.H>
#include <gtest.h>

const char* defaultjson;
const char* defaultRef;
const char* unifjson;
const char* unifRef;
const char* geomjson;
const char* geomRef;

int genTest(const char* jsonF, const char* ofname, const char* refname)
{
  std::string fname(jsonF);
  std::ifstream inputStream(fname);
  if(!inputStream.good())
  {
    std::cerr << "Error opening file " << defaultjson << std::endl;
    exit(1);
  }
  
  json inputjson;
  inputStream >> inputjson;
  for (const auto& prog : inputjson.array_range())
  {
    std::unique_ptr<NemDriver> nemdrvobj 
      = std::unique_ptr<NemDriver>(NemDriver::readJSON(prog));
  } 
  
  std::unique_ptr<meshBase> refMesh = meshBase::CreateUnique(refname);
  std::unique_ptr<meshBase> newMesh = meshBase::CreateUnique(ofname);

  return diffMesh(refMesh.get(),newMesh.get());
}

TEST(MeshGen, defaultGen)
{
  EXPECT_EQ(0,genTest(defaultjson,"hinge.vtu", defaultRef));
}

TEST(MeshGen, uniformRefinementGen)
{
  EXPECT_EQ(0,genTest(unifjson,"hingeUnif.vtu", unifRef));
}

TEST(MeshGen, geometryRefinementGen)
{
  EXPECT_EQ(0,genTest(geomjson,"hingeGeom.vtu",geomRef));
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 7);
  defaultjson = argv[1];
  defaultRef = argv[2];
  unifjson = argv[3];
  unifRef = argv[4];
  geomjson = argv[5];
  geomRef = argv[6];
  return RUN_ALL_TESTS();
}


