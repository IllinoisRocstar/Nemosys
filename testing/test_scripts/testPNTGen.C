#include <NemDriver.H>
#include <gtest.h>

const char* bench1_json;
const char* bench1_ref;
const char* bench5_json;
const char* bench5_ref;
const char* bench6_json;
const char* bench6_ref;

int genTest(const char* jsonF, const char* ofname, const char* refname)
{
  std::string fname(jsonF);
  std::ifstream inputStream(fname);
  if(!inputStream.good())
  {
    std::cerr << "Error opening file " << jsonF << std::endl;
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



TEST(PNTGen, Bench1)
{
  EXPECT_EQ(0,genTest(bench1_json,"bench1_conv.pntmesh", bench1_ref));
}

TEST(PNTGen, Bench5)
{
  EXPECT_EQ(0,genTest(bench5_json,"bench5_conv.pntmesh", bench5_ref));
}

TEST(PNTGen, Bench6)
{
  EXPECT_EQ(0,genTest(bench6_json,"bench6_conv.pntmesh", bench6_ref));
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 7);
  bench1_json = argv[1];
  bench1_ref = argv[2];
  bench5_json = argv[3];
  bench5_ref = argv[4];
  bench6_json = argv[5];
  bench6_ref = argv[6];
  return RUN_ALL_TESTS();
}


