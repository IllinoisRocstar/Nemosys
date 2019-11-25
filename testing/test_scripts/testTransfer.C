#include <gtest.h>

#include "TransferDriver.H"
#include "TransferBase.H"
#include "meshBase.H"

const char* pntSource;
const char* cellSource;
const char* targetF;
const char* pntRef;
const char* cellRef;

class TransferTest : public ::testing::Test 
{
  protected:
    // You can remove any or all of the following functions if its body
    // is empty.
  
    TransferTest()
    {
      target = meshBase::CreateShared(targetF);
    }
  
    virtual ~TransferTest()
    {}
  
    // Objects declared here can be used by all tests in the test case for orthoPoly.
    std::shared_ptr<meshBase> target;
};

TEST_F(TransferTest, pntDataTransfer)
{
  std::shared_ptr<meshBase> source = meshBase::CreateShared(pntSource);
  std::string method("Consistent Interpolation");
  // source.get()->transfer(target.get(),method);
  auto transfer = TransferDriver::CreateTransferObject(source.get(), target.get(), method);
  transfer->run();
  std::shared_ptr<meshBase> ref = meshBase::CreateShared(pntRef);
  target.get()->write("test.vtu");
  EXPECT_EQ(0,diffMesh(target.get(),ref.get()));
} 

TEST_F(TransferTest, cellDataTransfer)
{
  std::shared_ptr<meshBase> source = meshBase::CreateShared(cellSource);
  std::string method("Consistent Interpolation");
  // source.get()->transfer(target.get(),method);
  auto transfer = TransferDriver::CreateTransferObject(source.get(), target.get(), method);
  transfer->run();
  std::shared_ptr<meshBase> ref = meshBase::CreateShared(cellRef);
  EXPECT_EQ(0,diffMesh(target.get(),ref.get()));
} 

int main(int argc, char** argv) 
{
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 6);
  pntSource = argv[1];
  cellSource = argv[2];
  targetF = argv[3];
  pntRef = argv[4];
  cellRef = argv[5];
  return RUN_ALL_TESTS();
}

