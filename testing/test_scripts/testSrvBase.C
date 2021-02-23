#include <gtest.h>

#include "srvBase.H"

#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkSmartPointer.h>

#include "geoMeshBase.H"

// Tests for NEM::SRV::srvBase using minimal child class testSrvBase

class testGeoMeshBase : public NEM::MSH::geoMeshBase {
  using geoMeshBase::geoMeshBase;

 public:
  static testGeoMeshBase *New();
  vtkTypeMacro(testGeoMeshBase, geoMeshBase)

  void write(const std::string &) override {}
  void report(std::ostream &) const override {}

  void resetNative() override {}
};

vtkStandardNewMacro(testGeoMeshBase)

class testSrvBase : public NEM::SRV::srvBase {
 public:
  static testSrvBase *New();

 protected:
  testSrvBase() {
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
  }
  int RequestDataObject(vtkInformation *vtkNotUsed(request),
                        vtkInformationVector **vtkNotUsed(inputVector),
                        vtkInformationVector *outputVector) override {
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    auto *output = dynamic_cast<testGeoMeshBase *>(
        outInfo->Get(vtkDataObject::DATA_OBJECT()));

    if (!output) {
      output = testGeoMeshBase::New();
      outInfo->Set(vtkDataObject::DATA_OBJECT(), output);
      output->FastDelete();

      this->GetOutputPortInformation(0)->Set(vtkDataObject::DATA_EXTENT_TYPE(),
                                             output->GetExtentType());
    }

    return 1;
  }

  int RequestData(vtkInformation *vtkNotUsed(request),
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override {
    // get the info objects
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation *outInfo = outputVector->GetInformationObject(0);

    // get the input and output. Input may just have the geoMeshBase
    // interface, but output should be a concrete class.
    NEM::MSH::geoMeshBase *input = NEM::MSH::geoMeshBase::SafeDownCast(
        inInfo->Get(vtkDataObject::DATA_OBJECT()));
    testGeoMeshBase *output = testGeoMeshBase::SafeDownCast(
        outInfo->Get(vtkDataObject::DATA_OBJECT()));

    output->takeGeoMesh(input);

    return 1;
  }

  // see algorithm for more info
  int FillInputPortInformation(int vtkNotUsed(port),
                               vtkInformation *info) override {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "geoMeshBase");
    return 1;
  }
  int FillOutputPortInformation(int vtkNotUsed(port),
                                vtkInformation *info) override {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "testGeoMeshBase");
    return 1;
  }
};

vtkStandardNewMacro(testSrvBase)

TEST(srvBase, ExecuteVtkSmartPointer) {
  auto sb = vtkSmartPointer<testSrvBase>::New();
  auto in = vtkSmartPointer<testGeoMeshBase>::New();

  sb->SetInputDataObject(in);
  sb->Update();
  testGeoMeshBase::SafeDownCast(sb->GetOutputDataObject(0));
}

TEST(srvBase, ExecuteVtkNewVtkDelete) {
  auto *sb = testSrvBase::New();
  auto *in = testGeoMeshBase::New();

  sb->SetInputDataObject(in);
  sb->Update();
  testGeoMeshBase::SafeDownCast(sb->GetOutputDataObject(0));

  sb->Delete();
  in->Delete();
}

TEST(srvBase, ExecuteNewVtkDelete) {
  vtkNew<testSrvBase> sb;
  auto *in = new testGeoMeshBase;

  sb->SetInputDataObject(in);
  sb->Update();
  testGeoMeshBase::SafeDownCast(sb->GetOutputDataObject(0));

  in->Delete();
}

/*
// NOPE
TEST(srvBase, ExecuteConstructor) {
  testSrvBase sb{};
  testGeoMeshBase in{};

  sb.SetInputDataObject(&in);
  sb.Update();
  testGeoMeshBase *out =
      testGeoMeshBase::SafeDownCast(sb.GetOutputDataObject(0));
}

// NOPE
TEST(srvBase, ExecuteConstructorVtkDelete) {
  testSrvBase sb{};
  testGeoMeshBase in{};

  sb.SetInputDataObject(&in);
  sb.Update();
  testGeoMeshBase *out =
      testGeoMeshBase::SafeDownCast(sb.GetOutputDataObject(0));

  sb.Delete();
  in.Delete();
}

// NOPE
TEST(srvBase, ExecuteNewDelete) {
  auto *sb = new testSrvBase;
  auto *in = new testGeoMeshBase;

  sb->SetInputDataObject(in);
  sb->Update();
  testGeoMeshBase *out =
      testGeoMeshBase::SafeDownCast(sb->GetOutputDataObject(0));

  delete sb;
  delete in;
}
*/

int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);

  int res = RUN_ALL_TESTS();

  return res;
}
