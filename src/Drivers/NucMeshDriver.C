#include "Drivers/NucMeshDriver.H"

#include "Mesh/geoMeshFactory.H"
#include "Mesh/smeshGeoMesh.H"
#include "Services/NucMeshSrv.H"

namespace NEM {
namespace DRV {

NucMeshDriver::NucMeshDriver(Files file, Opts opts)
    : file_(std::move(file)), opts_(std::move(opts)) {}

NucMeshDriver::NucMeshDriver() : NucMeshDriver(Files{std::string{}}, Opts{}) {}

const NucMeshDriver::Files &NucMeshDriver::getFiles() const { return file_; }

void NucMeshDriver::setFiles(Files files) { file_ = std::move(files); }

const NucMeshDriver::Opts &NucMeshDriver::getOpts() const { return opts_; }

void NucMeshDriver::setOpts(Opts opts) { opts_ = std::move(opts); }

jsoncons::string_view NucMeshDriver::getProgramType() const {
  return programType;
}

vtkSmartPointer<NEM::MSH::geoMeshBase> NucMeshDriver::draw() const {
  vtkNew<NEM::SRV::NucMeshSrv> nucMeshRunner{};
  nucMeshRunner->SetConfiguration(opts_);
  nucMeshRunner->Update();
  return NEM::MSH::smeshGeoMesh::SafeDownCast(nucMeshRunner->GetOutput());
}

void NucMeshDriver::execute() const {
  auto outNative = this->draw();
  auto outType = vtkSmartPointer<NEM::MSH::geoMeshBase>::Take(
      NEM::MSH::New(file_.outputFile));
  outType->takeGeoMesh(outNative);
  outType->write(file_.outputFile);
}

}  // namespace DRV
}  // namespace NEM