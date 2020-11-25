#include "srvBase.H"

#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkStreamingDemandDrivenPipeline.h>

#include "vtkGeoMesh.H"

/* TODO:
 * nemInformation
 * https://gitlab.kitware.com/vtk/vtk/-/blob/master/Common/ExecutionModel/vtkAlgorithm.cxx
 * RequestData
 * https://gitlab.kitware.com/vtk/vtk/-/blob/master/Filters/Geometry/vtkGeometryFilter.cxx#L131-140
 * https://vtk.org/doc/nightly/html/classvtkAlgorithm.html
 */

namespace NEM {
namespace SRV {

srvBase::srvBase() {
  std::cout << "srvBase constructed" << std::endl;
}

srvBase::~srvBase() {
  this->ReferenceCount--;  // Prevents VTK warning regarding loose reference.
  std::cout << "srvBase destructed" << std::endl;
}

int srvBase::ProcessRequest(vtkInformation *request,
                            vtkInformationVector **inputVector,
                            vtkInformationVector *outputVector) {
  // Create an output object of the correct type.
  if (request->Has(vtkDemandDrivenPipeline::REQUEST_DATA_OBJECT())) {
    return this->RequestDataObject(request, inputVector, outputVector);
  }

  // generate the data
  if (request->Has(vtkDemandDrivenPipeline::REQUEST_DATA())) {
    return this->RequestData(request, inputVector, outputVector);
  }

  if (request->Has(vtkStreamingDemandDrivenPipeline::REQUEST_UPDATE_EXTENT())) {
    return this->RequestUpdateExtent(request, inputVector, outputVector);
  }

  // execute information
  if (request->Has(vtkDemandDrivenPipeline::REQUEST_INFORMATION())) {
    return this->RequestInformation(request, inputVector, outputVector);
  }

  return this->Superclass::ProcessRequest(request, inputVector, outputVector);
}

int srvBase::RequestInformation(
    vtkInformation *vtkNotUsed(request),
    vtkInformationVector **vtkNotUsed(inputVector),
    vtkInformationVector *vtkNotUsed(outputVector)) {
  // do nothing let subclasses handle it
  return 1;
}

NEM::MSH::geoMeshBase *srvBase::GetOutput() { return this->GetOutput(0); }

NEM::MSH::geoMeshBase *srvBase::GetOutput(int port) {
  return NEM::MSH::geoMeshBase::SafeDownCast(this->GetOutputDataObject(port));
}

int srvBase::RequestDataObject(vtkInformation *vtkNotUsed(request),
                               vtkInformationVector **vtkNotUsed(inputVector),
                               vtkInformationVector *vtkNotUsed(outputVector)) {
  return 0;
}

int srvBase::RequestUpdateExtent(
    vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector,
    vtkInformationVector *vtkNotUsed(outputVector)) {
  int numInputPorts = this->GetNumberOfInputPorts();
  for (int i = 0; i < numInputPorts; i++) {
    int numInputConnections = this->GetNumberOfInputConnections(i);
    for (int j = 0; j < numInputConnections; j++) {
      vtkInformation *inputInfo = inputVector[i]->GetInformationObject(j);
      inputInfo->Set(vtkStreamingDemandDrivenPipeline::EXACT_EXTENT(), 1);
    }
  }
  return 1;
}

int srvBase::RequestData(vtkInformation *vtkNotUsed(request),
                         vtkInformationVector **vtkNotUsed(inputVector),
                         vtkInformationVector *vtkNotUsed(outputVector)) {
  return 0;
}

}  // namespace SRV
}  // namespace NEM
