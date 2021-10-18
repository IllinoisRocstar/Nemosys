#include "Services/srvBase.H"

#include <vtkInformationExecutivePortKey.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkStreamingDemandDrivenPipeline.h>

#include "Mesh/vtkGeoMesh.H"
#include "Mesh/oshGeoMesh.H"
#include "Mesh/exoGeoMesh.H"
#include "Mesh/inpGeoMesh.H"

#ifdef HAVE_CFMSH
#  include "Mesh/foamGeoMesh.H"
#endif

#ifdef HAVE_GMSH
#  include "Mesh/gmshGeoMesh.H"
#endif

#ifdef HAVE_OCC
#  include "Mesh/smeshGeoMesh.H"
#endif

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

int srvBase::RequestInformation(vtkInformation *request,
                                vtkInformationVector **inputVector,
                                vtkInformationVector *outputVector) {
  // do nothing let subclasses handle it
  return 1;
}

NEM::MSH::geoMeshBase *srvBase::GetOutput() { return this->GetOutput(0); }

NEM::MSH::geoMeshBase *srvBase::GetOutput(int port) {
  return NEM::MSH::geoMeshBase::SafeDownCast(this->GetOutputDataObject(port));
}

int srvBase::RequestDataObject(vtkInformation *request,
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {
  for (int i = 0; i < outputVector->GetNumberOfInformationObjects(); ++i) {
    auto outInfo = outputVector->GetInformationObject(i);
    auto outAlgInfo = this->GetOutputPortInformation(i);
    if (!outAlgInfo->Has(vtkDataObject::DATA_TYPE_NAME())) {
      return 0;
    }
    std::string typeName = outAlgInfo->Get(vtkDataObject::DATA_TYPE_NAME());
    MSH::geoMeshBase *output = nullptr;
    if (typeName == "vtkGeoMesh") {
      output = MSH::vtkGeoMesh::New();
    } else if (typeName == "gmshGeoMesh") {
#ifdef HAVE_GMSH
      output = MSH::gmshGeoMesh::New();
#endif
    } else if (typeName == "oshGeoMesh") {
      output = MSH::oshGeoMesh::New();
    } else if (typeName == "exoGeoMesh") {
      output = MSH::exoGeoMesh::New();
    } else if (typeName == "inpGeoMesh") {
      output = MSH::inpGeoMesh::New();
    } else if (typeName == "foamGeoMesh") {
#ifdef HAVE_CFMSH
      output = MSH::foamGeoMesh::New();
#endif
    } else if (typeName == "smeshGeoMesh") {
#ifdef HAVE_OCC
      output = MSH::smeshGeoMesh::New();
#endif
    } else if (typeName == "geoMeshBase") {
      if (i < this->GetNumberOfInputPorts() &&
          inputVector[i]->GetNumberOfInformationObjects() > 0) {
        auto inObj = inputVector[i]
                ->GetInformationObject(0)
                ->Get(vtkDataObject::DATA_OBJECT());
        if (inObj) {
          output = MSH::geoMeshBase::SafeDownCast(inObj->NewInstance());
        }
      }
    }
    if (!output) {
      return 0;
    }
    outInfo->Set(vtkDataObject::DATA_OBJECT(), output);
    output->FastDelete();
    this->GetOutputPortInformation(i)->Set(
        vtkDataObject::DATA_EXTENT_TYPE(), output->GetExtentType());
  }
  return 1;
}

int srvBase::RequestUpdateExtent(vtkInformation *request,
                                 vtkInformationVector **inputVector,
                                 vtkInformationVector *outputVector) {
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

}  // namespace SRV
}  // namespace NEM
