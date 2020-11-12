#ifndef NEMOSYS_SRVBASE_H
#define NEMOSYS_SRVBASE_H

#include "nemosys_export.h"

#include <vtkAlgorithm.h>

#include "geoMeshBase.H"

namespace NEM {
namespace SRV {

/**
 * @class srvBase
 * @brief abstract class for services acting on @c geoMeshBase
 *
 * srvBase is an abstract class specifying the interface for all services acting
 * on a @c geoMeshBase mesh and geometry object. It derives from the
 * vtkAlgorithm, exposing the vtk robust algorithm pipeline for NEMoSys
 * services.
 */
class NEMOSYS_EXPORT srvBase : public vtkAlgorithm {
 public:
  vtkAbstractTypeMacro(srvBase, vtkAlgorithm)

  srvBase();
  ~srvBase() override;

 public:
  int ProcessRequest(vtkInformation *request, vtkInformationVector **inInfo,
                     vtkInformationVector *outInfo) override;

  //@{
  /**
   * Get the output GeoMeshBase for a port on this algorithm.
   */
  NEM::MSH::geoMeshBase *GetOutput();
  NEM::MSH::geoMeshBase *GetOutput(int);
  //@}

 protected:
  // convenience method
  virtual int RequestInformation(vtkInformation *request,
                                 vtkInformationVector **inputVector,
                                 vtkInformationVector *outputVector);

  /**
   * This is called by the superclass.
   * This is the method you should override.
   */
  virtual int RequestDataObject(vtkInformation *request,
                                vtkInformationVector **inputVector,
                                vtkInformationVector *outputVector);

  /**
   * This is called by the superclass.
   * This is the method you should override.
   */
  virtual int RequestData(vtkInformation *request,
                          vtkInformationVector **inputVector,
                          vtkInformationVector *outputVector);

  /**
   * This is called by the superclass.
   * This is the method you should override.
   */
  virtual int RequestUpdateExtent(vtkInformation *request,
                                  vtkInformationVector **inputVector,
                                  vtkInformationVector *outputVector);
};

}  // namespace SRV
}  // namespace NEM

#endif  // NEMOSYS_SRVBASE_H
