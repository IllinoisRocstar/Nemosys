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
#ifndef NEMOSYS_SRVBASE_H_
#define NEMOSYS_SRVBASE_H_

#include "nemosys_export.h"

#include <vtkAlgorithm.h>

namespace NEM {

namespace MSH {
class geoMeshBase;
}

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
  srvBase();
  ~srvBase() override;

  // convenience method
  virtual int RequestInformation(vtkInformation *request,
                                 vtkInformationVector **inputVector,
                                 vtkInformationVector *outputVector);

  /**
   * This is called by the superclass.
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
                          vtkInformationVector *outputVector) = 0;

  /**
   * This is called by the superclass.
   */
  virtual int RequestUpdateExtent(vtkInformation *request,
                                  vtkInformationVector **inputVector,
                                  vtkInformationVector *outputVector);
};

}  // namespace SRV
}  // namespace NEM

#endif  // NEMOSYS_SRVBASE_H_
