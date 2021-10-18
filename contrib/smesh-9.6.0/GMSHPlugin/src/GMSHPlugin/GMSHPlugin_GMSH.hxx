// Copyright (C) 2012-2015  ALNEOS
// Copyright (C) 2016-2020  EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
// See http://www.alneos.com/ or email : contact@alneos.fr
// See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//
#ifndef _GMSHPlugin_GMSH_HXX_
#define _GMSHPlugin_GMSH_HXX_

#include "GMSHPlugin_Defs.hxx"

#include "SMESH_Algo.hxx"
#include "SMESH_subMesh.hxx"
#include "SMESH_Mesh.hxx"
#include "StdMeshers_MaxElementVolume.hxx"
#include "Utils_SALOME_Exception.hxx"

//class GMSHPlugin_Hypothesis;

class GMSHPLUGIN_EXPORT GMSHPlugin_GMSH: public SMESH_3D_Algo
{
public:
  GMSHPlugin_GMSH(int hypId, SMESH_Gen* gen);
  virtual ~GMSHPlugin_GMSH();

  virtual bool CheckHypothesis(SMESH_Mesh& aMesh,
                               const TopoDS_Shape& aShape,
                               SMESH_Hypothesis::Hypothesis_Status& aStatus);

  virtual bool Compute(SMESH_Mesh& aMesh,
                       const TopoDS_Shape& aShape);

#ifdef WITH_SMESH_CANCEL_COMPUTE
  virtual void CancelCompute();
#endif

  virtual bool Evaluate(SMESH_Mesh& aMesh,
                        const TopoDS_Shape& aShape,
                        MapShapeNbElems& aResMap);

protected:
  const SMESHDS_Hypothesis* _hypothesis;
};

#endif
