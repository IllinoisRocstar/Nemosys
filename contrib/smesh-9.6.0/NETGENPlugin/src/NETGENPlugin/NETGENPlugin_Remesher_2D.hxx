// Copyright (C) 2007-2021  CEA/DEN, EDF R&D, OPEN CASCADE
//
// Copyright (C) 2003-2007  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
// CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
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
// See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//
// File      : NETGENPlugin_Remesher_2D.hxx
// Created   : Thu Sep 21 16:48:46 2017
// Author    : Edward AGAPOV (eap)


#ifndef __NETGENPlugin_Remesher_2D_HXX__
#define __NETGENPlugin_Remesher_2D_HXX__

#include "NETGENPlugin_Defs.hxx"

#include "SMESH_Algo.hxx"
#include "SMESH_Mesh.hxx"

class NETGENPLUGIN_EXPORT NETGENPlugin_Remesher_2D: public SMESH_2D_Algo
{
 public:
  NETGENPlugin_Remesher_2D(int hypId, SMESH_Gen* gen);

  virtual bool CheckHypothesis(SMESH_Mesh&                          theMesh,
                               const TopoDS_Shape&                  theShape,
                               SMESH_Hypothesis::Hypothesis_Status& theStatus);

  virtual bool Compute(SMESH_Mesh & theMesh, SMESH_MesherHelper* theHelper);

  virtual bool Compute(SMESH_Mesh& theMesh, const TopoDS_Shape& theShape);

  virtual void CancelCompute();

  virtual double GetProgress() const;


  virtual bool Evaluate(SMESH_Mesh&         theMesh,
                        const TopoDS_Shape& theShape,
                        MapShapeNbElems&    theResMap);

 protected:

  const SMESHDS_Hypothesis* _hypothesis;
};

#endif
