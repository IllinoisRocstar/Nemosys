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

#ifndef _GMSHPlugin_Mesher_HXX_
#define _GMSHPlugin_Mesher_HXX_

#include <stdio.h>
#include "GmshVersion.h"
// #if GMSH_MAJOR_VERSION >=4
// #include "gmsh.h"
// #else
// #include "Gmsh.h"
// #endif
#include "GmshConfig.h"
#include "GModelIO_OCC.h"
#include "GModelIO_GEO.h"
#include "Geo.h"
#if GMSH_MAJOR_VERSION >=4
#include "GEdge.h" 
#include "GFace.h"
#else
#include "GEdgeCompound.h"
#include "GFaceCompound.h"
#endif
#include "MElement.h"

#include "GMSHPlugin_Defs.hxx"
#include "StdMeshers_FaceSide.hxx"
#include "SMDS_MeshElement.hxx"
#include "SMESH_Algo.hxx"

#include <map>
#include <vector>
#include <set>

class SMESH_Mesh;
class SMESH_Comment;
class SMESHDS_Mesh;
class TopoDS_Shape;
//class TopTools_DataMapOfShapeShape;
//class TopTools_IndexedMapOfShape;
class GMSHPlugin_Hypothesis;


//=============================================================================
/*!
 * \brief This class calls the GMSH mesher of OCC geometry
 */
//=============================================================================

class GMSHPLUGIN_EXPORT GMSHPlugin_Mesher 
{
 public:
  // ---------- PUBLIC METHODS ----------

  GMSHPlugin_Mesher (SMESH_Mesh* mesh, const TopoDS_Shape& aShape);

  void SetParameters(const GMSHPlugin_Hypothesis*          hyp);

  bool Compute();

  bool Evaluate(MapShapeNbElems& aResMap);
  
  static float DistBoundingBox(const SBoundingBox3d& bounds, const SPoint3& point);

 private:
  SMESH_Mesh*          _mesh;
  const TopoDS_Shape&  _shape;
  int                  _algo2d;
  int                  _algo3d;
  int                  _recomb2DAlgo;
  bool                 _recombineAll;
  int                  _subdivAlgo;
  int                  _remeshAlgo;
  int                  _remeshPara;
  double               _smouthSteps;
  double               _sizeFactor;
  double               _minSize, _maxSize;
  bool                 _secondOrder, _useIncomplElem;
  bool                 _is2d;
  GModel*              _gModel;
  
  std::set<std::string> _compounds;
  
  void SetGmshOptions();
  void CreateGmshCompounds();
  void FillSMesh();

  class mymsg : public GmshMessage
  {
  private:
    GModel* _gModel;
  public:
    mymsg(GModel* _gModel) :  GmshMessage(), _gModel(_gModel) {}
    void operator()(std::string level, std::string msg);
  };
};

#endif
