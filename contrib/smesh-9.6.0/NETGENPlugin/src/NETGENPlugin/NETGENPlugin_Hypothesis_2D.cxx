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

//  NETGENPlugin : C++ implementation
// File      : NETGENPlugin_Hypothesis_2D.cxx
// Author    : Michael Sazonov (OCN)
// Date      : 28/03/2006
// Project   : SALOME
//=============================================================================
//
#include "NETGENPlugin_Hypothesis_2D.hxx"

#include <SMESH_Mesh.hxx>
#include <SMESH_Group.hxx>
#include <SMESHDS_GroupBase.hxx>

using namespace std;

//=============================================================================
/*!
 *  
 */
//=============================================================================
NETGENPlugin_Hypothesis_2D::NETGENPlugin_Hypothesis_2D (int hypId,
                                                        SMESH_Gen * gen)
  : NETGENPlugin_Hypothesis(hypId, gen)/*,
    _quadAllowed (GetDefaultQuadAllowed())*/
{
  _name = "NETGEN_Parameters_2D";
  _param_algo_dim = 2;
}

//=============================================================================
/*!
 *
 */
//=============================================================================
NETGENPlugin_RemesherHypothesis_2D::
NETGENPlugin_RemesherHypothesis_2D (int hypId, SMESH_Gen * gen)
  : NETGENPlugin_Hypothesis(hypId, gen),
    _ridgeAngle             ( DefaultRidgeAngle()             ),
    _edgeCornerAngle        ( DefaultEdgeCornerAngle()        ),
    _chartAngle             ( DefaultChartAngle()             ),
    _outerChartAngle        ( DefaultOuterChartAngle()        ),
    _restHChartDistFactor   ( DefaultRestHChartDistFactor()   ),
    _restHChartDistEnable   ( DefaultRestHChartDistEnable()   ),
    _restHLineLengthFactor  ( DefaultRestHLineLengthFactor()  ),
    _restHLineLengthEnable  ( DefaultRestHLineLengthEnable()  ),
    _restHCloseEdgeFactor   ( DefaultRestHCloseEdgeFactor()   ),
    _restHCloseEdgeEnable   ( DefaultRestHCloseEdgeEnable()   ),
    _restHSurfCurvFactor    ( DefaultRestHSurfCurvFactor()    ),
    _restHSurfCurvEnable    ( DefaultRestHSurfCurvEnable()    ),
    _restHEdgeAngleFactor   ( DefaultRestHEdgeAngleFactor()   ),
    _restHEdgeAngleEnable   ( DefaultRestHEdgeAngleEnable()   ),
    _restHSurfMeshCurvFactor( DefaultRestHSurfMeshCurvFactor()),
    _restHSurfMeshCurvEnable( DefaultRestHSurfMeshCurvEnable()),
    _keepExistingEdges      ( DefaultKeepExistingEdges()      ),
    _makeGroupsOfSurfaces   ( DefaultMakeGroupsOfSurfaces()   ),
    _fixedEdgeGroupID       ( -1                              ),
    _loadOnCancel           ( false                           )
{
  _name = "NETGEN_RemesherParameters_2D";
  _param_algo_dim = 2;
}

//=============================================================================
/*!
 *
 */
//=============================================================================

void NETGENPlugin_RemesherHypothesis_2D::SetRidgeAngle( double angle )
{
  if ( _ridgeAngle != angle )
  {
    _ridgeAngle = angle;
    NotifySubMeshesHypothesisModification();
  }
}

//=======================================================================
//function : SetEdgeCornerAngle
//purpose  : 
//=======================================================================

void NETGENPlugin_RemesherHypothesis_2D::SetEdgeCornerAngle( double angle )
{
  if ( _edgeCornerAngle != angle )
  {
    _edgeCornerAngle = angle;
    NotifySubMeshesHypothesisModification();
  }
}

//=======================================================================
//function : SetChartAngle
//purpose  : 
//=======================================================================

void NETGENPlugin_RemesherHypothesis_2D::SetChartAngle( double angle )
{
  if ( _chartAngle != angle )
  {
    _chartAngle = angle;
    NotifySubMeshesHypothesisModification();
  }
}

//=======================================================================
//function : SetOuterChartAngle
//purpose  : 
//=======================================================================

void NETGENPlugin_RemesherHypothesis_2D::SetOuterChartAngle( double angle )
{
  if ( _outerChartAngle != angle )
  {
    _outerChartAngle = angle;
    NotifySubMeshesHypothesisModification();
  }
}

//=======================================================================
//function : SetRestHChartDistFactor
//purpose  : 
//=======================================================================

void NETGENPlugin_RemesherHypothesis_2D::SetRestHChartDistFactor( double f )
{
  if ( _restHChartDistFactor != f )
  {
    _restHChartDistFactor = f;
    NotifySubMeshesHypothesisModification();
  }
}

//=======================================================================
//function : SetRestHChartDistEnable
//purpose  : 
//=======================================================================

void NETGENPlugin_RemesherHypothesis_2D::SetRestHChartDistEnable( bool enable )
{
  if ( _restHChartDistEnable != enable )
  {
    _restHChartDistEnable = enable;
    NotifySubMeshesHypothesisModification();
  }
}

//=======================================================================
//function : SetRestHLineLengthFactor
//purpose  : 
//=======================================================================

void NETGENPlugin_RemesherHypothesis_2D::SetRestHLineLengthFactor( double f )
{
  if ( _restHLineLengthFactor != f )
  {
    _restHLineLengthFactor = f;
    NotifySubMeshesHypothesisModification();
  }
}

//=======================================================================
//function : SetRestHLineLengthEnable
//purpose  : 
//=======================================================================

void NETGENPlugin_RemesherHypothesis_2D::SetRestHLineLengthEnable( bool enable )
{
  if ( _restHLineLengthEnable != enable )
  {
    _restHLineLengthEnable = enable;
    NotifySubMeshesHypothesisModification();
  }
}

//=======================================================================
//function : SetRestHCloseEdgeFactor
//purpose  : 
//=======================================================================

void NETGENPlugin_RemesherHypothesis_2D::SetRestHCloseEdgeFactor( double f )
{
  if ( _restHCloseEdgeFactor != f )
  {
    _restHCloseEdgeFactor = f;
    NotifySubMeshesHypothesisModification();
  }
}

//=======================================================================
//function : SetRestHCloseEdgeEnable
//purpose  : 
//=======================================================================

void NETGENPlugin_RemesherHypothesis_2D::SetRestHCloseEdgeEnable( bool enable )
{
  if ( _restHCloseEdgeEnable != enable )
  {
    _restHCloseEdgeEnable = enable;
    NotifySubMeshesHypothesisModification();
  }
}

//=======================================================================
//function : SetRestHSurfCurvFactor
//purpose  : 
//=======================================================================

void NETGENPlugin_RemesherHypothesis_2D::SetRestHSurfCurvFactor( double f )
{
  if ( _restHSurfCurvFactor != f )
  {
    _restHSurfCurvFactor = f;
    NotifySubMeshesHypothesisModification();
  }
}

//=======================================================================
//function : SetRestHSurfCurvEnable
//purpose  : 
//=======================================================================

void NETGENPlugin_RemesherHypothesis_2D::SetRestHSurfCurvEnable( bool enable )
{
  if ( _restHSurfCurvEnable != enable )
  {
    _restHSurfCurvEnable = enable;
    NotifySubMeshesHypothesisModification();
  }
}

//=======================================================================
//function : SetRestHEdgeAngleFactor
//purpose  : 
//=======================================================================

void NETGENPlugin_RemesherHypothesis_2D::SetRestHEdgeAngleFactor( double f )
{
  if ( _restHEdgeAngleFactor != f )
  {
    _restHEdgeAngleFactor = f;
    NotifySubMeshesHypothesisModification();
  }
}

//=======================================================================
//function : SetRestHEdgeAngleEnable
//purpose  : 
//=======================================================================

void NETGENPlugin_RemesherHypothesis_2D::SetRestHEdgeAngleEnable( bool enable )
{
  if ( _restHEdgeAngleEnable != enable )
  {
    _restHEdgeAngleEnable = enable;
    NotifySubMeshesHypothesisModification();
  }
}

//=======================================================================
//function : SetRestHSurfMeshCurvFactor
//purpose  :
//=======================================================================

void NETGENPlugin_RemesherHypothesis_2D::SetRestHSurfMeshCurvFactor( double f )
{
  if ( _restHSurfMeshCurvFactor != f )
  {
    _restHSurfMeshCurvFactor = f;
    NotifySubMeshesHypothesisModification();
  }
}

//=======================================================================
//function : SetRestHSurfMeshCurvEnable
//purpose  :
//=======================================================================

void NETGENPlugin_RemesherHypothesis_2D::SetRestHSurfMeshCurvEnable( bool enable )
{
  if ( _restHSurfMeshCurvEnable != enable )
  {
    _restHSurfMeshCurvEnable = enable;
    NotifySubMeshesHypothesisModification();
  }
}

//=======================================================================
//function : SetKeepExistingEdges
//purpose  :
//=======================================================================

void NETGENPlugin_RemesherHypothesis_2D::SetKeepExistingEdges( bool toKeep )
{
  if ( _keepExistingEdges != toKeep )
  {
    _keepExistingEdges = toKeep;
    NotifySubMeshesHypothesisModification();
  }
}

//=======================================================================
//function : SetMakeGroupsOfSurfaces
//purpose  :
//=======================================================================

void NETGENPlugin_RemesherHypothesis_2D::SetMakeGroupsOfSurfaces( bool toMake )
{
  if ( _makeGroupsOfSurfaces != toMake )
  {
    _makeGroupsOfSurfaces = toMake;
    NotifySubMeshesHypothesisModification();
  }
}

//=======================================================================
//function : SetFixedEdgeGroup
//purpose  : Set a group of edges whose nodes must not be moved
//=======================================================================

void NETGENPlugin_RemesherHypothesis_2D::SetFixedEdgeGroup( const SMESH_Group* edgeGroup )
{
  int id = edgeGroup ? edgeGroup->GetID() : -1;
  if ( id != _fixedEdgeGroupID )
  {
    _fixedEdgeGroupID = id;
    NotifySubMeshesHypothesisModification();
  }
}

//=======================================================================
//function : SetLoadMeshOnCancel
//purpose  : allow getting a current mesh existing upon CancelCompute()
//=======================================================================

void NETGENPlugin_RemesherHypothesis_2D::SetLoadMeshOnCancel( bool toLoad )
{
  if ( toLoad != _loadOnCancel )
  {
    _loadOnCancel = toLoad;
    NotifySubMeshesHypothesisModification();
  }
}

//=======================================================================
//function : GetFixedEdgeGroup
//purpose  : Return a group of edges whose nodes must not be moved
//=======================================================================

SMESH_Group*
NETGENPlugin_RemesherHypothesis_2D::GetFixedEdgeGroup( const SMESH_Mesh& mesh ) const
{
  SMESH_Group* group = mesh.GetGroup( _fixedEdgeGroupID );
  if ( group && group->GetGroupDS()->GetType() != SMDSAbs_Edge )
    group = NULL;

  return group;
}

//=============================================================================
/*!
 *
 */
//=============================================================================

std::ostream & NETGENPlugin_RemesherHypothesis_2D::SaveTo(std::ostream & save)
{
  NETGENPlugin_Hypothesis::SaveTo( save );
  save << " " << _ridgeAngle;

  save << " " << _edgeCornerAngle        ;
  save << " " << _chartAngle             ;
  save << " " << _outerChartAngle        ;
  save << " " << _restHChartDistFactor   ;
  save << " " << _restHChartDistEnable   ;
  save << " " << _restHLineLengthFactor  ;
  save << " " << _restHLineLengthEnable  ;
  save << " " << _restHCloseEdgeFactor   ;
  save << " " << _restHCloseEdgeEnable   ;
  save << " " << _restHSurfCurvFactor    ;
  save << " " << _restHSurfCurvEnable    ;
  save << " " << _restHEdgeAngleFactor   ;
  save << " " << _restHEdgeAngleEnable   ;
  save << " " << _restHSurfMeshCurvFactor;
  save << " " << _restHSurfMeshCurvEnable;
  save << " " << _keepExistingEdges      ;
  save << " " << _makeGroupsOfSurfaces   ;
  save << " " << _fixedEdgeGroupID       ;
  save << " " << _loadOnCancel           ;

  return save;
}

//=============================================================================
/*!
 *
 */
//=============================================================================

std::istream & NETGENPlugin_RemesherHypothesis_2D::LoadFrom(std::istream & load)
{
  NETGENPlugin_Hypothesis::LoadFrom( load );
  if ( !load )
    load.clear(ios::badbit | load.rdstate());

  load >> _ridgeAngle;

  if ( !load )
    _ridgeAngle = DefaultRidgeAngle();

  load >> _edgeCornerAngle;
  if ( !load )
    _edgeCornerAngle = DefaultEdgeCornerAngle();

  load >> _chartAngle;
  if ( !load )
    _chartAngle = DefaultChartAngle();

  load >> _outerChartAngle;
  if ( !load )
    _outerChartAngle = DefaultOuterChartAngle();

  load >> _restHChartDistFactor;
  if ( !load )
    _restHChartDistFactor = DefaultRestHChartDistFactor();

  load >> _restHChartDistEnable;
  if ( !load )
    _restHChartDistEnable = DefaultRestHChartDistEnable();

  load >> _restHLineLengthFactor;
  if ( !load )
    _restHLineLengthFactor = DefaultRestHLineLengthFactor();

  load >> _restHLineLengthEnable;
  if ( !load )
    _restHLineLengthEnable = DefaultRestHLineLengthEnable();

  load >> _restHCloseEdgeFactor;
  if ( !load )
    _restHCloseEdgeFactor = DefaultRestHCloseEdgeFactor();

  load >> _restHCloseEdgeEnable;
  if ( !load )
    _restHCloseEdgeEnable = DefaultRestHCloseEdgeEnable();

  load >> _restHSurfCurvFactor;
  if ( !load )
    _restHSurfCurvFactor = DefaultRestHSurfCurvFactor();

  load >> _restHSurfCurvEnable;
  if ( !load )
    _restHSurfCurvEnable = DefaultRestHSurfCurvEnable();

  load >> _restHEdgeAngleFactor;
  if ( !load )
    _restHEdgeAngleFactor = DefaultRestHEdgeAngleFactor();

  load >> _restHEdgeAngleEnable;
  if ( !load )
    _restHEdgeAngleEnable = DefaultRestHEdgeAngleEnable();

  load >> _restHSurfMeshCurvFactor;
  if ( !load )
    _restHSurfMeshCurvFactor = DefaultRestHSurfMeshCurvFactor();

  load >> _restHSurfMeshCurvEnable;
  if ( !load )
    _restHSurfMeshCurvEnable = DefaultRestHSurfMeshCurvEnable();

  load >> _keepExistingEdges;
  if ( !load )
    _keepExistingEdges = DefaultKeepExistingEdges();

  load >> _makeGroupsOfSurfaces;
  if ( !load )
    _makeGroupsOfSurfaces = DefaultMakeGroupsOfSurfaces();

  load >> _fixedEdgeGroupID;
  if ( !load )
    _fixedEdgeGroupID = -1;

  load >> _loadOnCancel;
  if ( !load )
    _loadOnCancel = false;

  return load;
}
