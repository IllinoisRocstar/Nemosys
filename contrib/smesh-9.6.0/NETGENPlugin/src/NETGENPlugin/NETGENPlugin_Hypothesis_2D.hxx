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
// File      : NETGENPlugin_Hypothesis_2D.hxx
// Author    : Michael Sazonov (OCN)
// Date      : 27/03/2006
// Project   : SALOME
//=============================================================================
//
#ifndef _NETGENPlugin_Hypothesis_2D_HXX_
#define _NETGENPlugin_Hypothesis_2D_HXX_

#include "NETGENPlugin_Defs.hxx"

#include "NETGENPlugin_Hypothesis.hxx"
#include "Utils_SALOME_Exception.hxx"

class SMESH_Group;

//  Parameters of NETGEN.
// This class is just to give 2D dimension, actually
// it inherits all behaviour of the parent 

class NETGENPLUGIN_EXPORT  NETGENPlugin_Hypothesis_2D: public NETGENPlugin_Hypothesis
{
public:

  NETGENPlugin_Hypothesis_2D(int hypId, SMESH_Gen * gen);

// private:
//   bool _quadAllowed;
};


//  Parameters of NETGEN remesher
//

class NETGENPLUGIN_EXPORT NETGENPlugin_RemesherHypothesis_2D: public NETGENPlugin_Hypothesis
{
 public:

  NETGENPlugin_RemesherHypothesis_2D(int hypId, SMESH_Gen * gen);

  void   SetRidgeAngle( double angle );
  double GetRidgeAngle() const{ return _ridgeAngle; }

  void   SetEdgeCornerAngle( double angle );
  double GetEdgeCornerAngle() const { return _edgeCornerAngle; }

  void   SetChartAngle( double angle );
  double GetChartAngle() const { return _chartAngle; }

  void   SetOuterChartAngle( double angle );
  double GetOuterChartAngle() const { return _outerChartAngle; }

  void   SetRestHChartDistFactor( double f );
  double GetRestHChartDistFactor() const { return _restHChartDistFactor; }

  void   SetRestHChartDistEnable( bool enable );
  bool   GetRestHChartDistEnable() const { return _restHChartDistEnable; }

  void   SetRestHLineLengthFactor( double f );
  double GetRestHLineLengthFactor() const { return _restHLineLengthFactor; }

  void   SetRestHLineLengthEnable( bool enable );
  bool   GetRestHLineLengthEnable() const { return _restHLineLengthEnable; }

  void   SetRestHCloseEdgeFactor( double f );
  double GetRestHCloseEdgeFactor() const { return _restHCloseEdgeFactor; }

  void   SetRestHCloseEdgeEnable( bool enable );
  bool   GetRestHCloseEdgeEnable() const { return _restHCloseEdgeEnable; }

  void   SetRestHSurfCurvFactor( double f );
  double GetRestHSurfCurvFactor() const { return _restHSurfCurvFactor; }

  void   SetRestHSurfCurvEnable( bool enable );
  bool   GetRestHSurfCurvEnable() const { return _restHSurfCurvEnable; }

  void   SetRestHEdgeAngleFactor( double f );
  double GetRestHEdgeAngleFactor() const { return _restHEdgeAngleFactor; }

  void   SetRestHEdgeAngleEnable( bool enable );
  bool   GetRestHEdgeAngleEnable() const { return _restHEdgeAngleEnable; }

  void   SetRestHSurfMeshCurvFactor( double f );
  double GetRestHSurfMeshCurvFactor() const { return _restHSurfMeshCurvFactor; }

  void   SetRestHSurfMeshCurvEnable( bool enable );
  bool   GetRestHSurfMeshCurvEnable() const { return _restHSurfMeshCurvEnable; }

  void   SetKeepExistingEdges( bool toKeep );
  bool   GetKeepExistingEdges() const { return _keepExistingEdges; }

  void   SetMakeGroupsOfSurfaces( bool toMake );
  bool   GetMakeGroupsOfSurfaces() const { return _makeGroupsOfSurfaces; }

  void   SetFixedEdgeGroup( const SMESH_Group* edgeGroup );
  int    GetFixedEdgeGroupID() const { return _fixedEdgeGroupID; }
  SMESH_Group* GetFixedEdgeGroup( const SMESH_Mesh& mesh ) const;

  void   SetLoadMeshOnCancel( bool toLoad );
  bool   GetLoadMeshOnCancel() const { return _loadOnCancel; }

  static double DefaultRidgeAngle()              { return 30.; }
  static double DefaultEdgeCornerAngle()         { return 60.; }
  static double DefaultChartAngle()              { return 15.; }
  static double DefaultOuterChartAngle()         { return 70.; }
  static double DefaultRestHChartDistFactor()    { return 1.2; }
  static bool   DefaultRestHChartDistEnable()    { return true; }
  static double DefaultRestHLineLengthFactor()   { return 0.5; }
  static bool   DefaultRestHLineLengthEnable()   { return true; }
  static double DefaultRestHCloseEdgeFactor()    { return 1.; }
  static bool   DefaultRestHCloseEdgeEnable()    { return true; }
  static double DefaultRestHSurfCurvFactor()     { return 1.; }
  static bool   DefaultRestHSurfCurvEnable()     { return false; }
  static double DefaultRestHEdgeAngleFactor()    { return 1.; }
  static bool   DefaultRestHEdgeAngleEnable()    { return false; }
  static double DefaultRestHSurfMeshCurvFactor() { return 1.; }
  static bool   DefaultRestHSurfMeshCurvEnable() { return false; }
  static bool   DefaultKeepExistingEdges()       { return false; }
  static bool   DefaultMakeGroupsOfSurfaces()    { return false; }

  virtual std::ostream & SaveTo(std::ostream & save);
  virtual std::istream & LoadFrom(std::istream & load);

 private:

  // STL charts
  double _ridgeAngle; // yellow edges angle (in degrees)
  double _edgeCornerAngle;
  double _chartAngle;
  double _outerChartAngle;

  // Mesh size: restrict h due to ...
  double _restHChartDistFactor;    // chart distance
  bool   _restHChartDistEnable;
  double _restHLineLengthFactor;   // line length
  bool   _restHLineLengthEnable;
  double _restHCloseEdgeFactor;    // close edges
  bool   _restHCloseEdgeEnable;
  double _restHSurfCurvFactor;     // surface curvature
  bool   _restHSurfCurvEnable;
  double _restHEdgeAngleFactor;    // edge angle
  bool   _restHEdgeAngleEnable;
  double _restHSurfMeshCurvFactor; // surface mesh curv
  bool   _restHSurfMeshCurvEnable;

  // SALOME features
  bool   _keepExistingEdges;
  bool   _makeGroupsOfSurfaces;
  int    _fixedEdgeGroupID;
  bool   _loadOnCancel;

};

#endif
