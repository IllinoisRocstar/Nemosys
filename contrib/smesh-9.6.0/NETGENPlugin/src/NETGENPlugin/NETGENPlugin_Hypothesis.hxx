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
// File      : NETGENPlugin_Hypothesis.hxx
// Author    : Michael Sazonov (OCN)
// Date      : 27/03/2006
// Project   : SALOME
//
#ifndef _NETGENPlugin_Hypothesis_HXX_
#define _NETGENPlugin_Hypothesis_HXX_

#include "NETGENPlugin_Defs.hxx"

#include "SMESH_Hypothesis.hxx"
#include "Utils_SALOME_Exception.hxx"

#include <map>

//  Parameters for work of NETGEN
//

class NETGENPLUGIN_EXPORT NETGENPlugin_Hypothesis: public SMESH_Hypothesis
{
public:

  NETGENPlugin_Hypothesis(int hypId, SMESH_Gen * gen);

  void   SetMaxSize(double theSize);
  double GetMaxSize() const { return _maxSize; }

  void   SetMinSize(double theSize);
  double GetMinSize() const { return _minSize; }

  void   SetSecondOrder(bool theVal);
  bool   GetSecondOrder() const { return _secondOrder; }

  void   SetOptimize(bool theVal);
  bool   GetOptimize() const { return _optimize; }

  enum Fineness
  {
    VeryCoarse,
    Coarse,
    Moderate,
    Fine,
    VeryFine,
    UserDefined
  };

  void   SetFineness(Fineness theFineness);
  Fineness GetFineness() const { return _fineness; }

  // the following 3 parameters are controlled by Fineness

  void   SetGrowthRate(double theRate);
  double GetGrowthRate() const { return _growthRate; }

  void   SetNbSegPerEdge(double theVal);
  double GetNbSegPerEdge() const { return _nbSegPerEdge; }

  void   SetNbSegPerRadius(double theVal);
  double GetNbSegPerRadius() const { return _nbSegPerRadius; }

  void   SetChordalErrorEnabled(bool value);
  double GetChordalErrorEnabled() const { return _chordalErrorEnabled; }
  void   SetChordalError(double value);
  double GetChordalError() const { return _chordalError; }

  typedef std::map<std::string, double> TLocalSize;
  void   SetLocalSizeOnEntry(const std::string& entry, double localSize);
  double GetLocalSizeOnEntry(const std::string& entry);
  const TLocalSize& GetLocalSizesAndEntries() const { return _localSize; }
  void   UnsetLocalSizeOnEntry(const std::string& entry);

  void   SetMeshSizeFile(const std::string& fileName);
  const std::string& GetMeshSizeFile() const { return _meshSizeFile; }

  void   SetQuadAllowed(bool theVal);
  bool   GetQuadAllowed() const { return _quadAllowed; }

  void   SetSurfaceCurvature(bool theVal);
  bool   GetSurfaceCurvature() const { return _surfaceCurvature; }

  void   SetFuseEdges(bool theVal);
  bool   GetFuseEdges() const { return _fuseEdges; }

  void   SetNbSurfOptSteps( int nb );
  int    GetNbSurfOptSteps() const { return _nbSurfOptSteps; }

  void   SetNbVolOptSteps( int nb );
  int    GetNbVolOptSteps() const { return _nbVolOptSteps; }

  void   SetElemSizeWeight( double size );
  double GetElemSizeWeight() const { return _elemSizeWeight; }

  void   SetWorstElemMeasure( int val );
  int    GetWorstElemMeasure() const { return _worstElemMeasure; }

  void   SetUseDelauney( bool toUse);
  bool   GetUseDelauney() const { return _useDelauney; }

  void   SetCheckOverlapping( bool toCheck );
  bool   GetCheckOverlapping() const { return _checkOverlapping; }

  void   SetCheckChartBoundary( bool toCheck );
  bool   GetCheckChartBoundary() const { return _checkChartBoundary; }

  // the default values (taken from NETGEN 4.5 sources)

  static Fineness GetDefaultFineness()          { return Moderate; }
  static bool     GetDefaultSecondOrder()       { return false; }
  static bool     GetDefaultQuadAllowed()       { return false; }
  static double   GetDefaultMaxSize()           { return 1000; }
  static double   GetDefaultGrowthRate()        { return 0.3; }
  static double   GetDefaultNbSegPerRadius()    { return 2; }
  static double   GetDefaultNbSegPerEdge()      { return 1; }
  static double   GetDefaultChordalError()      { return -1; } // disabled by default
  static bool     GetDefaultOptimize()          { return true; }
  static int      GetDefaultNbSurfOptSteps()    { return 3; }
  static int      GetDefaultNbVolOptSteps()     { return 3; }
  static double   GetDefaultElemSizeWeight()    { return 0.2; }
  static int      GetDefaultWorstElemMeasure()  { return 2; }
  static bool     GetDefaultSurfaceCurvature()  { return true; }
  static bool     GetDefaultUseDelauney()       { return true; }
  static bool     GetDefaultCheckOverlapping()  { return true; }
  static bool     GetDefaultCheckChartBoundary(){ return true; }
  static bool     GetDefaultFuseEdges()         { return true; }

  // Persistence
  virtual std::ostream & SaveTo  (std::ostream & save);
  virtual std::istream & LoadFrom(std::istream & load);

  /*!
   * \brief Does nothing
   * \param theMesh - the built mesh
   * \param theShape - the geometry of interest
   * \retval bool - always false
   */
  virtual bool SetParametersByMesh(const SMESH_Mesh* theMesh, const TopoDS_Shape& theShape);

  /*!
   * \brief Initialize my parameter values by default parameters.
   *  \retval bool - true if parameter values have been successfully defined
   */
  virtual bool SetParametersByDefaults(const TDefaults& dflts, const SMESH_Mesh* theMesh=0);

private:

  // General
  Fineness      _fineness;
  bool          _secondOrder;
  bool          _quadAllowed;

  // Mesh size
  double        _maxSize, _minSize;
  double        _growthRate;
  std::string   _meshSizeFile;
  double        _nbSegPerRadius;
  double        _nbSegPerEdge;
  // (SALOME additions)
  TLocalSize    _localSize;
  bool          _chordalErrorEnabled;
  double        _chordalError;

  // Optimizer
  bool          _optimize;
  int           _nbSurfOptSteps;
  int           _nbVolOptSteps;
  double        _elemSizeWeight;
  int           _worstElemMeasure;

  // Insider
  bool          _surfaceCurvature;
  bool          _useDelauney;
  bool          _checkOverlapping;
  bool          _checkChartBoundary;
  //bool          _blockFilling; -- not used by netgen
  // (SALOME additions)
  bool          _fuseEdges;
};

#endif
