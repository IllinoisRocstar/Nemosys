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
// File      : NETGENPlugin_Hypothesis.cxx
// Author    : Michael Sazonov (OCN)
// Date      : 28/03/2006
// Project   : SALOME
//
#include "NETGENPlugin_Hypothesis.hxx"

#include "NETGENPlugin_Mesher.hxx"
#include "SMESH_Mesh.hxx"

#include <utilities.h>

using namespace std;

//=============================================================================
/*!
 *
 */
//=============================================================================
NETGENPlugin_Hypothesis::NETGENPlugin_Hypothesis (int hypId, SMESH_Gen * gen)

  : SMESH_Hypothesis(hypId, gen),
    _fineness           (GetDefaultFineness()),
    _secondOrder        (GetDefaultSecondOrder()),
    _quadAllowed        (GetDefaultQuadAllowed()),
    _maxSize            (GetDefaultMaxSize()),
    _minSize            (0),
    _growthRate         (GetDefaultGrowthRate()),
    _nbSegPerRadius     (GetDefaultNbSegPerRadius()),
    _nbSegPerEdge       (GetDefaultNbSegPerEdge()),
    _chordalErrorEnabled(GetDefaultChordalError() > 0),
    _chordalError       (GetDefaultChordalError() ),
    _optimize           (GetDefaultOptimize()),
    _nbSurfOptSteps     (GetDefaultNbSurfOptSteps()),
    _nbVolOptSteps      (GetDefaultNbVolOptSteps()),
    _elemSizeWeight     (GetDefaultElemSizeWeight()),
    _worstElemMeasure   (GetDefaultWorstElemMeasure()),
    _surfaceCurvature   (GetDefaultSurfaceCurvature()),
    _useDelauney        (GetDefaultUseDelauney()),
    _checkOverlapping   (GetDefaultCheckOverlapping()),
    _checkChartBoundary (GetDefaultCheckChartBoundary()),
    _fuseEdges          (GetDefaultFuseEdges())
{
  _name = "NETGEN_Parameters";
  _param_algo_dim = 3;
}

//=============================================================================
/*!
 *
 */
//=============================================================================
void NETGENPlugin_Hypothesis::SetMaxSize(double theSize)
{
  if (theSize != _maxSize)
  {
    _maxSize = theSize;
    NotifySubMeshesHypothesisModification();
  }
}

//=============================================================================
/*!
 *  
 */
//=============================================================================
void NETGENPlugin_Hypothesis::SetMinSize(double theSize)
{
  if (theSize != _minSize)
  {
    _minSize = theSize;
    NotifySubMeshesHypothesisModification();
  }
}

//=============================================================================
/*!
 *  
 */
//=============================================================================
void NETGENPlugin_Hypothesis::SetSecondOrder(bool theVal)
{
  if (theVal != _secondOrder)
  {
    _secondOrder = theVal;
    NotifySubMeshesHypothesisModification();
  }
}

//=============================================================================
/*!
 *  
 */
//=============================================================================
void NETGENPlugin_Hypothesis::SetOptimize(bool theVal)
{
  if (theVal != _optimize)
  {
    _optimize = theVal;
    NotifySubMeshesHypothesisModification();
  }
}

//=============================================================================
/*!
 *  
 */
//=============================================================================
void NETGENPlugin_Hypothesis::SetFineness(Fineness theFineness)
{
  if (theFineness != _fineness)
  {
    _fineness = theFineness;
    // the predefined values are taken from NETGEN 4.5 sources
    switch (_fineness)
    {
    case VeryCoarse:
      _growthRate = 0.7;
      _nbSegPerEdge = 0.3;
      _nbSegPerRadius = 1;
      break;
    case Coarse:
      _growthRate = 0.5;
      _nbSegPerEdge = 0.5;
      _nbSegPerRadius = 1.5;
      break;
    case Fine:
      _growthRate = 0.2;
      _nbSegPerEdge = 2;
      _nbSegPerRadius = 3;
      break;
    case VeryFine:
      _growthRate = 0.1;
      _nbSegPerEdge = 3;
      _nbSegPerRadius = 5;
      break;
    case UserDefined:
      break;
    case Moderate:
    default:
      _growthRate = 0.3;
      _nbSegPerEdge = 1;
      _nbSegPerRadius = 2;
      break;
    }
    NotifySubMeshesHypothesisModification();
  }
}

//=============================================================================
/*!
 *  
 */
//=============================================================================
void NETGENPlugin_Hypothesis::SetGrowthRate(double theRate)
{
  if (theRate != _growthRate)
  {
    _growthRate = theRate;
    _fineness = UserDefined;
    NotifySubMeshesHypothesisModification();
  }
}

//=============================================================================
/*!
 *  
 */
//=============================================================================
void NETGENPlugin_Hypothesis::SetNbSegPerEdge(double theVal)
{
  if (theVal != _nbSegPerEdge)
  {
    _nbSegPerEdge = theVal;
    _fineness = UserDefined;
    NotifySubMeshesHypothesisModification();
  }
}

//=============================================================================
/*!
 *  
 */
//=============================================================================
void NETGENPlugin_Hypothesis::SetNbSegPerRadius(double theVal)
{
  if (theVal != _nbSegPerRadius)
  {
    _nbSegPerRadius = theVal;
    _fineness = UserDefined;
    NotifySubMeshesHypothesisModification();
  }
}

//=============================================================================
/*!
 *  
 */
//=============================================================================
void NETGENPlugin_Hypothesis::SetChordalErrorEnabled(bool theVal)
{
  if (theVal != _chordalErrorEnabled)
  {
    _chordalErrorEnabled = theVal;
    NotifySubMeshesHypothesisModification();
  }
}

//=============================================================================
/*!
 *  
 */
//=============================================================================
void NETGENPlugin_Hypothesis::SetChordalError(double theVal)
{
  if (theVal != _chordalError)
  {
    _chordalError = theVal;
    NotifySubMeshesHypothesisModification();
  }
}

//=============================================================================
/*!
 *  
 */
//=============================================================================
void NETGENPlugin_Hypothesis::SetLocalSizeOnEntry(const std::string& entry, double localSize)
{
  if(_localSize[entry] != localSize)
  {
    _localSize[entry] = localSize;
    NotifySubMeshesHypothesisModification();
  }
}

//=============================================================================
/*!
 *
 */
//=============================================================================
double NETGENPlugin_Hypothesis::GetLocalSizeOnEntry(const std::string& entry)
{
  TLocalSize::iterator it  = _localSize.find( entry );
  if ( it != _localSize.end() )
    return it->second;
  else
    return -1.0;
}

//=============================================================================
/*!
 *  
 */
//=============================================================================
void NETGENPlugin_Hypothesis::UnsetLocalSizeOnEntry(const std::string& entry)
{
  _localSize.erase(entry);
  NotifySubMeshesHypothesisModification();
}

//=============================================================================
/*!
 *  
 */
//=============================================================================
void NETGENPlugin_Hypothesis::SetMeshSizeFile(const std::string& fileName)
{
  if ( fileName != _meshSizeFile )
  {
    _meshSizeFile = fileName;
    NotifySubMeshesHypothesisModification();
  }
}

//=============================================================================
/*!
 *  
 */
//=============================================================================
void NETGENPlugin_Hypothesis::SetQuadAllowed(bool theVal)
{
  if (theVal != _quadAllowed)
  {
    _quadAllowed = theVal;
    NotifySubMeshesHypothesisModification();
  }
}

//=============================================================================
/*!
 *  
 */
//=============================================================================
void NETGENPlugin_Hypothesis::SetSurfaceCurvature(bool theVal)
{
  if (theVal != _surfaceCurvature)
  {
    _surfaceCurvature = theVal;
    NotifySubMeshesHypothesisModification();
  }
}

//=============================================================================
/*!
 *
 */
//=============================================================================
void NETGENPlugin_Hypothesis::SetFuseEdges(bool theVal)
{
  if (theVal != _fuseEdges)
  {
    _fuseEdges = theVal;
    NotifySubMeshesHypothesisModification();
  }
}

//=======================================================================
//function : SetNbSurfOptSteps
//purpose  : 
//=======================================================================

void NETGENPlugin_Hypothesis::SetNbSurfOptSteps( int theVal )
{
  if (theVal != _nbSurfOptSteps)
  {
    _nbSurfOptSteps = theVal;
    NotifySubMeshesHypothesisModification();
  }
}

//=======================================================================
//function : SetNbVolOptSteps
//purpose  : 
//=======================================================================

void NETGENPlugin_Hypothesis::SetNbVolOptSteps( int theVal )
{
  if (theVal != _nbVolOptSteps)
  {
    _nbVolOptSteps = theVal;
    NotifySubMeshesHypothesisModification();
  }
}

//=======================================================================
//function : SetElemSizeWeight
//purpose  : 
//=======================================================================

void NETGENPlugin_Hypothesis::SetElemSizeWeight( double theVal )
{
  if (theVal != _elemSizeWeight)
  {
    _elemSizeWeight = theVal;
    NotifySubMeshesHypothesisModification();
  }
}

//=======================================================================
//function : SetWorstElemMeasure
//purpose  : 
//=======================================================================

void NETGENPlugin_Hypothesis::SetWorstElemMeasure( int theVal )
{
  if (theVal != _worstElemMeasure)
  {
    _worstElemMeasure = theVal;
    NotifySubMeshesHypothesisModification();
  }
}

//=======================================================================
//function : SetUseDelauney
//purpose  : 
//=======================================================================

void NETGENPlugin_Hypothesis::SetUseDelauney( bool theVal )
{
  if (theVal != _useDelauney )
  {
    _useDelauney = theVal;
    NotifySubMeshesHypothesisModification();
  }
}

//=======================================================================
//function : SetCheckOverlapping
//purpose  : 
//=======================================================================

void NETGENPlugin_Hypothesis::SetCheckOverlapping( bool theVal )
{
  if (theVal != _checkOverlapping )
  {
    _checkOverlapping = theVal;
    NotifySubMeshesHypothesisModification();
  }
}

//=======================================================================
//function : SetCheckChartBoundary
//purpose  : 
//=======================================================================

void NETGENPlugin_Hypothesis::SetCheckChartBoundary( bool theVal )
{
  if (theVal != _checkChartBoundary)
  {
    _checkChartBoundary = theVal;
    NotifySubMeshesHypothesisModification();
  }
}

//=============================================================================
/*!
 *
 */
//=============================================================================
ostream & NETGENPlugin_Hypothesis::SaveTo(ostream & save)
{
  save << _maxSize << " " << _fineness;

  if (_fineness == UserDefined)
    save << " " << _growthRate << " " << _nbSegPerEdge << " " << _nbSegPerRadius;

  save << " " << (int)_secondOrder << " " << (int)_optimize;

  TLocalSize::iterator it_sm  = _localSize.begin();
  if (it_sm != _localSize.end()) {
    save << " " << "__LOCALSIZE_BEGIN__";
    for ( ; it_sm != _localSize.end(); ++it_sm ) {
      save << " " << it_sm->first
           << " " << it_sm->second << "%#"; // "%#" is a mark of value end
    }
    save << " " << "__LOCALSIZE_END__";
  }
  save << " " << _minSize;
  save << " " << _quadAllowed;
  save << " " << _surfaceCurvature;
  save << " " << _fuseEdges;

  save << " " << _meshSizeFile.size() << " " << _meshSizeFile;

  save << " " << ( _chordalErrorEnabled ? _chordalError : 0. );


  // added for option set completion

  save << " " << _nbSurfOptSteps;
  save << " " << _nbVolOptSteps;
  save << " " << _elemSizeWeight;
  save << " " << _worstElemMeasure;

  save << " " << _useDelauney;
  save << " " << _checkOverlapping;
  save << " " << _checkChartBoundary;

  return save;
}

//=============================================================================
/*!
 *  
 */
//=============================================================================
istream & NETGENPlugin_Hypothesis::LoadFrom(istream & load)
{
  bool isOK = true;
  int is;
  double val;

  isOK = static_cast<bool>(load >> val);
  if (isOK)
    _maxSize = val;
  else
    load.clear(ios::badbit | load.rdstate());

  isOK = static_cast<bool>(load >> is);
  if (isOK)
    SetFineness((Fineness) is);
  else
    load.clear(ios::badbit | load.rdstate());

  if (_fineness == UserDefined)
  {
    isOK = static_cast<bool>(load >> val);
    if (isOK)
      _growthRate = val;
    else
      load.clear(ios::badbit | load.rdstate());

    isOK = static_cast<bool>(load >> val);
    if (isOK)
      _nbSegPerEdge = val;
    else
      load.clear(ios::badbit | load.rdstate());

    isOK = static_cast<bool>(load >> val);
    if (isOK)
      _nbSegPerRadius = val;
    else
      load.clear(ios::badbit | load.rdstate());
  }

  isOK = static_cast<bool>(load >> is);
  if (isOK)
    _secondOrder = (bool) is;
  else
    load.clear(ios::badbit | load.rdstate());

  isOK = static_cast<bool>(load >> is);
  if (isOK)
    _optimize = (bool) is;
  else
    load.clear(ios::badbit | load.rdstate());

  std::string option_or_sm;
  bool hasLocalSize = false;

  isOK = static_cast<bool>(load >> option_or_sm);
  if (isOK)
    if (option_or_sm == "__LOCALSIZE_BEGIN__")
      hasLocalSize = true;

  std::string smEntry, smValue;
  while (isOK && hasLocalSize) {
    isOK = static_cast<bool>(load >> smEntry);
    if (isOK) {
      if (smEntry == "__LOCALSIZE_END__")
        break;
      isOK = static_cast<bool>(load >> smValue);
    }
    if (isOK) {
      std::istringstream tmp(smValue);
      double val;
      tmp >> val;
      _localSize[ smEntry ] = val;
    }
  }

  if ( !hasLocalSize && !option_or_sm.empty() )
    _minSize = atof( option_or_sm.c_str() );
  else
    load >> _minSize;

  isOK = static_cast<bool>( load >> is );
  if ( isOK )
    _quadAllowed = (bool) is;
  else
    _quadAllowed = GetDefaultQuadAllowed();

  isOK = static_cast<bool>( load >> is );
  if ( isOK )
    _surfaceCurvature = (bool) is;
  else
    _surfaceCurvature = GetDefaultSurfaceCurvature();

  isOK = static_cast<bool>( load >> is );
  if ( isOK )
    _fuseEdges = (bool) is;
  else
    _fuseEdges = GetDefaultFuseEdges();

  isOK = static_cast<bool>( load >> is >> std::ws ); // size of meshSizeFile
  if ( isOK && is > 0 )
  {
    _meshSizeFile.resize( is );
    load.get( &_meshSizeFile[0], is+1 );
  }

  isOK = static_cast<bool>(load >> val);
  if (isOK)
    _chordalError = val;
  else
    load.clear(ios::badbit | load.rdstate());
  _chordalErrorEnabled = ( _chordalError > 0 );


  // added for option set completion

  isOK = static_cast<bool>( load >> is );
  if ( isOK )
    _nbSurfOptSteps = is;

  isOK = static_cast<bool>( load >> is );
  if ( isOK )
    _nbVolOptSteps = is;

  isOK = static_cast<bool>( load >> val );
  if ( isOK )
    _elemSizeWeight =  val;

  isOK = static_cast<bool>( load >> is );
  if ( isOK )
    _worstElemMeasure = is;

  isOK = static_cast<bool>( load >> is );
  if ( isOK )
    _useDelauney = (bool) is;

  isOK = static_cast<bool>( load >> is );
  if ( isOK )
    _checkOverlapping = (bool) is;

  isOK = static_cast<bool>( load >> is );
  if ( isOK )
    _checkChartBoundary = (bool) is;

  return load;
}

//================================================================================
/*!
 * \brief Does nothing
 * \param theMesh - the built mesh
 * \param theShape - the geometry of interest
 * \retval bool - always false
 */
//================================================================================
bool NETGENPlugin_Hypothesis::SetParametersByMesh(const SMESH_Mesh*   /*theMesh*/,
                                                  const TopoDS_Shape& /*theShape*/)
{
  return false;
}

//================================================================================
/*!
 * \brief Initialize my parameter values by default parameters.
 *  \retval bool - true if parameter values have been successfully defined
 */
//================================================================================

bool NETGENPlugin_Hypothesis::SetParametersByDefaults(const TDefaults&  dflts,
                                                      const SMESH_Mesh* theMesh)
{
  _nbSegPerEdge = dflts._nbSegments;
  _maxSize      = dflts._elemLength;

  if ( dflts._shape && !dflts._shape->IsNull() )
    _minSize    = NETGENPlugin_Mesher::GetDefaultMinSize( *dflts._shape, _maxSize );
  else if ( theMesh && theMesh->HasShapeToMesh() )
    _minSize    = NETGENPlugin_Mesher::GetDefaultMinSize( theMesh->GetShapeToMesh(), _maxSize );

  if ( dflts._way == SMESH_Hypothesis::BY_AVERAGE_LENGTH )
  {
    _minSize      = dflts._elemLength / 100.;
    _nbSegPerEdge = 1;
    _chordalError = dflts._elemLength / 2.;
    _chordalErrorEnabled = true;
    _quadAllowed  = dflts._quadDominated;
  }

  return _nbSegPerEdge && _maxSize > 0;
}
