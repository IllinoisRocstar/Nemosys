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
#include "GMSHPlugin_Hypothesis.hxx"

#include "GMSHPlugin_Mesher.hxx"
#include "SMESH_Mesh.hxx"

#include <utilities.h>

using namespace std;


GMSHPlugin_Hypothesis::GMSHPlugin_Hypothesis (int hypId,
                                                  SMESH_Gen * gen)
  : SMESH_Hypothesis(hypId, gen),
    _algo2d         (automatic),
    _algo3d         (frontal3),
    _recomb2DAlgo   (standard),
    _recombineAll   (false),
    _subdivAlgo     (none),
    _remeshAlgo     (nosplit),
    _remeshPara     (harmonic),
    _smouthSteps    (1),
    _sizeFactor     (1),
    _minSize        (0),
    _maxSize        (1e22),
    _secondOrder    (false),
    _useIncomplElem (true)
{
  _name = "GMSH_Parameters";
  _param_algo_dim = 3;
}

void GMSHPlugin_Hypothesis::Set2DAlgo(Algo2D the2DAlgo)
{
  if (the2DAlgo != _algo2d)
  {
    _algo2d = the2DAlgo;
    NotifySubMeshesHypothesisModification();
  }
}

void GMSHPlugin_Hypothesis::Set3DAlgo(Algo3D the3DAlgo)
{
  if (the3DAlgo != _algo3d)
  {
    _algo3d = the3DAlgo;
    NotifySubMeshesHypothesisModification();
  }
}

void GMSHPlugin_Hypothesis::SetRecomb2DAlgo(Recomb2DAlgo theRecomb2DAlgo)
{
  if (theRecomb2DAlgo != _recomb2DAlgo)
  {
    _recomb2DAlgo = theRecomb2DAlgo;
    NotifySubMeshesHypothesisModification();
  }
}

void GMSHPlugin_Hypothesis::SetRecombineAll(bool theRecombineAll)
{
  if (theRecombineAll != _recombineAll)
  {
    _recombineAll = theRecombineAll;
    NotifySubMeshesHypothesisModification();
  }
}

void GMSHPlugin_Hypothesis::SetSubdivAlgo(SubdivAlgo theSubdivAlgo)
{
  if (theSubdivAlgo != _subdivAlgo)
  {
    _subdivAlgo = theSubdivAlgo;
    NotifySubMeshesHypothesisModification();
  }
}

void GMSHPlugin_Hypothesis::SetRemeshAlgo(RemeshAlgo theRemeshAlgo)
{
  if (theRemeshAlgo != _remeshAlgo)
  {
    _remeshAlgo = theRemeshAlgo;
    NotifySubMeshesHypothesisModification();
  }
}

void GMSHPlugin_Hypothesis::SetRemeshPara(RemeshPara theRemeshPara)
{
  if (theRemeshPara != _remeshPara)
  {
    _remeshPara = theRemeshPara;
    NotifySubMeshesHypothesisModification();
  }
}

void GMSHPlugin_Hypothesis::SetSmouthSteps(double theSmouthSteps)
{
  if (theSmouthSteps != _smouthSteps)
  {
    _smouthSteps = theSmouthSteps;
    NotifySubMeshesHypothesisModification();
  }
}

void GMSHPlugin_Hypothesis::SetSizeFactor(double theSizeFactor)
{
  if (theSizeFactor != _sizeFactor)
  {
    _sizeFactor = theSizeFactor;
    NotifySubMeshesHypothesisModification();
  }
}

void GMSHPlugin_Hypothesis::SetUseIncomplElem(bool theUseIncomplElem)
{
  if (theUseIncomplElem != _useIncomplElem)
  {
    _useIncomplElem = theUseIncomplElem;
    NotifySubMeshesHypothesisModification();
  }
}

void GMSHPlugin_Hypothesis::SetMaxSize(double theSize)
{
  if (theSize != _maxSize)
  {
    _maxSize = theSize;
    NotifySubMeshesHypothesisModification();
  }
}

void GMSHPlugin_Hypothesis::SetMinSize(double theSize)
{
  if (theSize != _minSize)
  {
    _minSize = theSize;
    NotifySubMeshesHypothesisModification();
  }
}

void GMSHPlugin_Hypothesis::SetSecondOrder(bool theVal)
{
  if (theVal != _secondOrder)
  {
    _secondOrder = theVal;
    NotifySubMeshesHypothesisModification();
  }
}

void GMSHPlugin_Hypothesis::SetIs2d(bool theIs2d)
{
  _is2d = theIs2d;
}


void GMSHPlugin_Hypothesis::SetCompoundOnEntry(const std::string& entry)
{
  if (_compounds.find(entry) == _compounds.end())
  {
    _compounds.insert(entry);
    NotifySubMeshesHypothesisModification();
  }
}

void GMSHPlugin_Hypothesis::UnsetCompoundOnEntry(const std::string& entry)
{
  if (_compounds.find(entry) != _compounds.end())
  {
    _compounds.erase(entry);
    NotifySubMeshesHypothesisModification();
  }
}

std::ostream & GMSHPlugin_Hypothesis::SaveTo(std::ostream & save)
{
  save << (int)_is2d << " " << _algo2d;
  if (!_is2d)
    save << " " << _algo3d;
  save << " " << _recomb2DAlgo         <<
          " " << (int)_recombineAll    <<
          " " << _subdivAlgo           <<
          " " << _remeshAlgo           <<
          " " << _remeshPara           <<
          " " << _smouthSteps          <<
          " " << _sizeFactor           <<
          " " << _maxSize              <<
          " " << _minSize              <<
          " " << (int)_secondOrder     <<
          " " << (int)_useIncomplElem  ;
  
  save << " " << "__COMPOUNDS_BEGIN__";
  for (TCompound::const_iterator it = _compounds.begin();  it != _compounds.end(); ++it )
    save << " " << *it << " ";
  save << " " << "__COMPOUNDS_END__";
  
  return save;
}

std::istream & GMSHPlugin_Hypothesis::LoadFrom(std::istream & load)
{
  bool isOK = true;
  int is;
  double val;

  isOK = static_cast<bool>(load >> is);
  if (isOK)
    _is2d = (bool)is;
  else
    load.clear(ios::badbit | load.rdstate());
  
  isOK = static_cast<bool>(load >> is);
  if (isOK)
    _algo2d = (Algo2D)is;
  else
    load.clear(ios::badbit | load.rdstate());
  
  if (!_is2d)
  {
    isOK = static_cast<bool>(load >> is);
    if (isOK)
      _algo3d = (Algo3D)is;
    else
      load.clear(ios::badbit | load.rdstate());
  }
  
  isOK = static_cast<bool>(load >> is);
  if (isOK)
    _recomb2DAlgo = (Recomb2DAlgo)is;
  else
    load.clear(ios::badbit | load.rdstate());
  
  isOK = static_cast<bool>(load >> is);
  if (isOK)
    _recombineAll = (bool)is;
  else
    load.clear(ios::badbit | load.rdstate());
  
  isOK = static_cast<bool>(load >> is);
  if (isOK)
    _subdivAlgo = (SubdivAlgo)is;
  else
    load.clear(ios::badbit | load.rdstate());
  
  isOK = static_cast<bool>(load >> is);
  if (isOK)
    _remeshAlgo = (RemeshAlgo)is;
  else
    load.clear(ios::badbit | load.rdstate());
  
  isOK = static_cast<bool>(load >> is);
  if (isOK)
    _remeshPara = (RemeshPara)is;
  else
    load.clear(ios::badbit | load.rdstate());
  
  isOK = static_cast<bool>(load >> val);
  if (isOK)
    _smouthSteps = val;
  else
    load.clear(ios::badbit | load.rdstate());
  
  isOK = static_cast<bool>(load >> val);
  if (isOK)
    _sizeFactor = val;
  else
    load.clear(ios::badbit | load.rdstate());
  
  isOK = static_cast<bool>(load >> val);
  if (isOK)
    _maxSize = val;
  else
    load.clear(ios::badbit | load.rdstate());
  
  isOK = static_cast<bool>(load >> val);
  if (isOK)
    _minSize = val;
  else
    load.clear(ios::badbit | load.rdstate());
  
  isOK = static_cast<bool>(load >> is);
  if (isOK)
    _secondOrder = (bool)is;
  else
    load.clear(ios::badbit | load.rdstate());
  
  isOK = static_cast<bool>(load >> is);
  if (isOK)
    _useIncomplElem = (bool)is;
  else
    load.clear(ios::badbit | load.rdstate());
  
  
  std::string entry;
  isOK = static_cast<bool>(load >> entry);
  if (isOK && entry == "__COMPOUNDS_BEGIN__")
  {
    while (isOK && entry != "__COMPOUNDS_END__")
    {
      isOK = static_cast<bool>(load >> entry);
      if (isOK && entry != "__COMPOUNDS_END__")
        _compounds.insert(entry);
    }
  }
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
bool GMSHPlugin_Hypothesis::SetParametersByMesh(const SMESH_Mesh*   theMesh,
                                                  const TopoDS_Shape& theShape)
{
  return false;
}

//================================================================================
/*!
 * \brief Initialize my parameter values by default parameters.
 *  \retval bool - true if parameter values have been successfully defined
 */
//================================================================================

bool GMSHPlugin_Hypothesis::SetParametersByDefaults(const TDefaults&  dflts,
                                                    const SMESH_Mesh* theMesh)
{
  //_nbSegPerEdge = dflts._nbSegments;
  //_maxSize      = dflts._elemLength;

  //if ( dflts._shape && !dflts._shape->IsNull() )
  //  _minSize    = GMSHPlugin_Mesher::GetDefaultMinSize( *dflts._shape, _maxSize );
  //else if ( theMesh && theMesh->HasShapeToMesh() )
  //  _minSize    = GMSHPlugin_Mesher::GetDefaultMinSize( theMesh->GetShapeToMesh(), _maxSize );

  //return _nbSegPerEdge && _maxSize > 0;
  return false;
}
