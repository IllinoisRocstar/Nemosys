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
#ifndef _GMSHPlugin_Hypothesis_HXX_
#define _GMSHPlugin_Hypothesis_HXX_

#include "GMSHPlugin_Defs.hxx"

#include "SMESH_Hypothesis.hxx"
#include "Utils_SALOME_Exception.hxx"

#include <set>

//  Parameters for work of GMSH
//

class GMSHPLUGIN_EXPORT GMSHPlugin_Hypothesis: public SMESH_Hypothesis
{
public:

  GMSHPlugin_Hypothesis(int hypId, SMESH_Gen * gen);

  enum Algo2D
  {
   automatic,
   meshadapt,
   delaunay,
   frontal,
   delaunayforquad,
   packingparallelograms
  };

  void Set2DAlgo(Algo2D the2DAlgo);
  Algo2D Get2DAlgo() const { return _algo2d; }
  
  enum Algo3D
  {
   frontal3,
   frontaldelaunay,
   fontalhex,
   mmg3d,
   rtree
  };

  void Set3DAlgo(Algo3D the3DAlgo);
  Algo3D Get3DAlgo() const { return _algo3d; }

  enum Recomb2DAlgo
  {
   standard,
   blossom
  };

  void SetRecomb2DAlgo(Recomb2DAlgo theRecomb2DAlgo);
  Recomb2DAlgo GetRecomb2DAlgo() const { return _recomb2DAlgo; }
  
  void SetRecombineAll(bool theRecombineAll);
  bool GetRecombineAll() const { return _recombineAll; }

  enum SubdivAlgo
  {
   none,
   allquads,
   allhexas
  };
  
  void SetSubdivAlgo(SubdivAlgo theSubdivAlgo);
  SubdivAlgo GetSubdivAlgo() const { return _subdivAlgo; }

  enum RemeshAlgo
  {
   nosplit,
   automaticR,
   automaticmetis
  };
  
  void SetRemeshAlgo(RemeshAlgo theRemeshAlgo);
  RemeshAlgo GetRemeshAlgo() const { return _remeshAlgo; }

  enum RemeshPara
  {
   harmonic,
   conformal,
   rbfharmonic
  };
  
  void SetRemeshPara(RemeshPara theRemeshPara);
  RemeshPara GetRemeshPara() const { return _remeshPara; }
  
  void SetSmouthSteps(double theSmouthSteps);
  double GetSmouthSteps() const { return _smouthSteps; }
  
  void SetSizeFactor(double theSizeFactor);
  double GetSizeFactor() const { return _sizeFactor; }
  
  void SetUseIncomplElem(bool theUseIncomplElem);
  bool GetUseIncomplElem() const { return _useIncomplElem; }
  
  void SetMaxSize(double theSize);
  double GetMaxSize() const { return _maxSize; }
  
  void SetMinSize(double theSize);
  double GetMinSize() const { return _minSize; }

  void SetSecondOrder(bool theVal);
  bool GetSecondOrder() const { return _secondOrder; }

  void SetIs2d(bool theIs2d);
  bool GetIs2d() const { return _is2d; }
  
  typedef std::set<std::string> TCompound;
  void SetCompoundOnEntry(const std::string& entry);
  const TCompound& GetCompoundOnEntries() const { return _compounds; }
  void UnsetCompoundOnEntry(const std::string& entry);
  
  // Persistence
  virtual std::ostream & SaveTo(std::ostream & save);
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
  Algo2D        _algo2d;
  Algo3D        _algo3d;
  Recomb2DAlgo  _recomb2DAlgo;
  bool          _recombineAll;
  SubdivAlgo    _subdivAlgo;
  RemeshAlgo    _remeshAlgo;
  RemeshPara    _remeshPara;
  double        _smouthSteps;
  double        _sizeFactor;
  double        _minSize, _maxSize;
  bool          _secondOrder, _useIncomplElem;
  bool          _is2d;
  TCompound     _compounds;
};

#endif
