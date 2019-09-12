// -*- C++ -*-
// -------------------------------------------------------------------
// MAdLib - Copyright (C) 2008-2009 Universite catholique de Louvain
//
// See the Copyright.txt and License.txt files for license information. 
// You should have received a copy of these files along with MAdLib. 
// If not, see <http://www.madlib.be/license/>
//
// Please report all bugs and problems to <contrib@madlib.be>
//
// Authors: Gaetan Compere, Jean-Francois Remacle
// -------------------------------------------------------------------

#ifndef _H_LOCALSIZEFIELD
#define _H_LOCALSIZEFIELD

#include "madlib_export.h"

#include "SizeFieldBase.h"
#include "DistanceFunction.h"

#include <set>
#include <utility>
#include <string>

// -------------------------------------------------------------------
/*
  This class stores a size field around an object. 
  * The object is defined by a set of geometrical entities.
  * The size field is defined by a radius and strings giving the 
    size in function of the distance to the object.
  * Returns:
    - If the distance is superior to the radius, the class 
      returns a huge size.
    - Otherwise, it can return one of these sizes:
      - an isotropic size (evaluated from 'sizeN')
      - an anisotropic size with a main axis in the normal 
        direction to the object and a corresponding value ('sizeN'), 
        and two axis in the tangent directions with another value 
        ('sizeT').
*/
// -------------------------------------------------------------------

namespace MAd {

  class MAdStringFieldEvaluator;

  // -------------------------------------------------------------------
  class MADLIB_EXPORT LocalSizeField : public SizeFieldBase {

  public:

    LocalSizeField(pMesh m, std::string name="", bool _distToFaces=false);
    LocalSizeField(const LocalSizeField& _lsf);
    ~LocalSizeField();

    sFieldType getType() const { return LOCALSFIELD; }

    void scale(double);

    void addGeometricEntity(int type, int tag);
    void updateTree();

    void setIsoSize(double _radius, std::string _sizeN);
    void setAnisoSize(double _radius, std::string _sizeN, std::string _sizeT);
    void setCurvatureLimiter(double onTgSize, double _maxCurv=1.e12);

    // get the size at a location (allocate space!)
    pMSize getSize(const pVertex) const;
    pMSize getSizeOnEntity(const pEntity, const double[3]) const;

    // edge length (squared)
    double SF_VV_lengthSq(const pVertex, const pVertex) const;
    double SF_XYZ_lengthSq(const double[3], const double[3],
                           const pMSize, const pMSize=NULL) const;

    // face area (squared)
    double SF_F_areaSq(const pFace) const;
    double SF_XYZ_areaSq(const double[3][3], const pMSize,
                         const double[3]) const;

    // region volume
    double SF_R_volume(const pRegion) const;
    double SF_XYZ_volume(const double[4][3], const pMSize) const;

    // center and its associated size
    double SF_E_center(const pEdge, double[3], double * reducSq, pMSize *) const;
    double SF_VV_center(const pVertex, const pVertex,
                        double[3], double * reducSq, pMSize *) const;

    // returns the distance to the walls
    double getDistance(const pVertex pv) const;

    // divergence of the curvature
    bool getCurvature(const pVertex pv, double *c) const;
    
  private:

    void collectEntities(std::set<pVertex> * verts, 
                         std::set<pEntity> * ents) const;

  private:

    pMesh mesh;

    std::set<pGEntity> geomEntities;
    int geoDim;

    bool isotropic;
    double radius;
    std::string sizeN, sizeT;
    MAdStringFieldEvaluator * sEvalN, * sEvalT;

    bool distToFaces; // true:  distance computed to the faces of the wall (exact)
                      // false: distance computed to the vertices of the wall
    distanceFunction dFct; // tool to compute distance and its derivatives

    // limiters based on curvature
    bool limit;
    double tgSizeLimitCoef;
    double maxCurv;
  };

}

// -------------------------------------------------------------------
#endif
