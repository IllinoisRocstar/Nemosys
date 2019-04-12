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
// Authors: Jean-Francois Remacle, Gaetan Compere
// -------------------------------------------------------------------

#ifndef _H_MAD_GMSHENTITIES
#define _H_MAD_GMSHENTITIES

#ifdef _HAVE_GMSH_

#include "Physical.h"

#include "gmsh/GEntity.h"
#include "gmsh/discreteRegion.h"
#include "gmsh/discreteFace.h"
#include "gmsh/discreteEdge.h"
#include "gmsh/discreteVertex.h"

namespace MAd {

  class GmshModel;
  class GmshGEntity;
  class GmshGRegion;
  class GmshGFace;
  class GmshGEdge;
  class GmshGVertex;

  // -------------------------------------------------------------------
  class GmshGEntity {

  public:
    GmshGEntity(GmshModel *m, Physical *p): model(m), phys(p)
    {}
    GmshGEntity(const GmshGEntity& ge): model(ge.model),phys(ge.phys)
    {}
    virtual ~GmshGEntity()
    {}

    void setPhysical(int d, int t);

    virtual int dim() const=0; // to make it pure virtual (GC)
    virtual int tag() const=0;
    int pDim() const { if(phys) {return phys->dim();} else {return -1;} }
    int pTag() const { if(phys) {return phys->tag();} else {return  0;} }

  protected:
    GmshModel * model;
    Physical * phys;
  };

  // -------------------------------------------------------------------
  class GmshGEntityLessThan {
  public:
    bool operator()(const GmshGEntity *ent1, const GmshGEntity *ent2) const
    { return ent1->tag() < ent2->tag(); }
  };

  // -------------------------------------------------------------------
  class GmshGRegion : public GmshGEntity {

  private:
    GRegion * gr;

  public:
    GmshGRegion(GmshModel& model, int tag, Physical *_p=NULL);
    GmshGRegion(const GmshGRegion& ggr);
    virtual ~GmshGRegion();
    int dim() const { return 3; }
    virtual int tag() const { return gr->tag(); }

    std::list<GmshGFace *> faces() const;
  };

  // -------------------------------------------------------------------
  class GmshGFace : public GmshGEntity {

  private:
    GFace * gf;

  public:
    GmshGFace(GmshModel& model, int tag, Physical *_p=NULL);
    GmshGFace(const GmshGFace& ggf);
    virtual ~GmshGFace();
    int dim() const { return 2; }
    virtual int tag() const { return gf->tag(); }
    
    GFace * getGFace() const { return gf; }

    int numRegions() const { return gf->numRegions(); }
    std::list<GmshGEdge*> edges() const; 

    // compute the parameters UV from a point XYZ
    void XYZtoUV(const double X, const double Y, const double Z,
                 double &U, double &V, const double relax,
                 const bool onSurface=true) const 
    { gf->XYZtoUV(X,Y,Z,U,V,relax,onSurface); }

    // return the point on the face corresponding to the given parameter
    GPoint point(double par1, double par2) const;

    // compute the min and max curvatures and the corresponding directions
    // return the max curvature
    // outputs have to be allocated before calling this function
    double curvatures(const SPoint2 &param, SVector3 *dirMax, SVector3 *dirMin,
                      double *curvMax, double *curvMin) const;

  // return the curvature computed as the divergence of the normal
    double curvatureDiv(const SPoint2 &param) const;
    
    // compute, in parametric space, the interpolation from pt1 to pt2
    // along a geodesic of the surface
    SPoint2 geodesic(const SPoint2 &pt1, const SPoint2 &pt2, double t);
  };

  // -------------------------------------------------------------------
  class GmshGEdge : public GmshGEntity {

  private:
    GEdge * ge;

  public:
    GmshGEdge(GmshModel& model, int tag, Physical *_p=NULL);
    GmshGEdge(const GmshGEdge& gge);
    virtual ~GmshGEdge();
    int dim() const { return 1; }
    virtual int tag() const { return ge->tag(); }

    std::list<GmshGVertex *> vertices() const;
    GmshGVertex *getBeginVertex() const;
    GmshGVertex *getEndVertex() const;

    // reparamaterize the point onto the given face
    SPoint2 reparamOnFace(const GmshGFace *face, double epar, int dir) const;

    // get bounds of parametric coordinate 
    Range<double> parBounds(int i) const;
  
    // get the curvature
    double curvature(double par) const;

    // true if the edge is a seam for the given face.
    bool isSeam(const GmshGFace *face) const;

    // get the point for the given parameter location
    GPoint point(double p) const;
  };

  // -------------------------------------------------------------------
  class GmshGVertex : public GmshGEntity {

  private:
    GVertex * gv;

  public:
    GmshGVertex(GmshModel& model, int tag, Physical *_p=NULL);
    GmshGVertex(const GmshGVertex& ggv);
    virtual ~GmshGVertex();
    int dim() const { return 0; }
    virtual int tag() const { return gv->tag(); }

    // get the edges that this vertex bounds
    std::list<GmshGEdge*> edges() const;

    SPoint2 reparamOnFace(const GmshGFace *gf, int i) const;

    // return true if this vertex is on a seam of the given face
    bool isOnSeam(const GmshGFace *gf) const;
  };

  // -------------------------------------------------------------------
  typedef class GmshGEntity MAdGEntity;
  typedef class GmshGEntityLessThan MAdGEntityLessThan;
  typedef class GmshGRegion MAdGRegion;
  typedef class GmshGFace   MAdGFace;
  typedef class GmshGEdge   MAdGEdge;
  typedef class GmshGVertex MAdGVertex;

  // -------------------------------------------------------------------
}

#endif

#endif
