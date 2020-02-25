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

#ifdef _HAVE_GMSH_

#include "GmshEntities.h"
#include "GmshModel.h"
#include "gmsh/Context.h"

namespace MAd {

  // -------------------------------------------------------------------
  // -------------------------------------------------------------------
  void GmshGEntity::setPhysical(int d, int t)
  {
//     if ( phys ) phys.delGE(this);
    phys = model->getPhysical(d,t);
//     phys.addGE(this);
  }

  // -------------------------------------------------------------------
  // -------------------------------------------------------------------
  GmshGRegion::GmshGRegion(GmshModel& m, int tag, Physical *_p) : 
    GmshGEntity(&m,_p) 
  {
    gr = model->getModel()->getRegionByTag(tag);
    if ( !gr ) {
      gr = new discreteRegion(model->getModel(),tag);
      model->getModel()->add(gr);
    }
  }

  // -------------------------------------------------------------------
  GmshGRegion::GmshGRegion(const GmshGRegion& ggr) : 
    GmshGEntity(ggr),gr(ggr.gr) {}

  // -------------------------------------------------------------------
   GmshGRegion::~GmshGRegion() {}

  // -------------------------------------------------------------------
  std::vector<GmshGFace *> GmshGRegion::faces() const
  {
    std::vector<GmshGFace*> res;
    std::vector<GFace*> lst = gr->faces();
    std::vector<GFace*>::const_iterator it = lst.begin();
    for (; it != lst.end(); it++ ){
      res.push_back( model->getFaceByTag( (*it)->tag() ) );
    }
    return res;
  }

  // -------------------------------------------------------------------
  // -------------------------------------------------------------------
  GmshGFace::GmshGFace(GmshModel& m, int tag, Physical *_p) : 
    GmshGEntity(&m, _p)
  {
    gf = model->getModel()->getFaceByTag(tag);
    if ( !gf ) {
      gf = new discreteFace(model->getModel(),tag);
      model->getModel()->add(gf);
    }
  }

  // -------------------------------------------------------------------
  GmshGFace::GmshGFace(const GmshGFace& ggf) : 
    GmshGEntity(ggf),gf(ggf.gf) {}

  // -------------------------------------------------------------------
  GmshGFace::~GmshGFace() {}

  // -------------------------------------------------------------------
  std::vector<GmshGEdge*> GmshGFace::edges() const
  {
    std::cout << "1" << std::endl;
    std::vector<GmshGEdge*> res;
    //std::cout << "2" << std::endl;
    std::vector<GEdge*> lst = gf->edges();
    std::vector<GEdge*>::const_iterator it = lst.begin();
    for (; it != lst.end(); it++ ){
      res.push_back( model->getEdgeByTag( (*it)->tag() ) );
    }
    return res;
  }

  // -------------------------------------------------------------------
  GPoint GmshGFace::point(double par1, double par2) const
  {
    return gf->point(par1,par2);
  }

  // -------------------------------------------------------------------
  double GmshGFace::curvatures(const SPoint2 &param, SVector3 *dirMax, 
                               SVector3 *dirMin, double *curvMax, 
                               double *curvMin) const
  {
    return gf->curvatures(param,*dirMax,*dirMin,*curvMax,*curvMin);
  }

  // -------------------------------------------------------------------
  double GmshGFace::curvatureDiv(const SPoint2 &param) const
  {
    return gf->curvatureDiv(param);
  }

  // -------------------------------------------------------------------
  SPoint2 GmshGFace::geodesic(const SPoint2 &pt1, const SPoint2 &pt2, double t)
  {
//    return gf->geodesic(pt1,pt2,t);

  // by default we assume that straight lines are geodesics
  if(CTX::instance()->mesh.secondOrderExperimental && gf->geomType() != GEntity::Plane ){
      // FIXME: this is buggy -- remove the CTX option once we do it in
      // a robust manner
      GPoint gp1 = point(pt1.x(), pt1.y());
      GPoint gp2 = point(pt2.x(), pt2.y());
      SPoint2 guess = pt1 + (pt2 - pt1) * t;
      GPoint gp = gf->closestPoint(SPoint3(gp1.x() + t * (gp2.x() - gp1.x()),
                                           gp1.y() + t * (gp2.y() - gp1.y()),
                                           gp1.z() + t * (gp2.z() - gp1.z())),
                                   (double*)guess);
      if (gp.g())
          return SPoint2(gp.u(), gp.v());
      else
          return pt1 + (pt2 - pt1) * t;
  }
  else{
      return pt1 + (pt2 - pt1) * t;
  }

  }

  // -------------------------------------------------------------------
  // -------------------------------------------------------------------
  GmshGEdge::GmshGEdge(GmshModel& m, int tag, Physical *_p) : 
    GmshGEntity(&m,_p) 
  {
    ge = model->getModel()->getEdgeByTag(tag);
    if ( !ge ) {
      ge = new discreteEdge(model->getModel(),tag,NULL,NULL);
      model->getModel()->add(ge);
    }
  }

  // -------------------------------------------------------------------
  GmshGEdge::GmshGEdge(const GmshGEdge& gge) : 
    GmshGEntity(gge),ge(gge.ge) {}

  // -------------------------------------------------------------------
  GmshGEdge::~GmshGEdge() {}

  // -------------------------------------------------------------------
  std::vector<GmshGVertex *> GmshGEdge::vertices() const
  {
    std::vector<GmshGVertex *> res;
    res.push_back(getBeginVertex());
    res.push_back(getEndVertex());
    return res;
  }

  // -------------------------------------------------------------------
  GmshGVertex *GmshGEdge::getBeginVertex() const
  {
    GVertex * v = ge->getBeginVertex();
    return model->getVertexByTag(v->tag());
  }

  // -------------------------------------------------------------------
  GmshGVertex *GmshGEdge::getEndVertex() const
  {
    GVertex * v = ge->getEndVertex();
    return model->getVertexByTag(v->tag());
  }

  // -------------------------------------------------------------------
  SPoint2 GmshGEdge::reparamOnFace(const GmshGFace *face, double epar, int dir) const
  {
    return ge->reparamOnFace(face->getGFace(),epar,dir);
  }

  // -------------------------------------------------------------------
  Range<double> GmshGEdge::parBounds(int i) const
  {
    return ge->parBounds(i);
  }

  // -------------------------------------------------------------------
  double GmshGEdge::curvature(double par) const
  {
    return ge->curvature(par);
  }

  // -------------------------------------------------------------------
  bool GmshGEdge::isSeam(const GmshGFace *face) const
  {
    return ge->isSeam(face->getGFace());
  }

  // -------------------------------------------------------------------
  GPoint GmshGEdge::point(double p) const
  {
    return ge->point(p);
  }

  // -------------------------------------------------------------------
  // -------------------------------------------------------------------
  GmshGVertex::GmshGVertex(GmshModel& m, int tag, Physical *_p) : 
    GmshGEntity(&m,_p) 
  {
    gv = model->getModel()->getVertexByTag(tag);
    if ( !gv ) {
      gv = new discreteVertex(model->getModel(),tag);
      model->getModel()->add(gv);
    }
  }
  // -------------------------------------------------------------------
  GmshGVertex::GmshGVertex(const GmshGVertex& ggv) : 
    GmshGEntity(ggv),gv(ggv.gv) {}

  // -------------------------------------------------------------------
  GmshGVertex::~GmshGVertex() {}

  // -------------------------------------------------------------------
  std::vector<GmshGEdge*> GmshGVertex::edges() const
  {
    std::vector<GmshGEdge*> res;
    std::vector<GEdge*> lst = gv->edges();
    std::vector<GEdge*>::const_iterator it = lst.begin();
    for (; it != lst.end(); it++ ){
      res.push_back( model->getEdgeByTag( (*it)->tag() ) );
    }
    return res;
  }

  // -------------------------------------------------------------------
  SPoint2 GmshGVertex::reparamOnFace(const GmshGFace *gf, int i) const
  {
    return gv->reparamOnFace(gf->getGFace(),i);
  }

  // -------------------------------------------------------------------
  bool GmshGVertex::isOnSeam(const GmshGFace *gf) const
  {
    return gv->isOnSeam(gf->getGFace());
  }

  // -------------------------------------------------------------------
}

#endif
