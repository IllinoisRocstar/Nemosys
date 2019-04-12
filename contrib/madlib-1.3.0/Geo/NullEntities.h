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

#ifndef _H_NULLENTITIES
#define _H_NULLENTITIES

#ifndef _HAVE_GMSH_

#include "Physical.h"
#include <stdlib.h>

namespace MAd {

  class NullModel;

  // -------------------------------------------------------------------
  class NullGEntity {

  public:
    NullGEntity(int tag, NullModel *m, Physical *_p=NULL): _tag(tag),model(m),phys(_p) {}
    NullGEntity(const NullGEntity& ge): _tag(ge._tag),model(ge.model),phys(ge.phys) {}
    virtual ~NullGEntity() {}

    void setPhysical(int d, int t);

    virtual int dim() const=0;
    virtual int tag() const { return _tag; }
    int pDim() const { if(phys) {return phys->dim();} else {return -1;} }
    int pTag() const { if(phys) {return phys->tag();} else {return  0;} }

  private:
    int _tag;
    Physical * phys;
    NullModel * model;
  };

  // -------------------------------------------------------------------
  class NullGEntityLessThan {
  public:
    bool operator()(const NullGEntity *ent1, const NullGEntity *ent2) const
    { return ent1->tag() < ent2->tag(); }
  };

  // -------------------------------------------------------------------
  class NullGRegion : public NullGEntity {

  public:
    NullGRegion(int tag, NullModel *m, Physical *_p=NULL) : NullGEntity(tag, m, _p) {}
    NullGRegion(const NullGRegion& gr) : NullGEntity(gr) {}
    virtual ~NullGRegion() {}
    int dim() const { return 3; }
  };

  // -------------------------------------------------------------------
  class NullGFace : public NullGEntity {

  public:
    NullGFace(int tag, NullModel *m, Physical *_p=NULL) : NullGEntity(tag, m, _p) {}
    NullGFace(const NullGFace& gf) : NullGEntity(gf) {}
    virtual ~NullGFace() {}
    int dim() const { return 2; }
  };

  // -------------------------------------------------------------------
  class NullGEdge : public NullGEntity {

  public:
    NullGEdge(int tag, NullModel *m, Physical *_p=NULL) : NullGEntity(tag, m, _p) {}
    NullGEdge(const NullGEdge& ge) : NullGEntity(ge) {}
    virtual ~NullGEdge() {}
    int dim() const { return 1; }
  };

  // -------------------------------------------------------------------
  class NullGVertex : public NullGEntity {

  public:
    NullGVertex(int tag, NullModel *m, Physical *_p=NULL) : NullGEntity(tag, m, _p) {}
    NullGVertex(const NullGVertex& gv) : NullGEntity(gv) {}
    virtual ~NullGVertex() {}
    int dim() const { return 0; }
  };

  typedef class NullGEntity MAdGEntity;
  typedef class NullGEntityLessThan MAdGEntityLessThan;
  typedef class NullGRegion MAdGRegion;
  typedef class NullGFace   MAdGFace;
  typedef class NullGEdge   MAdGEdge;
  typedef class NullGVertex MAdGVertex;
  // -------------------------------------------------------------------
}

#endif

#endif

