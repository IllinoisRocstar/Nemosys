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

#ifndef _H_GM_ITERATORS
#define _H_GM_ITERATORS

#include "ModelInterface.h"

namespace MAd {

  // -------------------------------------------------------------------
  class GM_RegionIterator {

  private:
    MAdModel * model;
    MAdModel::riter iter;

  public:
    GM_RegionIterator(MAdModel * _model):
      model(_model)
    { iter = model->firstRegion(); }
    GM_RegionIterator(const GM_RegionIterator& it):
      model(it.model),iter(it.iter) {}
    ~GM_RegionIterator() {}

    void reset() { iter = model->firstRegion(); }
    MAdGRegion * next()
    {
      if ( iter == model->lastRegion() ) return NULL;
      MAdGRegion * res = *iter;
      iter++;
      return res;
    }
  };

  // -------------------------------------------------------------------
  class GM_FaceIterator {

  private:
    MAdModel * model;
    MAdModel::fiter iter;

  public:
    GM_FaceIterator(MAdModel * _model):
      model(_model)
    { iter = model->firstFace(); }
    GM_FaceIterator(const GM_FaceIterator& it):
      model(it.model),iter(it.iter) {}
    ~GM_FaceIterator() {}

    void reset() { iter = model->firstFace(); }
    MAdGFace * next()
    {
      if ( iter == model->lastFace() ) return NULL;
      MAdGFace * res = *iter;
      iter++;
      return res;
    }
  };

  // -------------------------------------------------------------------
  class GM_EdgeIterator {

  private:
    MAdModel * model;
    MAdModel::eiter iter;

  public:
    GM_EdgeIterator(MAdModel * _model):
      model(_model)
    { iter = model->firstEdge(); }
    GM_EdgeIterator(const GM_EdgeIterator& it):
      model(it.model),iter(it.iter) {}
    ~GM_EdgeIterator() {}

    void reset() { iter = model->firstEdge(); }
    MAdGEdge * next()
    {
      if ( iter == model->lastEdge() ) return NULL;
      MAdGEdge * res = *iter;
      iter++;
      return res;
    }
  };

  // -------------------------------------------------------------------
  class GM_VertexIterator {

  private:
    MAdModel * model;
    MAdModel::viter iter;

  public:
    GM_VertexIterator(MAdModel * _model):
      model(_model)
    { iter = model->firstVertex(); }
    GM_VertexIterator(const GM_VertexIterator& it):
      model(it.model),iter(it.iter) {}
    ~GM_VertexIterator() {}

    void reset() { iter = model->firstVertex(); }
    MAdGVertex * next()
    {
      if ( iter == model->lastVertex() ) return NULL;
      MAdGVertex * res = *iter;
      iter++;
      return res;
    }
  };

  // -------------------------------------------------------------------

}

#endif
