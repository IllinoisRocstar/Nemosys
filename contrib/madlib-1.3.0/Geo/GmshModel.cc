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

#include "GmshModel.h"
using std::set;

namespace MAd {

  // -------------------------------------------------------------------
  GmshModel::~GmshModel()
  {
    if (model) delete model;
    clean();
  }

  // -------------------------------------------------------------------
  void GmshModel::clean()
  {
    riter itR = regions.begin();
    for (; itR != regions.end(); itR++) delete (*itR);
    regions.clear();
    fiter itF = faces.begin();
    for (; itF != faces.end(); itF++) delete (*itF);
    faces.clear();
    eiter itE = edges.begin();
    for (; itE != edges.end(); itE++) delete (*itE);
    edges.clear();
    viter itV = vertices.begin();
    for (; itV != vertices.end(); itV++) delete (*itV);
    vertices.clear();
  }

  // -------------------------------------------------------------------
  int GmshModel::readMSH(const std::string &filename)
  {
    clean();
    int flag = model->readMSH(filename);
    importFromModel();
    return flag;
  }

  // MS
  // -------------------------------------------------------------------
  int GmshModel::readVTK(const std::string &filename)
  {
    clean();
    int flag = model->readVTK(filename);
    importFromModel();
    return flag;
  }
  // MS END

  // -------------------------------------------------------------------
  int GmshModel::readGEO(const std::string &filename)
  {
    clean();
    int flag = model->readGEO(filename);
    importFromModel();
    return flag;
  }
  
  // -------------------------------------------------------------------
  int GmshModel::readSTEP(const std::string &filename)
  {
    clean();
    int flag = model->readOCCSTEP(filename);
    importFromModel();
    return flag;
  }

  // -------------------------------------------------------------------
  int GmshModel::readBREP(const std::string &filename)
  {
    clean();
    int flag =  model->readOCCBREP(filename);
    importFromModel();
    return flag;
  }

  // -------------------------------------------------------------------
  int GmshModel::readIGES(const std::string &filename)
  { 
    clean();
    int flag = model->readOCCIGES(filename);
    importFromModel();
    return flag;
  }
  
  // -------------------------------------------------------------------
  void GmshModel::importFromModel()
  {
    std::set<GRegion*, GEntityLessThan>::iterator rit = model->firstRegion();
    for (; rit != model->lastRegion(); rit++) getRegionByTag((*rit)->tag());
    std::set<GFace*, GEntityLessThan>::iterator fit = model->firstFace();
    for (; fit != model->lastFace(); fit++) getFaceByTag((*fit)->tag());
    std::set<GEdge*, GEntityLessThan>::iterator eit = model->firstEdge();
    for (; eit != model->lastEdge(); eit++) getEdgeByTag((*eit)->tag());
    std::set<GVertex*, GEntityLessThan>::iterator vit = model->firstVertex();
    for (; vit != model->lastVertex(); vit++) getVertexByTag((*vit)->tag());
  }

  // -------------------------------------------------------------------
  GmshGEntity * GmshModel::getEntityByTag(int dim, int tag)
  {
    switch(dim) {
    case 3: return getRegionByTag(tag);
    case 2: return getFaceByTag(tag);
    case 1: return getEdgeByTag(tag);
    case 0: return getVertexByTag(tag);
    }
    return NULL;
  }

  // -------------------------------------------------------------------
  GmshGRegion * GmshModel::getRegionByTag(int tag)
  {
    std::set<GmshGRegion *,GmshGEntityLessThan>::const_iterator it = regions.begin();
    for (; it != regions.end(); it++) {
      if ( (*it)->tag() == tag ) return *it;
    }
    GmshGRegion * newRegion = new GmshGRegion(*this,tag);
    regions.insert(newRegion);
    return newRegion;
  }

  // -------------------------------------------------------------------
  GmshGFace * GmshModel::getFaceByTag(int tag)
  {
    set<GmshGFace *,GmshGEntityLessThan>::const_iterator it = faces.begin();
    for (; it != faces.end(); it++) {
      if ( (*it)->tag() == tag ) return *it;
    }
    GmshGFace * newFace = new GmshGFace(*this,tag);
    faces.insert(newFace);
    return newFace;
  }

  // -------------------------------------------------------------------
  GmshGEdge * GmshModel::getEdgeByTag(int tag)
  {
    set<GmshGEdge *,GmshGEntityLessThan>::const_iterator it = edges.begin();
    for (; it != edges.end(); it++) {
      if ( (*it)->tag() == tag ) return *it;
    }
    GmshGEdge * newEdge = new GmshGEdge(*this,tag);
    edges.insert(newEdge);
    return newEdge;
  }

  // -------------------------------------------------------------------
  GmshGVertex * GmshModel::getVertexByTag(int tag)
  {
    set<GmshGVertex *,GmshGEntityLessThan>::const_iterator it = vertices.begin();
    for (; it != vertices.end(); it++) {
      if ( (*it)->tag() == tag ) return *it;
    }
    GmshGVertex * newVertex = new GmshGVertex(*this,tag);
    vertices.insert(newVertex);
    return newVertex;
  }

  // -------------------------------------------------------------------
}

#endif
