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

#ifndef _H_GMSHMODEL
#define _H_GMSHMODEL

#ifdef _HAVE_GMSH_

#include "MAdModel.h"

#include "gmsh/GModel.h"
#include "GmshEntities.h"
#include <list>
#include <set>

namespace MAd {

  // -------------------------------------------------------------------
  class GmshModel : public MAdModel {

  private:
    GModel * model;

  public:

    GmshModel(std::string name=""): 
      MAdModel(name)
    { 
      model = new GModel(_name);
      GModel::current(GModel::list.size() - 1);
    }
    ~GmshModel();

    ModelType type() const { return GMSHMODEL; }

    // Gmsh mesh file format
    int readMSH(const std::string &filename);
  
    // Gmsh native CAD format
    int readGEO(const std::string &filename);
  
    // OCC formats
    int readSTEP(const std::string &filename);
    int readBREP(const std::string &filename);
    int readIGES(const std::string &filename);

    GModel * getModel() { return model; }

    // get the number of entities in this model
    int getNumRegions()  const { return regions.size();  }
    int getNumFaces()    const { return faces.size();    }
    int getNumEdges()    const { return edges.size();    }
    int getNumVertices() const { return vertices.size(); }

    // find the entity with the given tag and create it if not found
    GmshGEntity * getEntityByTag(int dim, int tag);
    GmshGRegion * getRegionByTag(int);
    GmshGFace   * getFaceByTag  (int);
    GmshGEdge   * getEdgeByTag  (int);
    GmshGVertex * getVertexByTag(int);

    // get an iterator initialized to the first/last entity in this model.
    riter firstRegion() { return regions.begin(); }
    fiter firstFace()   { return faces.begin(); }
    eiter firstEdge()   { return edges.begin(); }
    viter firstVertex() { return vertices.begin(); }
    riter lastRegion()  { return regions.end(); }
    fiter lastFace()    { return faces.end(); }
    eiter lastEdge()    { return edges.end(); }
    viter lastVertex()  { return vertices.end(); }

  private:

    void importFromModel();

    void clean();

    std::set<GmshGRegion *,GmshGEntityLessThan> regions;
    std::set<GmshGFace   *,GmshGEntityLessThan> faces;
    std::set<GmshGEdge   *,GmshGEntityLessThan> edges;
    std::set<GmshGVertex *,GmshGEntityLessThan> vertices;

  };

}

// -------------------------------------------------------------------

#endif

#endif
