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

#ifndef _H_NULLMODEL
#define _H_NULLMODEL

#ifndef _HAVE_GMSH_

#include "MAdModel.h"
#include "NullEntities.h"
#include <set>

namespace MAd {

  // -------------------------------------------------------------------
  class NullModel : public MAdModel {

  public:

    NullModel(std::string name=""): MAdModel(name) {}
    ~NullModel();

    ModelType type() const { return NULLMODEL; }

    // Gmsh mesh file format
    int readMSH(const std::string &filename);
  
    // Gmsh native CAD format
    int readGEO(const std::string &filename);
  
    // OCC formats
    int readSTEP(const std::string &filename);
    int readBREP(const std::string &filename);
    int readIGES(const std::string &filename);

    // get the number of entities in this model
    int getNumRegions()  const { return regions.size();  }
    int getNumFaces()    const { return faces.size();    }
    int getNumEdges()    const { return edges.size();    }
    int getNumVertices() const { return vertices.size(); }

    // find the entity with the given tag and create it if not found
    MAdGEntity * getEntityByTag(int dim, int tag);
    MAdGRegion * getRegionByTag(int);
    MAdGFace   * getFaceByTag  (int);
    MAdGEdge   * getEdgeByTag  (int);
    MAdGVertex * getVertexByTag(int);

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

    // delete all geometric entities
    void clean();

    std::set<MAdGRegion *,MAdGEntityLessThan> regions;
    std::set<MAdGFace   *,MAdGEntityLessThan> faces;
    std::set<MAdGEdge   *,MAdGEntityLessThan> edges;
    std::set<MAdGVertex *,MAdGEntityLessThan> vertices;

  };

}

// -------------------------------------------------------------------

#endif

#endif
