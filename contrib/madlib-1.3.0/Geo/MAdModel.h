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

#ifndef _H_MADMODEL
#define _H_MADMODEL

#include "Physical.h"

#ifdef _HAVE_GMSH_
#include "GmshEntities.h"
#else
#include "NullEntities.h"
#endif

#include <string>
#include <set>

namespace MAd {

  // -------------------------------------------------------------------
  enum ModelType {
    NULLMODEL,
    GMSHMODEL
  };

  // -------------------------------------------------------------------
  class MAdModel {

  public:

    MAdModel(std::string name="") { _name = name; }
    virtual ~MAdModel() {};

    virtual ModelType type() const = 0;

    // Gmsh mesh file format
    virtual int readMSH(const std::string &name) = 0;
  
    // Gmsh native CAD format
    virtual int readGEO(const std::string &name) = 0;
  
    // OCC formats
    virtual int readSTEP(const std::string &name) = 0;
    virtual int readBREP(const std::string &name) = 0;
    virtual int readIGES(const std::string &name) = 0;

    // get the number of entities in this model
    virtual int getNumRegions()  const = 0;
    virtual int getNumFaces()    const = 0;
    virtual int getNumEdges()    const = 0;
    virtual int getNumVertices() const = 0;

    // find the entity with the given tag and create it if not found
    virtual MAdGEntity * getEntityByTag(int dim, int tag) = 0;
    virtual MAdGRegion * getRegionByTag(int) = 0;
    virtual MAdGFace   * getFaceByTag  (int) = 0;
    virtual MAdGEdge   * getEdgeByTag  (int) = 0;
    virtual MAdGVertex * getVertexByTag(int) = 0;

    typedef std::set<MAdGRegion*, MAdGEntityLessThan>::iterator riter;
    typedef std::set<MAdGFace*,   MAdGEntityLessThan>::iterator fiter;
    typedef std::set<MAdGEdge*,   MAdGEntityLessThan>::iterator eiter;
    typedef std::set<MAdGVertex*, MAdGEntityLessThan>::iterator viter;

    // get an iterator initialized to the first/last entity in this model.
    virtual riter firstRegion() = 0;
    virtual fiter firstFace() = 0;
    virtual eiter firstEdge() = 0;
    virtual viter firstVertex() = 0;
    virtual riter lastRegion() = 0;
    virtual fiter lastFace() = 0;
    virtual eiter lastEdge() = 0;
    virtual viter lastVertex() = 0;

    Physical * getPhysical(int d, int t);
    bool physical() const { if (physicals.size()){return true;} return false; }

  protected:
    std::string _name;
    std::set<Physical *> physicals;
  };

}

// -------------------------------------------------------------------

#endif
