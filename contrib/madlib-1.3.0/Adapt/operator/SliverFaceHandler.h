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

#ifndef _H_SLIVERFACEHANDLER
#define _H_SLIVERFACEHANDLER

#include "MAdOperatorBase.h"

#include <string>

namespace MAd {

  // -------------------------------------------------------------------
  class sliverFaceHandler {

  public:

    sliverFaceHandler(pMesh, DiscreteSF *);
    ~sliverFaceHandler() {}

  public:

    void removeSliverFaces(int * nbIn=NULL, int * nbOut=NULL);

    void setSizeField(DiscreteSF *);

    // whether we can use operations on boundaries
    void collapseOnBoundary(bool, double);

    // whether we can use operations that add nodes
    void newVertexPermission(bool);
    bool getNewVertexPermission() const { return addNodes; }

    void enableReport(std::string);

  private:

    int findOperation(pFace, pMAdOperator *);

    bool F_isSliver(pFace);
    bool F_isSliver(pFace, double *);

    void reportSliver(pFace);

  private:

    pMesh mesh;
    DiscreteSF * sizeField;
    MeshQualityManager& mqm;

    bool collapseOnBdry;
    double dATol;

    bool addNodes; // whether or not we can use operators that add nodes 
                   // (edge split)

    bool   reportFailures;
    int    reportId;
    std::string reportPrefix;

  };

  // -------------------------------------------------------------------
  inline sliverFaceHandler::sliverFaceHandler(pMesh m, DiscreteSF * sf):
    mesh(m),sizeField(sf), mqm(MeshQualityManagerSgl::instance()), 
    collapseOnBdry(false), dATol(MAdTOL), addNodes(true), 
    reportFailures(false),reportId(0),reportPrefix("")
  {}

  // -------------------------------------------------------------------
  inline void sliverFaceHandler::setSizeField(DiscreteSF * sf)
  {
    sizeField = sf;
  }

  // -------------------------------------------------------------------
  inline void sliverFaceHandler::collapseOnBoundary(bool cob, double tolerance)
  {
    collapseOnBdry = cob;
    dATol = tolerance;
  }

  // -------------------------------------------------------------------
  inline void sliverFaceHandler::newVertexPermission(bool allow)
  {
    addNodes = allow;
  }

  // -------------------------------------------------------------------
  inline void sliverFaceHandler::enableReport(std::string prefix)
  {
    reportFailures = true; 
    reportPrefix   = prefix;
  }

  // -------------------------------------------------------------------

}

#endif

