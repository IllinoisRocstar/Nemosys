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

#ifndef _H_SLIVERREGIONHANDLER
#define _H_SLIVERREGIONHANDLER

#include "MAdOperatorBase.h"

#include <set>
#include <iostream>
#include <string>

#define MAKESTATS 1

namespace MAd {

  // -------------------------------------------------------------------
  class sliverRegionHandler {

  public:

    sliverRegionHandler(pMesh, DiscreteSF *);
    ~sliverRegionHandler() {}

  public:

    void removeSliverRegions(int * nbIn=NULL, int * nbOut=NULL);

    void setSizeField(DiscreteSF *);

    // whether we can use operations on boundaries
    void collapseOnBoundary(bool, double);
    void swapOnBoundary(bool, double);

    // whether we can use operations that add nodes
    void newVertexPermission(bool);
    bool getNewVertexPermission() const { return addNodes; }

    void enableReport(std::string);
    void setTestAllOperators(bool);

  private:

    int findOperation(pRegion, pMAdOperator *);
    bool R_isSliver(pRegion);
    bool R_isSliver(pRegion, double *);

    int getSliverType(const pRegion, pEntity[4]);

    void reportSliver(pRegion);

  private:

    pMesh mesh;
    DiscreteSF * sizeField;

    MeshQualityManager& mqm;

    bool collapseOnBdry; // whether or not we can collapse boundary edges
    double dVTolClps;
    bool swapOnBdry; // whether or not we can swap boundary edges
    double dVTolSwap;

    bool addNodes; // whether or not we can use operators that add nodes 
                   // (edge split, face collapse, DESC)

    bool   reportFailures;
    int    reportId;
    std::string reportPrefix;

    bool testOperators;

#ifdef MAKESTATS

  private:

    // statistics:

    int numTypeI_In, numTypeI_Out, numTypeII_In, numTypeII_Out;
    int numSwapTested_I, numSwapAvailable_I, numSwapOperated_I, numSwapToSliver_I;
    int numEClpsTested_I, numEClpsAvailable_I, numEClpsOperated_I, numEClpsToSliver_I;
    int numFClpsTested_I, numFClpsAvailable_I, numFClpsOperated_I, numFClpsToSliver_I;
    int numDESCTested_I, numDESCAvailable_I, numDESCOperated_I, numDESCToSliver_I;
    int numSplitTested_I, numSplitAvailable_I, numSplitOperated_I, numSplitToSliver_I;
    int numVMoveTested_I, numVMoveAvailable_I, numVMoveOperated_I, numVMoveToSliver_I;

    int numSwapTested_II, numSwapAvailable_II, numSwapOperated_II, numSwapToSliver_II;
    int numEClpsTested_II, numEClpsAvailable_II, numEClpsOperated_II, numEClpsToSliver_II;
    int numFSwapTested_II, numFSwapAvailable_II, numFSwapOperated_II, numFSwapToSliver_II;
    int numFClpsTested_II, numFClpsAvailable_II, numFClpsOperated_II, numFClpsToSliver_II;
    int numVMoveTested_II, numVMoveAvailable_II, numVMoveOperated_II, numVMoveToSliver_II;

    // remaining slivers
    std::set<pRegion> sliversOut;

#endif

    // test operators
    std::set<pRegion> sliversIn;
    int tested_I, tested_II, notSolved_I, notSolved_II, solvedWithSliver_I, solvedWithSliver_II;
    int notSolvedDuplicated_I, notSolvedDuplicated_II, solvedWithSliverDuplicated_I, solvedWithSliverDuplicated_II;

    int swapSolves_I,clpsSolves_I,fclpsSolves_I,descSolves_I,splitSolves_I,nodeSolves_I;
    int swapSolvesToSliver_I,clpsSolvesToSliver_I,fclpsSolvesToSliver_I,descSolvesToSliver_I,splitSolvesToSliver_I,nodeSolvesToSliver_I;

    int swapSolves_II, fswapSolves_II, clpsSolves_II, fclpsSolves_II, nodeSolves_II;
    int swapSolvesToSliver_II, fswapSolvesToSliver_II, clpsSolvesToSliver_II, fclpsSolvesToSliver_II, nodeSolvesToSliver_II;

    // remaining slivers
    std::set<pRegion> notSolved;

  public:

    void testOperations(pRegion);

    void printStats(std::ostream&) const;
    void printOperatorsTest(std::ostream&) const;

  };

  // -------------------------------------------------------------------
  inline sliverRegionHandler::sliverRegionHandler(pMesh m, DiscreteSF * sf):
    mesh(m),sizeField(sf), mqm(MeshQualityManagerSgl::instance()), 
    collapseOnBdry(false),dVTolClps(MAdTOL),swapOnBdry(false),dVTolSwap(MAdTOL), 
    addNodes(true), 
    reportFailures(false),reportId(0),reportPrefix(""),testOperators(false)
#ifdef MAKESTATS
    , numTypeI_In(0), numTypeI_Out(0), numTypeII_In(0), numTypeII_Out(0),
    numSwapTested_I(0), numSwapAvailable_I(0), numSwapOperated_I(0), numSwapToSliver_I(0),
    numEClpsTested_I(0), numEClpsAvailable_I(0), numEClpsOperated_I(0), numEClpsToSliver_I(0),
    numFClpsTested_I(0), numFClpsAvailable_I(0), numFClpsOperated_I(0), numFClpsToSliver_I(0),
    numDESCTested_I(0), numDESCAvailable_I(0), numDESCOperated_I(0), numDESCToSliver_I(0),
    numSplitTested_I(0), numSplitAvailable_I(0), numSplitOperated_I(0), numSplitToSliver_I(0),
    numVMoveTested_I(0), numVMoveAvailable_I(0), numVMoveOperated_I(0), numVMoveToSliver_I(0),
    numSwapTested_II(0), numSwapAvailable_II(0), numSwapOperated_II(0), numSwapToSliver_II(0),
    numEClpsTested_II(0), numEClpsAvailable_II(0), numEClpsOperated_II(0), numEClpsToSliver_II(0),
    numFSwapTested_II(0), numFSwapAvailable_II(0), numFSwapOperated_II(0), numFSwapToSliver_II(0),
    numFClpsTested_II(0), numFClpsAvailable_II(0), numFClpsOperated_II(0), numFClpsToSliver_II(0),
    numVMoveTested_II(0), numVMoveAvailable_II(0), numVMoveOperated_II(0), numVMoveToSliver_II(0),
#endif
    tested_I(0), tested_II(0), notSolved_I(0), notSolved_II(0), solvedWithSliver_I(0), solvedWithSliver_II(0), 
    notSolvedDuplicated_I(0), notSolvedDuplicated_II(0), solvedWithSliverDuplicated_I(0), solvedWithSliverDuplicated_II(0),
    swapSolves_I(0),clpsSolves_I(0),fclpsSolves_I(0),descSolves_I(0),splitSolves_I(0),nodeSolves_I(0), swapSolvesToSliver_I(0),clpsSolvesToSliver_I(0),fclpsSolvesToSliver_I(0),descSolvesToSliver_I(0),splitSolvesToSliver_I(0),nodeSolvesToSliver_I(0), 
    swapSolves_II(0), fswapSolves_II(0), clpsSolves_II(0), fclpsSolves_II(0), nodeSolves_II(0), swapSolvesToSliver_II(0), fswapSolvesToSliver_II(0), clpsSolvesToSliver_II(0), fclpsSolvesToSliver_II(0), nodeSolvesToSliver_II(0)
  {}

  // -------------------------------------------------------------------
  inline void sliverRegionHandler::setSizeField(DiscreteSF * sf)
  {
    sizeField = sf;
  }

  // -------------------------------------------------------------------
  inline void sliverRegionHandler::collapseOnBoundary(bool cob, double tolerance)
  {
    collapseOnBdry = cob;
    dVTolClps = tolerance;
  }

  // -------------------------------------------------------------------
  inline void sliverRegionHandler::swapOnBoundary(bool sob, double tolerance)
  {
    swapOnBdry = sob;
    dVTolSwap = tolerance;
  }

  // -------------------------------------------------------------------
  inline void sliverRegionHandler::newVertexPermission(bool allow)
  {
    addNodes = allow;
  }

  // -------------------------------------------------------------------
  inline void sliverRegionHandler::enableReport(std::string prefix)
  {
    reportFailures = true; 
    reportPrefix   = prefix;
  }

  // -------------------------------------------------------------------
  inline void sliverRegionHandler::setTestAllOperators(bool test)
  {
    testOperators = test;
  }

  // -------------------------------------------------------------------

}

#endif
