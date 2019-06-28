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

#ifndef _H_ADAPTINTERFACE
#define _H_ADAPTINTERFACE

#include "madlib_export.h"

// from Geo
#include "ModelInterface.h"

// from Mesh
#include "MeshDataBaseInterface.h"
#include "CheckMesh.h"
#ifdef PARALLEL
#include "MeshDataBaseComm.h"
#endif

// from Adapt
#include "MAdOutput.h"
#include "CallbackDefinition.h"

#include <set>
#include <string>
#include <vector>

namespace MAd {

  class SizeFieldManager;
  class geoMatcher;
  class mobileObjectSet;
  class MeshParametersManager;
  
  // Local mesh modification operators
  class edgeSplitOp;
  class edgeCollapseOp;
  class faceCollapseOp;
  class DESCOp;
  class edgeSwapOp;
  class faceSwapOp;
  class vertexMoveOp;
  class regionRemoveOp;
  class sliverFaceHandler;
  class sliverRegionHandler;
  
  // -------------------------------------------------------------------

  enum algorithmDefinition {
    SPT_SWP_SLV_CPS,
    CPS_SWP_SLV_SPT,
    SLV_CPS_SWP_SPT
  };
#ifdef PARALLEL
  enum loadBalanceAlgorithm {
    DEFAULT_ALGORITHM,
    METIS_ALGORITHM
  };
#endif
  // -------------------------------------------------------------------
  // -------------------------------------------------------------------
  class MADLIB_EXPORT MeshAdapter {

  public:
  
    MeshAdapter(pMesh m, pSField sf=NULL);
    ~MeshAdapter();

    // --------------------------------------------

    // gives an overview of the current parameters
    void printParameters() const;

    void addCallback(CBFunction CB, void* userData, 
                     CBFunc_move CB_move=0, void* userData_move=0);
    void addSizeField(pSField sf);

    // whether or not the global size field is smoothed and maximum gradient
    void setSizeFieldSmoothing(bool, double maxGrad=1.);

    // set the maximum number of iterations available to reach the 'convergence' 
    // of the global mesh adaptation procedure
    void setMaxIterationsNumber(int max);

#ifdef PARALLEL
    // set load balancing algorithm
    void setLoadBalancingAlgorithm( loadBalanceAlgorithm lbAlgo );
    // set the data Exchanger
//     void setDataExchanger( MDB_DataExchanger* dataExch);
#endif

    // impose the interval governing edges length in the transformed space
    void setEdgeLenSqBounds(double lower, double upper);

    // set permission for operators to modify slightly boundaries 
    // (default: false) and tolerance for the relative volume/area 
    // modifications.
    // Edge collapse and face collapse:
    void setCollapseOnBoundary(bool accept=true, double tolerance=1.e-6);
    // Edge swap:
    void setSwapOnBoundary(bool accept=true, double tolerance=1.e-6);

    // impose the element quality from which a swap is not required
    void setNoSwapQuality(double noSwapQuality);

    // impose a threshold rate of improvement for edge and face swaps
    void setSwapMinImproveRatio(double ratio);

    // impose the maximum quality of a sliver
    void setSliverQuality(double sliverQuality);

    // Allow or forbid operations creating a sliver. 
    // If allowed and a bound is prescribed, it is allowed only if the edge is 
    // longer/shorter than the bound (bound is expressed as the square of 
    // adimensional edge length)
    void setSliverPermissionInESplit   (bool perm, double lenSqBound=-1.);
    void setSliverPermissionInECollapse(bool perm, double lenSqBound=-1.);

    // tell if you want to project new vertices on geometric entities
    void setGeoTracking(bool enable, bool cavityEqualMesh=false,
                        int cavityThickness=2, double chi=1.0,
                        bool strictChecking=false, bool forceRelocation=false);

    // set a value for a very huge edge length (default=1.e14)
    void setInfiniteLength(double length);

    // frequency at which the size fields are updated inside the global loop
    // (useful for local and analytical SF)
    void setSFUpdateFrequency(int freq);

    // Set verbosity:
    //  < 0 : no detail
    //    1 : global procedure details
    //  > 2 : iterations details
    void setVerbosity(int _verbosity);

    // constraint a mesh/geometric entity and all its downward entities
    void clearConstraints() const;
    void setConstraint(pEntity e)        const;
    void setConstraint(pGEntity ge)      const;
    void setConstraint(int type, int id) const;
    // unconstrain a geometric entity
    void removeConstraint(int type, int id) const;
    void removeConstraint(pGEntity ge)      const;

    // manage the physical time
    void   incrementTime(double dt);
    void   setTime(double t);
    double getTime() const;
    void   updateSizeField();

    // functions to keep track of the initial coordinates
    void storeInitialCoordinates();
    void removeStoredCoordinates();

    // add predefined mobile objects
    void registerObjects(mobileObjectSet* objs);

    // will attach/get the datas to the nodes of the mesh
    // order in vector is the node iterator order
    void registerData  (std::string name, const std::vector<double>)               const;
    void registerVData (std::string name, const std::vector<std::vector<double> >) const;
    void getMeshData   (std::string name, std::vector<double> *)                   const;
    void getMeshVData  (std::string name, std::vector<std::vector<double> > *)     const;
    void removeData    (std::string name) const;
    void removeVData   (std::string name) const;

  public:

    // ---------------- OPERATIONS ----------------

    // ----- Level 1 -----
    // ---> Elementary operations

    // Split the edge and if 'checkSize' is true, check that the two resulting 
    // edges are not short edges.
    bool splitEdge (pEdge e, bool checkSize=false);
    bool collapseEdge (pEdge e);
    void collapseEdgeBrute (pEdge e); // use it only if you know what your are doing
    bool collapseFace(pFace f, pEdge e);
    bool DSplitCollapseEdge(pRegion pr, pEdge edge1, pEdge edge2);
    bool swapEdge (pEdge e);
    bool swapFace (pFace f);
    bool removeRegion(pRegion region);
    bool moveVertex (pVertex v, double dxyz[3]);
    bool putVertex (pVertex v, double xyz[3]);
    //    bool moveVertices (std::multiset<vDisplacement,vDisplacementLess>& vDisps);

    // ----- Level 2 -----
    // ---> Loops on one elementary operation

    // node repositioning
    double LaplaceSmoothing();

    // topology operations
    int eSplitLoop();
    int eCollapseLoop();
    int eSplitCollapseLoop();
    int edgeSwapLoop();
    int faceSwapLoop();
    int splitEveryEdgeOnce();

    // slivers handling
    int removeSlivers();

    // geometry matching
    void snapVertices();

    // ----- Level 3 -----
    // ---> Procedures with a global objective
    int optimiseEdgeLength();
    int optimiseElementShape();
    int splitLongestEdges();
    int runOneIter();
    void uglyTheMesh(double avgQualThresh, int maxIt);
    int removeNegativeElements();

    // ----- Level 4 -----
    // ---> Global procedure
    void run();

    // ----- objects motion -----
    // move boundaries without repositioning nodes in the volume
    int partlyMoveObjects        (double t, double dt, double* part);
    // move boundaries and reposition all nodes with an elasticity analogy
    //   subAdaptation:    allow sliver elimination and mesh optimization if necessary
    //   qualityThreshold: at which quality do we optimize the mesh (if allowed)
    //   chi:              stiffness alteration coefficient (-1 = none)
    //   meshIsCavity:     true = elastic computation on whole mesh
    //   cavityThickness:  nb layers of elements if mesh is not the cavity
    void moveObjectsAndReposition (double t, double dt, 
                                   bool subAdaptation=true,
                                   double qualityThreshold=0.,
                                   double chi=-1., 
                                   bool meshIsCavity=true, 
                                   int cavityThickness=3);
    
    // --------------------------------------------

  public:

    // ------ Diagnostics ------

    // get informations on mesh quality
    void getStatistics(double * meanQuality, double * worstQuality) const;

    // get information on applied mesh modifications
    void getModificationsInfo(int * nSplit, int * nColl, 
                              int * nSwap, double * cpuSplit,
                              double * cpuColl, double * cpuSwap,
                              double * cpuSliv) const;
    
    // about all datas attached to the nodes
    void nodalDataDiagnostics(std::ostream& out) const;
  
    // journals listing all operations tested or applied
    void setDebugLevel(int debug) { debugLevel = debug; }
    void openJournal() const;
    void setReferenceJournal(std::string& name) const;
    void flushJournal(std::ostream& out) const;

    // sliver outputs
    void enableSliverReports();
    void testSliverOperators(bool test);

    // performs several checks to check the validity of the mesh
    bool checkTheMesh(int verbose=1, 
                      std::ostream& out=std::cout, 
                      MeshStatus * status=NULL) const;

    // get infos about mobile objects
    void infoMobileObjects(std::ostream& out=std::cout) const;

  public:

    // ------ Outputs ------

    // set the path to output directory
    void setOutputPrefix(std::string prefix);

    // write mesh with required postpro data in 'pos' format (Gmsh)
    void writePos(std::string fn, MAdOutputData type=OD_CONSTANT) const;

    // write mesh in 'msh' format (Gmsh)
    void writeMsh(std::string fn) const;

    // write a .pos file with the distance to walls for every local size field
    void writeDistanceToWalls(std::string fnBase) const;
    // write a .pos file with the 'volumic' curvature for every local size field
    void writeVolumicCurvature(std::string fnBase) const;

    // get global data over the mesh
    void printStatistics(std::ostream& out) const;
    void printSliverRegionStatistics(std::ostream& out) const;

  public:

    // save all available informations to output directory and abort
    void abort(int line=-1, const char* file=NULL) const;

    // --------------------------------------------

  private:

    void setDefaultValues();
    void buildOperators();
    void removeOperators();

  private:

    pMesh mesh;
    SizeFieldManager * SFManager;

    mobileObjectSet * objects;

    // ----- Local mesh modification operators -----
    edgeSplitOp          *  eSplitOp;
    edgeCollapseOp       *  eCollapseOp;
    faceCollapseOp       *  fCollapseOp;
    DESCOp               *  descOp;
    edgeSwapOp           *  eSwapOp;
    faceSwapOp           *  fSwapOp;
    vertexMoveOp         *  vMoveOp;
    regionRemoveOp       *  rRegionOp;
    sliverFaceHandler    *  sliverFOp;
    sliverRegionHandler  *  sliverROp;

    // ----- Geometry related -----
    geoMatcher * geoTracker;

    // ----- Adaptation parameters -----
    algorithmDefinition algorithm;
    int maxIterationsNumber;
    MeshParametersManager& mpm;
#ifdef PARALLEL
    loadBalanceAlgorithm load_balance_algorithm;
    MDB_DataExchanger*   dataExchanger;
#endif
    int updateSFFrequency;

    // ----- Output parameters -----
    int verbosity;
    std::string outPrefix;
    int debugLevel;
  };

  // -------------------------------------------------------------------

} // End of namespace MAd

#endif

