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

// from Common
#include "MAdMessage.h"

#include "AdaptInterface.h"
#include "MAdResourceManager.h"
#include "MAdTimeManager.h"
#include "NodalDataManager.h"
#include "MAdStatistics.h"
#include "ModelConstraintManager.h"
#include "History.h"
#include "SizeFieldManager.h"
#include "DistanceFunction.h"
#include "MeshParametersManager.h"
#include "CallbackManager.h"

// operators
#include "MAdOperatorBase.h"
#include "EdgeSplitOp.h"
#include "EdgeCollapseOp.h"
#include "FaceCollapseOp.h"
#include "DESCOp.h"
#include "EdgeSwapOp.h"
#include "FaceSwapOp.h"
#include "VertexMoveOp.h"
#include "RegionRemoveOp.h"
#include "SliverFaceHandler.h"
#include "SliverRegionHandler.h"

// repositioning
#include "GeoMatcher.h"
#include "NodesRepositioningOp.h"
#include "LaplaceSmoothingOp.h"
#include "MobileObject.h"

// standard C/C++
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <sstream>

#ifdef PARALLEL
#include "mpi.h"
#include "MeshDataBaseParallelInterface.h"
#endif
using std::cerr;
using std::cout;
using std::endl;
using std::set;
using std::vector;
using std::multiset;
using std::string;
using std::stringstream;

namespace MAd {

  // -------------------------------------------------------------------
  void MeshAdapterAbortFct(void * data) {
    ((MeshAdapter *) data) -> abort();
  }

  // -------------------------------------------------------------------
  MeshAdapter::MeshAdapter(pMesh m, pSField sf): 
    mesh(m),
    eSplitOp(0), eCollapseOp(0), fCollapseOp(0), descOp(0), 
    eSwapOp(0), fSwapOp(0), vMoveOp(0), rRegionOp(0), sliverFOp(0), 
    sliverROp(0), geoTracker(0), 
    mpm(MeshParametersManagerSgl::instance()), 
    updateSFFrequency(0), verbosity(2), outPrefix(""), debugLevel(0)
  {
    SFManager = new SizeFieldManager(mesh,sf);

    MAdMsgSgl            ::instance().registerAbortFct(MeshAdapterAbortFct,this);
    MAdTimeManagerSgl    ::instance().initialize();
    MAdStatisticsSgl     ::instance().initialize();
    CallBackManagerSgl   ::instance().initialize();
    NodalDataManagerSgl  ::instance().initialize(mesh);
    MeshQualityManagerSgl::instance().initialize(mesh,SFManager->getSizeField(),MEANRATIO);
    HistorySgl           ::instance().initialize();
    HistorySgl           ::instance().closeJournal();
    ModelConstraintManagerSgl::instance().initialize(M_model(mesh));
    MeshParametersManagerSgl::instance().initialize();

    if (mesh) buildOperators();

    setDefaultValues();

    if ( debugLevel >= 1 && !checkTheMesh() ) abort(__LINE__,__FILE__);
  }

  // -------------------------------------------------------------------
  MeshAdapter::~MeshAdapter()
  {
    removeOperators();
    if (geoTracker) delete geoTracker;
    if (SFManager) delete SFManager;

    MAdTimeManagerSgl    ::instance().finalize();
    MAdStatisticsSgl     ::instance().finalize();
    CallBackManagerSgl   ::instance().finalize();
    NodalDataManagerSgl  ::instance().finalize();
    MeshQualityManagerSgl::instance().finalize();
    HistorySgl           ::instance().finalize();
    ModelConstraintManagerSgl::instance().finalize();
    MeshParametersManagerSgl::instance().finalize();
  }

  // -------------------------------------------------------------------
  // -------------------------------------------------------------------
  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  void MeshAdapter::setDefaultValues()
  {
    algorithm = CPS_SWP_SLV_SPT;
    maxIterationsNumber = 10;

    setEdgeLenSqBounds             ( 1./3. , 3. );

    if ( M_dim(mesh) == 2 ) {
      setNoSwapQuality             ( 0.2 );
    }
    else { //  to be checked
      setNoSwapQuality             ( 0.1 );
    }
    mpm.setSwapMinImproveRatio     ( 1.0 );

    mpm.setSliverTriBound          ( 0.01 );
    mpm.setSliverTetBound          ( 0.02 );

    // infinite loops control
    setSliverPermissionInESplit    ( true, 10. );
    setSliverPermissionInECollapse ( true, 0.1 );

    setCollapseOnBoundary          ( true, 1.e-6 );
    setSwapOnBoundary              ( true, 1.e-6 );

#ifdef PARALLEL
    load_balance_algorithm = DEFAULT_ALGORITHM;
    dataExchanger = new PWLSFieldDE((PWLSField*) SFManager->getSizeField());
#endif

    updateSFFrequency = 0;
  }

  // -------------------------------------------------------------------
  void MeshAdapter::printParameters() const
  {
    mpm.diagnostics();
  }

  // -------------------------------------------------------------------
  void MeshAdapter::setEdgeLenSqBounds(double lower, double upper)
  {
    if ( sqrt(lower) > 0.5*sqrt(upper) ) {
      cout << "Warning: edge length interval is too small ("
           << sqrt(lower) << ", " << sqrt(upper) <<")\n";
    }

    mpm.setLowerLengthSqBound(lower);
    mpm.setUpperLengthSqBound(upper);
  }

  // -------------------------------------------------------------------
  void MeshAdapter::setCollapseOnBoundary(bool accept, double tolerance) 
  {
    eCollapseOp->collapseOnBoundary(accept,tolerance);
    fCollapseOp->collapseOnBoundary(accept,tolerance);
    sliverROp  ->collapseOnBoundary(accept,tolerance);
    sliverFOp  ->collapseOnBoundary(accept,tolerance);
  }

  // -------------------------------------------------------------------
  void MeshAdapter::setSwapOnBoundary(bool accept, double tolerance) 
  {
    eSwapOp  ->swapOnBoundary(accept,tolerance);
    sliverROp->swapOnBoundary(accept,tolerance);
  }

  // -------------------------------------------------------------------
  void MeshAdapter::setMaxIterationsNumber(int max)
  {
    maxIterationsNumber = max;
  }

  // -------------------------------------------------------------------
  void MeshAdapter::setNoSwapQuality(double noSwapQuality)
  {
    mpm.setNoSwapQuality(noSwapQuality);
  }

  // -------------------------------------------------------------------
  void MeshAdapter::setSwapMinImproveRatio(double ratio)
  {
    mpm.setSwapMinImproveRatio(ratio);
  }

  // -------------------------------------------------------------------
  void MeshAdapter::setSliverQuality(double sliverQuality)
  {
    mpm.setSliverTriBound(sliverQuality);
    mpm.setSliverTetBound(sliverQuality);
  }

  // -------------------------------------------------------------------
  void MeshAdapter::setSliverPermissionInESplit(bool perm, double bound)
  {
    mpm.setSliverPermissionInESplit(perm, bound);
  }

  // -------------------------------------------------------------------
  void MeshAdapter::setSliverPermissionInECollapse(bool perm, double bound)
  {
    mpm.setSliverPermissionInECollapse(perm, bound);
  }

  // -------------------------------------------------------------------
  void MeshAdapter::setGeoTracking(bool track, bool cavityEqualMesh, 
                                   int cavityThickness, double chi,
                                   bool strictChecking, bool forceRelocation)
  {
    if ( track && !geoTracker) {
      geoTracker = new geoMatcher(mesh,this);
      geoTracker->setCavityEqualMesh(cavityEqualMesh,cavityThickness);
      geoTracker->setStiffnessAlterationCoef(chi);
      geoTracker->setStrictChecking(strictChecking);
      geoTracker->setForceRelocation(forceRelocation);
      if ( !sliverFOp && !sliverROp ) {
        MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                    "Build sliver handler before setting up geo tracker");
      }
      geoTracker->setSliverFaceHandler(sliverFOp);
      geoTracker->setSliverRegionHandler(sliverROp);
    }
    else if ( !track && geoTracker) { 
      delete geoTracker;
      geoTracker = NULL;
    }
  }

  // -------------------------------------------------------------------
  void MeshAdapter::setInfiniteLength(double length)
  {
    mpm.setBigLength(length);
  }

  // -------------------------------------------------------------------
  void MeshAdapter::setSFUpdateFrequency(int freq)
  {
    updateSFFrequency = freq;
  }

  // -------------------------------------------------------------------
  void MeshAdapter::setVerbosity(int _verbosity)
  {
    verbosity = _verbosity;
  }

  // -------------------------------------------------------------------
  void MeshAdapter::clearConstraints() const
  {
    DeleteConstraint(mesh);
    ModelConstraintManagerSgl::instance().clear();
  }

  // -------------------------------------------------------------------
  void MeshAdapter::setConstraint(pEntity e) const
  {
    EN_constrain(e);

    switch (EN_type(e)) {
    case 1:
      EN_constrain( (pEntity)E_vertex((pEdge)e,0) );
      EN_constrain( (pEntity)E_vertex((pEdge)e,1) );
      break;
    case 2:
      setConstraint( (pEntity)F_edge((pFace)e,0) );
      setConstraint( (pEntity)F_edge((pFace)e,1) );
      setConstraint( (pEntity)F_edge((pFace)e,2) );
      break;
    case 3:
      setConstraint( (pEntity)R_face((pRegion)e,0) );
      setConstraint( (pEntity)R_face((pRegion)e,1) );
      setConstraint( (pEntity)R_face((pRegion)e,2) );
      setConstraint( (pEntity)R_face((pRegion)e,3) );
      break;    
    }
  }

  // -------------------------------------------------------------------
  void MeshAdapter::setConstraint(int type, int id) const
  {
    ModelConstraintManagerSgl::instance().constrain(type,id);
  }

  // -------------------------------------------------------------------
  void MeshAdapter::removeConstraint(int type, int id) const
  {
    ModelConstraintManagerSgl::instance().unconstrain(type,id);
  }

  // -------------------------------------------------------------------
  void MeshAdapter::setConstraint(pGEntity e) const
  {
    ModelConstraintManagerSgl::instance().constrain(e);
  }

  // -------------------------------------------------------------------
  void MeshAdapter::removeConstraint(pGEntity e) const
  {
    ModelConstraintManagerSgl::instance().unconstrain(e);
  }

  // -------------------------------------------------------------------
  void MeshAdapter::removeOperators()
  {
    if (eSplitOp)    { delete eSplitOp;     eSplitOp=0;    }
    if (eCollapseOp) { delete eCollapseOp;  eCollapseOp=0; }
    if (fCollapseOp) { delete fCollapseOp;  fCollapseOp=0; }
    if (descOp)      { delete descOp;       descOp=0;      }
    if (eSwapOp)     { delete eSwapOp;      eSwapOp=0;     }
    if (fSwapOp)     { delete fSwapOp;      fSwapOp=0;     }
    if (vMoveOp)     { delete vMoveOp;      vMoveOp=0;     }
    if (rRegionOp)   { delete rRegionOp;    rRegionOp=0;   }
    if (sliverFOp)   { delete sliverFOp;    sliverFOp=0;   }
    if (sliverROp)   { delete sliverROp;    sliverROp=0;   }
  }

  // -------------------------------------------------------------------
  void MeshAdapter::buildOperators()
  {
    removeOperators();

    eSplitOp    = new edgeSplitOp         (mesh,SFManager->getSizeField());
    eCollapseOp = new edgeCollapseOp      (mesh,SFManager->getSizeField());
    fCollapseOp = new faceCollapseOp      (mesh,SFManager->getSizeField());
    descOp      = new DESCOp              (mesh,SFManager->getSizeField());
    eSwapOp     = new edgeSwapOp          (mesh,SFManager->getSizeField());
    fSwapOp     = new faceSwapOp          (mesh,SFManager->getSizeField());
    vMoveOp     = new vertexMoveOp        (mesh,SFManager->getSizeField(),false);
    rRegionOp   = new regionRemoveOp      (mesh,SFManager->getSizeField());
    sliverROp   = new sliverRegionHandler (mesh,SFManager->getSizeField());
    sliverFOp   = new sliverFaceHandler   (mesh,SFManager->getSizeField());
  }

  // -------------------------------------------------------------------
  void MeshAdapter::addCallback(CBFunction CB, void* userData, 
                                CBFunc_move CB_move, void* userData_move)
  { 
    CallBackManagerSgl::instance().registerCallBack(CB,userData);
    if (CB_move) {
      CallBackManagerSgl::instance().registerCallBackMove(CB_move,userData_move);
    }
  }

  // -------------------------------------------------------------------
  void MeshAdapter::addSizeField(pSField sf)
  {
    SFManager->addSizeField(sf);
  }

  // -------------------------------------------------------------------
  void MeshAdapter::setSizeFieldSmoothing(bool enable, double maxGrad)
  {
    SFManager->setSmoothing(enable,maxGrad);
  }

  // -------------------------------------------------------------------
  void MeshAdapter::incrementTime(double dt) 
  {
    MAdTimeManagerSgl::instance().incrementTime(dt);
    updateSizeField();
  }

  // -------------------------------------------------------------------
  void MeshAdapter::setTime(double t) 
  {
    MAdTimeManagerSgl::instance().setTime(t);
    updateSizeField();
  }

  // -------------------------------------------------------------------
  double MeshAdapter::getTime() const {
    return MAdTimeManagerSgl::instance().getTime();
  }

  // -------------------------------------------------------------------
  void MeshAdapter::updateSizeField()
  {
    SFManager->update();
  }

  // -------------------------------------------------------------------
  void MeshAdapter::storeInitialCoordinates()
  {
    if ( getTime() != 0. ) {
      cout << "Warning: storing coordinates at time " << getTime() << endl;
    }
    NodalDataManagerSgl::instance().storeCoordinates();
  }

  // -------------------------------------------------------------------
  void MeshAdapter::removeStoredCoordinates()
  {
    NodalDataManagerSgl::instance().removeCoordinates();
  }

  // -------------------------------------------------------------------
  void MeshAdapter::registerObjects(mobileObjectSet* objs)
  {
    set<mobileObject*> allObj = objs->getObjects();
  
    set<mobileObject*>::const_iterator it = allObj.begin();
    for (; it != allObj.end(); it++) {
      set<LocalSizeField* > sizes = (*it)->getSizes();
      set<LocalSizeField* >::iterator it = sizes.begin();
      for (; it != sizes.end(); it++) {
        SFManager->addSizeField(*it);
      }
    }

    objects = objs;
  }

  // -------------------------------------------------------------------
  void MeshAdapter::registerData  (string name, 
                                   const vector<double> data) const
  {
    NodalDataManagerSgl::instance().registerData(name,data);
  }

  // -------------------------------------------------------------------
  void MeshAdapter::registerVData (string name, 
                                   const vector<vector<double> > data) const
  {
    NodalDataManagerSgl::instance().registerVData(name,data);
  }

  // -------------------------------------------------------------------
  void MeshAdapter::getMeshData   (string name, 
                                   vector<double> * data) const
  {
    NodalDataManagerSgl::instance().getMeshData(name,data);
  }

  // -------------------------------------------------------------------
  void MeshAdapter::getMeshVData  (string name, 
                                   vector<vector<double> > * data) const
  {
    NodalDataManagerSgl::instance().getMeshVData(name,data);
  }

  // -------------------------------------------------------------------
  void MeshAdapter::removeData    (string name) const
  {
    NodalDataManagerSgl::instance().removeData(name);
  }

  // -------------------------------------------------------------------
  void MeshAdapter::removeVData   (string name) const
  {
    NodalDataManagerSgl::instance().removeVData(name);
  }

  // -------------------------------------------------------------------
  // -------------------------------------------------------------------
  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  // LEVEL 1 OPERATIONS
  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  bool MeshAdapter::splitEdge(pEdge edge, bool checkSize)
  {
    if (debug)
    {
    std::cout << "======================================" << std::endl;
    std::cout << "SPLITEDGE" << std::endl;
    std::cout << "Edge points: " << std::endl;
    std::cout << "p1: x,y,z = " << edge->p1->X << ", " << edge->p1->Y << ", " << edge->p1->Z << std::endl;
    std::cout << "p2: x,y,z = " << edge->p2->X << ", " << edge->p2->Y << ", " << edge->p2->Z << std::endl;
    }


    // set edge to split
    double reducSq = eSplitOp->setSplitEdge(edge);

    // check that the edge is not splitted in short edges
    if ( checkSize )
      {
        double lenSq = SFManager->getSizeField()->SF_E_lengthSq(edge);
        if ( ( reducSq * lenSq ) <= mpm.getLowerLengthSqBound() ) return false;
      }

    // do the checks
    double worstShape;
    if ( !eSplitOp->evaluate(&worstShape) ) return false;

    // check that we do not create a worse sliver ( if asked for )
    // if we can create new slivers, check that the edge is long enough to allow it
    if ( worstShape < mpm.getSliverTetBound() ) 
      {
        double oriWorst;
        MeshQualityManagerSgl::instance().E_worstShape(edge, &oriWorst);
        if ( worstShape < oriWorst ) {
          if ( !mpm.getSliverPermissionInESplit() ) return false;
          if ( mpm.getSliverUpperLengthSqBound() > 0. &&
               mpm.getSliverUpperLengthSqBound() > SFManager->getSizeField()->SF_E_lengthSq(edge) )
            {
              return false;
            }
        }
      }

    // do the job
    eSplitOp->apply();
    return true;
  }

  // -------------------------------------------------------------------
  void MeshAdapter::collapseEdgeBrute(pEdge edge)
  {
    if (debug)
    {
    // tw note: here is important
    std::cout << "======================================" << std::endl;
    std::cout << "COLLAPSEEDGEBRUTE" << std::endl;
    std::cout << "Edge points: " << std::endl;
    std::cout << "p1: x,y,z = " << edge->p1->X << ", " << edge->p1->Y << ", " << edge->p1->Z << std::endl;
    std::cout << "p2: x,y,z = " << edge->p2->X << ", " << edge->p2->Y << ", " << edge->p2->Z << std::endl;
  }

      pVertex vDel = E_vertex(edge,0);
      pVertex vTgt = E_vertex(edge,1);
      if ( V_whatInType(vTgt) != 0 ) {
        vDel = E_vertex(edge,1);
        vTgt = E_vertex(edge,0);
      }
      assert(V_whatInType(vTgt)==0);
      eCollapseOp->setCollapseEdge(edge,vDel,vTgt);
      double shape;
      eCollapseOp->evaluate(&shape);
      eCollapseOp->apply();
  }

  // -------------------------------------------------------------------
  bool MeshAdapter::collapseEdge(pEdge edge)
  {
    if (debug)
    {
    // tw note: here is important
    std::cout << "======================================" << std::endl;
    std::cout << "COLLAPSEEDGE" << std::endl;
    std::cout << "Edge points: " << std::endl;
    std::cout << "p1: x,y,z = " << edge->p1->X << ", " << edge->p1->Y << ", " << edge->p1->Z << std::endl;
    std::cout << "p2: x,y,z = " << edge->p2->X << ", " << edge->p2->Y << ", " << edge->p2->Z << std::endl;
    }

    double shapes[2] = {-1.,-1.};
    int ok[2] = {0,0};

    for (int iDir=0; iDir<2; iDir++) {

      pVertex vDel = E_vertex(edge,iDir);
      pVertex vTgt = E_vertex(edge,1-iDir);
      eCollapseOp->setCollapseEdge(edge,vDel,vTgt);

      // do the checks
      ok[iDir] = eCollapseOp->evaluate(&(shapes[iDir]));

      // check that we do not create a worse sliver ( if asked for )
      // if we can create new slivers, check that the edge is short enough to allow it
      if ( shapes[iDir] < mpm.getSliverTetBound() ) 
        {
          double oriWorst;
          MeshQualityManagerSgl::instance().E_worstShape(edge, &oriWorst);
          if ( shapes[iDir] < oriWorst ) {
            if ( !mpm.getSliverPermissionInECollapse() ) ok[iDir] = 0;
            if ( mpm.getSliverLowerLengthSqBound() > 0. &&
                 mpm.getSliverLowerLengthSqBound() < SFManager->getSizeField()->SF_E_lengthSq(edge) )
              {
                ok[iDir] = 0;
              }
          }
        }
    }
    if (debug)
    {
    std::cout << "----------------" << std::endl;  
    std::cout << "ok[0] = " << ok[0] << std::endl;
    std::cout << "ok[1] = " << ok[1] << std::endl;
    std::cout << "shapes[1] = " << shapes[1] << std::endl;
  }

    // find the best direction
    int best = -1; double worst = -1.;
    if ( ok[0] ) { best = 0; worst = shapes[0]; }
    if ( ok[1] && shapes[1] > worst ) {
      best = 1;
      worst = shapes[1];
    }
  
    // if possible apply an edge collapse
    if ( best >= 0 ) {
      pVertex vDel = E_vertex(edge,best);
      pVertex vTgt = E_vertex(edge,1-best);
      eCollapseOp->setCollapseEdge(edge,vDel,vTgt);
      eCollapseOp->apply();
      if (debug){
      std::cout << "performing collapse..." << std::endl;
      }
      return true;
    }
    if (debug)
    {
    std::cout << "not performing collapse..." << std::endl;
    }
    return false;
  }

  // -------------------------------------------------------------------
  bool MeshAdapter::collapseFace(pFace face, pEdge edge)
  {
    if (debug)
    {
    // tw note: here is important
    std::cout << "======================================" << std::endl;
    std::cout << "COLLAPSEFACE" << std::endl;
    std::cout << "Edge points: " << std::endl;
    std::cout << "p1: x,y,z = " << edge->p1->X << ", " << edge->p1->Y << ", " << edge->p1->Z << std::endl;
    std::cout << "p2: x,y,z = " << edge->p2->X << ", " << edge->p2->Y << ", " << edge->p2->Z << std::endl;

    std::cout << "face points: " << std::endl;
    for (int iNode = 0; iNode < face->getNbNodes(); iNode++)
    {
      std::cout << "Coors = " << face->getNode(iNode)->X << ", " << face->getNode(iNode)->Y << ", " << face->getNode(iNode)->Z << std::endl;
    }
    }

    for( int ClpsOnvt = 1; ClpsOnvt >=0; ClpsOnvt-- )
      {
        fCollapseOp->reset(face, edge, ClpsOnvt);
        
        // do the checks
        double worstShape;
        if ( !fCollapseOp->evaluate(&worstShape) ) continue;
        
        // check that we do not create a sliver
        if ( worstShape < mpm.getSliverTetBound() ) continue;
        
        // do the job
        fCollapseOp->apply();
        return true;
      }
    return false;
  }

  // -------------------------------------------------------------------
  bool MeshAdapter::DSplitCollapseEdge(pRegion pr, pEdge edge1, pEdge edge2)
  {
    if (debug)
    {
    // tw note: here is important
    std::cout << "======================================" << std::endl;
    std::cout << "DSPLITCOLLAPSEEDGE" << std::endl;
    std::cout << "Edge1 points: " << std::endl;
    std::cout << "p1: x,y,z = " << edge1->p1->X << ", " << edge1->p1->Y << ", " << edge1->p1->Z << std::endl;
    std::cout << "p2: x,y,z = " << edge1->p2->X << ", " << edge1->p2->Y << ", " << edge1->p2->Z << std::endl;
    std::cout << "Edge2 points: " << std::endl;
    std::cout << "p1: x,y,z = " << edge2->p1->X << ", " << edge2->p1->Y << ", " << edge2->p1->Z << std::endl;
    std::cout << "p2: x,y,z = " << edge2->p2->X << ", " << edge2->p2->Y << ", " << edge2->p2->Z << std::endl;
  }


    descOp->setDESC(pr,edge1,edge2);

    // do the checks
    double worstShape;
    if ( !descOp->evaluate(&worstShape) ) return false;

    // check that we do not create a sliver
    if ( worstShape < mpm.getSliverTetBound() ) return false;

    // do the job
    descOp->apply();
    return true;
  }

  // -------------------------------------------------------------------
  bool MeshAdapter::swapEdge(pEdge edge)
  {
    if (debug)
    {
    // tw note: here is important
    std::cout << "======================================" << std::endl;
    std::cout << "SWAPEDGE" << std::endl;
    std::cout << "Edge points: " << std::endl;
    std::cout << "p1: x,y,z = " << edge->p1->X << ", " << edge->p1->Y << ", " << edge->p1->Z << std::endl;
    std::cout << "p2: x,y,z = " << edge->p2->X << ", " << edge->p2->Y << ", " << edge->p2->Z << std::endl;
  }

    eSwapOp->setSwapEdge(edge);

    // do the checks
    double worstShape;
    if ( !eSwapOp->evaluate(&worstShape) ) return false;

    // check that we do not create a sliver
    if ( worstShape < mpm.getSliverTetBound() ) return false;

    // do the job
    eSwapOp->apply();
    return true;
  }

  // -------------------------------------------------------------------
  bool MeshAdapter::swapFace(pFace face)
  {
    if (debug)
    {
        // tw note: here is important
    std::cout << "======================================" << std::endl;
    std::cout << "SWAPFACE" << std::endl;
    std::cout << "face points: " << std::endl;
    for (int iNode = 0; iNode < face->getNbNodes(); iNode++)
    {
      std::cout << "Coors = " << face->getNode(iNode)->X << ", " << face->getNode(iNode)->Y << ", " << face->getNode(iNode)->Z << std::endl;
    }
    }

    fSwapOp->setSwapFace(face);
  
    // do the checks
    double worstShape;
    if ( !fSwapOp->evaluate(&worstShape) ) return false;

    // check that we do not create a sliver
    if ( worstShape < mpm.getSliverTetBound() )  return false;

    // do the job
    fSwapOp->apply();
    return true;
  }

  // -------------------------------------------------------------------
  bool MeshAdapter::moveVertex (pVertex v, double dxyz[3])
  {
    if (debug)
    {
    std::cout << "MOVEVERTEX" << std::endl;
    std::cout << "v: x,y,z = " << v->X << ", " << v->Y << ", " << v->Z << std::endl;
  }
    return vMoveOp->move(v,dxyz);
  }

  // -------------------------------------------------------------------
  bool MeshAdapter::removeRegion(pRegion region)
  {
    if (debug)
    {
    // tw note: here is important
    std::cout << "======================================" << std::endl;
    std::cout << "REMOVEREGION" << std::endl;
  }

    rRegionOp->setRegion(region);
  
    // do the checks
    double worstShape;
    if ( !rRegionOp->evaluate(&worstShape) ) return false;

    // check that we do not create a sliver
    if ( worstShape < mpm.getSliverTetBound() )  return false;

    // do the job
    rRegionOp->apply();
    return true;
  }

  // -------------------------------------------------------------------
  bool MeshAdapter::putVertex (pVertex v, double xyz[3])
  {
    if (debug)
    {
    std::cout << "PUTVERTEX" << std::endl;
    std::cout << "v: x,y,z = " << v->X << ", " << v->Y << ", " << v->Z << std::endl;
  }
    double oriPos[3], dxyz[3];
    V_coord(v,oriPos);
    for (int i=0; i<3; i++) dxyz[i] = xyz[i] - oriPos[i];
    return moveVertex(v,dxyz);
  }

  // -------------------------------------------------------------------
//   bool MeshAdapter::moveVertices (multiset<vDisplacement,vDisplacementLess>& vDisps)
//   {
//     return vMoveOp->move(vDisps);
//   }

  // -------------------------------------------------------------------
  // LEVEL 2 OPERATIONS
  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  double MeshAdapter::LaplaceSmoothing()
  {
    MAdResourceManager& tm = MAdResourceManagerSgl::instance();

    double t0 = tm.getTime();

    double convergenceCriterion = 1e-2;
    double L2Disp = 0.;
    LaplaceSmoothingOp* laplOp = new LaplaceSmoothingOp(mesh,SFManager->getSizeField());
    laplOp->run(&L2Disp);
    double L2Disp0 = L2Disp; int i = 0;
    double L2DispTot = 0.;
    while ( ( L2Disp > convergenceCriterion * L2Disp0 ) && ( L2Disp > MAdTOL ) ) {
      laplOp->run(&L2Disp);
      i++;
      L2DispTot += L2Disp;
    }
    if (laplOp) delete laplOp;

    double dt = tm.getTime() - t0;
    if ( verbosity >= 2 ) {
      cout << "Performed a Laplace smoothing: total L2 norm of the displacements: "<<L2DispTot<<" in "<<dt<<" seconds\n";
    }

    if ( debugLevel >= 3 && !checkTheMesh() ) abort(__LINE__,__FILE__);

    return L2DispTot;
  }

  // -------------------------------------------------------------------
  int MeshAdapter::eSplitLoop() 
  {
    MAdResourceManager& tm = MAdResourceManagerSgl::instance();
    MAdStatistics&  stat = MAdStatisticsSgl ::instance();

    double t0 = tm.getTime();

    int nSplit = 0;
    int countE = 0; int oriNbEdges = M_numEdges(mesh);
    EIter ei = M_edgeIter(mesh);
    pEdge edge;
    while ( (countE < oriNbEdges) && ( edge = EIter_next(ei) ) ) {
      double lengthSq = SFManager->getSizeField()->SF_E_lengthSq(edge);
      if ( lengthSq > mpm.getUpperLengthSqBound() ) if (splitEdge(edge,true)) nSplit++;
      countE++;
    }
    EIter_delete(ei);
        
    double dt = tm.getTime() - t0;
    stat.addCPUESplits(dt);
    stat.addNumESplits(nSplit);
    if ( verbosity >= 2 ) {
      cout << "Performed "<< nSplit<<" edge splits in "<<dt<<" seconds\n";
    }

    if ( debugLevel >= 3 && !checkTheMesh() ) abort(__LINE__,__FILE__);

    // --- Geometry tracking ---
    if (nSplit) snapVertices();
    
    return nSplit;
  }

  // -------------------------------------------------------------------
  int MeshAdapter::eCollapseLoop() 
  {
    MAdResourceManager& tm = MAdResourceManagerSgl::instance();
    MAdStatistics&  stat = MAdStatisticsSgl ::instance();

    double t0 = tm.getTime();

    int nCollapse = 0;
    int countE = 0; int oriNbEdges = M_numEdges(mesh);
    EIter ei = M_edgeIter(mesh);
    pEdge edge;
    while ( (countE < oriNbEdges) && ( edge = EIter_next(ei) ) ) {
      double lengthSq = SFManager->getSizeField()->SF_E_lengthSq(edge);
      if ( lengthSq < mpm.getLowerLengthSqBound() ) if (collapseEdge(edge)) nCollapse++;
      countE++;
    }
    EIter_delete(ei);

    double dt = tm.getTime() - t0;
    stat.addCPUECollapses(dt);
    stat.addNumECollapses(nCollapse);
    if ( verbosity >= 2 ) {
      cout << "Performed "<< nCollapse<<" edge collapses in "<<dt<<" seconds\n";
    }

    if ( debugLevel >= 3 && !checkTheMesh() ) abort(__LINE__,__FILE__);

    return nCollapse;
  }

  // -------------------------------------------------------------------
  int MeshAdapter::eSplitCollapseLoop() 
  {
    int nSplit = 0;
    int nCollapse = 0;
    int countE = 0; int oriNbEdges = M_numEdges(mesh);
    EIter ei = M_edgeIter(mesh);
    pEdge edge;
    while ( (countE < oriNbEdges) && ( edge = EIter_next(ei) ) ) {
      double lengthSq = SFManager->getSizeField()->SF_E_lengthSq(edge);
      if ( lengthSq > mpm.getUpperLengthSqBound() ) if (splitEdge(edge,true)) nSplit++;
      if ( lengthSq < mpm.getLowerLengthSqBound() ) if (collapseEdge(edge)) nCollapse++;
      countE++;
    }
    EIter_delete(ei);

    if ( debugLevel >= 3 && !checkTheMesh() ) abort(__LINE__,__FILE__);

    return nSplit+nCollapse;
  }

  // -------------------------------------------------------------------
  int MeshAdapter::edgeSwapLoop()
  {
    if (debug)
    {
    // tw note: here is important
    std::cout << "======================================" << std::endl;
    std::cout << "EDGESWAPLOOP" << std::endl;
  }

    MAdResourceManager& tm = MAdResourceManagerSgl::instance();
    MAdStatistics&  stat = MAdStatisticsSgl ::instance();

    double t0 = tm.getTime();

    int nb_eswap = 0;
    pEdge edge;
    EIter ei = M_edgeIter(mesh);
    int ne = M_numEdges(mesh);
    int i = 0;
    while( ( edge = EIter_next(ei) ) )
    {

    if (debug)
    {
    std::cout << "Edge points: " << std::endl;
    std::cout << "p1: x,y,z = " << edge->p1->X << ", " << edge->p1->Y << ", " << edge->p1->Z << std::endl;
    std::cout << "p2: x,y,z = " << edge->p2->X << ", " << edge->p2->Y << ", " << edge->p2->Z << std::endl;
  }

      if( i > ne ) break;
      i++;
  
      if( EN_constrained((pEntity)edge) )  continue;
      if( E_numRegions(edge) > eSwapOp->getMaxNumRgns() )  continue;
  
      double worstShape;
      MeshQualityManagerSgl::instance().E_worstShape(edge, &worstShape);
  
      if( worstShape > mpm.getNoSwapQuality() ) continue;
  
      eSwapOp->setSwapEdge(edge);

      double newWorstShape;
      if( eSwapOp->evaluate(&newWorstShape) )
      {
        if( newWorstShape > mpm.getSwapMinImproveRatio()*worstShape )
        {
          if( ( eSwapOp->getMaxLenSq() < mpm.getUpperLengthSqBound() ) &&
            ( eSwapOp->getMinLenSq() > mpm.getLowerLengthSqBound() ) )
          {
            eSwapOp->apply();
            nb_eswap++;
          }
        }
      }
    }
    EIter_delete(ei);

    double dt = tm.getTime() - t0;
    stat.addCPUESwaps(dt);
    stat.addNumESwaps(nb_eswap);
    if ( verbosity >= 2 ) {
      cout << "Performed "<< nb_eswap<<" edge swaps in "<<dt<<" seconds\n";
    }

    if ( debugLevel >= 3 && !checkTheMesh() ) abort(__LINE__,__FILE__);

    return nb_eswap;
  }

  // -------------------------------------------------------------------
  int MeshAdapter::faceSwapLoop() 
  {

    if (debug)
    {
    std::cout << "======================================" << std::endl;
    std::cout << "FACESWAPLOOP" << std::endl;
  }

    int nfswap = 0;
    pFace face;
    FIter fi = M_faceIter(mesh);
    int nf = M_numFaces(mesh);
    int i = 0;
    while ( ( face = FIter_next(fi) ) ) {

    if (debug)
    {
        // tw note: here is important
    std::cout << "face points: " << std::endl;
    for (int iNode = 0; iNode < face->getNbNodes(); iNode++)
    {
      std::cout << "Coors = " << face->getNode(iNode)->X << ", " << face->getNode(iNode)->Y << ", " << face->getNode(iNode)->Z << std::endl;
    }
  }

      if( i > nf ) break;
      i++;

      if( EN_constrained((pEntity)face) )  continue;

      double oriWorst;
      MeshQualityManagerSgl::instance().F_worstShape(face, &oriWorst);

      if ( oriWorst > mpm.getNoSwapQuality() ) continue;

      fSwapOp->setSwapFace(face);
      double newWorst;
      if( fSwapOp->evaluate(&newWorst) ) {
        if ( newWorst > mpm.getSwapMinImproveRatio()*oriWorst ) {
          fSwapOp->apply();
          nfswap++;
        }
      }
    }
    FIter_delete(fi);

    if ( debugLevel >= 3 && !checkTheMesh() ) abort(__LINE__,__FILE__);

    return nfswap;
  }

  // -------------------------------------------------------------------
  int MeshAdapter::splitEveryEdgeOnce()
  {
    MAdResourceManager& tm = MAdResourceManagerSgl::instance();
    double t0 = tm.getTime();

    // list initial edges
    std::set<pEdge> initEdges;
    EIter ei = M_edgeIter(mesh);
    pEdge edge;
    while ( ( edge = EIter_next(ei) ) ) initEdges.insert(edge);
    EIter_delete(ei);
    
    // split every edge
    int numOp = 0;
    std::set<pEdge>::const_iterator it = initEdges.begin();
    for (; it != initEdges.end(); it++ ) {
      if ( splitEdge( *it ) ) numOp++;
    }

    double dt = tm.getTime() - t0;
    if ( verbosity >= 2 ) {
      cout << "Performed "<< numOp<<" edge splits in "<<dt<<" seconds\n";
    }

    if ( debugLevel >= 3 && !checkTheMesh() ) abort(__LINE__,__FILE__);

    // --- Geometry tracking ---
    if (numOp) snapVertices();
    
    return numOp;
  }

  // -------------------------------------------------------------------
  int MeshAdapter::removeSlivers()
  {
    MAdResourceManager& tm = MAdResourceManagerSgl::instance();
    MAdStatistics&  stat = MAdStatisticsSgl ::instance();

    double t0 = tm.getTime();

    int nbBefore, nbAfter;
    if (M_dim(mesh) == 3) sliverROp->removeSliverRegions(&nbBefore,&nbAfter);
    else                  sliverFOp->removeSliverFaces  (&nbBefore,&nbAfter);

    double dt = tm.getTime() - t0;
    stat.addCPURSlivers(dt);
    if ( verbosity >= 2 ) {
      cout << "Removed slivers " << nbBefore << " -> " << nbAfter << " in "<<dt<<" seconds\n";
    }

    if ( debugLevel >= 3 && !checkTheMesh() ) abort(__LINE__,__FILE__);

    return ( nbBefore - nbAfter );
  }

  // -------------------------------------------------------------------
  void MeshAdapter::snapVertices()
  {
    if (geoTracker) {
      MAdResourceManager& tm = MAdResourceManagerSgl::instance();
      double t0 = tm.getTime();

      if ( !(geoTracker->snap()) ) {
        MAdMsgSgl::instance().error(__LINE__,__FILE__,"Snapping failed");
      }

      double dt = tm.getTime() - t0;
      MAdMsgSgl::instance().info(__LINE__,__FILE__,
                                 "Performed vertex snapping in %f seconds",dt);
    }

    if ( debugLevel >= 3 && !checkTheMesh() ) abort(__LINE__,__FILE__);
  }

  // -------------------------------------------------------------------
  // LEVEL 3 OPERATIONS
  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  int MeshAdapter::optimiseEdgeLength()
  {
    int nbMaxItOptiLen = 10;
    int nSplColl = 0;
    for (int iter=0; iter < nbMaxItOptiLen; iter++ ) {
      int nSC = eSplitCollapseLoop();
      cout << "Applied "<<nSC<<" split or collapses\n";
      if (!nSC) break;
      nSplColl += nSC;
    }
    return nSplColl;
  }

  // -------------------------------------------------------------------
  int MeshAdapter::optimiseElementShape()
  {
    int nbMaxItOptiShp = 10;
    int nESwapTot = 0, nFSwapTot = 0;
    int nESwap = 0, nFSwap = 0;
    for (int iter=0; iter < nbMaxItOptiShp; iter++ ) {
      // find ill-shaped elements around each edge
      nESwap = edgeSwapLoop();
      cout << "Applied " << nESwap<< " edge swaps"<<endl;
      if ( M_dim(mesh) == 3 ) {
        // find ill-shaped elements around each face
        nFSwap = faceSwapLoop();
        cout << "Applied " << nFSwap<< " face swaps"<<endl;
      }
    
      if (nESwap+nFSwap == 0) break;
      nESwapTot += nESwap;
      nFSwapTot += nFSwap;
    }
    return ( nESwapTot + nFSwapTot );
  }

  // -------------------------------------------------------------------
  int MeshAdapter::splitLongestEdges() 
  {
    MAdResourceManager& tm = MAdResourceManagerSgl::instance();
    MAdStatistics&  stat = MAdStatisticsSgl ::instance();

    double t0 = tm.getTime();

    int nSplitTot = 0;
    double boundSq = 2.*mpm.getUpperLengthSqBound();

    while (1)
      {
        int nSplit = 0;
        
        EIter ei = M_edgeIter(mesh);
        pEdge edge;
        double lengthSq;
        while ( ( edge = EIter_next(ei) ) ) {
          lengthSq = SFManager->getSizeField()->SF_E_lengthSq(edge);
          if ( lengthSq > boundSq ) if (splitEdge(edge,true)) nSplit++;
        }
        EIter_delete(ei);

        nSplitTot += nSplit;
        if ( debugLevel >= 3 && !checkTheMesh() ) abort(__LINE__,__FILE__);

        if ( nSplit == 0 ) break;
      }
    
    double dt = tm.getTime() - t0;
    stat.addCPUESplits(dt);
    stat.addNumESplits(nSplitTot);
    if ( verbosity >= 2 ) {
      cout << "Performed "<< nSplitTot<<" (very long) edge splits in "<<dt<<" seconds\n";
    }

    // --- Geometry tracking ---
    if (nSplitTot) snapVertices();

    return nSplitTot;
  }

  // -------------------------------------------------------------------
  int MeshAdapter::runOneIter()
  {
    MAdResourceManager& tm = MAdResourceManagerSgl::instance();

    // --- constrain parallel interfaces ---
#ifdef PARALLEL
    { 
      double t0= tm.getTime();
      UpdateParallelConstraint(mesh);
      double dt = tm.getTime() - t0;
      MAdMsgSgl::instance().info(__LINE__,__FILE__,
                                 "Performed constraint of parallel interfaces in %f seconds",dt);
      MPI_Barrier(MPI_COMM_WORLD);
    }
#endif

    // --- perform local mesh modifications ---
    int numTopoChg = 0;

    // first, split very long edges (avoid infinite loops split/collapse)
//     numTopoChg += splitLongestEdges();
    
    if ( algorithm == CPS_SWP_SLV_SPT )
      {
        numTopoChg += eCollapseLoop(); // Collapse short edges
        numTopoChg += edgeSwapLoop();  // Swap edges to improve elements quality
        numTopoChg += removeSlivers(); // Remove slivers
        numTopoChg += eSplitLoop();    // Split long edges
      }

    if ( algorithm == SPT_SWP_SLV_CPS )
      {
        numTopoChg += eSplitLoop();    // Split long edges
        numTopoChg += edgeSwapLoop();  // Swap edges to improve elements quality
        numTopoChg += removeSlivers(); // Remove slivers
        numTopoChg += eCollapseLoop(); // Collapse short edges
      }

    if ( algorithm == SLV_CPS_SWP_SPT )
      {
        numTopoChg += removeSlivers(); // Remove slivers
        numTopoChg += eCollapseLoop(); // Collapse short edges
        numTopoChg += edgeSwapLoop();  // Swap edges to improve elements quality
        numTopoChg += eSplitLoop();    // Split long edges
      }

    // --- Reposition the vertices ---
    //LaplaceSmoothing();
//     LaplaceSmoothing(FAST);

#ifdef PARALLEL
    // --- move parallel interfaces ---
    { 
      MPI_Barrier(MPI_COMM_WORLD);
      double t0 = tm.getTime();
      DeleteParallelConstraint(mesh);
      double dt = tm.getTime() - t0;
      if( verbosity >= 2 ) {
        MAdMsgSgl::instance().info(__LINE__,__FILE__,
                                   "Performed unconstraint of parallel interface in %f seconds",dt);
      }

      t0 = tm.getTime();
      switch( load_balance_algorithm  ) {      
      case DEFAULT_ALGORITHM:
        Balance2( mesh, *dataExchanger );
        break;
      case METIS_ALGORITHM:
#ifdef _HAVE_PARMETIS_
// #warning "Temporary removed calls to load balancing with Metis (bug in it)"
//         BalanceMetis2( mesh, *dataExchanger );
#else
        MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                    "Could not balance load: Parmetis not enabled");
#endif
        break;
      default:
        MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                    "Unknown load balancing algorithm: %d",
                                    load_balance_algorithm);        
      }
      dt = tm.getTime() - t0;
      if( verbosity >= 2 ) {
        MAdMsgSgl::instance().info(__LINE__,__FILE__,
                                   "Performed load balancing in %f seconds",dt);
      }
    }
#endif

    // --- clean up the mesh by deleting dead elements ---
    M_clean(mesh);

// #warning "debug"
//     M_writeMsh(mesh,"test.msh",2);

    // --- check the mesh (debug) ---
    if ( debugLevel >= 2 && !checkTheMesh() )  abort(__LINE__,__FILE__);

    return numTopoChg;
  }

  // -------------------------------------------------------------------
  void MeshAdapter::uglyTheMesh(double avgQualThresh, int maxIt)
  {
        if (debug)
    {
    std::cout << "======================================" << std::endl;
    std::cout << "UGLYTHEMESH" << std::endl;
  }

    MeshQualityManagerSgl::instance().evaluateStatistics();
    double meanQuality  = MeshQualityManagerSgl::instance().getMeanShape();
  
    int iter = 0;
    while ( (meanQuality > avgQualThresh) && (iter < maxIt) ) {

      int numOp = 0;
      int oriNumFaces = M_numFaces(mesh);
      int count = 0;
      FIter fi = M_faceIter(mesh);
      pFace pf;
      while ( ( pf = FIter_next(fi) ) && ( count < oriNumFaces ) ) {
        int flag = 0;
        for (int j=0; j<3; j++) {
          pEdge pe = F_edge(pf,j);
          for( int ClpsOnvt = 0; ClpsOnvt < 2; ClpsOnvt++ ) {
            fCollapseOp->reset(pf, pe, ClpsOnvt);
            double worstShape;
            if ( fCollapseOp->evaluate(&worstShape) ) {
              fCollapseOp->apply();
              flag = 1;
              break;
            }
          }
          if ( flag ) {
            numOp++;
            break;
          }
        }
        count++;
      }
      FIter_delete(fi);
      cout << "Num face collapses: "<<numOp<<endl;
    
      meanQuality  = MeshQualityManagerSgl::instance().getMeanShape();
      iter++;
    }
  }

  // -------------------------------------------------------------------
  //! Removes elements with negative volume and returns \ingroup adaptation
  //! number of remaining negative volumes
  int MeshAdapter::removeNegativeElements()
  {
    int nNeg = 0;
    int nRem = 0;
    RIter rIt= M_regionIter(mesh);
    while ( pRegion region = RIter_next(rIt) ) {
      if (R_volume(region) < 0.) {
        nNeg++;
        if ( removeRegion(region) ) nRem++;
      }
    }
    RIter_delete(rIt);

    cout << "Negative elements removal: "
         << nNeg << " detected, "
         << nRem << " removed\n";

    return (nNeg - nRem);
  }

  // -------------------------------------------------------------------
  // LEVEL 4 OPERATIONS
  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  void MeshAdapter::run()
  {
    // --- Adaptation procedure ---

    MAdResourceManager& tm = MAdResourceManagerSgl::instance();
    double t0= tm.getTime();

    cout << "\n +++ Starting a MDB mesh adaptation procedure +++\n";
    if ( verbosity >= 2 ) cout<<"\n";

    for (int iter=1; iter <= maxIterationsNumber; iter++) {
    
      if ( verbosity >= 1 ) {
        cout<<"   --- Mesh adaptation iteration "<<iter<<" ---\n"<<endl;
      }
      int num = runOneIter();

//       // for debugging
//       stringstream ss;
//       string iterStr;  ss << iter;  ss >> iterStr;
//       string debugName = outPrefix + "iter" + iterStr + ".msh";
//       M_writeMsh(mesh,debugName.c_str(),2);

#ifdef PARALLEL
      int numGlob = 0;
      MPI_Allreduce(&num,&numGlob,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif

      if ( updateSFFrequency && (iter%updateSFFrequency) == 0 ) updateSizeField();

      if ( verbosity >= 2 ) cout<<"\n";
#ifdef PARALLEL
      if ( numGlob == 0 ) break;
#else
      if ( num == 0 ) break;
#endif

      if ( iter == maxIterationsNumber  && num != 0 )
        {
          MAdStatisticsSgl::instance().addInfiniteLoops(1);
          if ( verbosity >= 1 ) {
            MAdMsgSgl::instance().warning(__LINE__,__FILE__,
                                          "Infinite loop detected");
          }
          break;
        }
    }

    // --- Load balancing ---
#ifdef PARALLEL
    MPI_Barrier(MPI_COMM_WORLD);
    if( load_balance_algorithm != METIS_ALGORITHM ) {   
#ifdef _HAVE_PARMETIS_ 
// #warning "Temporary removed calls to load balancing with Metis (bug in it)"
//       BalanceMetis2( mesh, *dataExchanger );
#endif
    }
#endif
    
    if ( debugLevel >= 1 && !checkTheMesh() )  abort(__LINE__,__FILE__);
    cout << "\n +++ Ending a MDB mesh adaptation procedure (" 
         << tm.getTime() - t0 << " seconds) +++\n";
    if ( verbosity >= 2 ) cout<<"\n";
  }

#ifdef PARALLEL
  // -------------------------------------------------------------------
  void MeshAdapter::setLoadBalancingAlgorithm( loadBalanceAlgorithm lbAlgo )
  {
    load_balance_algorithm = lbAlgo;
  }
  // -------------------------------------------------------------------
//   void MeshAdapter::setDataExchanger( MDB_DataExchanger* dataExch )
//   {
//     dataExchanger = dataExch;
//   }
#endif

  // -------------------------------------------------------------------
  int MeshAdapter::partlyMoveObjects(double t, double dt, double* part)
  {
    return objects->partlyMove(*vMoveOp,t,dt,part);
  }

  // -------------------------------------------------------------------
  void MeshAdapter::moveObjectsAndReposition(double t, double dt,
                                             bool subAdaptation,
                                             double qualityThreshold,
                                             double chi, bool meshIsCavity,
                                             int cavityThickness)
  {
    objects->moveAndReposition(mesh, t, dt, subAdaptation, qualityThreshold, chi, 
                               meshIsCavity, cavityThickness, this);

    /*
    objects->setupElasticRepositioning(mesh, t, dt, qualityThreshold, chi, 
                                       meshIsCavity, cavityThickness);
    int subIter = 0;
    double ratio = 0.;
    int achieved = -1;
    while ( achieved != 2 ) {
      achieved = objects->reposition(&ratio);
      MAdMsgSgl::instance().info(-1,__FILE__,
                                 "Advanced repositioning, achieved: %d, ratio: %f",achieved,ratio);
// #warning "output paper"
//       string name = "subIter";
//       stringstream ss;
//       string iterStr;  ss << subIter;  ss >> iterStr;
//       name = name + iterStr + "_1.pos";
//       writePos(name.c_str(),OD_MEANRATIO);

      if ( debugLevel >= 3 && !checkTheMesh() )  abort(__LINE__,__FILE__);
      if ( achieved <= 1 ) {
        if ( !subAdaptation ) {
          MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                      "Could not advance objects, ratio reached: %f",ratio);
        }
        else if ( !removeSlivers() && !optimiseElementShape() ) {
          MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                      "Could not advance objects, ratio reached: %f",ratio);
        }
      }
// #warning "output paper"
//       name = "subIter";
//       stringstream ss2;
//       iterStr ="";  ss2 << subIter;  ss2 >> iterStr;
//       name = name + iterStr + "_2.pos";
//       writePos(name.c_str(),OD_MEANRATIO);

      subIter++;
    }
    */

    if ( debugLevel >= 1 && !checkTheMesh() )  abort(__LINE__,__FILE__);
  }

  // -------------------------------------------------------------------
  // -------------------------------------------------------------------
  // -------------------------------------------------------------------
  void MeshAdapter::nodalDataDiagnostics(std::ostream& out) const
  {
    NodalDataManagerSgl::instance().diagnostics(out);
  }

  // -------------------------------------------------------------------
  void MeshAdapter::openJournal() const
  {
    HistorySgl::instance().openJournal();
  }

  // -------------------------------------------------------------------
  void MeshAdapter::setReferenceJournal(std::string& name) const
  {
    HistorySgl::instance().loadJournal(name);
  }

  // -------------------------------------------------------------------
  void MeshAdapter::flushJournal(std::ostream& out) const
  {
    HistorySgl::instance().flushJournal(out);
  }

  // -------------------------------------------------------------------
  void MeshAdapter::enableSliverReports()
  {
    if ( !sliverROp ) {
      std::cerr<<"Error: no sliver region handler\n";
      exit(0);
    }

    sliverROp->enableReport(outPrefix);

    if ( !sliverFOp ) {
      std::cerr<<"Error: no sliver face handler\n";
      exit(0);
    }

    sliverFOp->enableReport(outPrefix);
  }

  // -------------------------------------------------------------------
  void MeshAdapter::testSliverOperators(bool test)
  {
    sliverROp->setTestAllOperators(test);
  }

  // -------------------------------------------------------------------
  bool MeshAdapter::checkTheMesh(int verbosity, std::ostream& out, 
                                 MeshStatus * status) const
  {
    MAdMsgSgl::instance().info(-1,__FILE__,
                               "Checking mesh validity (debug level = %d)",
                               debugLevel);
    return checkMesh(mesh,CHECK_ALL,verbosity,out,status);
  }

  // -------------------------------------------------------------------
  void MeshAdapter::infoMobileObjects(std::ostream& out) const
  {
    objects->describe(out);
  }

  // -------------------------------------------------------------------
  // -------------------------------------------------------------------
  // -------------------------------------------------------------------

  void MeshAdapter::setOutputPrefix(string prefix)
  {
    outPrefix = prefix;
  }

  // -------------------------------------------------------------------
  void MeshAdapter::writePos(string fn, MAdOutputData type) const
  {
    string fullName = outPrefix + fn;
    MAdGmshOutput(mesh, (const pSField) SFManager->getSizeField(), 
                  fullName.c_str(), type);
  }

  // -------------------------------------------------------------------
  void MeshAdapter::writeMsh(string fn) const
  {
    string fullName = outPrefix + fn;
    M_writeMsh(mesh,fullName.c_str(),2,NULL);
  }

  // -------------------------------------------------------------------
  void MeshAdapter::writeDistanceToWalls(string fnb) const
  {
    std::set<LocalSizeField*> localSF = SFManager->getLocalSizeFields();
    std::set<LocalSizeField*>::const_iterator it = localSF.begin();
    for (; it != localSF.end(); it++)
      {
        string fullName = outPrefix + fnb + "_" + (*it)->getName() + ".pos";
        MAdGmshOutput(mesh, (const pSField) *it, fullName.c_str(),
                      OD_DISTANCE);
      }
  }

  // -------------------------------------------------------------------
  void MeshAdapter::writeVolumicCurvature(string fnb) const
  {
    std::cout << __FILE__ << __LINE__ << std::endl;
    std::set<LocalSizeField*> localSF = SFManager->getLocalSizeFields();
    std::set<LocalSizeField*>::const_iterator it = localSF.begin();
    for (; it != localSF.end(); it++)
      {
        std::cout << __FILE__ << __LINE__ << std::endl;
        string fullName = outPrefix + fnb + "_" + (*it)->getName() + ".pos";
        MAdGmshOutput(mesh, (const pSField) *it, fullName.c_str(),
                      OD_CURVATURE_DIV);
      }
  }

  // -------------------------------------------------------------------
  void MeshAdapter::getModificationsInfo(int * nSplit, int * nColl, 
                                         int * nSwap, double * cpuSplit,
                                         double * cpuColl, double * cpuSwap,
                                         double * cpuSliv) const
  {
    *nSplit   = MAdStatisticsSgl::instance().getNumESplits();
    *nColl    = MAdStatisticsSgl::instance().getNumECollapses();
    *nSwap    = MAdStatisticsSgl::instance().getNumESwaps();
    *cpuSplit = MAdStatisticsSgl::instance().getCPUESplits();
    *cpuColl  = MAdStatisticsSgl::instance().getCPUECollapses();
    *cpuSwap  = MAdStatisticsSgl::instance().getCPUESwaps();
    *cpuSliv  = MAdStatisticsSgl::instance().getCPUSlivers();
  }

  // -------------------------------------------------------------------
  void MeshAdapter::getStatistics(double * meanQuality, 
                                  double * worstQuality) const
  {
    MeshQualityManagerSgl::instance().evaluateStatistics();
    *meanQuality  = MeshQualityManagerSgl::instance().getMeanShape();
    *worstQuality = MeshQualityManagerSgl::instance().getWorstShape();
  }

  // -------------------------------------------------------------------
  void MeshAdapter::printStatistics(std::ostream& out) const
  {
    MeshQualityManagerSgl::instance().evaluateStatistics();
    MeshQualityManagerSgl::instance().printStatistics(out);
    out << "\n\n";
    MAdStatisticsSgl::instance().print(out);
    out << "\n\n";
    if (geoTracker) geoTracker->printFailures(out);
  }

  // -------------------------------------------------------------------
  void MeshAdapter::printSliverRegionStatistics(std::ostream& out) const 
  {
    if ( sliverROp ) sliverROp->printStats(out);
    //   if ( sliverROp ) sliverROp->printOperatorsTest(out);
  }
  // -------------------------------------------------------------------
  void MeshAdapter::abort(int line, const char* file) const
  {
    cerr << "\n"
         << " ****************************\n"
         << "  MeshAdapter is aborting... \n"
         << " ****************************\n\n";    

    string stName = outPrefix + "aborted_statistics";
    std::ofstream stOut(stName.c_str());
    printStatistics(stOut);

    string slName = outPrefix + "aborted_slivers";
    std::ofstream slOut(slName.c_str());
    printSliverRegionStatistics(slOut);

    string jName = outPrefix + "aborted_journal";
    std::ofstream jOut(jName.c_str());
    flushJournal(jOut);

    writePos("aborted.pos",OD_MEANRATIO);
    writeMsh("aborted.msh");

    if ( file && line >=0 ) {
      MAdMsgSgl::instance().info(-1,__FILE__,
                                 "Abort function called from file %s (line %d)",
                                 file,line);
    }
    else {
      MAdMsgSgl::instance().info(-1,__FILE__,"Abort!");
    }

    flush(cout);
    flush(cerr);

    exit(EXIT_FAILURE);
  }

  // -------------------------------------------------------------------

} // End of namespace MAd

