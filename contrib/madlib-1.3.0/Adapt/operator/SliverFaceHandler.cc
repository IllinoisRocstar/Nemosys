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

#include "SliverFaceHandler.h"
#include "EdgeSwapOp.h"
#include "EdgeCollapseOp.h"
#include "FaceCollapseOp.h"
#include "MeshParametersManager.h"
#include "MAdOutput.h"


// standard C/C++
#include <set>
#include <stdio.h>
using std::set;
#include <sstream>
using std::stringstream;
using std::string;

namespace MAd {

  // -------------------------------------------------------------------
  void sliverFaceHandler::removeSliverFaces(int* nbSliverIn, int* nbSliverOut)
  {
    // --- list all slivers ---

    set<pFace> slivers;
    if (nbSliverIn) (*nbSliverIn) = 0;
    FIter fIter = M_faceIter(mesh);
    while( pFace face = FIter_next(fIter) ) {
      if ( F_isSliver(face) ) {
        slivers.insert(face);
        if (nbSliverIn) (*nbSliverIn)++;
      }
    }
    FIter_delete(fIter);

    // --- iterate on slivers elimination attempts ---

    int iter = 0;
    int maxIter = 5;
    bool eliminated = true;
    *nbSliverOut = *nbSliverIn;
    while ( ( eliminated       ) && 
            ( *nbSliverOut > 0 ) && 
            ( iter < maxIter   ) )  {

      eliminated = false;

      set<pFace> oldSlivers;

      // --- process the slivers list ---
      set<pFace>::iterator sIter = slivers.begin();
      set<pFace>::iterator sLast = slivers.end();
      for (; sIter != sLast; sIter++) {

        pFace sliver = *sIter;

        // --- check if the sliver still exists ---
        if ( oldSlivers.find(sliver) != oldSlivers.end() ) continue;

        // --- check if the sliver is still a sliver ---
        if ( !F_isSliver(sliver) ) {
          oldSlivers.insert(sliver);
          continue;
        }

        // --- find an appropriate modification ---
        pMAdOperator modif = NULL;
        int result = findOperation(sliver,&modif);

        if ( result <= 0 && reportFailures ) reportSliver(sliver);
      
        if ( result >= 1 && modif ) {

          // --- prepare to remove deleted faces from sliver list ---
          pPList cavity;
          modif->getCavity(&cavity);
          void * tmp = 0;
          while( pFace pf = (pFace)PList_next(cavity,&tmp) ) {
            oldSlivers.insert(pf);
          }
          PList_delete(cavity);
        
          // --- eliminate the sliver ---
          modif->apply();
          eliminated = true;
          oldSlivers.insert(sliver);
        }

        if (modif) delete modif;
      }

      // --- count remaining and list slivers ---
      slivers.clear();
      if (nbSliverOut) (*nbSliverOut) = 0;
      fIter = M_faceIter(mesh);
      while( pFace face = FIter_next(fIter) ) {
        if ( F_isSliver(face) ) {
          slivers.insert(face);
          if (nbSliverOut) (*nbSliverOut)++;
        }
      }
      FIter_delete(fIter);

      iter++;
    }
  }

  // -------------------------------------------------------------------
  //  Try to find a modification that eliminates the sliver triangle.
  //  Return:
  //    0: no modification found
  //    1: a modification is found but with a resulting sliver (with better shape)
  //    2: a modification is found with no remaining sliver
  int sliverFaceHandler::findOperation(pFace face, 
                                       pMAdOperator* modification)
  {  
    MeshParametersManager& mpm = MeshParametersManagerSgl::instance();

    *modification = NULL;

    double bestShape = 0.;
    mqm.getShape(face,0,&bestShape);

    if ( bestShape >= mpm.getSliverTriBound() ) {
      printf ("Warning: processing a face which is not a sliver in findOperation\n");
      return 0;
    }

    // check edge swaps
    // -----------------

    edgeSwapOp eSwap(mesh,sizeField);
//     if( collapseOnBdry )  eSwap.setCheckVolume(0.,dATol);
  
    for (int i=0; i<3; i++) {
  
      eSwap.setSwapEdge(F_edge(face,i));
    
      double worst;
      if ( eSwap.evaluate(&worst) && worst > bestShape ) {
        bestShape = worst;
        if (*modification) delete *modification;
        *modification = new edgeSwapOp(eSwap);
        if (bestShape > mpm.getSliverTriBound()) return 2;
      }
    }

    // check edge collapses
    // ---------------------

    edgeCollapseOp eCollapse(mesh,sizeField);
    eCollapse.collapseOnBoundary(collapseOnBdry,dATol);

    for(int i=0; i<3; i++) {
    
      pEdge edge = F_edge(face,i);
    
      for (int j=0; j<2; j++) {
      
        eCollapse.setCollapseEdge(edge,E_vertex(edge,j),E_vertex(edge,(j+1)%2));
 
        double worst;
        if (eCollapse.evaluate(&worst) && worst > bestShape ) {
          bestShape = worst;
          if (*modification) delete *modification;
          *modification = new edgeCollapseOp(eCollapse);
          if (bestShape > mpm.getSliverTriBound()) return 2;
          //         if (bestShape > mpm.getSliverTriBound()) { printf("Applied edge collapse\n"); return 2;}
        }
      }
    }

    // check face collapse
    // --------------------

    if ( addNodes )
      {
        faceCollapseOp fCollapse(mesh,sizeField);
        fCollapse.collapseOnBoundary(collapseOnBdry,dATol);
        
        for(int iE=0; iE<3; iE++) {
          
          pEdge edge = F_edge(face,iE);
          
          for (int iDir=0; iDir<2; iDir++) {
            
            fCollapse.reset(face,edge,(bool)iDir);
            
            double worst;
            if (fCollapse.evaluate(&worst) && worst > bestShape ) {
              bestShape = worst;
              if (*modification) delete *modification;
              *modification = new faceCollapseOp(fCollapse);
              if (bestShape > mpm.getSliverTriBound()) return 2;
              //         if (bestShape > mpm.getSliverTriBound()) { printf("Applied face collapse\n"); return 2;}
            }
          }
        }
      }

    if( *modification )  return 1;

    return 0;
  }

  // -------------------------------------------------------------------
  bool sliverFaceHandler::F_isSliver(pFace face)
  {
    MeshParametersManager& mpm = MeshParametersManagerSgl::instance();

    double shape;
    mqm.getShape(face,0,&shape);

    if( shape < mpm.getSliverTriBound() ) 
      return true;

    return false;
  }

  // -------------------------------------------------------------------
  bool sliverFaceHandler::F_isSliver(pFace face, double * shape)
  {
    MeshParametersManager& mpm = MeshParametersManagerSgl::instance();

    mqm.getShape(face,0,shape);

    if( *shape < mpm.getSliverTriBound() ) 
      return true;

    return false;
  }

  // -------------------------------------------------------------------
  void sliverFaceHandler::reportSliver(pFace face)
  {
    stringstream ss;
    string idStr;  ss << reportId;  ss >> idStr;

    string name = reportPrefix + "sliver" + idStr + ".pos";
  
    pPList cavity = PList_new();

    PList_append(cavity,(pEntity)face);
  
    for (int iV=0; iV<F_numVertices(face); iV++) {

      pVertex vertex = F_vertex(face, iV);
    
      pPList vFaces = V_faces(vertex);
      void * temp = 0;
      while ( pFace pf = (pFace)PList_next(vFaces,&temp) ) {
        PList_appUnique(cavity,(pEntity)pf);
      }
      PList_delete(vFaces);
    }

    printPosEntities(cavity,name.c_str(),OD_MEANRATIO,sizeField);

    PList_delete(cavity);

    reportId++;
  }

  // -------------------------------------------------------------------

}
