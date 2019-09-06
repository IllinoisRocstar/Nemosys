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

#include "SliverRegionHandler.h"
#include "EdgeSplitOp.h"
#include "EdgeCollapseOp.h"
#include "FaceCollapseOp.h"
#include "DESCOp.h"
#include "EdgeSwapOp.h"
#include "FaceSwapOp.h"
#include "VertexMoveOp.h"
#include "MeshParametersManager.h"
#include "MAdOutput.h"
#include "MathUtils.h"
#include "DistanceFunction.h"

#include <sstream>
using std::stringstream;
using std::set;
using std::ostream;
using std::string;
using std::endl;

namespace MAd {

  // -------------------------------------------------------------------
  void sliverRegionHandler::removeSliverRegions(int* nbSliverIn, int* nbSliverOut)
  {
    // --- list all slivers ---

    set<pRegion> slivers;
    if (nbSliverIn) (*nbSliverIn) = 0;
    RIter rIter = M_regionIter(mesh);
    while( pRegion region = RIter_next(rIter) ) {
      if ( R_isSliver(region) ) {
        slivers.insert(region);
        if (nbSliverIn) (*nbSliverIn)++;
      }
    }
    RIter_delete(rIter);

    // --- iterate on slivers elimination attempts ---

    int iter = 0;
    int maxIter = 5;
    bool eliminated = true;
    *nbSliverOut = *nbSliverIn;
    while ( ( eliminated       ) && 
            ( *nbSliverOut > 0 ) && 
            ( iter < maxIter   ) )  {

      eliminated = false;

      set<pRegion> oldSlivers;

      // --- process the slivers list ---
      set<pRegion>::iterator sIter = slivers.begin();
      set<pRegion>::iterator sLast = slivers.end();
      for (; sIter != sLast; sIter++) {

        pRegion sliver = *sIter;

        // --- check if the sliver still exists ---
        if ( oldSlivers.find(sliver) != oldSlivers.end() ) continue;

        // --- check if the sliver is still a sliver ---
        if ( !R_isSliver(sliver) ) {
          oldSlivers.insert(sliver);
          continue;
        }

        // --- find an appropriate modification ---
        pMAdOperator modif = NULL;
        int result = findOperation(sliver,&modif);

        if ( testOperators ) testOperations(sliver);

        if ( result <= 0 && reportFailures ) reportSliver(sliver);
      
        if ( result >= 1 ) {

          assert(modif);

          // --- prepare to remove deleted regions from sliver list ---
          pPList affReg;
          modif->getCavity(&affReg);
          void * tmp = 0;
          while( pRegion pr = (pRegion)PList_next(affReg,&tmp) ) {
            oldSlivers.insert(pr);
          }
          PList_delete(affReg);
        
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
      rIter = M_regionIter(mesh);
      while( pRegion region = RIter_next(rIter) ) {
        if ( R_isSliver(region) ) {
          slivers.insert(region);
          if (nbSliverOut) (*nbSliverOut)++;
        }
      }
      RIter_delete(rIter);

      iter++;
    }
  }

  // -------------------------------------------------------------------
  //  Try to find an operation that eliminates the sliver.
  //  Returns:
  //    0: no modification found
  //    1: a modification is found but with a resulting sliver (with better shape)
  //    2: a modification is found with no remaining sliver
  int sliverRegionHandler::findOperation(pRegion region, 
                                         pMAdOperator* operation)
  {  
    MeshParametersManager& mpm = MeshParametersManagerSgl::instance();

    *operation = NULL;

    double bestShape = 0.;
    mqm.getShape(region,&bestShape);

    if ( bestShape >= mpm.getSliverTetBound() ) {
      printf ("Warning: processing a tet which is not a sliver in handleSliverTet\n");
      return 0;
    }

    edgeCollapseOp eCollOp (mesh,sizeField);
    eCollOp.collapseOnBoundary(collapseOnBdry,dVTolClps);

    faceCollapseOp fCollOp (mesh,sizeField);
    fCollOp.collapseOnBoundary(collapseOnBdry,dVTolClps);

    edgeSwapOp     eSwapOp (mesh,sizeField);
    eSwapOp.swapOnBoundary(swapOnBdry,dVTolSwap);

    DESCOp         descOp  (mesh,sizeField);
    edgeSplitOp    eSplitOp(mesh,sizeField);
    faceSwapOp     fSwapOp (mesh,sizeField);
    vertexMoveOp   vMoveOp (mesh,sizeField,false);

    pEntity keyEnts[4];
    int sType = getSliverType(region, keyEnts);
    if( sType == 0 ) return 0;

    double worst;

    // --- TYPE I slivers ---
    if( sType == 1 ) 
      {
#ifdef MAKESTATS
        numTypeI_In++;
#endif

        // check split 
        // ------------
        if ( addNodes ) {
          for(int i=0; i<2; i++) {
#ifdef MAKESTATS
            numSplitTested_I++;
#endif
            eSplitOp.setSplitEdge((pEdge)keyEnts[i]);
        
            if( E_whatInType((pEdge)keyEnts[i]) != 3 ) continue;

            // determine the position of the new vertex
//             double newPos[3];
//             E_cavityCenter( (pEdge)keyEnts[i], newPos );
//             double t = E_linearParams((pEdge)keyEnts[i],newPos);
//             if( t>1. || t<0. ) continue;
//             eSplitOp.setSplitPos(newPos);

            worst = -1.;
            if( eSplitOp.evaluate(&worst) && worst > bestShape &&
                eSplitOp.getMinLenSq() > mpm.getSliverLowerLengthSqBound() &&
                eSplitOp.getMaxLenSq() < mpm.getSliverUpperLengthSqBound() ) {
              bestShape = worst;
              if(*operation) delete *operation;
              *operation = new edgeSplitOp(eSplitOp);
#ifdef MAKESTATS
              numSplitAvailable_I++;
              if (bestShape < mpm.getSliverTetBound()) numSplitToSliver_I++;
#endif
              if (bestShape > mpm.getSliverTetBound()) break;
            }
          }
          if( *operation && bestShape > mpm.getSliverTetBound() ) {
#ifdef MAKESTATS
            numSplitOperated_I++;
#endif
            return 2;
          }
        }

        // check edge collapses
        // ---------------------
        for(int i=0; i<6; i++) {
          pEdge edge = R_edge(region,i);
          for (int j=0; j<2; j++) {
#ifdef MAKESTATS
            numEClpsTested_I++;
#endif
            eCollOp.setCollapseEdge(edge,E_vertex(edge,j),E_vertex(edge,1-j));

            worst = -1.;
            if ( eCollOp.evaluate(&worst) && worst > bestShape &&
                 eCollOp.getMinLenSq() > mpm.getSliverLowerLengthSqBound() &&
                 eCollOp.getMaxLenSq() < mpm.getSliverUpperLengthSqBound() ) {
              bestShape = worst;
              if (*operation) delete *operation;
              *operation = new edgeCollapseOp(eCollOp);
#ifdef MAKESTATS
              numEClpsAvailable_I++;
              if (worst < mpm.getSliverTetBound()) numEClpsToSliver_I++;
#endif 
              if (bestShape > mpm.getSliverTetBound()) break;
            }
          }
          if (bestShape > mpm.getSliverTetBound()) break;
        }
        if( *operation && bestShape > mpm.getSliverTetBound() ) {
#ifdef MAKESTATS
          numEClpsOperated_I++;
#endif
          return 2;
        }

        // check double edge split collapse
        // ---------------------------------
        if ( addNodes ) {
          for(int i=0; i<2; i++) {
#ifdef MAKESTATS
            numDESCTested_I++;
#endif
            descOp.setDESC(region,(pEdge)keyEnts[i%2],(pEdge)keyEnts[(i+1)%2]);

            worst = -1.;
            if( descOp.evaluate(&worst) && worst > bestShape &&
                descOp.getMinLenSq() > mpm.getSliverLowerLengthSqBound() &&
                descOp.getMaxLenSq() < mpm.getSliverUpperLengthSqBound() ) {
              bestShape = worst;
              if(*operation) delete *operation;
              *operation = new DESCOp(descOp);
#ifdef MAKESTATS
              numDESCAvailable_I++;
              if (worst < mpm.getSliverTetBound()) numDESCToSliver_I++;
#endif
            }
            if( *operation && bestShape > mpm.getSliverTetBound() ) {
#ifdef MAKESTATS
              numDESCOperated_I++;
#endif
              return 2;          
            }
          }
        }

        // check face collapse
        // --------------------
        if ( addNodes ) {
          for(int i=0; i<2; i++) {
            for( int ClpsOnvt=1; ClpsOnvt>=0; ClpsOnvt-- ) {
#ifdef MAKESTATS
              numFClpsTested_I++;
#endif
              fCollOp.reset((pFace)keyEnts[i+2],(pEdge)keyEnts[i], ClpsOnvt);
  
              worst = -1.;
              if( fCollOp.evaluate(&worst) && worst > bestShape &&
                  fCollOp.getMinLenSq() > mpm.getSliverLowerLengthSqBound() &&
                  fCollOp.getMaxLenSq() < mpm.getSliverUpperLengthSqBound() ) {
                bestShape = worst;
                if (*operation) delete *operation;
                *operation = new faceCollapseOp(fCollOp); 
#ifdef MAKESTATS
                numFClpsAvailable_I++;
                if (worst < mpm.getSliverTetBound()) {
                  numFClpsToSliver_I++;
                }
#endif
                if( *operation && bestShape > mpm.getSliverTetBound() ) {
#ifdef MAKESTATS
                  numFClpsOperated_I++;
#endif
                  return 2;
                }
              }
            }
          }
        }
      
        // check edge swaps
        // -----------------
        for(int i=0; i<2; i++) {
#ifdef MAKESTATS
          numSwapTested_I++;
#endif
          eSwapOp.setSwapEdge((pEdge)keyEnts[i]);

          worst = -1.;
          if( eSwapOp.evaluate(&worst) && worst > bestShape &&
              eSwapOp.getMinLenSq() > mpm.getSliverLowerLengthSqBound() &&
              eSwapOp.getMaxLenSq() < mpm.getSliverUpperLengthSqBound() ) {
            bestShape = worst;
            if (*operation) delete *operation;
            *operation = new edgeSwapOp(eSwapOp); 
#ifdef MAKESTATS
            numSwapAvailable_I++;
            if (worst < mpm.getSliverTetBound()) numSwapToSliver_I++;
#endif
            if (bestShape > mpm.getSliverTetBound()) break;
          }
        }
        if( *operation && bestShape > mpm.getSliverTetBound() ) {
#ifdef MAKESTATS
          numSwapOperated_I++;
#endif
          return 2;
        }

        // check vertex move 
        // ------------------
        pPList verts = R_vertices(region);
        void * temp = NULL;
        while ( pVertex pv = (pVertex)PList_next(verts,&temp) ) {
#ifdef MAKESTATS
          numVMoveTested_I++;
#endif
          if( V_whatInType(pv) != 3 ) continue;
        
          vMoveOp.setPositionToOptimal(pv);

          worst = -1.;
          if( vMoveOp.evaluate(&worst) && worst > bestShape &&
              vMoveOp.getMinLenSq() > mpm.getSliverLowerLengthSqBound() &&
              vMoveOp.getMaxLenSq() < mpm.getSliverUpperLengthSqBound() ) {
            bestShape = worst;
            if(*operation) delete *operation;
            *operation = new vertexMoveOp(vMoveOp);
#ifdef MAKESTATS
            numVMoveAvailable_I++;
            if (bestShape < mpm.getSliverTetBound()) numVMoveToSliver_I++;
#endif
          }
        }
        PList_delete(verts);
        if( *operation && bestShape > mpm.getSliverTetBound() ) {
#ifdef MAKESTATS
          numVMoveOperated_I++;
#endif
          return 2;
        }
        

#ifdef MAKESTATS
        if( *operation ) {
          switch ((*operation)->type()) {
          case MAd_ECOLLAPSE:  { numEClpsOperated_I++; break; }
          case MAd_ESWAP:      { numSwapOperated_I++;  break; }
          case MAd_FCOLLAPSE:  { numFClpsOperated_I++; break; }
          case MAd_DESPLTCLPS: { numDESCOperated_I++;  break; }
          case MAd_ESPLIT:     { numSplitOperated_I++; break; }
          case MAd_VERTEXMOVE: { numVMoveOperated_I++; break; }
          default: throw;
          }
        }
#endif

        if( *operation )  return 1;

#ifdef MAKESTATS
        set<pRegion>::const_iterator sIter = sliversOut.find(region);
        if (sIter == sliversOut.end() ) {
          sliversOut.insert(region);
          numTypeI_Out++;
        }
#endif
      }

    // --- TYPE II slivers ---
    else 
      {
#ifdef MAKESTATS
        numTypeII_In++;
#endif

        // check edge collapses
        // ---------------------
        for(int i=0; i<6; i++) {
#ifdef MAKESTATS
          numEClpsTested_II++;
#endif
          pEdge edge = R_edge(region,i);
          for (int j=0; j<2; j++) {
            eCollOp.setCollapseEdge(edge,E_vertex(edge,j),E_vertex(edge,(j+1)%2));

            worst = -1.;
            if( eCollOp.evaluate(&worst) && worst > bestShape &&
                eCollOp.getMinLenSq() > mpm.getSliverLowerLengthSqBound() &&
                eCollOp.getMaxLenSq() < mpm.getSliverUpperLengthSqBound() ) {
              bestShape = worst;
              if (*operation) delete *operation;
              *operation = new edgeCollapseOp(eCollOp);
#ifdef MAKESTATS
              numEClpsAvailable_II++;
              if (worst < mpm.getSliverTetBound()) numEClpsToSliver_II++;
#endif 
              if (bestShape > mpm.getSliverTetBound()) break;
            }
          }
          if (bestShape > mpm.getSliverTetBound()) break;
        }
        if( *operation && bestShape > mpm.getSliverTetBound() ) {
#ifdef MAKESTATS
          numEClpsOperated_II++;
#endif
          return 2;
        }

        // check face collapse
        // --------------------
        if ( addNodes ) {
          for( int ClpsOnvt=1; ClpsOnvt>=0; ClpsOnvt-- ) {
#ifdef MAKESTATS
            numFClpsTested_II++;
#endif
            fCollOp.reset((pFace)keyEnts[2],(pEdge)keyEnts[1], ClpsOnvt);
  
            worst = -1.;
            if( fCollOp.evaluate(&worst) && worst > bestShape &&
                fCollOp.getMinLenSq() > mpm.getSliverLowerLengthSqBound() &&
                fCollOp.getMaxLenSq() < mpm.getSliverUpperLengthSqBound() ) {
              bestShape = worst;
              if (*operation) delete *operation;
              *operation = new faceCollapseOp(fCollOp);
#ifdef MAKESTATS
              numFClpsAvailable_II++;
              if (bestShape < mpm.getSliverTetBound()) numFClpsToSliver_II++;
#endif
            }
            if( *operation && bestShape > mpm.getSliverTetBound() ) {
#ifdef MAKESTATS
              numFClpsOperated_II++;
#endif
              return 2;
            }
          }
        }

        // check edge swaps
        // -----------------
        for(int i=0; i<3; i++) {
#ifdef MAKESTATS
          numSwapTested_II++;
#endif
          pEdge edge = F_edge((pFace)keyEnts[0],i);
          eSwapOp.setSwapEdge(edge); 

          worst=-1.;
          if( eSwapOp.evaluate(&worst) && worst > bestShape &&
              eSwapOp.getMinLenSq() > mpm.getSliverLowerLengthSqBound() &&
              eSwapOp.getMaxLenSq() < mpm.getSliverUpperLengthSqBound() ) {
            bestShape = worst;
            if (*operation) delete *operation;
            *operation = new edgeSwapOp(eSwapOp); 
#ifdef MAKESTATS
            numSwapAvailable_II++;
            if (worst < mpm.getSliverTetBound()) numSwapToSliver_II++;
#endif
            if (bestShape > mpm.getSliverTetBound()) break;
          }
        }
        if( *operation && bestShape > mpm.getSliverTetBound() ) {
#ifdef MAKESTATS
          numSwapOperated_II++;
#endif
          return 2;
        }
      
        // check face swap on the key face
        // --------------------------------
#ifdef MAKESTATS
        numFSwapTested_II++;
#endif
        fSwapOp.setSwapFace((pFace)keyEnts[0]);

        worst = -1.;
        if( fSwapOp.evaluate(&worst) && worst > bestShape &&
            fSwapOp.getMinLenSq() > mpm.getSliverLowerLengthSqBound() &&
            fSwapOp.getMaxLenSq() < mpm.getSliverUpperLengthSqBound() ) {
          bestShape = worst;
          if(*operation) delete *operation;
          *operation = new faceSwapOp(fSwapOp);
#ifdef MAKESTATS
          numFSwapAvailable_II++;
          if (bestShape < mpm.getSliverTetBound()) numFSwapToSliver_II++;
#endif
        } 
        if( *operation && bestShape > mpm.getSliverTetBound() ) {
#ifdef MAKESTATS
          numFSwapOperated_II++;
#endif
          return 2;
        }

        // check vertex move 
        // ------------------
        pPList verts = R_vertices(region);
        void * temp = NULL;
        while ( pVertex pv = (pVertex)PList_next(verts,&temp) ) {
#ifdef MAKESTATS
          numVMoveTested_II++;
#endif
          if( V_whatInType(pv) != 3 ) continue;
        
          vMoveOp.setPositionToOptimal(pv);

          worst = -1.;
          if( vMoveOp.evaluate(&worst) && worst > bestShape &&
              vMoveOp.getMinLenSq() > mpm.getSliverLowerLengthSqBound() &&
              vMoveOp.getMaxLenSq() < mpm.getSliverUpperLengthSqBound() ) {
            bestShape = worst;
            if(*operation) delete *operation;
            *operation = new vertexMoveOp(vMoveOp);
#ifdef MAKESTATS
            numVMoveAvailable_II++;
            if (bestShape < mpm.getSliverTetBound()) numVMoveToSliver_II++;
#endif
          }
        }
        PList_delete(verts);
        if( *operation && bestShape > mpm.getSliverTetBound() ) {
#ifdef MAKESTATS
          numVMoveOperated_II++;
#endif
          return 2;
        }

#ifdef MAKESTATS
        if( *operation ) {
          switch ((*operation)->type()) {
          case MAd_ECOLLAPSE:  { numEClpsOperated_II++;    break; }
          case MAd_ESWAP:      { numSwapOperated_II++;     break; }
          case MAd_FSWAP:      { numFSwapOperated_II++;    break; }
          case MAd_FCOLLAPSE:  { numFClpsOperated_II++;    break; }
          case MAd_VERTEXMOVE: { numVMoveOperated_II++;    break; }
          default: throw;
          }
        }
#endif

        if( *operation )  return 1;

#ifdef MAKESTATS
        set<pRegion>::const_iterator sIter = sliversOut.find(region);
        if (sIter == sliversOut.end() ) {
          sliversOut.insert(region);
          numTypeII_Out++;
        }
#endif
      }
  
    return 0;
  }

  // -------------------------------------------------------------------
  bool sliverRegionHandler::R_isSliver(pRegion region, double* shape)
  {
    MeshParametersManager& mpm = MeshParametersManagerSgl::instance();

    mqm.getShape(region,shape);

    if( *shape < mpm.getSliverTetBound() ) return true;
    return false;
  }

  // -------------------------------------------------------------------
  bool sliverRegionHandler::R_isSliver(pRegion region)
  {
    MeshParametersManager& mpm = MeshParametersManagerSgl::instance();

    double shape;
    mqm.getShape(region,&shape);

    if( shape < mpm.getSliverTetBound() ) return true;
    return false;
  }

  // -------------------------------------------------------------------
  // Gets the type of sliver and the related key entities.
  // Returns 0 : failure
  //         1 : type I
  //         2 : type II
  int sliverRegionHandler::getSliverType(const pRegion region,
                                         pEntity keyEnts[4])
  {
    pFace   base = R_face(region,0);
    pVertex oppV = R_fcOpVt(region,base);

    pEdge edges[6];    // the three ordered edges of the tet with base to be base
    pVertex verts[3];  // the three ordered vertices of the face
    double fxyz[3][3]; // coordinates of the three vertices of the face
    double pxyz[3];    // coordinates of the opposite vertex
    pFace faces[4];    // ordered faces of the tet
    double area[4];    // areas of the faces
  
    for(int i=0; i<3; i++) {
      edges[i] = F_edge(base,i);
      if(F_edgeDir(base,i)) verts[i] = E_vertex(edges[i],0);
      else                  verts[i] = E_vertex(edges[i],1);
      V_coord(verts[i],fxyz[i]);
    }
  
    faces[0] = base;
    for(int i=0; i<4; i++) {
      pFace face = R_face(region,i);
      if( face == base ) continue;
      if( F_inClosure(face,(pEntity)edges[0]) ) { faces[1]=face; continue; }
      if( F_inClosure(face,(pEntity)edges[1]) ) { faces[2]=face; continue; }
      if( F_inClosure(face,(pEntity)edges[2]) ) faces[3]=face;
    }
    
    for(int i=1; i<3; i++) {
      for(int j=0; j<3; j++) {
        pEdge edge = F_edge(faces[i],j);
        if( edge == edges[i-1] ) continue;
        if( E_inClosure(edge,(pEntity)verts[0]) ) { edges[3]=edge; continue; }
        if( E_inClosure(edge,(pEntity)verts[1]) ) { edges[4]=edge; continue; }
        if( E_inClosure(edge,(pEntity)verts[2]) ) edges[5]=edge;
      }
    }
    
    V_coord(oppV,pxyz);
    double proj[3];
    int bit;
    pointToTriangle(fxyz,pxyz,proj,&bit,area);

    //            010=2   | 011=3  /  001=1
    //                    |       /
    //        ------------+--e2--+-----------
    //                  v0|     /v2
    //                    | 7  /
    //            110=6   e0  e1
    //                    |  /
    //                    | /    101=5
    //                    |/
    //                  v1+
    //                   /|
    //                  / |
    //                   4
    

    // determine the key entities
    switch( bit ) {
    case 1:{
      int Emap[] = {0,4,3};
      int Fmap[] = {0,2,3};
      keyEnts[0] = (pEntity)faces[1];
      int index = indexOfMin(area[0],area[2],area[3]);
      keyEnts[1] = (pEntity)edges[Emap[index]];
      keyEnts[2] = (pEntity)faces[Fmap[index]];
      break;
    }
    case 2: {
      int Emap[] = {1,4,5};
      int Fmap[] = {0,1,3};
      keyEnts[0] = (pEntity)faces[2];
      int index = indexOfMin(area[0],area[1],area[3]);
      keyEnts[1] = (pEntity)edges[Emap[index]];
      keyEnts[2] = (pEntity)faces[Fmap[index]];
      break;   
    }
    case 3: {
      keyEnts[0] = (pEntity)edges[2];
      keyEnts[1] = (pEntity)edges[4];
      keyEnts[2] = (pEntity)( (area[0]<area[3]) ? faces[0]:faces[3] );
      keyEnts[3] = (pEntity)( (area[1]<area[2]) ? faces[1]:faces[2] );
      break; 
    }   
    case 4: {
      int Emap[] = {2,3,5};
      int Fmap[] = {0,1,2};
      keyEnts[0] = (pEntity)faces[3];
      int index = indexOfMin(area[0],area[1],area[2]);
      keyEnts[1] = (pEntity)edges[Emap[index]];
      keyEnts[2] = (pEntity)faces[Fmap[index]];
      break; 
    }  
    case 5: {
      keyEnts[0] = (pEntity)edges[1];
      keyEnts[1] = (pEntity)edges[3];
      keyEnts[2] = (pEntity)( (area[0]<area[2]) ? faces[0]:faces[2] );
      keyEnts[3] = (pEntity)( (area[1]<area[3]) ? faces[1]:faces[3] );
      break;
    }
    case 6: {
      keyEnts[0] = (pEntity)edges[0];
      keyEnts[1] = (pEntity)edges[5];
      keyEnts[2] = (pEntity)( (area[0]<area[1]) ? faces[0]:faces[1] );
      keyEnts[3] = (pEntity)( (area[2]<area[3]) ? faces[2]:faces[3] );
      break;
    }
    case 7: {
      int Emap[] = {0,1,2};
      int Fmap[] = {1,2,3};
      keyEnts[0] = (pEntity)faces[0];
      int index = indexOfMin(area[1],area[2],area[3]);
      keyEnts[1] = (pEntity)edges[Emap[index]];
      keyEnts[2] = (pEntity)faces[Fmap[index]];
      break; 
    }
    case 0:
    case 8:
    case 9:
    case 10: {
      break;
    }
    default:
      MAdMsgSgl::instance().error(__LINE__,__FILE__,"not a valid value %d",bit);
    }
    
    if ( bit==0 || bit > 7 ) return 0;
    if ( bit==3 || bit==5 || bit==6 ) return 1;
    return 2;
  }

  // -------------------------------------------------------------------
  void sliverRegionHandler::printStats(ostream& out) const {

#ifdef MAKESTATS
    out << "\n*** Statistics about the sliverRegionHandler ***\n\n";
    out << "Slivers considered: " << numTypeI_In+numTypeII_In<<endl;
    out << "  - Type I: \t"<< numTypeI_In  << "  ( " << numTypeI_Out  << " not resolved )\n";
    out << "  - Type II:\t"<< numTypeII_In << "  ( " << numTypeII_Out << " not resolved )\n\n";

    out << " \n--- Type I slivers ---\n\n";
    out << "\t     nbTestedOp\t nbAvailableOp \t nbOperated \t nbLeadToASliver\n";
    out << "ESplit:\t\t"<<numSplitTested_I << "\t\t" << numSplitAvailable_I<< "\t\t" << numSplitOperated_I<< "\t\t" <<numSplitToSliver_I << "\n";
    out << "EClps:\t\t"<< numEClpsTested_I << "\t\t" << numEClpsAvailable_I << "\t\t" << numEClpsOperated_I<< "\t\t" <<numEClpsToSliver_I << "\n";
    out << "DESpltClps:\t"<<numDESCTested_I << "\t\t" << numDESCAvailable_I << "\t\t" <<numDESCOperated_I << "\t\t" <<numDESCToSliver_I << "\n";
    out << "FClps:\t\t"<<numFClpsTested_I << "\t\t" << numFClpsAvailable_I<< "\t\t" << numFClpsOperated_I<< "\t\t" <<numFClpsToSliver_I << "\n";
    out << "ESwap:\t\t"<< numSwapTested_I << "\t\t" << numSwapAvailable_I << "\t\t" << numSwapOperated_I<< "\t\t" <<numSwapToSliver_I << "\n";
    out << "VMove:\t\t"<<numVMoveTested_I << "\t\t" << numVMoveAvailable_I<< "\t\t" << numVMoveOperated_I<< "\t\t" <<numVMoveToSliver_I << "\n";

    out << " \n--- Type II slivers ---\n\n";
    out << "\t     nbTestedOp\t nbAvailableOp \t nbOperated \t nbLeadToASliver\n";
    out << "EClps:\t\t"<<numEClpsTested_II << "\t\t" << numEClpsAvailable_II<< "\t\t" << numEClpsOperated_II<< "\t\t" <<numEClpsToSliver_II << "\n";
    out << "FClps:\t\t"<<numFClpsTested_II << "\t\t" << numFClpsAvailable_II << "\t\t" <<numFClpsOperated_II << "\t\t" << numFClpsToSliver_II<< "\n";
    out << "ESwap:\t\t"<<numSwapTested_II << "\t\t" << numSwapAvailable_II<< "\t\t" << numSwapOperated_II<< "\t\t" <<numSwapToSliver_II << "\n";
    out << "FSwap:\t\t"<<numFSwapTested_II << "\t\t" << numFSwapAvailable_II << "\t\t" <<numFSwapOperated_II << "\t\t" << numFSwapToSliver_II<< "\n";
    out << "VMove:\t\t"<<numVMoveTested_II << "\t\t" << numVMoveAvailable_II<< "\t\t" << numVMoveOperated_II<< "\t\t" <<numVMoveToSliver_II << "\n";

    out << "\n\n";

#else

    cout << "Could not write slivers statistics because SliverRegionHandler was compiled without the flag 'MAKESTATS'\n";

#endif

  }


  // -------------------------------------------------------------------
  void sliverRegionHandler::testOperations(pRegion region)
  {
    if (!R_isSliver(region)) return;

    bool solvedSlivGlob = false;
    bool solvedGlob = false;

    bool solvedSliv = false;
    bool solved = false;

    double bestShape = 0.;
    mqm.getShape(region,&bestShape);

    MeshParametersManager& mpm = MeshParametersManagerSgl::instance();

    edgeCollapseOp eCollOp (mesh,sizeField);
    eCollOp.collapseOnBoundary(collapseOnBdry,dVTolClps);

    faceCollapseOp fCollOp (mesh,sizeField);
    fCollOp.collapseOnBoundary(collapseOnBdry,dVTolClps);

    edgeSwapOp     eSwapOp (mesh,sizeField);
    eSwapOp.swapOnBoundary(swapOnBdry,dVTolSwap);

    DESCOp         descOp(mesh,sizeField);
    edgeSplitOp    eSplitOp(mesh,sizeField);
    faceSwapOp     fSwapOp(mesh,sizeField);
    vertexMoveOp   vMoveOp(mesh,sizeField,false);

    pEntity keyEnts[4];
    int sType = getSliverType(region, keyEnts);
    if( sType == 0 ) return;

    double worst;

    // --- TYPE I slivers ---
    if( sType == 1 ) 
      {
        tested_I++;

        // check split 
        // ------------
        solvedSliv = false;
        solved = false;
        for(int i=0; i<2; i++) {

          eSplitOp.setSplitEdge((pEdge)keyEnts[i]);
        
          if( E_whatInType((pEdge)keyEnts[i]) != 3 ) continue;

          // determine the position of the new vertex
//           double newPos[3];
//           E_cavityCenter( (pEdge)keyEnts[i], newPos );
//           double t = E_linearParams((pEdge)keyEnts[i],newPos);
//           if( t>1. || t<0. ) continue;
//           eSplitOp.setSplitPos(newPos);

          worst = -1.;
          if( eSplitOp.evaluate(&worst) && worst > bestShape &&
              eSplitOp.getMinLenSq() > mpm.getSliverLowerLengthSqBound() &&
              eSplitOp.getMaxLenSq() < mpm.getSliverUpperLengthSqBound() ) {
            if (worst < mpm.getSliverTetBound()) {
              solvedSliv = true;
            }
            else {
              solved = true;
              break;
            }
          }
        }
        if ( solved ) { solvedGlob=true; splitSolves_I++; }
        else if ( solvedSliv ) { solvedSlivGlob=true; splitSolvesToSliver_I++; }

        // check edge collapses
        // ---------------------
        solvedSliv = false;
        solved = false;
        for(int i=0; i<6; i++) {
          pEdge edge = R_edge(region,i);
          for (int j=0; j<2; j++) {
            eCollOp.setCollapseEdge(edge,E_vertex(edge,j),E_vertex(edge,1-j));

            worst = -1.;
            if ( eCollOp.evaluate(&worst) && worst > bestShape &&
                 eCollOp.getMinLenSq() > mpm.getSliverLowerLengthSqBound() &&
                 eCollOp.getMaxLenSq() < mpm.getSliverUpperLengthSqBound() ) {
              if (worst < mpm.getSliverTetBound()) {
                solvedSliv = true;
              }
              else {
                solved = true;
                break;
              }
            }
          }
          if (worst > mpm.getSliverTetBound()) break;
        }
        if ( solved ) { solvedGlob=true; clpsSolves_I++; }
        else if ( solvedSliv ) { solvedSlivGlob=true; clpsSolvesToSliver_I++; }

        // check double edge split collapse
        // ---------------------------------
        solvedSliv = false;
        solved = false;
        for(int i=0; i<2; i++) {
          descOp.setDESC(region,(pEdge)keyEnts[i%2],(pEdge)keyEnts[(i+1)%2]);
  
          worst = -1.;
          if( descOp.evaluate(&worst) && worst > bestShape &&
              descOp.getMinLenSq() > mpm.getSliverLowerLengthSqBound() &&
              descOp.getMaxLenSq() < mpm.getSliverUpperLengthSqBound() ) {
            if (worst < mpm.getSliverTetBound()) {
              solvedSliv = true;
            }
            else {
              solved = true;
            }
          }
        }
        if ( solved ) { solvedGlob=true; descSolves_I++; }
        else if ( solvedSliv ) { solvedSlivGlob=true; descSolvesToSliver_I++; }

        // check face collapse
        // --------------------
        solvedSliv = false;
        solved = false;
        for(int i=0; i<2; i++) {
          for( int ClpsOnvt=1; ClpsOnvt>=0; ClpsOnvt-- ) {
            fCollOp.reset((pFace)keyEnts[i+2],(pEdge)keyEnts[i], ClpsOnvt);

            worst = -1.;
            if( fCollOp.evaluate(&worst) && worst > bestShape &&
                fCollOp.getMinLenSq() > mpm.getSliverLowerLengthSqBound() &&
                fCollOp.getMaxLenSq() < mpm.getSliverUpperLengthSqBound() ) {
              if (worst < mpm.getSliverTetBound()) {
                solvedSliv = true;
              }
              else {
                solved = true;
                break;
              }
            }
          }
          if( worst > mpm.getSliverTetBound() ) break;
        }
        if ( solved ) { solvedGlob=true; fclpsSolves_I++; }
        else if ( solvedSliv ) { solvedSlivGlob=true; fclpsSolvesToSliver_I++; }
 
        // check edge swaps (2 edges)
        // ---------------------------
        solvedSliv = false;
        solved = false;
        for(int i=0; i<2; i++) {
          eSwapOp.setSwapEdge((pEdge)keyEnts[i]);

          worst = -1.;
          if( eSwapOp.evaluate(&worst) && worst > bestShape &&
              eSwapOp.getMinLenSq() > mpm.getSliverLowerLengthSqBound() &&
              eSwapOp.getMaxLenSq() < mpm.getSliverUpperLengthSqBound() ) {
            if (worst < mpm.getSliverTetBound()) {
              solvedSliv = true;
            }
            else {
              solved = true;
              break;
            }
          }
        }
        if ( solved ) { solvedGlob=true; swapSolves_I++; }
        else if ( solvedSliv ) { solvedSlivGlob=true; swapSolvesToSliver_I++; }

        // check vertex move 
        // ------------------
        solvedSliv = false;
        solved = false;
        pPList verts = R_vertices(region);
        void * temp = NULL;
        while ( pVertex pv = (pVertex)PList_next(verts,&temp) ) {
          if( V_whatInType(pv) != 3 ) continue;
        
          vMoveOp.setPositionToOptimal(pv);

          worst = -1.;
          if( vMoveOp.evaluate(&worst) && worst > bestShape &&
              vMoveOp.getMinLenSq() > mpm.getSliverLowerLengthSqBound() &&
              vMoveOp.getMaxLenSq() < mpm.getSliverUpperLengthSqBound() ) {
            if (worst < mpm.getSliverTetBound()) {
              solvedSliv = true;
            }
            else {
              solved = true;
              break;
            }
          }
        }
        PList_delete(verts);
        if ( solved ) { solvedGlob=true; nodeSolves_I++; }
        else if ( solvedSliv ) { solvedSlivGlob=true; nodeSolvesToSliver_I++; }

        if (!solvedGlob) {
          set<pRegion>::const_iterator sIter = notSolved.find(region);
          if (sIter == notSolved.end() ) {
            notSolved.insert(region);
            notSolved_I++;
            if (solvedSlivGlob) solvedWithSliver_I++;
          }
          else {
            notSolvedDuplicated_I++;
            if (solvedSlivGlob) solvedWithSliverDuplicated_I++;          
          }
        }

      }

    // --- TYPE II slivers ---
    else 
      {
        tested_II++;

        // check edge collapses
        // ---------------------
        solvedSliv = false;
        solved = false;
        for(int i=0; i<6; i++) {
          pEdge edge = R_edge(region,i);
          for (int j=0; j<2; j++) {
            eCollOp.setCollapseEdge(edge,E_vertex(edge,j),E_vertex(edge,(j+1)%2));

            worst = -1.;
            if ( eCollOp.evaluate(&worst) && worst > bestShape &&
                 eCollOp.getMinLenSq() > mpm.getSliverLowerLengthSqBound() &&
                 eCollOp.getMaxLenSq() < mpm.getSliverUpperLengthSqBound() ) {
              if (worst < mpm.getSliverTetBound()) {
                solvedSliv = true;
              }
              else {
                solved = true;
                break;
              }
            }
          }
          if (worst > mpm.getSliverTetBound()) break;
        }
        if ( solved ) { solvedGlob=true; clpsSolves_II++; }
        else if ( solvedSliv ) { solvedSlivGlob=true; clpsSolvesToSliver_II++; }

        // check face collapse
        // --------------------

        solvedSliv = false;
        solved = false;
        for( int ClpsOnvt=1; ClpsOnvt>=0; ClpsOnvt-- ) {
          fCollOp.reset((pFace)keyEnts[2],(pEdge)keyEnts[1], ClpsOnvt);
  
          worst = -1.;
          if( fCollOp.evaluate(&worst) && worst > bestShape &&
              fCollOp.getMinLenSq() > mpm.getSliverLowerLengthSqBound() &&
              fCollOp.getMaxLenSq() < mpm.getSliverUpperLengthSqBound() ) {
            if (worst < mpm.getSliverTetBound()) {
              solvedSliv = true;
              break;
            }
            else {
              solved = true;
              break;
            }
          }
        }
        if ( solved ) { solvedGlob=true; fclpsSolves_II++; }
        else if ( solvedSliv ) { solvedSlivGlob=true; fclpsSolvesToSliver_II++; }
      
        // check edge swaps
        // -----------------
        solvedSliv = false;
        solved = false;
        for(int i=0; i<3; i++) {
          pEdge edge = F_edge((pFace)keyEnts[0],i);
          eSwapOp.setSwapEdge(edge); 

          worst = -1.;
          if( eSwapOp.evaluate(&worst) && worst > bestShape &&
              eSwapOp.getMinLenSq() > mpm.getSliverLowerLengthSqBound() &&
              eSwapOp.getMaxLenSq() < mpm.getSliverUpperLengthSqBound() ) {
            if (worst < mpm.getSliverTetBound()) {
              solvedSliv = true;
            }
            else {
              solved = true;
              break;
            }
          }
        }
        if ( solved ) { solvedGlob=true; swapSolves_II++; }
        else if ( solvedSliv ) { solvedSlivGlob=true; swapSolvesToSliver_II++; }
      
        // check face swap
        // ----------------
        solvedSliv = false;
        solved = false;
        fSwapOp.setSwapFace((pFace)keyEnts[0]);

        worst = -1.;
        if( fSwapOp.evaluate(&worst) && worst > bestShape &&
            fSwapOp.getMinLenSq() > mpm.getSliverLowerLengthSqBound() &&
            fSwapOp.getMaxLenSq() < mpm.getSliverUpperLengthSqBound() ) {
          if (worst < mpm.getSliverTetBound()) solvedSliv = true;
          else                                 solved = true;
        } 
        if ( solved ) { solvedGlob=true; fswapSolves_II++; }
        else if ( solvedSliv ) { solvedSlivGlob=true; fswapSolvesToSliver_II++; }

        // check vertex move 
        // ------------------
        solvedSliv = false;
        solved = false;
        pPList verts = R_vertices(region);
        void * temp = NULL;
        while ( pVertex pv = (pVertex)PList_next(verts,&temp) ) {
          if( V_whatInType(pv) != 3 ) continue;
        
          vMoveOp.setPositionToOptimal(pv);

          worst = -1.;
          if( vMoveOp.evaluate(&worst) && worst > bestShape &&
              vMoveOp.getMinLenSq() > mpm.getSliverLowerLengthSqBound() &&
              vMoveOp.getMaxLenSq() < mpm.getSliverUpperLengthSqBound() ) {
            if (worst < mpm.getSliverTetBound()) {
              solvedSliv = true;
            }
            else {
              solved = true;
              break;
            }
          }
        }
        PList_delete(verts);
        if ( solved ) { solvedGlob=true; nodeSolves_II++; }
        else if ( solvedSliv ) { solvedSlivGlob=true; nodeSolvesToSliver_II++; }

        if (!solvedGlob) {
          set<pRegion>::const_iterator sIter = notSolved.find(region);
          if (sIter == notSolved.end() ) {
            notSolved.insert(region);
            notSolved_II++;
            if (solvedSlivGlob) solvedWithSliver_II++;
          }
          else {
            notSolvedDuplicated_II++;
            if (solvedSlivGlob) solvedWithSliverDuplicated_II++;          
          }
        }
      }
  
  }

  // -------------------------------------------------------------------
  void sliverRegionHandler::printOperatorsTest(ostream& out) const
  {
    int testedUnique_I = tested_I-solvedWithSliverDuplicated_I-notSolvedDuplicated_I;
    int solved_I = tested_I-solvedWithSliver_I-notSolved_I-solvedWithSliverDuplicated_I-notSolvedDuplicated_I;
    int testedUnique_II = tested_II-solvedWithSliverDuplicated_II-notSolvedDuplicated_II;
    int solved_II = tested_II-solvedWithSliver_II-notSolved_II-solvedWithSliverDuplicated_II-notSolvedDuplicated_II;

    out << "\n*** SliverRegionHandler statistics (every operation is tested) ***\n\n";
    out << "Slivers considered: " << testedUnique_I+testedUnique_II<<endl;
    out << "\t    Tested\t Tested dupl.\t   Solved  \t New sliver\t New sl. dupl.\t Not solved\t Not solved dupl.\n";
    out << "  - Type I: \t"<< testedUnique_I  
        << "\t\t" << tested_I-testedUnique_I
        << "\t\t" << solved_I
        << "\t\t" << solvedWithSliver_I 
        << "\t\t" << solvedWithSliverDuplicated_I 
        << "\t\t" << notSolved_I 
        << "\t\t" << notSolvedDuplicated_I 
        << "\n";
    out << "  - Type II: \t"<< testedUnique_II
        << "\t\t" << tested_II-testedUnique_II
        << "\t\t" << solved_II
        << "\t\t" << solvedWithSliver_II 
        << "\t\t" << solvedWithSliverDuplicated_II 
        << "\t\t" << notSolved_II 
        << "\t\t" << notSolvedDuplicated_II
        << "\n";

    out << " \n--- Type I slivers ---\n\n";
    out << "\t     Solved\t Solved with a sliver\n";
    out << "ESwap:\t\t"<< swapSolves_I << "\t\t" << swapSolvesToSliver_I << "\n";
    out << "EClps:\t\t"<< clpsSolves_I << "\t\t" << clpsSolvesToSliver_I << "\n";
    out << "FClps:\t\t"<< fclpsSolves_I << "\t\t" << fclpsSolvesToSliver_I << "\n";
    out << "DESpltClps:\t"<< descSolves_I << "\t\t" << descSolvesToSliver_I << "\n";
    out << "ESplt:\t\t"<< splitSolves_I << "\t\t" << splitSolvesToSliver_I << "\n";
    out << "Reloc:\t\t"<< nodeSolves_I << "\t\t" << nodeSolvesToSliver_I << "\n";

    out << " \n--- Type II slivers ---\n\n";
    out << "\t     Solved\t Solved with a sliver\n";
    out << "ESwap:\t\t"<< swapSolves_II << "\t\t" << swapSolvesToSliver_II << "\n";
    out << "FSwap:\t\t"<< fswapSolves_II << "\t\t" << fswapSolvesToSliver_II << "\n";
    out << "EClps:\t\t"<< clpsSolves_II << "\t\t" << clpsSolvesToSliver_II << "\n";
    out << "FClps:\t\t"<< fclpsSolves_II << "\t\t" << fclpsSolvesToSliver_II << "\n";
    out << "Reloc:\t\t"<< nodeSolves_II << "\t\t" << nodeSolvesToSliver_II << "\n";

    out << "\n\n";

  }

  // -------------------------------------------------------------------
  void sliverRegionHandler::reportSliver(pRegion region)
  {
    stringstream ss;
    string idStr;  ss << reportId;  ss >> idStr;

    string name = reportPrefix + "sliver" + idStr + ".pos";
  
    pPList cavity = PList_new();

    PList_append(cavity, (pEntity)region);
  
    for (int iV=0; iV<R_numVertices(region); iV++) {

      pVertex vertex = R_vertex(region, iV);
    
      pPList vRegs = V_regions(vertex);
      void * temp = 0;
      while ( pRegion pr = (pRegion)PList_next(vRegs,&temp) ) {
        PList_appUnique(cavity,(pEntity)pr);
      }
      PList_delete(vRegs);
    }

    printPosEntities(cavity,name.c_str(),OD_MEANRATIO,sizeField);

    PList_delete(cavity);

    reportId++;
  }

  // -------------------------------------------------------------------

}
