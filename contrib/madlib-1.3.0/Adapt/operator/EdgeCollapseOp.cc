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

#include "EdgeCollapseOp.h"
#include "OperatorTools.h"
#include <map>
#include <cmath>

using std::map;

namespace MAd {

  // -------------------------------------------------------------------
  bool edgeCollapseOp::checkConstraints() const
  {
    if ( EN_constrained((pEntity)edgeDel) ) return false;
    if ( EN_constrained((pEntity)vDel   ) ) return false;
    return true;
  }

  // -------------------------------------------------------------------
  bool edgeCollapseOp::checkGeometry()
  {
    int typeD = V_whatInType(vDel);
    if ( typeD == 0 ) return false;

    if ( constrainBoundary ) {
      switch( typeD ) {
      case 3: break;
      case 2:
        {
          bool ok = true;
          pPList eRgns = E_regions(edgeDel);
          if( PList_size(eRgns) != 0 )  ok = false;
          PList_delete(eRgns);
          if (!ok) return false;
          break;
        }
      case 1: return false;
      case 0: return false;
      }
    }
    else {
      // check in edgeDel
      if ( V_whatIn(vDel) != E_whatIn(edgeDel) ) return false;
      
      if ( V_whatInType(vDel) <= 2 )
        {
          // check in every tet around edgeDel (= to be collapsed):
          //   - no dimension reduction surface on surface
          pEdge oppE;
          pFace newSurfaceFace;
          pPList vDelFaces = V_faces(vDel);
          void * temp = NULL;
          pFace pf;
          while ( ( pf = (pFace)PList_next(vDelFaces,&temp) ) ) {
            if ( F_inClosure(pf,(pEntity)edgeDel) ) continue;
            if ( F_whatInType(pf) == 2 ) {
              oppE           = F_vtOpEd(pf,vDel);
              newSurfaceFace = F_exist(oppE,vTgt);
              if ( newSurfaceFace && F_whatInType(newSurfaceFace) == 2 ) {
                PList_delete(vDelFaces);
                return false;
              }
            }
          }
          PList_delete(vDelFaces);
      
          
          // check in every face around edge del (= to be collapsed):
          //   - the face exists (that's why we loop over edges, not faces)
          //   - no new contact between surface and line
          //   - no new contact between surface and surface
          //   - no dimension reduction line on line
          //
          //                     vTgt 
          //
          //                   /  |
          //                  /   |
          //                 /    |
          //                /     |  oppE
          //     edgeDel   /      |
          //              /  pf   |
          //             /        |
          //            /         |
          //
          //          vDel ----- oppV
          //
          //                 pe
          //

          pVertex oppV;
          pEdge pe;
          pPList vDelEdges = V_edges(vDel);
          temp = NULL;
          while ( ( pe = (pEdge)PList_next(vDelEdges,&temp) ) ) {
            if ( pe == edgeDel ) continue;
            int d1 = E_whatInType(pe);
            if ( d1 <= 2 ) {
              oppV = E_otherVertex(pe,vDel);
              oppE = E_exist(oppV,vTgt);
              if ( oppE ) {
                int d2 = E_whatInType(oppE);
                if ( d2 <= 2 ) {
                  pf = F_exist(oppE,vDel);

                  // check if the face exists
                  if ( !pf ) {
                    PList_delete(vDelEdges); return false;
                  }

                  // check for line on line dimension reduction
                  if ( d1 == 1 && d2 == 1 ) {
                    PList_delete(vDelEdges); return false;
                  }

                  pGEntity fge = F_whatIn(pf);

                  // check for new contacts between a line and a surface
                  if ( d1 == 1 && d2 == 2 && fge != E_whatIn(oppE) ) {
                    PList_delete(vDelEdges); return false;
                  }
                  if ( d1 == 2 && d2 == 1 && fge != E_whatIn(pe) ) {
                    PList_delete(vDelEdges); return false;
                  }

                  // check for new contacts between 2 surfaces
                  if ( d1 == 2 && d2 == 2 ) {
                    if ( fge != E_whatIn(pe) || fge != E_whatIn(oppE) ) {
                      PList_delete(vDelEdges); return false;
                    }
                  }
                }
              }
            }
          }
          PList_delete(vDelEdges);
          
        }
    }

    return true;
  }

  // -------------------------------------------------------------------
  bool edgeCollapseOp::evaluateShapes2D()
  {
    pPList vDelFaces = V_faces(vDel);

    double worstShape = mqm.getElementEvaluator()->bestShapeEver();

    map<pGEntity,double> area_after;

    void * temp1 = 0;
    while ( pFace face = (pFace)PList_next(vDelFaces,&temp1) ) {

      // Faces surrounding edgeDel will be deleted
      if (F_inClosure(face,(pEntity)edgeDel)) continue;

      // Gather the coordinates and sizes of the new face
      double xyz[3][3];
      pMSize fSizes[3] = {NULL,NULL,NULL};

      pPList fVerts = F_vertices(face,1);
      void * temp2 = 0;
      int iV = 0;
      while ( pVertex vertex = (pVertex)PList_next(fVerts, &temp2) ) {
        if (vertex == vDel) {
          fSizes[iV] = sizeField->findSize(vTgt);
          V_coord(vTgt,xyz[iV]);
        }
        else {
          fSizes[iV] = sizeField->findSize(vertex);
          V_coord(vertex,xyz[iV]);
        }
        iV++;
      }
      PList_delete(fVerts);

      // Check if the shape is acceptable and compare it to the worst shape
      double shape = 0.;
      double normalBefore[3];
      F_normal(face,normalBefore);
      if ( !mqm.getElementEvaluator()->XYZ_F_shape(xyz, fSizes, normalBefore, &shape) ) {
        PList_delete(vDelFaces); 
        return false;
      }
      else {
        if ( shape < worstShape ) worstShape = shape;
      }

      // Compute the area of the new face
      if( !constrainBoundary && V_whatInType(vDel) < 2 ) {

        pGEntity g = F_whatIn(face);
        std::map<pGEntity,double>::iterator it = area_after.find( g );
        if( it != area_after.end() )
         (*it).second += XYZ_F_area(xyz,normalBefore);
        else
         area_after[g] = XYZ_F_area(xyz,normalBefore);
      }
    }

    // If the edge is on a boundary, check that the area of
    // the old and the new cavities are the same for each Geo entities
    if( !constrainBoundary && V_whatInType(vDel) < 2 ) {
      map< pGEntity,double>  area_before;
      map<pGEntity,double>::iterator it, it2;

      void * temp3 = 0;
      while ( pFace face = (pFace)PList_next(vDelFaces,&temp3) ) {
        double normal[3];
        F_normal(face,normal);
        pGEntity g = F_whatIn(face);
        it = area_before.find( g );
        if( it != area_before.end() )
          (*it).second += F_area(face,normal);
        else
          area_before[g] = F_area(face,normal);
      }

      for( it=area_before.begin(); it!=area_before.end(); ++it  )
      {
        double before = (*it).second;
        double after( 0.0 );
        it2 = area_after.find( (*it).first );
        if( it2 != area_after.end() )
          after = (*it2).second;

        // Check area change
        double dA = fabs( (before-after) / before );
        if( dA > dATol ) {
          PList_delete(vDelFaces);
          return false;
        }
      }
    }
    results->setWorstShape(worstShape);

    PList_delete(vDelFaces);

    return true;
  }

  // -------------------------------------------------------------------
  bool edgeCollapseOp::evaluateShapes()
  {
    pPList vDelRgns = V_regions(vDel);

    // 2D case
    if( PList_size(vDelRgns)==0 ) {
      bool flag = evaluateShapes2D();
      PList_delete(vDelRgns);
      return flag;
    }

    // 3D case
    else {

      double worstShape = mqm.getElementEvaluator()->bestShapeEver();
      double volAfter = 0.;

      void * temp1 = 0;
      while ( pRegion region = (pRegion)PList_next(vDelRgns,&temp1) ) {
      
        // Regions surrounding edgeDel will be deleted
        if (R_inClosure(region,(pEntity)edgeDel)) continue;

        // Gather the coordinates and sizes of the new region
        double xyz[4][3];
        pMSize rSizes[4] = {NULL,NULL,NULL,NULL};

        pPList rVerts = R_vertices(region);
        void * temp2 = 0;
        int iV = 0;
        while ( pVertex vertex = (pVertex)PList_next(rVerts, &temp2) ) {
          if (vertex == vDel) {
            rSizes[iV] = sizeField->findSize(vTgt);
            V_coord(vTgt,xyz[iV]);
          }
          else {
            rSizes[iV] = sizeField->findSize(vertex);
            V_coord(vertex,xyz[iV]);
          }
          iV++;
        }
        PList_delete(rVerts);

        // Check if the shape is acceptable and compare it to the worst shape
        double shape = 0.;
        if ( !mqm.getElementEvaluator()->XYZ_R_shape(xyz, rSizes, &shape) ) {
          PList_delete(vDelRgns); 
          return false;
        }
        else {
          if ( shape < worstShape ) worstShape = shape;
        }

        // Compute the volume of the new region
        if( !constrainBoundary && V_whatInType(vDel)!=3 ) {
          volAfter += R_XYZ_volume (xyz);
        }
      }

      // If the edge is on a boundary, check that the volume of 
      // the old and the new cavities are the same
      if( !constrainBoundary && V_whatInType(vDel)!=3 ) {

        double volBefore = 0.;
      
        void * temp3 = 0;
        while ( pRegion region = (pRegion)PList_next(vDelRgns,&temp3) ) {
          volBefore += R_volume (region);
        }

        double dV = fabs( (volBefore-volAfter) / volBefore );
        if( dV > dVTol ) {
          PList_delete(vDelRgns);
          return false;
        }
      }
    
      results->setWorstShape(worstShape);
    }
  
    PList_delete(vDelRgns); 
  
    return true;
  }

  // -------------------------------------------------------------------
  void edgeCollapseOp::evaluateLengths() const
  {
    double minSq = MAdBIG, maxSq = 0.;

    double xyzTgt[3];
    V_coord(vTgt,xyzTgt);
    pMSize sizeTgt = sizeField->findSize(vTgt);

    // list nodes linked to target vertex by an edge
    std::set<pVertex> vTgtLinked;
    pPList vTgtEdges = V_edges(vTgt);
    void * temp1 = 0;
    while ( pEdge edge = (pEdge) PList_next(vTgtEdges,&temp1) ) {
      vTgtLinked.insert(E_vertex(edge,0));
      vTgtLinked.insert(E_vertex(edge,1));
    }
    PList_delete(vTgtEdges);

    // find nodes linked to vDel and not linked to vTgt
    pPList vDelEdges = V_edges(vDel);
    void * temp2 = 0;
    while ( pEdge edge = (pEdge) PList_next(vDelEdges,&temp2) ) {
      for (int iV=0; iV<2; iV++) {
        pVertex pV = E_vertex(edge,iV);
        if ( pV == vDel ) continue;
        if ( vTgtLinked.find(pV) != vTgtLinked.end() ) continue;
      
        double xyzV[3];
        V_coord(pV,xyzV);
        pMSize sizeV = sizeField->findSize(pV);
      
        double newSq = sizeField->SF_XYZ_lengthSq(xyzTgt,xyzV,sizeTgt,sizeV);
        if( newSq < minSq ) { minSq = newSq; }
        if( newSq > maxSq ) { maxSq = newSq; }
      }
    }
    PList_delete(vDelEdges);

    results->setMaxLenSq(maxSq);
    results->setMinLenSq(minSq);
  }

  // -------------------------------------------------------------------
  void edgeCollapseOp::getCavity(pPList * cavity) const
  {
    if ( dim == 3 ) *cavity = V_regions(vDel);
    else            *cavity = V_faces(vDel);
  }

  // -------------------------------------------------------------------
  void edgeCollapseOp::apply()
  {
    E_collapse(mesh, edgeDel, vDel, vTgt);

    HistorySgl::instance().add((int)type(),OPERATOR_APPLY,1);
  }

  // -------------------------------------------------------------------

}

