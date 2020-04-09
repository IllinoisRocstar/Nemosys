// -------------------------------------------------------------------
// MAdLib - Copyright (C) 2008-2009 Universite catholique de Louvain
//
// See the Copyright.txt and License.txt files for license information. 
// You should have received a copy of these files along with MAdLib. 
// If not, see <http://www.madlib.be/license/>
//
// Please report all bugs and problems to <contrib@madlib.be>
//
// Authors: Arnaud Francois, Gaetan Compere, Jean-Francois Remacle
// -------------------------------------------------------------------

#include "EdgeSwapOp.h"
#include "MAdDefines.h"
#include "CallbackManager.h"
#include "MAdOutput.h"
#include "MathUtils.h"
#include "MAdMessage.h"

#include <iostream>
#include <sstream>
using std::cout;
using std::endl;
using std::cerr;
using std::stringstream;
#include <math.h>
#include <assert.h>

namespace MAd
{
  // -------------------------------------------------------------------
  bool edgeSwapOp::checkConstraints() const
  {
    if( EN_constrained((pEntity)edge) ) return false;
    return true;
  }

  // -------------------------------------------------------------------
  bool edgeSwapOp::checkGeometry()
  {
    // set the available swap configurations for a given edge 
    // regarding geometry
    if( M_numRegions(mesh)==0 ) return checkGeometry2D();  // 2D mesh
    else                        return checkGeometry3D();  // 3D mesh
    return false;
  }

  // -------------------------------------------------------------------
  bool edgeSwapOp::checkGeometry2D()
  {
    if ( E_numFaces(edge)   != 2 ) return false;
    if ( E_whatInType(edge) <= 1 ) return false;
    
    // check for dimension reduction
    pFace face0 = E_face(edge,0);
    pVertex pv0 = F_edOpVt(face0,edge);
    pFace face1 = E_face(edge,1);
    pVertex pv1 = F_edOpVt(face1,edge);
    if ( E_exist(pv0,pv1) ) return false;

    return true;
  }

  // -------------------------------------------------------------------
  // Set the swap configuration to be used. Return false if the edge swap
  // operation is not allowed
  bool edgeSwapOp::checkGeometry3D()
  {
    int nbf = E_numFaces(edge);  // nb of faces attached to the edge
    if( nbf < 3 || nbf >7 )      // check if swap configuration is implemented
      return false;

    // check the dimension of the geometrical entity on which edge is classified
    bool boundarySwap = false;
    int gDim = E_whatInType(edge);
    if( gDim == 1 ) return false;
    else if ( gDim == 2 )
      {
        if ( constrainBoundary ) return false;
        else {
          checkVol = true;
          boundarySwap = true;
        }
      }

    // get 'face' that connects to edge. if one classified on a surface, take that
    // one, otherwise take any (e.g. last one)
    int nbBFaces = 0;
    pFace face;
    pFace bFace=NULL;
    for( int i=0; i < nbf; i++ )
      {
        face = E_face(edge,i);
        if( F_whatInType(face) == 2 ) {
          nbBFaces++;
          bFace = face;
        }
      }
    if (bFace) face = bFace;

    // reject swaps on non-manifold cavities
    if ( nbBFaces > 2 || nbBFaces == 1 ) {
      exportCavity("swap_nonmanif1.pos");
      MAdMsgSgl::instance().warning(__LINE__,__FILE__,
                                    "Found a non-manifold domain, rejecting edge swap");
      return false;
    }

    // get region connected to face
    pRegion region = F_region(face,0);
    if( !region ) region = F_region(face,1);

    // Get vertices of the polyhedron around 'edge'
    vertices.clear();
    vertices.resize(nbf);
    pFace current_face = face;
    pRegion current_region = region;
    for( int i=0; i < nbf; i++ )
      {
        vertices[i] = F_edOpVt(current_face, edge);
        
        if (i == (nbf-1)) {
          if ( boundarySwap && current_region != NULL ) {
            // non-manifold-domain: a geometrical face with 2 regions
            exportCavity("swap_nonmanif2.pos");
            MAdMsgSgl::instance().warning(__LINE__,__FILE__,
                                          "Found a non-manifold domain, rejecting edge swap");
            return false;
          }
        }
        else {
          pFace next_face = E_otherFace(edge, current_face, current_region);
          pRegion next_region = F_region(next_face, 0);
          if( next_region == current_region ) {
            next_region = F_region(next_face, 1);
          }
          current_face   = next_face;
          current_region = next_region;
        }
      }

    // Check dimension reduction: if swap on boundary, check that no edge 
    // links the first and last points of the crown. 
    if ( boundarySwap ) {
      if ( E_exist(vertices[0],vertices[nbf-1]) ) return false;
    }

    // Find top and bottom vertices from orientation of the polyhedron
    // v2 is top vertex if ( V01 x V02 ) . V32 > 0
    pVertex  v0 = vertices[0];
    pVertex  v1 = vertices[1];
    pVertex  v2 = E_vertex(edge, 0);
    pVertex  v3 = E_vertex(edge, 1);

    double p[4][3];
    V_coord( v0, p[0] );
    V_coord( v1, p[1] );
    V_coord( v2, p[2] );
    V_coord( v3, p[3] );
    double a[3] = { p[1][0]-p[0][0], p[1][1]-p[0][1], p[1][2]-p[0][2] };
    double b[3] = { p[2][0]-p[0][0], p[2][1]-p[0][1], p[2][2]-p[0][2] };
    double c[3] = { p[2][0]-p[3][0], p[2][1]-p[3][1], p[2][2]-p[3][2] };
    double ab[3] = { a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0] };  // a x b
    if( ab[0]*c[0] + ab[1]*c[1] + ab[2]*c[2] > 0.0 )
      {
        vt = v2;
        vb = v3;
      }
    else
      {
        vt = v3;
        vb = v2;
      }

    swap_config.set( nbf );
    return true;
  }

  // -------------------------------------------------------------------
  bool edgeSwapOp::evaluateShapes()
  {
    // select the best triangulation and evaluate worst shape
    if( M_numRegions(mesh)==0 ) return evaluateShapes2D();
    else                        return evaluateShapes3D();
    return false;
  }

  // -------------------------------------------------------------------
  // loop over all possible swap configurations and select the best triangulation
  bool edgeSwapOp::evaluateShapes3D()
  {
    // if edge is classified on face type need to check volume
    double volume = 0.0;
    if( checkVol )
      {
        pRegion region;
        pPList regions = E_regions (edge);
        void *temp=0;
        while( ( region = (pRegion)PList_next(regions, &temp) ) )
          volume += R_volume (region);
        PList_delete (regions);
      }
    int nb_triangles = swap_config.nb_triangles();
    // here nb_triangles is maximum 35 (for config = 7), 
    // ansi c++ does not support variable lengths for arrays
    double shapes[35];     // worst shape for each configuration
    bool   valid_shapes[35];
    double volumes[35];    // volumes of the new swap configurations

    // -- loop on all possible triangulations and compute associated element shapes --
    for( int i=0; i< nb_triangles; i++ )
      {
        double xyz[4][3];
        pMSize s[4];

        pVertex v0 = vertices[ swap_config.triangle(i,0) ];
        pVertex v1 = vertices[ swap_config.triangle(i,1) ];
        pVertex v2 = vertices[ swap_config.triangle(i,2) ];

        V_coord(v0, xyz[0]);
        V_coord(v1, xyz[1]);
        V_coord(v2, xyz[2]);
        V_coord(vt, xyz[3]);

        s[0] = sizeField->findSize(v0);
        s[1] = sizeField->findSize(v1);
        s[2] = sizeField->findSize(v2);
        s[3] = sizeField->findSize(vt);

        // evaluate 'top' Tet [v0,v1,v2,vt]
        double shape;
        bool shape_ok = mqm.getElementEvaluator()->XYZ_R_shape( xyz, s , &shape );
        valid_shapes[i] = shape_ok;
        shapes[i] = shape;

        if( !shape_ok )    // no need to compute shape of bottom tet since
          {                  // configuration will be rejected
            valid_shapes[i] = false;
            continue;
          }

        // evaluate 'bottom' Tet [v0,v2,v1,vb]
        V_coord(v2, xyz[1]);
        V_coord(v1, xyz[2]);
        V_coord(vb, xyz[3]);
        s[1] = sizeField->findSize(v2);
        s[2] = sizeField->findSize(v1);
        s[3] = sizeField->findSize(vb);

        valid_shapes[i] = mqm.getElementEvaluator()->XYZ_R_shape( xyz, s , &shape );

        // keep worst value of top-bottom tets for swap configuration i in shapes[i]
        if( shape < shapes[i] ) shapes[i] = shape;

        // volume computation
        if( checkVol )
          {
            volumes[i]  = R_XYZ_volume( xyz ); // bottom tet
            V_coord(vt,xyz[3]);
            volumes[i] += R_XYZ_volume( xyz ); // top tet
          }
      }

    // -- Select best valid edge swap configuration --
    double optimum = 0.0;
    for( int i=0; i < swap_config.nb_triangulations(); i++ )
      {
        // find worst element shape of this configuration
        double worst_shape = MAdBIG;    // assume the shape is always between 0.0 and 1.0
        bool   valid = false;
        for( int j=0; j < swap_config.nb_tri_triangulation(); j++ )
          {
            int t = swap_config.triangulation(i,j);  // test triangle t
            valid = valid_shapes[t];
            if( !valid )
              break;

            if( shapes[t] < worst_shape )
              worst_shape = shapes[t];
          }

        // if shape is not valid, check next swap configuration
        if( !valid )
          continue;

        // check volume change
        if( checkVol )
          {
            double new_volume = 0.0;
            for( int j=0; j < swap_config.nb_tri_triangulation(); j++ )
              new_volume += volumes[ swap_config.triangulation(i,j) ];

            if( fabs(volume-new_volume)/volume > dVTol )
              continue;
          }

        if( worst_shape > optimum )
          {
            optimum = worst_shape;  // worst shape of the configuration
            conf =  i;              // best swap configuration
          }
      }

    // no valid swap configuration was found
    if( optimum <= MAdTOL ) return false;

    results->setWorstShape( optimum );

    return true;
  }

  // -------------------------------------------------------------------
  bool edgeSwapOp::evaluateShapes2D()
  {
    double worstShape = MAdBIG;
  
    // evaluate shape of new triangles
    pFace f0 = E_face(edge,0);
    pFace f1 = E_face(edge,1);
    pVertex v0 = E_vertex(edge,0);
    pVertex v1 = E_vertex(edge,1);
    pVertex v2 = F_edOpVt(f0, edge);
    pVertex v3 = F_edOpVt(f1, edge);
  
    double xyz[3][3];
    V_coord(v0, xyz[0]);
    V_coord(v2, xyz[1]);
    V_coord(v3, xyz[2]);
  
    pMSize s[3];
    s[0] = sizeField->findSize(v0);
    s[1] = sizeField->findSize(v2);
    s[2] = sizeField->findSize(v3);
  
    // evaluate shape of new triangles [v0,v2,v3] && [v1,v2,v3]
    double normal[3];
    double e02[3],e03[3];
    diffVec(xyz[1],xyz[0],e02);
    diffVec(xyz[2],xyz[0],e03);
    crossProd(e02,e03,normal);
    double shape;
    // [v0,v2,v3]
    if( !mqm.getElementEvaluator()->XYZ_F_shape(xyz, s, normal, &shape) ) return false;
    worstShape = shape;
  
    V_coord(v1, xyz[0]);
    s[0] = sizeField->findSize(v1);
  
    normal[0] = -normal[0];
    normal[1] = -normal[1];
    normal[2] = -normal[2];
  
    if( !mqm.getElementEvaluator()->XYZ_F_shape(xyz, s, normal, &shape) ) return false;
  
    if( worstShape > shape ) worstShape = shape;
  
    results->setWorstShape(worstShape);
    // check the area of the new faces in case of '3D' surf mesh
  
    return true;
  }

  // -------------------------------------------------------------------
  //  calculate the max/min size after the modification is applied
  void edgeSwapOp::evaluateLengths() const
  {
    double maxSwap, minSwap;

    // consider the new edges created after swap
    if( M_numRegions(mesh) == 0 )   // 2D mesh
      {
        pFace face0 = E_face(edge, 0);
        pVertex v0  = F_edOpVt(face0,edge);
        pFace face1 = E_face(edge, 1);
        pVertex v1  = F_edOpVt(face1,edge);
        maxSwap = sizeField->SF_VV_lengthSq(v0,v1);
        minSwap = maxSwap;
      }
    else
      {
        // 3D mesh
        minSwap = MAdBIG;
        maxSwap = 0.;

        // evaluate the triangles of the c'th configuration
        pVertex v[3];
        for( int j=0; j < swap_config.nb_tri_triangulation(); j++ )
          {
            int t = swap_config.triangulation(conf,j);
            for( int i=0; i<3; i++ )
              {
                v[i] = vertices[ swap_config.triangle(t,i) ];
              }

            for( int i=0; i<3; i++ )
              {
                // (AF) we could also pre-define the new edges to be tested in the EdgeSwapConfig
                // to avoid the E_exist loop
                if( !E_exist(v[i],v[(i+1)%3]) )   /// 0-1,1-2,2-0
                  {
                    double lSq = sizeField->SF_VV_lengthSq(v[i], v[(i+1)%3]);
                    if( lSq > maxSwap ) maxSwap = lSq;
                    if( lSq < minSwap ) minSwap = lSq;
                  }
              }
          }
      }

    // initialization
    double lMinSq = MAdBIG;
    double lMaxSq = 0.0;
    if( maxSwap > lMaxSq ) lMaxSq = maxSwap;
    if( minSwap < lMinSq ) lMinSq = minSwap;

    results->setMaxLenSq(lMaxSq);
    results->setMinLenSq(lMinSq);
  }

  // -------------------------------------------------------------------
  void edgeSwapOp::apply()
  {
    if( M_numRegions(mesh) !=0 ) apply3D();
    else                         apply2D();
    HistorySgl::instance().add((int)type(),OPERATOR_APPLY,1);
  }

  // -------------------------------------------------------------------
  int edgeSwapOp::apply2D()
  {
    // 2D mesh case
    pPList old_faces = PList_new();
    pPList new_faces = PList_new();

    pFace f0 = E_face(edge,0);
    pFace f1 = E_face(edge,1);
    pVertex v0 = E_vertex(edge,0);
    pVertex v1 = E_vertex(edge,1);
    pVertex v2 = F_edOpVt(f0, edge);
    pVertex v3 = F_edOpVt(f1, edge);

    pGEntity gface = F_whatIn(f0);  // classification

    // new triangles [v0,v2,v3] && [v1,v2,v3]
    pEdge newe  = M_createE(mesh, v2, v3, gface);
    pEdge e02 = F_findEdge(f0,v0,v2);
    pEdge e03 = F_findEdge(f1,v0,v3);
    pEdge e12 = F_findEdge(f0,v1,v2);
    pEdge e13 = F_findEdge(f1,v1,v3);
    pEdge fEdges[3];
    fEdges[0] = e02; fEdges[1] = newe; fEdges[2] = e03;
    pFace newf0 = M_createF( mesh, 3, fEdges, gface);
    fEdges[0] = e12; fEdges[1] = newe; fEdges[2] = e13;
    pFace newf1 = M_createF( mesh, 3, fEdges, gface);

    PList_append(new_faces, (pEntity)newf0);
    PList_append(new_faces, (pEntity)newf1);

    // call callback
    CallBackManagerSgl::instance().callCallBacks(old_faces, new_faces, MAd_ESWAP, (pEntity)edge);
    PList_delete(old_faces);
    PList_delete(new_faces);

    // deletion
    M_removeFace(mesh,f0);
    M_removeFace(mesh,f1);
    M_removeEdge(mesh,edge);

    edge = 0;
    conf = -1;

    return 1;

  }

  // -------------------------------------------------------------------
  int edgeSwapOp::apply3D()
  {
    pPList eregs = E_regions(edge);   // cavity before swap
    pGRegion greg = R_whatIn( (pRegion)PList_item(eregs, 0) );

    pPList newRegions = PList_new();

    if( conf < 0 ) {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                  "Edge swap configuration not defined, call evaluate() first");
    }

    // create new regions for the selected swap configuration
    pVertex v1[4], v2[4];
    for( int i=0; i < swap_config.nb_tri_triangulation(); i++ )
      {
        int t = swap_config.triangulation(conf,i);
        for( int j=0; j<3; j++ ) {
          v1[j] = vertices[ swap_config.triangle(t,j) ];
        }
        v1[3] = vt;
        v2[0] = v1[0]; v2[1] = v1[2]; v2[2] = v1[1]; v2[3] = vb;

        //create new regions
        pRegion r1 = M_createR( mesh, 4, v1, (pGEntity)greg );
        pRegion r2 = M_createR( mesh, 4, v2, (pGEntity)greg );
        PList_append( newRegions, (pEntity)r1 );
        PList_append( newRegions, (pEntity)r2 );
      }

    // If swap on a boundary, classify new boundary entities (1 edge, 2 faces)
    if( E_whatInType(edge) == 2 )
      {
        pGEntity bgent = E_whatIn(edge);

        // the new edge to be classified is supposed to be between the first 
        // and the last vertices of the configuration (we started to list nodes
        // at a classified face in checkGeometry)
        pEdge bEdge = E_exist(vertices[0],vertices[swap_config.get()-1]);
        E_setWhatIn(bEdge,bgent);

        // now classify the faces
        pPList eFaces = E_faces(bEdge);
        void * temp = NULL;
        pFace face;
        pVertex oppV;
        int count = 0;
        while ( ( face = (pFace)PList_next(eFaces,&temp) ) ) {
          oppV = F_edOpVt(face,bEdge);
          if ( oppV == vt || oppV == vb ) {
            assert ( F_numRegions(face) == 1 );
            F_setWhatIn(face,bgent);
            count++;
          }
        }
        PList_delete(eFaces);
        assert(count == 2);
      }

    // call callback
    CallBackManagerSgl::instance().callCallBacks(eregs, newRegions, MAd_ESWAP, (pEntity)edge );
    PList_delete(newRegions);

    // remove old regions connected to edge
    pRegion frgn;
    void   *temp = 0;
    while( ( frgn = (pRegion)PList_next(eregs, &temp) ) )
      M_removeRegion(mesh,frgn);
    PList_delete(eregs);

    // remove old faces connected to edge
    pFace face;
    pPList faces = E_faces(edge);
    void   *temp1 = 0;
    while( ( face = (pFace)PList_next(faces, &temp1) ) )
      M_removeFace(mesh, face);
    PList_delete(faces);

    // remove edge
    M_removeEdge(mesh, edge);

    edge = 0;
    conf = -1;

    return 1;
  }

  // -------------------------------------------------------------------
  void edgeSwapOp::getCavity(pPList * cavity) const
  {
    if ( dim == 3 )   *cavity = E_regions(edge);
    else              *cavity = E_faces(edge);
  }

  // -------------------------------------------------------------------
  void edgeSwapOp::reportSwap( )
  {
    stringstream ss;
    std::string idStr;  ss << reportId;  ss >> idStr;
    std::string name = reportPrefix + "swap" + idStr + ".pos";
    pPList cavity = E_regions( edge );
    printPosEntities(cavity,name.c_str(),OD_CONSTANT, 0 );
    PList_delete(cavity);
    reportId++;
  }

  // -------------------------------------------------------------------
}
