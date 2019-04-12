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

#include "CheckMesh.h"
#include "MAdMessage.h"

#include <stdio.h>

namespace MAd {

  // -------------------------------------------------------------------
  // checkGeometricalCompatibility:

  // *** For 2D meshes ***

  // --- Faces ---
  //   - check every face is classified on a surface

  // --- Edges ---
  // For every edge:
  //   - if is classified on a surface: the edge is used by exactly 
  //     2 faces with the same classification as the edge,
  //   - if classified on a line: the edge is used by exactly 1 face.
  //   - is not classified on a point
  
  // --- Vertices ---
  // For every vertex:
  //   - if classified on a surface: all edges using it have to be 
  //     classified on the same surface,
  //   - if classified on a line: exactly two edges classified on the 
  //     same line use it. No other edge classified on a line use it.
  //   - if classified on a point, at least two edges classified 
  //     on different lines use it.

  // *** For 3D meshes ***

  // --- Regions ---
  //  - check every mesh region is classified on a model region

  // --- Faces ---
  // For every face:
  //   - if classified on a model region: the face is used by exactly 
  //     2 mesh regions with the same classification as the face,
  //   - if classified on a surface: the face is used by exactly 1 
  //     mesh region.
  //   - is not classified on a line
  //   - is not classified on a point

  // --- Edges ---
  // For every edge:
  //   - if classified on a model region: check that all faces using 
  //     it are classified on the same model region,
  //   - if classified on a surface: check that exactly 2 faces using it 
  //     are classified on the same surface, and that no face is classified 
  //     on a different surface,
  //   - if classified on a line: check that at least two faces using it 
  //     are classified on surfaces (it can be the same surface: revoluted 
  //     cylinder).
  //   - is not classified on a point

  // --- Vertices ---
  // For every vertex:
  //   - if classified on a model region: all edges using it have to be 
  //     classified on the same model region,
  //   - if classified on a surface: at least 3 edges using it are 
  //     classified on the same surface, and no edge can be classified 
  //     on a line or on a different surface.
  //   - if classified on a line: exactly 2 edges using it are 
  //     classified on the same line and no other edge is classified 
  //     on a line,
  //   - if classified on a point: at least 1 edge using it is 
  //     classified on a line (only 1: for revoluted cones).
  
  int checkGeomCompatibility(MDB_Mesh * mesh, int verbose, std::ostream& out)
  {
#ifdef PARALLEL
    // In a mesh partitioned with MDB, the parallel boundary 
    // faces are not classified on a model face and the following 
    // check fails
    return 1;
#endif

// #warning "check geometric compatibility not implemented for 2D meshes"
    if ( M_dim(mesh) < 3 ) return 1;

    int gDim;

    int flag = 1;

    // --- Regions ---

    pRegion pr;
    RIter rIt= M_regionIter(mesh);
    while ( ( pr = RIter_next(rIt) ) )
      {
        gDim = R_whatInType(pr);
        if ( gDim != 3 ) {
          flag = 0;
          if (verbose) {
            out <<  "Mesh region not classified on a model region\n";
            R_info_topology(pr,out);
          }
          break;
        }
      }
    RIter_delete(rIt);
    if ( !flag ) return 0;

    // -- Faces ---

    pFace pf;
    FIter fIt= M_faceIter(mesh);
    while ( flag && ( pf = FIter_next(fIt) ) )
      {
        gDim = F_whatInType(pf);
        switch (gDim) {
        case 3: {
          int nRgn = F_numRegions(pf); 
          if ( nRgn != 2 ) {
            flag = 0;
            if (verbose) {
              out << "Face classif on dim 3 with not exactly 2 regions\n";
              F_info(pf,"",out);
            }
          }
          if ( ( F_whatIn(pf) != (pGEntity)R_whatIn( F_region(pf,0) ) ) || 
               ( F_whatIn(pf) != (pGEntity)R_whatIn( F_region(pf,1) ) ) ) {
            flag = 0;
            if (verbose) {
              out << "Face classif on dim 3 used by regions classif on different entities\n";
              F_info(pf,"",out);
            }
          }
          break;
        }
        case 2: {
          int nRgn = F_numRegions(pf); 
          if ( nRgn != 1 ) {
            flag = 0;
            if (verbose) {
              out << "Face classif on dim 2 with not exactly 1 region\n";
              F_info(pf,"",out);
            }
          }
          break;
        }
        default: {
          flag = 0;
          if (verbose) {
            out << "Face classif on wrong dimension\n";
            F_info(pf,"",out);
          }
        }
        }
      }
    FIter_delete(fIt);
    if ( !flag ) return 0;

    // -- Edges ---

    pEdge pe;
    pGEntity pge;
    pPList eFaces;
    void * temp;
    EIter eIt = M_edgeIter(mesh);
    while ( flag && ( pe = EIter_next(eIt) ) )
      {
        pge  = E_whatIn(pe);
        gDim = GEN_type(pge);
        switch (gDim) {
        case 3: {
          eFaces = E_faces(pe);
          temp = NULL;
          while ( ( pf = (pFace) PList_next(eFaces,&temp) ) ) {
            if ( F_whatIn(pf) != pge ) {
              flag = 0;
              if (verbose) {
                out << "Edge classif on a region but used by faces "
                    << "not classif on the same region\n";
                E_info(pe,"",out);
              }
            }
          }
          PList_delete(eFaces);
          break;
        }
        case 2: {
          int nbBF = 0;
          eFaces = E_faces(pe);
          temp = NULL;
          while ( ( pf = (pFace) PList_next(eFaces,&temp) ) ) {
            if ( F_whatInType(pf) == 2 ) {
              if ( F_whatIn(pf) != pge ) {
                flag = 0;
                if (verbose) {
                  out << "Edge classif on a surface used by a face classif "
                      << "on another surface with tag " << GEN_tag(F_whatIn(pf)) <<"\n";
                  E_info(pe,"",out);
                }
              }
              nbBF++;
            }
          }
          PList_delete(eFaces);
          if ( nbBF != 2 ) {
            flag = 0;
            if (verbose) {
              out << "Edge classif on a surface and not used by exactly "
                  << "2 faces classif on surfaces\n";
              E_info(pe,"",out);
            }
          }
          break;
        }
        case 1: {
          int nbBF = 0;
          eFaces = E_faces(pe);
          temp = NULL;
          while ( ( pf = (pFace) PList_next(eFaces,&temp) ) ) {
            if ( F_whatInType(pf) == 2 ) nbBF++;
          }
          PList_delete(eFaces);
          if ( nbBF < 2 ) {
            flag = 0;
            if (verbose) {
              out << "Edge classif on a line and not used by at least"
                  << " 2 faces classif on surfaces\n";
              E_info(pe,"",out);
            }
          }
          break;
        }
        default: {
          flag = 0;
          if (verbose) {
            out << "Edge classif on wrong dimension\n";
            E_info(pe,"",out);
          }
        }
        }
      }
    EIter_delete(eIt);
    if ( !flag ) return 0;

    // -- Vertices ---

    pPList vEdges;
    pVertex pv;
    VIter vIt = M_vertexIter(mesh);
    while ( flag && ( pv = VIter_next(vIt) ) )
      {
        pge  = V_whatIn(pv);
        gDim = GEN_type(pge);
        switch (gDim) {
        case 3: {
          vEdges = V_edges(pv);
          temp = NULL;
          while ( ( pe = (pEdge) PList_next(vEdges,&temp) ) ) {
            if ( E_whatIn(pe) != pge ) {
              flag = 0;
              if (verbose) {
                out << "Vertex classif on a region but used by edges "
                    << "not classif on the same region\n";
                V_info(pv,"",out);
              }
            }
          }
          PList_delete(vEdges);
          break;
        }
        case 2: {
          int nbE = 0;
          vEdges = V_edges(pv);
          temp = NULL;
          while ( ( pe = (pEdge) PList_next(vEdges,&temp) ) ) {
            if ( E_whatInType(pe) == 3 ) continue;
            if ( E_whatIn(pe) != pge ) {
              flag = 0;
              if (verbose) {
                out << "Vertex classif on a surface but used by edges "
                    << "not classif on the same surface or on a region\n";
                V_info(pv,"",out);
              }
            }
            nbE++;
          }
          PList_delete(vEdges);
          if ( nbE < 3 ) {
            flag = 0;
            if (verbose) {
              out << "Vertex classif on a surface and not used by at least"
                  << " 3 edges classif on surfaces\n";
              V_info(pv,"",out);
            }
          }
          break;
        }
        case 1: {
          int nbE = 0;
          vEdges = V_edges(pv);
          temp = NULL;
          while ( ( pe = (pEdge) PList_next(vEdges,&temp) ) ) {
            if ( E_whatInType(pe) > 1 ) continue;
            if ( E_whatIn(pe) != pge ) {
              flag = 0;
              if (verbose) {
                out << "Vertex classif on a line and used by an edge "
                    << "classif on another line\n";
                V_info(pv,"",out);
              }
            }
            nbE++;
          }
          PList_delete(vEdges);
          if ( nbE != 2 ) {
            flag = 0;
            if (verbose) {
              out << "Vertex classif on a line and not used by exactly"
                  << " 2 edges classif on the same line\n";
              V_info(pv,"",out);
            }
          }
          break;
        }
        case 0: {
          int nbE = 0;
          std::set<pGEntity> glines;
          vEdges = V_edges(pv);
          temp = NULL;
          while ( ( pe = (pEdge) PList_next(vEdges,&temp) ) ) {
            if ( E_whatInType(pe) > 1 ) continue;
            nbE++;
            glines.insert(E_whatIn(pe));
          }
          PList_delete(vEdges);
          if ( nbE == 0 ) {
            flag = 0;
            if (verbose) {
              out << "Vertex classif on a point and not used by an"
                  << " edge classif on a line\n";
              V_info(pv,"",out);
            }
          }
//           if ( nbE < 2 ) {
//             flag = 0;
//             if (verbose) {
//               out << "Vertex classif on a point and used by less than"
//                   << " 2 edges classif on lines\n";
//               V_info(pv,"",out);
//             }
//           }
//           if ( glines.size() < 2 ) {
//             flag = 0;
//             if (verbose) {
//               out << "Vertex classif on a point and not used by edges"
//                   << " classified on different lines\n";
//               V_info(pv,"",out);
//             }
//           }
          break;
        }
        default: {
          flag = 0;
          if (verbose) {
            out << "Vertex classif on non-existing dimension\n";
            V_info(pv,"",out);
          }
        }
        }
      }
    VIter_delete(vIt);
    if ( !flag ) return 0;

    return flag;
    
  }

  // -------------------------------------------------------------------
  int checkRegionsVolume(MDB_Mesh * mesh, int verbose, std::ostream& out)
  {
    int flag = 1;
    pRegion r;
    RIter rIt= M_regionIter(mesh);
    while ( ( r = RIter_next(rIt) ) ) {
      if (R_volume(r) < 0.) {
        flag = 0;
        if (verbose) {
          out << "Negative volume found !\n";
          R_info_topology(r,out);
        }
        break;
      }  
    }
    RIter_delete(rIt);
    return flag;
  }

  // -------------------------------------------------------------------
  int checkEdgeToRegionConnectivity(MDB_Mesh * mesh, int verbose, 
                                    std::ostream& out)
  {
#ifdef PARALLEL
    // In a mesh partitioned with MDB, the parallel boundary 
    // faces are not classified on a model face and the following 
    // check fails
    return 1;
#endif

    if ( M_dim(mesh) <= 2 ) return 1;

    int flag = 1;
    pEdge edge;
    EIter eIt = M_edgeIter(mesh);
    while ( ( edge = EIter_next(eIt) ) ) {

      if ( E_whatInType(edge) != 3 ) continue;

      pFace face = E_face(edge,0);
      pRegion start_region = F_region(face,0);
      
      pFace current_face = face;
      pRegion current_region = start_region;
      pFace next_face;
      pRegion next_region = NULL;
      while ( next_region != start_region )
      {
        next_face = E_otherFace(edge, current_face, current_region);
        next_region = F_region(next_face, 0);
        if( next_region == current_region ) {
          next_region = F_region(next_face, 1);
        }
        current_face   = next_face;
        current_region = next_region;
        
        if ( !next_region ) {
          flag = 0;
          if (verbose) {
            out << "Found a wrong edge to region connectivity !\n";
            E_info(edge,"",out);
          }
          break;
        }
      }
    }
    EIter_delete(eIt);
    return flag;
  }

  // -------------------------------------------------------------------
  // Obsolete: included in checkGeomCompatibility()
  int checkFaceToRegionConnectivity(MDB_Mesh * mesh, int verbose, std::ostream& out) {

#ifdef PARALLEL
    // In a mesh partitioned with MDB, the parallel boundary 
    // faces are not classified on a model face and the following 
    // check fails
    return 1;
#endif

    if ( M_dim(mesh) <= 2 ) return 1;

    int flag = 1;
    pFace f;
    FIter fIt= M_faceIter(mesh);
    while ( ( f = FIter_next(fIt) ) ) {
      int ok = 1;
      int nRgn = F_numRegions(f); 
      if (nRgn==0) out << "nRgn: " << nRgn << "\n";
      int type = F_whatInType(f);
      if (              nRgn == 0 ) ok = 0;
      if ( type == 2 && nRgn != 1 ) ok = 0;
      if ( type == 3 && nRgn != 2 ) ok = 0;
      if (!ok) {
        flag = 0;
        if (verbose) {
          out << "Found a wrong face to region connectivity !\n";
          F_info(f,"",out);
        }
        break;
      }  
    }
    FIter_delete(fIt);
    return flag;
  }

  // -------------------------------------------------------------------
  int checkFaceToVertexConnectivity(MDB_Mesh * mesh, int verbose, std::ostream& out)
  {
    int flag = 1;
    pFace f;
    pVertex v0,v1,v2;
    FIter fIt = M_faceIter(mesh);
    while ( ( f = FIter_next(fIt) ) ) {
      v0 = F_vertex(f,0);
      v1 = F_vertex(f,1);
      v2 = F_vertex(f,2);
      if ( v0 == v1 || v0 == v2 || v1 == v2 ) {
        flag = 0;
        if (verbose) {
          out << "Found a wrong face to vertex connectivity !\n";
          F_info(f,"",out);
        }
        break;
      }
    }
    FIter_delete(fIt);
    return flag;
  }

  // -------------------------------------------------------------------
  // Check that there is at least one face attached to each edge in 2D, 
  // one region and two faces in 3D.
  int checkEdgeConnectivity(MDB_Mesh * mesh, int verbose, std::ostream& out) {

#ifdef PARALLEL
    // In a mesh partitioned with MDB, the parallel boundary 
    // edges are not classified on a model line and the following 
    // check fails
    return 1;
#endif

    int dim = M_dim(mesh);
    if ( dim <= 1 ) return 1;

    int flag = 1;
    pEdge e;
    EIter eIt= M_edgeIter(mesh);
    while ( ( e = EIter_next(eIt) ) ) {
      int nFace = E_numFaces(e);
      int nRgn  = E_numRegions(e);
      if ( dim == 2 ) {
        if ( (nRgn != 0) || (nFace < 1) ) {
          flag = 0; break;
        }
      }
      else if ( dim == 3 ) {
        if ( (nRgn < 1) || (nFace < 2) ) {
          flag = 0; break;
        }
      }
    }
    EIter_delete(eIt);
    if ( (!flag) && verbose ) {
      out << "Found an ill-connected edge !\n";
      E_info(e,"",out);
    }
    return flag;
  }

  // -------------------------------------------------------------------
  int checkEntityPointers(MDB_Mesh * mesh, int verbose, std::ostream& out) {

    int flag = 1;

    VIter vIt = M_vertexIter(mesh);
    pVertex pv;
    while ( ( pv = VIter_next(vIt) ) ) {
      if ( !pv ) {
        flag = 0;
        if (verbose) out << "Found a null pointer when iterating on vertices\n";
      } 
    }
    VIter_delete(vIt);

    EIter eIt = M_edgeIter(mesh);
    pEdge pe;
    while ( ( pe = EIter_next(eIt) ) ) {
      if ( !pe ) {
        flag = 0;
        if (verbose) out << "Found a null pointer when iterating on edges\n";
      } 
    }
    EIter_delete(eIt);

    FIter fIt = M_faceIter(mesh);
    pFace pf;
    while ( ( pf = FIter_next(fIt) ) ) {
      if ( !pf ) {
        flag = 0;
        if (verbose) out << "Found a null pointer when iterating on faces\n";
      } 
    }
    FIter_delete(fIt);

    RIter rIt = M_regionIter(mesh);
    pRegion pr;
    while ( ( pr = RIter_next(rIt) ) ) {
      if ( !pr ) {
        flag = 0;
        if (verbose) out << "Found a null pointer when iterating on regions\n";
      } 
    }
    RIter_delete(rIt);
  
    return flag;
  }

  // -------------------------------------------------------------------
  int checkIterators(MDB_Mesh * mesh, int verbose, std::ostream& out) {

    int flag = 1;

    int nbMeshV = M_numVertices(mesh);
    int nbIterV = 0;
    VIter vIt= M_vertexIter(mesh);
    while ( VIter_next(vIt) )  nbIterV++;
    VIter_delete(vIt);
    if ( nbMeshV != nbIterV ) {
      flag = 0;
      if (verbose) {
        out << "Incoherent number of vertices: " 
            << nbMeshV << " in mesh, " 
            << nbIterV << " in iterator\n";
      }
    }

    int nbMeshE = M_numEdges(mesh);
    int nbIterE = 0;
    EIter eIt= M_edgeIter(mesh);
    while ( EIter_next(eIt) )  nbIterE++;
    EIter_delete(eIt);
    if ( nbMeshE != nbIterE ) {
      flag = 0;
      if (verbose) {
        out << "Incoherent number of edges: "
            << nbMeshE << " in mesh, "
            << nbIterE << " in iterator\n";
      }
    }

    int nbMeshF = M_numFaces(mesh);
    int nbIterF = 0;
    FIter fIt= M_faceIter(mesh);
    while ( FIter_next(fIt) )  nbIterF++;
    FIter_delete(fIt);
    if ( nbMeshF != nbIterF ) {
      flag = 0;
      if (verbose) {
        out << "Incoherent number of faces: "
            << nbMeshF << " in mesh, "
            << nbIterF << " in iterator\n";
      }
    }

    int nbMeshR = M_numRegions(mesh);
    int nbIterR = 0;
    RIter rIt= M_regionIter(mesh);
    while ( RIter_next(rIt) )  nbIterR++;
    RIter_delete(rIt);
    if ( nbMeshR != nbIterR ) {
      flag = 0;
      if (verbose) {
        out << "Incoherent number of regions: "
            << nbMeshR << " in mesh, "
            << nbIterR << " in iterator\n";
      }
    }

    return flag;
  }

  // -------------------------------------------------------------------
  int checkParameters(MDB_Mesh * mesh, int verbose, std::ostream& out)
  {
    if ( !M_isParametric(mesh) ) return 1;

    int flag = 1;

    pVertex pv;
    double tmp0,tmp1;
    VIter vIt= M_vertexIter(mesh);
    while ( ( pv = VIter_next(vIt) ) ) {
      int gDim = V_whatInType(pv);
      bool param = V_params(pv,&tmp0,&tmp1);
      if ( ( gDim==0 || gDim==1 || gDim==2 ) && !param ) {
        flag = 0;
        if (verbose) {
          out << "Non parametric point classified on geo entity with dim " <<gDim<<"\n";
          V_info(pv,"",out);
        }
      }
      if ( ( gDim==3 ) && param ) {
        flag = 0;
        if (verbose) {
          out << "Parametric point classified on geo entity with dim " <<gDim<<"\n";
          V_info(pv,"",out);
        }
      }
    }
    VIter_delete(vIt);

    return flag;
  }

  // -------------------------------------------------------------------
  bool checkMesh(MDB_Mesh * mesh, checkType type, int verbose, 
                 std::ostream& out, MeshStatus * status) {

    switch (type) {
    case CHECK_ALL: {
      if ( !checkRegionsVolume(mesh, verbose, out) ) {
        if(status) *status = NEGATIVE_VOLUME; 
        return false;
      }
      if ( !checkGeomCompatibility(mesh, verbose, out) ) {
        if(status) *status = GEOM_INCOMPATIBILITY;
        return false;
      }
      if ( !checkEdgeToRegionConnectivity(mesh, verbose, out) ) {
        if(status) *status = WRONG_EDGE_TO_RGN_CONN;
        return false;
      }
      if ( !checkFaceToRegionConnectivity(mesh, verbose, out) ) {
        if(status) *status = WRONG_FACE_TO_RGN_CONN;
        return false;
      }
      if ( !checkFaceToVertexConnectivity(mesh, verbose, out) ) {
        if(status) *status = WRONG_FACE_TO_VTX_CONN;
        return false;
      }
      if ( !checkEdgeConnectivity(mesh, verbose, out) ) {
        if(status) *status = WRONG_EDGE_CONN;
        return false;
      }
      if ( !checkEntityPointers(mesh, verbose, out) ) {
        if(status) *status = WRONG_ENTITY_POINTERS;
        return false;
      }
      if ( !checkIterators(mesh, verbose, out) ) {
        if(status) *status = WRONG_ITERATORS;
        return false;
      }
      if ( !checkParameters(mesh, verbose, out) ) {
        if(status) *status = WRONG_PARAMETERS;
        return false;
      }
      break;
    }
    case CHECK_VOLUME: {
      if ( !checkRegionsVolume(mesh, verbose, out) ) {
        if(status) *status = NEGATIVE_VOLUME;
        return false;
      }
      break;
    }
    case CHECK_GEOM_COMPATIBILITY: {
      if ( !checkGeomCompatibility(mesh, verbose, out) ) {
        if(status) *status = GEOM_INCOMPATIBILITY;
        return false;
      }
      break;
    }
    case CHECK_EDGE_TO_RGN_CONN: {
      if ( !checkEdgeToRegionConnectivity(mesh, verbose, out) ) {
        if(status) *status = WRONG_EDGE_TO_RGN_CONN;
        return false;
      }
      break;
    }
    case CHECK_FACE_TO_RGN_CONN: {
      if ( !checkFaceToRegionConnectivity(mesh, verbose, out) ) {
        if(status) *status = WRONG_FACE_TO_RGN_CONN;
        return false;
      }
      break;
    }
    case CHECK_FACE_TO_VTX_CONN: {
      if ( !checkFaceToVertexConnectivity(mesh, verbose, out) ) {
        if(status) *status = WRONG_FACE_TO_VTX_CONN;
        return false;
      }
      break;
    }
    case CHECK_EDGE_CONN: {
      if ( !checkEdgeConnectivity(mesh, verbose, out) ) {
        if(status) *status = WRONG_EDGE_CONN;
        return false;
      }
      break;
    }
    case CHECK_ENTITY_POINTERS: {
      if ( !checkEntityPointers(mesh, verbose, out) ) {
        if(status) *status = WRONG_ENTITY_POINTERS;
        return false;
      }
      break;
    }
    case CHECK_ITERATORS: {
      if ( !checkIterators(mesh, verbose, out) ) {
        if(status) *status = WRONG_ITERATORS;
        return false;
      }
      break;
    }
    case CHECK_PARAMETERS: {
      if ( !checkParameters(mesh, verbose, out) ) {
        if(status) *status = WRONG_PARAMETERS;
        return false;
      }
      break;
    }
    default: {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                  "Not a valid checkType: %d",type);
      return false;
    }
    }

    if(status) *status = VALID;
    return true;
  }

  // -------------------------------------------------------------------

}
