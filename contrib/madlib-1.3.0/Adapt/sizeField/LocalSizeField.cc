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

#include "LocalSizeField.h"
#include "IsoMeshSize.h"
#include "AnisoMeshSize.h"
#include "MathUtils.h"
#include "MAdMessage.h"
#include "MeshParametersManager.h"
#include "MAdTimeManager.h"
#include "MAdResourceManager.h"
#include "DistanceFunction.h"
#include "MAdStringFieldEvaluator.h"

#include <iostream>
#include <stdlib.h>
using std::cerr;
using std::cout;
using std::endl;
using std::set;
using std::string;

namespace MAd {

  // -------------------------------------------------------------------
  LocalSizeField::LocalSizeField(pMesh m, string name, bool _distToFaces):
    SizeFieldBase(name), mesh(m), geoDim(-1),
    isotropic(true), radius(-1.), sizeN(""), sizeT(""),
    sEvalN(NULL), sEvalT(NULL), 
    distToFaces(_distToFaces), dFct(mesh, distToFaces), 
    limit(false), tgSizeLimitCoef(MAdBIG), maxCurv(MAdBIG)
  {
#ifdef PARALLEL
    MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                "Local size fields not supported in parallel");
#endif
  }
  
  // -------------------------------------------------------------------
  LocalSizeField::LocalSizeField(const LocalSizeField& _lsf): dFct(NULL,false)
  {
    throw;
  }
  
  // -------------------------------------------------------------------
  LocalSizeField::~LocalSizeField()  
  {
    if ( sEvalN ) delete sEvalN;
    if ( sEvalT ) delete sEvalT;
  }

  // -------------------------------------------------------------------
  // -------------------------------------------------------------------
  void LocalSizeField::addGeometricEntity(int type, int tag)
  {
    // check that the new entity has the same dimension as the previous ones
    if ( geoDim < 0 ) geoDim = type;
    else if ( geoDim != type ) {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                  "Trying to insert a geo entity with dim %d while another entity with dim %d was previously inserted",
                                  type, geoDim);
    }

    geomEntities.insert( GM_entityByTag(M_model(mesh),type,tag) );
  }

  // -------------------------------------------------------------------
  void LocalSizeField::setIsoSize(double _radius, string _sizeN)
  {
    isotropic = true;
    radius = _radius;
    sizeN = _sizeN;
    sizeT = "";
    if( sEvalN ) delete sEvalN;
    sEvalN = new MAdStringFieldEvaluator(1,sizeN.c_str());
  }

  // -------------------------------------------------------------------
  void LocalSizeField::setAnisoSize(double _radius, string _sizeN, 
                                    string _sizeT)
  {
    isotropic = false;
    radius = _radius;
    sizeN = _sizeN;
    sizeT = _sizeT;
    if( sEvalN ) delete sEvalN;
    if( sEvalT ) delete sEvalT;
    sEvalN = new MAdStringFieldEvaluator(1,sizeN.c_str());
    sEvalT = new MAdStringFieldEvaluator(1,sizeT.c_str());
  }

  // -------------------------------------------------------------------
  void LocalSizeField::setCurvatureLimiter(double onTgSize, double _maxCurv)
  {
    limit = true;
    tgSizeLimitCoef = onTgSize;
    maxCurv = _maxCurv;
  }

  // -------------------------------------------------------------------
  pMSize LocalSizeField::getSize(const pVertex pv) const
  {
    // get the distance
    double dist = dFct.getDistance(pv);

    double dxyz[3] = {0., 0., 0.};  dxyz[0] = dist;

    // get the time
    double time = MAdTimeManagerSgl::instance().getTime();

    // get the sizes for this distance at this time
    double sizeN_d, sizeT_d;
    sEvalN->eval(dxyz,time,&sizeN_d);
    if ( !isotropic ) sEvalT->eval(dxyz,time,&sizeT_d);

    pMSize theSize;

    if (dist >= radius) {
      theSize = new IsoMeshSize (MeshParametersManagerSgl::instance().getBigLength());
    }
    else {
      if ( isotropic ) {
        theSize = new IsoMeshSize( sizeN_d );
      }
      else {
        double nor[3];
        if ( !dFct.getGradient(pv,nor) ) {
          MAdMsgSgl::instance().error(__LINE__,__FILE__,"Gradient not available for vertex %p",pv);
        }

        if ( fabs(nor[0]) <= MAdTOL && fabs(nor[1]) <= MAdTOL && fabs(nor[2]) <= MAdTOL ) {
          theSize = new IsoMeshSize( sizeN_d );
        }
        else {
          if ( limit )
            {
              double curvR;
              if ( !dFct.getCurvature(pv,&curvR) ) {
                MAdMsgSgl::instance().error(__LINE__,__FILE__,"Curvature not available for vertex %p",pv);
              }
              if ( curvR > MAdTOL ) {
                curvR = std::min(curvR,maxCurv);
                curvR = 1. / curvR;
                // limit the tangent size regarding the curvature radius of the wall
                sizeT_d = std::min(sizeT_d, tgSizeLimitCoef * curvR);
                sizeT_d = std::max(sizeT_d,sizeN_d);
              }
            }

          theSize = new AnisoMeshSize( nor, sizeN_d, sizeT_d );
        }
      }
    }
    return theSize;
  }

  // -------------------------------------------------------------------
  pMSize LocalSizeField::getSizeOnEntity(const pEntity pe, 
                                         const double xyz[3]) const
  {
    // get the distance
    double dist = dFct.computeDistance(xyz);

    double dxyz[3] = {0., 0., 0.};  dxyz[0] = dist;

    // get the time
    double time = MAdTimeManagerSgl::instance().getTime();

    // get the sizes for this distance at this time
    double sizeN_d, sizeT_d;
    sEvalN->eval(dxyz,time,&sizeN_d);
    if ( !isotropic ) sEvalT->eval(dxyz,time,&sizeT_d);

    pMSize theSize;

    if (dist > radius) {
      theSize = new IsoMeshSize (MeshParametersManagerSgl::instance().getBigLength());
    }
    else {
      if ( isotropic ) {
        theSize = new IsoMeshSize( sizeN_d );
      }
      else {
        double nor[3];
        if ( !dFct.getGradientOnEntity(pe,xyz,nor) ) {
          MAdMsgSgl::instance().error(__LINE__,__FILE__,"Gradient not available for entity %p",pe);
        }

        if ( fabs(nor[0]) <= MAdTOL && fabs(nor[1]) <= MAdTOL && fabs(nor[2]) <= MAdTOL ) {
          theSize = new IsoMeshSize( sizeN_d );
        }
        else {
          if ( limit )
            {
              double curvR;
              if ( !dFct.getCurvatureOnEntity(pe,xyz,&curvR) ) {
                MAdMsgSgl::instance().error(__LINE__,__FILE__,"Curvature not available for entity %p",pe);
              }
              if ( curvR > MAdTOL ) {
                curvR = std::min(curvR,maxCurv);
                curvR = 1. / curvR;
                // limit the tangent size regarding the curvature radius of the wall
                sizeT_d = std::min(sizeT_d, tgSizeLimitCoef * curvR);
                sizeT_d = std::max(sizeT_d,sizeN_d);
              }
            }

          theSize = new AnisoMeshSize( nor, sizeN_d, sizeT_d );
        }
      }
    }
    return theSize;
  }

  // -------------------------------------------------------------------
  // Length squared computation
  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  double LocalSizeField::SF_VV_lengthSq(const pVertex pV0, const pVertex pV1) const
  {
    double xyz[2][3];
    V_coord(pV0,xyz[0]);
    V_coord(pV1,xyz[1]);

    pMSize pS[2];
    pS[0] = getSize(pV0);
    pS[1] = getSize(pV1);

    double lSq = SF_XYZ_lengthSq(xyz[0],xyz[1],pS[0],pS[1]);
  
    if ( pS[0] ) delete pS[0];
    if ( pS[1] ) delete pS[1];

    return lSq;
  }

  // -------------------------------------------------------------------
  double LocalSizeField::SF_XYZ_lengthSq(const double xyz0[3], 
                                         const double xyz1[3], 
                                         const pMSize pS0, 
                                         const pMSize pS1) const
  {
    if( pS0 )
      {
        double e[3];
        diffVec(xyz0,xyz1,e);
        double lenSq0 = pS0->normSq(e);
        if ( pS1 )
          {
            double lenSq1 = pS1->normSq(e);
            return sqrt(lenSq0*lenSq1);
          }
        else return lenSq0;
      }
    else {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,"No size defined");
    }
    return 0.;
  }

  // -------------------------------------------------------------------
  // Area squared computation
  // -------------------------------------------------------------------

  double LocalSizeField::SF_F_areaSq(const pFace face) const
  {
    double area = 0.;

    double xyz[3][3];
    F_coordP1(face,xyz);
  
    void * temp = 0;
    pPList fVerts = F_vertices(face,1);
    while( pVertex pV = (pVertex)PList_next(fVerts,&temp) )
      {
        pMSize pS = getSize(pV);
        area += SF_XYZ_areaSq(xyz,pS,0);
        if (pS) delete pS;
      }
    PList_delete(fVerts);
  
    area /= F_numVertices(face);

    return area;
  }

  // -------------------------------------------------------------------
  double LocalSizeField::SF_XYZ_areaSq(const double fxyz[3][3], 
                                       const pMSize pS, 
                                       const double norDir[3]) const
  {
    if( !pS ) {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,"No size defined");
    }

    // get the two first edges
    double e01[3],e02[3];
    diffVec(fxyz[1],fxyz[0],e01);
    diffVec(fxyz[2],fxyz[0],e02);
    
    double nor[3];
    crossProd(e01,e02,nor);

    double l1SqInv = 1. / pS->lengthSqInDir(e01);
    double l2SqInv = 1. / pS->lengthSqInDir(e02);

    if( norDir && dotProd(norDir,nor) < MAdTOL ) return 0.;

    double areaSq = 0.25 * dotProd(nor,nor) * l1SqInv * l2SqInv;
    if( areaSq < MAdTOL ) return 0.;

    return areaSq;
  }

  // -------------------------------------------------------------------
  // Volume computation
  // -------------------------------------------------------------------

  double LocalSizeField::SF_R_volume(const pRegion region) const
  {
    double vol = 0.;

    double xyz[4][3];
    R_coordP1(region,xyz);

    pPList rVerts = R_vertices(region);
    void * temp = 0;
    while( pVertex pV = (pVertex)PList_next(rVerts,&temp) )
      {
        pMSize pS = getSize(pV);
        vol += SF_XYZ_volume(xyz,pS);
        if (pS) delete pS;
      }
    PList_delete(rVerts);

    vol /= R_numVertices(region);

    return vol;
  }

  // -------------------------------------------------------------------
  double LocalSizeField::SF_XYZ_volume(const double xyz[4][3], const pMSize pS) const
  {
    if( !pS ) {
      printf("Error in LocalSizeField::volume: no size given\n");
      throw;
    }

    double physVol = R_XYZ_volume(xyz);
  
    return ( physVol / (pS->size(0)*pS->size(1)*pS->size(2)) );
  }

  // -------------------------------------------------------------------
  // Center of edge computation
  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  double LocalSizeField::SF_E_center(const pEdge edge, double center[3], 
                                     double * reducSq, pMSize * cSize) const
  {
    return SF_VV_center(E_vertex(edge,0),E_vertex(edge,1),center,reducSq,cSize);
  }

  // -------------------------------------------------------------------
  double LocalSizeField::SF_VV_center(const pVertex v0, const pVertex v1,
                                      double center[3], double * reducSq, 
                                      pMSize * cSize) const
  {
    double xyz[2][3];
    V_coord(v0,xyz[0]);
    V_coord(v1,xyz[1]);

    pMSize pS[2];
    pS[0] = getSize(v0);
    pS[1] = getSize(v1);

    double cParam = SF_XYZ_center(xyz,pS,center,reducSq,cSize);

    if ( pS[0] ) delete pS[0];
    if ( pS[1] ) delete pS[1];

    return cParam;
  }

  // -------------------------------------------------------------------
  bool LocalSizeField::getCurvature(const pVertex pv, double *c) const
  {
    return dFct.getCurvature(pv,c);
  }

  // -------------------------------------------------------------------
  double LocalSizeField::getDistance(const pVertex pv) const
  {
    return dFct.getDistance(pv);
  }

  // ----------------------------------------------------------------------
  // -------------------------------------------------------------------
  void LocalSizeField::scale(double fact)
  {
    printf("Not implemented (LocalSizeField::scale)\n");
    throw;  
  }

  // ----------------------------------------------------------------------
  // ----------------------------------------------------------------------
  void LocalSizeField::updateTree()
  {
#if 0
    dFct.computeAllDistancesEDP(geomEntities);
    dFct.outputDistance("distance.pos");

          // gather all regions with a vertex in the scope of the size field
          set<pRegion> allR;
          RIter rit = M_regionIter(mesh);
          while (pRegion pr = RIter_next(rit)) {
            allR.insert(pr);
          }
          RIter_delete(rit);

#warning "changed gradient computation ton constant gradients"
          dFct.computeGradientAndCurvature(allR);

//           // gather all faces 
//           set<pFace> allF;
//           FIter fit = M_faceIter(mesh);
//           while (pFace pf = FIter_next(fit)) {
//             allF.insert(pf);
//           }
//           FIter_delete(fit);
//     dFct.computeGradientAndCurvature2D(allF);
    dFct.outputGradAtVertices("grad.pos");
    if ( limit ) dFct.outputCurvature("curv.pos");
    exit(0);
#endif

    MAdResourceManager& tm = MAdResourceManagerSgl::instance();
    double t0 = tm.getTime();

    set<pVertex> vertices;
    set<pEntity> entities;
    collectEntities(&vertices, &entities);
    dFct.computeTree(vertices,entities);

    printf("Computed tree in %f sec\n",tm.getTime()-t0);

    if ( !isotropic )
      {
        // curvatures are not implemented in 2D
        if ( M_dim(mesh) < 3 ) dFct.computeGradientAtVertices();
        else {
          
// #warning "changed gradient computation ton constant gradients"
//           dFct.computeAllDistances();
          double t1 = tm.getTime();
          dFct.computeAllDistAndGrad();
          printf("Computed dist and grad in %f sec\n",tm.getTime()-t1);

// #warning "hack for exact dist"
//           dFct.computeAllDistances();

          double t2 = tm.getTime();

          // gather all regions with a vertex in the scope of the size field
          set<pRegion> allR;
          pPList vRegs;
          VIter vit = M_vertexIter(mesh);
          while (pVertex pv = VIter_next(vit)) {
            if ( dFct.getDistance(pv) <= radius ) {
              vRegs = V_regions(pv);
              for (int i=0; i< PList_size(vRegs); i++) allR.insert((pRegion)PList_item(vRegs,i));
              PList_delete(vRegs);
            }
          }
          VIter_delete(vit);

// #warning "changed gradient computation ton constant gradients"
//           dFct.computeGradientAndCurvature(allR);
          if ( limit ) {
            dFct.computeCurvature(allR);
// #warning "smoothing the curvature (3)"
//             dFct.limitCurvature(maxCurv);
//             dFct.smoothCurvature(50.);
          }
          printf("Computed curvatures in %f sec\n",tm.getTime()-t2);

// #warning "debug: output curvatures in volume"
//           double t3 = tm.getTime();
//           dFct.outputDistance("dist.pos");
//           dFct.outputGradAtVertices("grad.pos");
//           if ( limit ) dFct.outputCurvature("curv.pos");
//           printf("Outputs in %f sec\n",tm.getTime()-t3);
//           exit(0);
        }
      }
    printf("Updated tree in %f sec\n",tm.getTime()-t0);
  }

  // ----------------------------------------------------------------------
  void LocalSizeField::collectEntities(set<pVertex> * verts, 
                                       set<pEntity> * ents) const
  {
    verts->clear();
    ents->clear();

    set<pGEntity>::const_iterator it    = geomEntities.begin();
    set<pGEntity>::const_iterator itEnd = geomEntities.end();
    for (; it != itEnd; it++) 
      {
        switch ( GEN_type(*it) ) {
        case 0:
          {
            VIter vit = M_vertexIter(mesh);
            while (pVertex pv = VIter_next(vit))
              {
                if ( EN_whatIn((pEntity)pv) == *it ) {
                  ents->insert((pEntity)pv);
                  verts->insert(pv);
                }
              }
            VIter_delete(vit);
            break;
          }
        case 1:
          {
            EIter eit = M_edgeIter(mesh);
            while (pEdge pe = EIter_next(eit))
              {
                if ( EN_whatIn((pEntity)pe) == *it )
                  {
                    ents->insert((pEntity)pe);
                    verts->insert(E_vertex (pe,0));
                    verts->insert(E_vertex (pe,1));
                  }
              }
            EIter_delete(eit);
          }
        case 2:
          {
            FIter fit = M_faceIter(mesh);
            while (pFace pf = FIter_next(fit))
              {
                if ( EN_whatIn((pEntity)pf) == *it )
                  {
                    ents->insert((pEntity)pf);
                    verts->insert(F_vertex (pf,0));
                    verts->insert(F_vertex (pf,1));
                    verts->insert(F_vertex (pf,2));
                  }
              }
            FIter_delete(fit);
          }
        case 3:
          {
            RIter rit = M_regionIter(mesh);
            while (pRegion pr = RIter_next(rit))
              {
                if ( EN_whatIn((pEntity)pr) == *it )
                  {
                    ents->insert((pEntity)pr);
                    verts->insert(R_vertex (pr,0));
                    verts->insert(R_vertex (pr,1));
                    verts->insert(R_vertex (pr,2));
                    verts->insert(R_vertex (pr,3));
                  }
              }
            RIter_delete(rit);
          }
        }
      }
  }

  // ----------------------------------------------------------------------

}
