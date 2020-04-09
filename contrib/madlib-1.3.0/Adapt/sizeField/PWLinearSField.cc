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

#include "PWLinearSField.h"
#include "CallbackManager.h"
#include "IsoMeshSize.h"
#include "AnisoMeshSize.h"
#include "MathUtils.h"
#include "MAdMessage.h"
#include "MeshParametersManager.h"

#include <stdio.h>
#include <math.h>
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <queue>
#include <stdlib.h>
#include <string.h>

using namespace MAd;

#ifdef PARALLEL
static int PWLSFDE_tag = 86586745;
#endif

// -------------------------------------------------------------------
void PWLinearSFCBFunction (pPList before, pPList after, void * data,
                           operationType type , pEntity ppp) 
{
  PWLSField * pwl = (PWLSField *)(data);

  switch (type) {
  case MAd_ESPLIT: {

    double xyz[3];
    V_coord((pVertex)ppp,xyz);

    // find the old edge
    void *tmp=0;
    pEntity pE = PList_next(before,&tmp);

    // interpolate size at new location
    pMSize newSize = pwl->getSizeOnEntity(pE,xyz);

    pwl->setSize(ppp,newSize);    
    break;
  } 
  case MAd_ECOLLAPSE: {
    pwl->deleteSize(ppp);
    break;
  }
  case MAd_FSWAP:
  case MAd_ESWAP: {
    break;
  }
  case MAd_RREMOVE: {
    void * temp = NULL;
    while ( pEntity pE = PList_next(before,&temp) ) {
      if ( EN_type(pE) == 0 ) {
        pwl->deleteSize( pE );
      }
    }
    break;
  }
  default: {
    MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                "Not implemented for mesh modification %d",
                                type);
  }
  }
}

namespace MAd {

  // -------------------------------------------------------------------
  PWLSField::PWLSField(pMesh m, std::string name): DiscreteSF(m, name)
  {
    CallBackManagerSgl::instance().registerCallBack(PWLinearSFCBFunction,this);
  }
 
  // -------------------------------------------------------------------
  PWLSField::~PWLSField()
  {
    cleanUp();
    CallBackManagerSgl::instance().unregisterCallBack(PWLinearSFCBFunction,this);
  }

  // -------------------------------------------------------------------
  void PWLSField::cleanUp() 
  {
    VIter iter = M_vertexIter(mesh);
    while( pVertex pV = VIter_next(iter) ) deleteSize((pEntity)pV);
    VIter_delete(iter);
  }

  // -------------------------------------------------------------------
  void PWLSField::intersect(const pSField extField) 
  {
    VIter iter = M_vertexIter(mesh);
    while( pVertex pV = VIter_next(iter) )
      {
        // get the size for the extern size field
        pMSize extS = extField->getSize(pV);
        if (!extS) continue;

        // get the size for this size field
        pMSize intS = this->findSize(pV);
        if (!intS) {
          setSize((pEntity)pV,extS);
        }
        else {
          pMSize newS = MS_intersect(intS,extS);
          delete extS;
          setSize((pEntity)pV,newS);
        }
      }
    VIter_delete(iter);
  }

  // -------------------------------------------------------------------
  // Smooth the size field by limiting the size gradient along an edge 
  // to maxGrad
  void PWLSField::smooth(double maxGrad)
  {
    // edges still to be checked
    std::queue<pEdge> toCheck;

    // check every edge at least once
    EIter eIt = M_edgeIter(mesh);
    while( pEdge edge = EIter_next(eIt) ) toCheck.push(edge);
    EIter_delete(eIt);
  
    while ( !toCheck.empty() ) {

      std::set<pEdge> toAdd;

      // check the content of toCheck
      while ( !toCheck.empty() ) {
        smoothOnEdge(toCheck.front(),maxGrad,&toAdd);
        toCheck.pop();
      }

      // move content of toAdd to toCheck
      std::set<pEdge>::const_iterator eIter = toAdd.begin();
      for (; eIter != toAdd.end(); eIter++ ) toCheck.push(*eIter);
      toAdd.clear();
    }
  }

  // -------------------------------------------------------------------
  void PWLSField::smoothOnEdge(const pEdge edge, double maxGrad, 
                               std::set<pEdge>* toAdd)
  {
    pVertex pV[2];
    pMSize pMS[2];
    double h[2];
    for (int iV=0; iV<2; iV++) {
      pV[iV] = E_vertex(edge,iV);
      pMS[iV] = findSize(pV[iV]);
      if ( pMS[iV]->getType() != ISOTROPIC ) {
        printf("Error: PWLSField::smoothOnEdge not implemented for anisotropic sizes\n");
        exit(1);
      }
      h[iV] = ( (IsoMeshSize*)(pMS[iV]) )->size();
    }
  
    double physLength = E_length(edge);

    double grad = fabs( h[1] - h[0] ) / physLength;
    if ( grad > maxGrad ) {
    
      int toScale = 1;
      if ( h[0] > h[1] ) toScale = 0;
    
      double maxSize = maxGrad * physLength + h[1-toScale];
      ( (IsoMeshSize*)(pMS[toScale]) )->setSize(maxSize-MAdTOL);

      for( int iE=0; iE<V_numEdges(pV[toScale]); iE++ ) {
        pEdge newEdge = V_edge(pV[toScale],iE);
        if ( newEdge == edge ) continue;
        toAdd->insert(newEdge);
      }
    }
  }

  // -------------------------------------------------------------------
  void PWLSField::setCurrentSize() 
  {
    VIter vi = M_vertexIter(mesh);
    while ( pVertex vert = VIter_next(vi) ) {
      double lenSq = V_meanEdgeLenSq(vert);
      pMSize pSV = new IsoMeshSize(sqrt(lenSq));
      setSize((pEntity)vert,pSV);
    }
    VIter_delete(vi);
  }

  // -------------------------------------------------------------------
  void PWLSField::setCurvatureSize(bool aniso, double alpha, double hMin) 
  {
    if ( !(M_isParametric(mesh)) ) {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                  "Curvature not available without a geometry (mesh is not parametric");
    }

#ifdef _HAVE_GMSH_
    double bigLen = MeshParametersManagerSgl::instance().getBigLength();
    double curvMaxBound = 10. / hMin;
    pVertex vert;
    VIter vi = M_vertexIter(mesh);
    while ( ( vert = VIter_next(vi) ) )
      {
        pGEntity pge = V_whatIn(vert);
        int gdim = GEN_type(pge);
        pMSize pSV;
        if ( aniso && gdim==2 )
          {
            double u[2];
            V_params(vert,&(u[0]),&(u[1]));
            double dir[3][3], curv[3];
            GF_curvatures((pGFace)pge, u, dir[0], dir[1], &(curv[0]), &(curv[1]), curvMaxBound);
            crossProd(dir[0],dir[1],dir[2]);
            curv[2] = -1.;
            
            double h[3] = { -1., -1., -1. };
            for (int iD=0; iD<3; iD++) {
              if ( curv[iD] <= MAdTOL ) h[iD] = bigLen;
              else {
                h[iD] = 1. / ( curv[iD] * alpha );
                h[iD] = std::min(std::max(h[iD],hMin),bigLen);
              }
            }

            pSV = new AnisoMeshSize(dir,h);
          }
        else 
          {
            double curv = -1.;
            switch(gdim) {
            case 3: break;
            case 2: {
              double u[2];
              V_params(vert,&(u[0]),&(u[1]));
              curv = GF_curvatureDiv((pGFace)pge, u, curvMaxBound);
              break;
            }
            case 1: {
              double u, tmp;
              V_params(vert,&u,&tmp);
              curv = GE_curvature((pGEdge)pge, u, curvMaxBound);
              break;
            }
            case 0: {
              curv = -1.;
              std::vector<pGEdge> gEdges = GV_edges((pGVertex)pge);
              std::vector<pGEdge>::const_iterator eIter = gEdges.begin();
              for (; eIter != gEdges.end(); eIter++) {
                pGEdge pGEd = *eIter;
                double u;
                GV_reparamOnEdge((pGVertex)pge, pGEd, &u);
                double tmpcurv = GE_curvature(pGEd, u, curvMaxBound);
                if ( tmpcurv > curv ) curv = tmpcurv;
              }
              break;
            }
            }
        
            double h = -1.;
            if ( curv <= MAdTOL ) h = bigLen;
            else {
              h = 1. / ( curv * alpha );
              h = std::min(std::max(h,hMin),bigLen);
            }
            pSV = new IsoMeshSize(h);
          }
        setSize((pEntity)vert,pSV);
      }
    VIter_delete(vi);
#else
    MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                "Curvature not available without Gmsh");
#endif
  }

  // -------------------------------------------------------------------
  void PWLSField::setAllVSizes(pMSize pS)
  {
    VIter vi = M_vertexIter(mesh);
    while ( pVertex vert = VIter_next(vi) ) {
      pMSize pSCopy = MS_copy(pS);
      setSize((pEntity)vert,pSCopy);
    }
    VIter_delete(vi);
  }

  // -------------------------------------------------------------------
  void PWLSField::setAllVSizes(double dirs[3][3], // three unit vectors
                               double h[3])
  {
    pMSize pS = new AnisoMeshSize(dirs,h);
    setAllVSizes(pS);
  }

  // -------------------------------------------------------------------
  void PWLSField::setAllVSizes(double h)
  {
    pMSize pS = new IsoMeshSize(h);
    setAllVSizes(pS);
  }

  // -------------------------------------------------------------------
  void PWLSField::scale(double fact)
  {
    VIter vIter = M_vertexIter(mesh);
    while( pVertex pV = VIter_next(vIter) )
      {
        pMSize pS = findSize(pV);
        pS->scale(fact);
      }
    VIter_delete(vIter);
  }

  // -------------------------------------------------------------------
  // Length squared computation
  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  double PWLSField::SF_VV_lengthSq(const pVertex v0, const pVertex v1) const
  {
    double xyz[2][3];
    V_coord(v0,xyz[0]);
    V_coord(v1,xyz[1]);

    pMSize pS[2];
    pS[0] = findSize(v0);
    pS[1] = findSize(v1);

    return SF_XYZ_lengthSq(xyz[0],xyz[1],pS[0],pS[1]);
  }

  // -------------------------------------------------------------------
  double PWLSField::SF_XYZ_lengthSq(const double xyz0[3], const double xyz1[3], 
                                    const pMSize pS0, const pMSize pS1) const
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

  // -------------------------------------------------------------------
  double PWLSField::SF_F_areaSq(const pFace face) const
  {
    double area = 0.;

    double xyz[3][3];
    F_coordP1(face,xyz);
  
    void * temp = 0;
    pPList fVerts = F_vertices(face,1);
    while( pVertex pV = (pVertex)PList_next(fVerts,&temp) )
      {
        pMSize pS = findSize(pV);
        area += SF_XYZ_areaSq(xyz,pS,0);
      }
    PList_delete(fVerts);
  
    area /= F_numVertices(face);

    return area;
  }

  // -------------------------------------------------------------------
  double PWLSField::SF_XYZ_areaSq(const double fxyz[3][3], const pMSize pS, 
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

  // -------------------------------------------------------------------
  double PWLSField::SF_R_volume(const pRegion region) const
  {
    double vol = 0.;

    double xyz[4][3];
    R_coordP1(region,xyz);

    pPList rVerts = R_vertices(region);
    void * temp = 0;
    while( pVertex pV = (pVertex)PList_next(rVerts,&temp) )
      {
        pMSize pS = findSize(pV);
        vol += SF_XYZ_volume(xyz,pS);
      }
    PList_delete(rVerts);

    vol /= R_numVertices(region);

    return vol;
  }

  // -------------------------------------------------------------------
  double PWLSField::SF_XYZ_volume(const double xyz[4][3], const pMSize pS) const
  {
    if( !pS ) {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,"No size defined");
    }

    double physVol = R_XYZ_volume(xyz);
    double h[3];
    pS->sizes(h);
    return ( physVol / ( h[0] * h[1] * h[2]) );
  }

  // -------------------------------------------------------------------
  // Center of edge computation
  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  double PWLSField::SF_E_center(const pEdge edge, double center[3], 
                                double * reducSq, pMSize * cSize) const
  {
    return SF_VV_center(E_vertex(edge,0),E_vertex(edge,1),center,reducSq,cSize);
  }

  // -------------------------------------------------------------------
  double PWLSField::SF_VV_center(const pVertex v0, const pVertex v1,
                                 double center[3], double * reducSq, 
                                 pMSize * cSize) const
  {
    double xyz[2][3];
    V_coord(v0,xyz[0]);
    V_coord(v1,xyz[1]);

    pMSize pS[2];
    pS[0] = findSize(v0);
    pS[1] = findSize(v1);
    
    return SF_XYZ_center(xyz,pS,center,reducSq,cSize);
  }

  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  const pMSize PWLSField::findSize(const pVertex pV) const
  {
    void * size;
    if( EN_getDataPtr( (pEntity)pV, pMSizeFieldId, &size) ) return (pMSize)size;
    return NULL;
  }

  // -------------------------------------------------------------------
  pMSize PWLSField::findSize(const pVertex pV)
  {
    void * size;
    if( EN_getDataPtr( (pEntity)pV, pMSizeFieldId, &size) ) return (pMSize)size;
    return NULL;
  }

  // -------------------------------------------------------------------
  pMSize PWLSField::getSize(const pVertex pV) const
  {
    void * temp;
    if( EN_getDataPtr((pEntity)pV,pMSizeFieldId,&temp) ) {
      return MS_copy((pMSize)temp);
    }
    return NULL;
  }

  // -------------------------------------------------------------------
  pMSize PWLSField::getSizeOnEntity(const pEntity entity, 
                                    const double xyz[3]) const
  {
    int type = EN_type(entity);
    switch(type) {
    case 0: return getSize((pVertex)entity);
    case 1: return getSizeOnEdge((pEdge)entity,xyz);
    case 2: return getSizeOnFace((pFace)entity,xyz);
    case 3: return getSizeOnRegion((pRegion)entity,xyz);
    }
    return NULL;
  }

  // -------------------------------------------------------------------
  pMSize PWLSField::getSizeOnEdge(const pEdge edge, 
                                  const double xyz[3]) const
  {
    double u = E_linearParams(edge,xyz);
    return getSizeOnEdgeParam(edge,u);
  }

  // -------------------------------------------------------------------
  pMSize PWLSField::getSizeOnEdgeParam(const pEdge pE, 
                                       const double u) const
  {
    pMSize pS0 = findSize( E_vertex(pE,0) );
    pMSize pS1 = findSize( E_vertex(pE,1) );
  
    return MS_interpolate( pS0, pS1, u );
  }

  // -------------------------------------------------------------------
  pMSize PWLSField::getSizeOnFace(const pFace face, 
                                  const double xyz[3]) const
  {
    double u[2];
    F_linearParams(face,xyz,u);
    return getSizeOnFaceParam(face,u);
  }

  // -------------------------------------------------------------------
  pMSize PWLSField::getSizeOnFaceParam(const pFace face, 
                                       const double u[2]) const
  {
    pMSize pS0 = findSize( F_vertex(face,0) );
    pMSize pS1 = findSize( F_vertex(face,1) );
    pMSize pS2 = findSize( F_vertex(face,2) );

    return MS_interpolate( pS0, pS1, pS2, u[0], u[1] );
  }

  // -------------------------------------------------------------------
  pMSize PWLSField::getSizeOnRegion(const pRegion region, 
                                    const double xyz[3]) const
  {
    double u[3];
    R_linearParams(region, xyz, u);
    return getSizeOnRegionParam(region, u);
  }

  // -------------------------------------------------------------------
  pMSize PWLSField::getSizeOnRegionParam(const pRegion region, 
                                         const double u[3]) const
  {
    pMSize pS0 = findSize( R_vertex(region,0) );
    pMSize pS1 = findSize( R_vertex(region,1) );
    pMSize pS2 = findSize( R_vertex(region,2) );
    pMSize pS3 = findSize( R_vertex(region,3) );

    return MS_interpolate( pS0, pS1, pS2, pS3, u[0], u[1], u[2] );
  }

  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
#ifdef PARALLEL
  PWLSFieldDE::PWLSFieldDE(PWLSField *f) : 
    MDB_DataExchanger(PWLSFDE_tag), field(f) {}
  
  PWLSFieldDE::~PWLSFieldDE() {}
  
  // user allocates sends a message of _size size related to mesh entity pe to proc iProc
  void * PWLSFieldDE::sendData (pEntity pe,    // in
                                int iProcDest, // in
                                int &_size ) {
    if(EN_type(pe)==0){
      pMSize pS = field->findSize((pVertex) pe);
      if ( pS->getType() == ANISOTROPIC ) {
        MAdMsgSgl::instance().warning(__LINE__,__FILE__,
                                      "Anisotropic sizes exchange not implemented");
      }
      _size = sizeof(double);
      double *msg = (double * )malloc(_size);
      double size = pS->size(0);
      msg[0] = size;
      return msg;
    } 
    else {
      _size = 0;
      return 0;
    }
  }
  
  // mesh entity pe recieves data *buf form proc iProc.
  // The user shall NOT delete the message !!
  void PWLSFieldDE::receiveData (pEntity pe,      //in
                                 int iProcSender, //in
                                 void *buf ) {
    if(EN_type(pe)==0){
      double *msg = (double *) buf;
      if(buf){
        field->setSize((pEntity)pe,msg[0]);
      }
    }
    //pEntity * msg = (pEntity *) buf;
    //assert(pe==*msg);
  }
  
  // After migration  and if the entity is deleted on the proc,
  // the related data should be removed to avoid memory leak
  void PWLSFieldDE::deleteExternalData(pEntity pe) const {
    if(EN_type(pe)==0){ 
      pMSize pS = field->findSize((pVertex) pe);
      if( pS ) delete pS;
    }
  }
#endif
  
}
