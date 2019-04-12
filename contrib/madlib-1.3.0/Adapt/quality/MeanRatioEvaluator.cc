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

#include "MeanRatioEvaluator.h"
#include "MeshSizeBase.h"
#include "MAdDefines.h"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
using std::string;

namespace MAd {

  // -------------------------------------------------------------------
  meanRatioEvaluator::meanRatioEvaluator(const DiscreteSF * f): 
    elementEvaluatorBase(f)
  {}

  // -------------------------------------------------------------------
  int meanRatioEvaluator::F_shape(const pFace face, const double normal[3],
                                  double * shape) const
  {
    double fxyz[3][3];
    pMSize pS[3];

    void *iter=0;
    int i = 0;
    pPList fvts=F_vertices(face,1);
    while( pVertex vt=(pVertex)PList_next(fvts,&iter) ) {
      pS[i]=sizeField->findSize(vt);
      V_coord(vt,fxyz[i++]);  
    }
    PList_delete(fvts);

    return XYZ_F_shape(fxyz,pS,normal,shape);  
  }

  // -------------------------------------------------------------------
  int meanRatioEvaluator::F_shapeWithDisp(const pFace face, const double normal[3],
                                          double disp[3][3], double * shape) const
  {
    double fxyz[3][3];
    pMSize pS[3];

    void *iter=0;
    int i = 0;
    pPList fvts=F_vertices(face,1);
    while( pVertex vt=(pVertex)PList_next(fvts,&iter) ) {
      pS[i]=sizeField->findSize(vt);
      V_coord(vt,fxyz[i++]);  
    }
    PList_delete(fvts);

    for (int j=0; j<3; j++) for (int k=0; k<3; k++) fxyz[j][k] += disp[j][k];

    return XYZ_F_shape(fxyz,pS,normal,shape);  
  }

  // -------------------------------------------------------------------
  int meanRatioEvaluator::XYZ_F_shape(const double xyz[3][3], const pMSize pS[3], 
                                      const double nor[3], double * shape) const
  {
    if( sizeField->getType() == NULLSFIELD ) {
      if( ! XYZ_F_shape(xyz,pS[0],nor,shape) )
        { *shape=0; return 0; }
      return 1;
    }

    // use central size to compute shape
    pMSize pSZ = MS_interpolate(pS[0], pS[1], pS[2], 0.33, 0.33);
    int res = XYZ_F_shape(xyz, pSZ, nor, shape);
    if ( pSZ ) delete pSZ;
    return res;

//     // compute shapes using sizes at every vertex and retain the worst one
//     double shp[3];
//     for(int i=0; i<3; i++) {
//       if( ! XYZ_F_shape(xyz, pS[i], nor, &shp[i]) )
//         { *shape = 0.; return 0; }
//     }
//     *shape = std::min(std::min(shp[0],shp[1]),shp[2]);
//     return 1;
  }

  // -------------------------------------------------------------------
  int meanRatioEvaluator::XYZ_F_shape(const double xyz[3][3], const pMSize pS, 
                                      const double nor[3], double * shape) const
  {
    double aSq,sumSq;
  
    aSq = sizeField->SF_XYZ_areaSq(xyz,pS,nor);

    if( aSq <= 0. ) { *shape = 0.; return 0; } 
  
    sumSq  = sizeField->SF_XYZ_lengthSq(xyz[1],xyz[0],pS);
    sumSq += sizeField->SF_XYZ_lengthSq(xyz[2],xyz[0],pS);
    sumSq += sizeField->SF_XYZ_lengthSq(xyz[1],xyz[2],pS);
    
    *shape = 48.*aSq/(sumSq*sumSq);

    if(*shape < MAdTOL) return 0;
    else if (*shape > (1.+MAdTOL) ) *shape = 1.;

    return 1;
  }

  // -------------------------------------------------------------------
  int meanRatioEvaluator::R_shape(const pRegion rgn, double * shape) const
  {
    double rxyz[4][3];
    pMSize pS[4];
  
    void *iter=0;
    int i=0;
    pPList verts=R_vertices(rgn);
    while( pVertex vt=(pVertex)PList_next(verts,&iter) ) {
      pS[i]=sizeField->findSize(vt);
      V_coord(vt,rxyz[i++]);
    }  
    PList_delete(verts);

    return XYZ_R_shape(rxyz,pS,shape);
  }

  // -------------------------------------------------------------------
  int meanRatioEvaluator::R_shapeWithDisp(const pRegion rgn, double disp[4][3], 
                                          double * shape) const
  {
    double rxyz[4][3];
    pMSize pS[4];
  
    void *iter=0;
    int i=0;
    pPList verts=R_vertices(rgn);
    while( pVertex vt=(pVertex)PList_next(verts,&iter) ) {
      pS[i]=sizeField->findSize(vt);
      V_coord(vt,rxyz[i++]);
    }  
    PList_delete(verts);

    for (int j=0; j<4; j++) for (int k=0; k<3; k++) rxyz[j][k] += disp[j][k];

    return XYZ_R_shape(rxyz,pS,shape);
  }

  // -------------------------------------------------------------------
  int meanRatioEvaluator::XYZ_R_shape(const double xyz[4][3], const pMSize pS[4], 
                                      double * shape) const
  {
    if( sizeField->getType() == NULLSFIELD ) {
      if( ! XYZ_R_shape(xyz,pS[0],shape) )
        { *shape=0.; return 0; }
      return 1;
    }

    // use central size
    pMSize pSZ = MS_interpolate(pS[0], pS[1], pS[2], pS[3], 0.25, 0.25, 0.25);

    // use minimal size
//     pMSize pSZ;
//     double hmin = MAdBIG;
//     for(int i=0; i<4; i++) {
//       if( pS[i] ) {
//         double s = pS[i]->minSize();
//         if( s <  hmin ) { hmin = s; pSZ = pS[i]; }
//       }
//     }
//     assert( hmin != MAdBIG );

    int res = XYZ_R_shape(xyz,pSZ,shape);
    if ( pSZ ) delete pSZ;
    return res;
  }

  // -------------------------------------------------------------------
  int meanRatioEvaluator::XYZ_R_shape(const double xyz[4][3], const pMSize pS, 
                                      double * shape) const
  {
    double vol,sumSq;
  
    vol = sizeField->SF_XYZ_volume(xyz,pS);

    if ( vol < 0. ) { *shape = 0.; return 0; }
  
    sumSq  = sizeField->SF_XYZ_lengthSq(xyz[1],xyz[0],pS);
    sumSq += sizeField->SF_XYZ_lengthSq(xyz[2],xyz[0],pS);
    sumSq += sizeField->SF_XYZ_lengthSq(xyz[3],xyz[0],pS);
    sumSq += sizeField->SF_XYZ_lengthSq(xyz[3],xyz[1],pS);
    sumSq += sizeField->SF_XYZ_lengthSq(xyz[3],xyz[2],pS);
    sumSq += sizeField->SF_XYZ_lengthSq(xyz[1],xyz[2],pS);

    *shape = 15552.*vol*vol/(sumSq*sumSq*sumSq);

    if(*shape < MAdTOL) return 0;
    else if (*shape > (1.+MAdTOL) ) *shape = 1.;

    return 1;
  }

  // -------------------------------------------------------------------

}
