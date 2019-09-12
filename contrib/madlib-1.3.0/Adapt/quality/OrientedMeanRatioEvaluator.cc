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

#include "OrientedMeanRatioEvaluator.h"
#include "AnisoMeshSize.h"
#include "MAdDefines.h"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
using std::string;

namespace MAd {

  // -------------------------------------------------------------------
  orientedMeanRatioEvaluator::orientedMeanRatioEvaluator(const DiscreteSF * f): 
    meanRatioEvaluator(f)
  {}

  // -------------------------------------------------------------------
  int orientedMeanRatioEvaluator::XYZ_F_shape(const double xyz[3][3], const pMSize pS, 
                                              const double nor[3], double * shape) const
  {
    double aSq,sumSq;
  
    aSq = sizeField->SF_XYZ_areaSq(xyz,pS,nor);

    if( aSq <= 0. ) { *shape = 0.; return 0; } 
  
    sumSq  = sizeField->SF_XYZ_lengthSq(xyz[1],xyz[0],pS);
    sumSq += sizeField->SF_XYZ_lengthSq(xyz[2],xyz[0],pS);
    sumSq += sizeField->SF_XYZ_lengthSq(xyz[1],xyz[2],pS);
    
    *shape = 48.*aSq/(sumSq*sumSq);

    if ( pS->getType() == ANISOTROPIC )
      {
        double tmp[3]; diffVec(xyz[1],xyz[0],tmp);
        double cosa = fabs( ((AnisoMeshSize*)pS)->angleWithDir0(tmp) );
        *shape *= ( fabs( cosa - 0.5 ) + 0.5 );
        diffVec(xyz[2],xyz[0],tmp);
        cosa = fabs( ((AnisoMeshSize*)pS)->angleWithDir0(tmp) );
        *shape *= ( fabs( cosa - 0.5 ) + 0.5 );
        diffVec(xyz[2],xyz[1],tmp);
        cosa = fabs( ((AnisoMeshSize*)pS)->angleWithDir0(tmp) );
        *shape *= ( fabs( cosa - 0.5 ) + 0.5 );
      }

    if(*shape < MAdTOL) return 0;
    else if (*shape > (1.+MAdTOL) ) *shape = 1.;

    return 1;
  }

  // -------------------------------------------------------------------
  int orientedMeanRatioEvaluator::XYZ_R_shape(const double xyz[4][3], const pMSize pS, 
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

    if ( pS->getType() == ANISOTROPIC )
      {
        double dir0[3];  pS->direction(0,dir0); normalizeVec(dir0,dir0);
        double orientFactor = 1.;
        
        double tmp[3]; diffVec(xyz[1],xyz[0],tmp); normalizeVec(tmp,tmp);
        double cosa = fabs( dotProd(dir0,tmp) );
        orientFactor *= ( fabs( cosa - 0.5 ) + 0.5 );
        
        diffVec(xyz[2],xyz[0],tmp); normalizeVec(tmp,tmp);
        cosa = fabs( dotProd(dir0,tmp) );
        orientFactor *= ( fabs( cosa - 0.5 ) + 0.5 );
        
        diffVec(xyz[3],xyz[0],tmp); normalizeVec(tmp,tmp);
        cosa = fabs( dotProd(dir0,tmp) );
        orientFactor *= ( fabs( cosa - 0.5 ) + 0.5 );
        
        diffVec(xyz[2],xyz[1],tmp); normalizeVec(tmp,tmp);
        cosa = fabs( dotProd(dir0,tmp) );
        orientFactor *= ( fabs( cosa - 0.5 ) + 0.5 );
        
        diffVec(xyz[3],xyz[1],tmp); normalizeVec(tmp,tmp);
        cosa = fabs( dotProd(dir0,tmp) );
        orientFactor *= ( fabs( cosa - 0.5 ) + 0.5 );
        
        diffVec(xyz[3],xyz[2],tmp); normalizeVec(tmp,tmp);
        cosa = fabs( dotProd(dir0,tmp) );
        orientFactor *= ( fabs( cosa - 0.5 ) + 0.5 );

        *shape *= orientFactor;
      }

    if(*shape < MAdTOL) return 0;
    else if (*shape > (1.+MAdTOL) ) *shape = 1.;

    return 1;
  }

  // -------------------------------------------------------------------

}
