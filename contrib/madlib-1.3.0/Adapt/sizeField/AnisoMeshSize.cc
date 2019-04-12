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

#include "AnisoMeshSize.h"
#include "MathUtils.h"
#include "MAdMessage.h"
#include "MAdDefines.h"
#include "MeshParametersManager.h"

#include <math.h>
#include <iostream>
using std::string;

namespace MAd {

  // -------------------------------------------------------------------
  AnisoMeshSize::AnisoMeshSize(double dirs[3][3], double _h[3]):
    MeshSizeBase()
  {
    if( _h[0] <= 0. || _h[1] <= 0. || _h[2] <= 0. ) {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                  "Negative size(s): %f %f %f",
                                  _h[0], _h[1], _h[2]);
    }

    double t1[3], t2[3], t3[3];
    double l1,l2,l3;

    // sort the sizes from lowest to highest and fill first direction (GCRemark: not useful)
    int hMin = indexOfMin(_h[0], _h[1], _h[2]);
    l1 = 1. / (_h[hMin] * _h[hMin]);
    double normInv = 1. / sqrt ( dirs[hMin][0] * dirs[hMin][0] + 
                                 dirs[hMin][1] * dirs[hMin][1] + 
                                 dirs[hMin][2] * dirs[hMin][2]   );
    for (int i=0; i<3; i++) t1[i] = dirs[hMin][i] * normInv;

    hMin++;
    if( _h[(hMin+1)%3] > _h[(hMin)%3] )
      {
        l2 = 1. / (_h[hMin%3] * _h[hMin%3]);
        l3 = 1. / (_h[(hMin+1)%3] * _h[(hMin+1)%3]);
        hMin  = hMin%3;
      }
    else
      {
        l2 = 1. / (_h[(hMin+1)%3] * _h[(hMin+1)%3]);
        l3 = 1. / (_h[hMin%3] * _h[hMin%3]);
        hMin = (hMin+1)%3;
      }

    // project dir 1 so that it is perpendicular to dir 0
    double cosa = dotProd(t1,dirs[hMin]);
    double vec[3];
    for( int i=0; i<3; i++ ) vec[i] = dirs[hMin][i] - cosa * t1[i];
    normInv = 1. / sqrt ( vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2] );
    for (int i=0; i<3; i++) t2[i] = vec[i] * normInv;
    
    // find last dir as the direction orthogonal to the two previous ones
    crossProd(t1,t2,t3);

    // build the metric
    M = MAdMetric(l1,l2,l3,t1,t2,t3);
  }

  // -------------------------------------------------------------------
  AnisoMeshSize::AnisoMeshSize(double dir[3], double hDir, double hTg):
    MeshSizeBase()
  {
    if( hDir <= 0. || hTg <= 0. ) {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,"Negative size(s): %f %f",
                                  hDir, hTg);
    }
    
    double e[3][3];
    double len[3];

    // see if 'dir' is the direction of the minimal length
    int iDir = 0;
    if ( hDir > hTg ) iDir = 2;

    // set the 'dir' direction and its length
    normalizeVec(dir,e[iDir]);
    len[iDir] = 1. / ( hDir * hDir );

    // set other sizes in the right order
    len[(iDir+1)%3] = 1. / ( hTg * hTg );
    len[(iDir+2)%3] = len[(iDir+1)%3];

    // find a perpendicular direction to dir
    double dummy[3];
    dummy[0] = 2. * ( e[iDir][0] + 1.3654364 );
    dummy[1] = 3. * ( e[iDir][1] + 3.1368136 );
    dummy[2] = 5. * ( e[iDir][2] + 7.3683686 );
    normalizeVec(dummy,dummy);
    double cosa = dotProd(e[iDir],dummy);
    for( int i=0; i<3; i++ ) dummy[i] -= cosa * e[iDir][i];
    normalizeVec(dummy,e[(iDir+1)%3]);
    
    // find last dir as the direction orthogonal to the two previous ones
    crossProd(e[iDir],e[(iDir+1)%3],e[(iDir+2)%3]);

    // build the metric
    M = MAdMetric(len[0],len[1],len[2],e[0],e[1],e[2]);
  }

  // -------------------------------------------------------------------
  AnisoMeshSize::AnisoMeshSize(double h):
    MeshSizeBase()
  {
    M = MAdMetric( 1./(h*h) );
  }

  // -------------------------------------------------------------------
  AnisoMeshSize::AnisoMeshSize(const AnisoMeshSize &pm):
    MeshSizeBase()
  {
    M = pm.M;
  }

  // -------------------------------------------------------------------
  MeshSizeType AnisoMeshSize::getType() const 
  {
    return ANISOTROPIC;
  }

  // -------------------------------------------------------------------
  void AnisoMeshSize::intersect(const pMSize pMS1, const pMSize pMS2)
  {
    M = intersection(pMS1->getMetric(),pMS2->getMetric());
  }

  // -------------------------------------------------------------------
  void AnisoMeshSize::interpolate(const pMSize pMS0, const pMSize pMS1, 
                                  double param)
  {
    M = interpolation(pMS0->getMetric(),pMS1->getMetric(),param);
  }

  // -------------------------------------------------------------------
  void AnisoMeshSize::interpolate(const pMSize pMS0, const pMSize pMS1, 
                                  const pMSize pMS2, double u, double v)
  {
    M = interpolation(pMS0->getMetric(), pMS1->getMetric(), 
                      pMS2->getMetric(), u, v);
  }

  // -------------------------------------------------------------------
  void AnisoMeshSize::interpolate(const pMSize pMS0, const pMSize pMS1, 
                                  const pMSize pMS2, const pMSize pMS3, 
                                  double u, double v, double w)
  {
    M = interpolation(pMS0->getMetric(), pMS1->getMetric(), 
                      pMS2->getMetric(), pMS3->getMetric(), 
                      u, v, w);
  }

  // -------------------------------------------------------------------
  double AnisoMeshSize::size(int i) const
  {
    doubleMatrix V = doubleMatrix(3,3);
    doubleVector S = doubleVector(3);
    M.eig(V,S,true);
    double s = 1. / sqrt( S(i) );
    if ( isnan(s) ) return MeshParametersManagerSgl::instance().getBigLength();
    return s;
  }

  // -------------------------------------------------------------------
  void AnisoMeshSize::sizes(double _h[3]) const
  {
    doubleMatrix V = doubleMatrix(3,3);
    doubleVector S = doubleVector(3);
    M.eig(V,S,true);
    for (int i=0; i<3; i++) {
      _h[i] = ( 1. / sqrt( S(i) ) );
      if ( isnan(_h[i]) ) _h[i] = MeshParametersManagerSgl::instance().getBigLength();
    }
  }

  // -------------------------------------------------------------------
  double AnisoMeshSize::direction(int i, double dir[3]) const 
  {
    doubleMatrix V = doubleMatrix(3,3);
    doubleVector S = doubleVector(3);
    M.eig(V,S,true); // each column of V is a direction, S gives corresponding 1/h^2
    for (int iC=0; iC<3; iC++) dir[iC] = V(iC,i);
    double s = 1. / sqrt( S(i) );
    if ( isnan(s) ) return MeshParametersManagerSgl::instance().getBigLength();
    return s;
  }

  // -------------------------------------------------------------------
  void AnisoMeshSize::scale(int dir, double factor)
  {
    MAdMsgSgl::instance().error(__LINE__,__FILE__,"Not implemented");
  }

  // -------------------------------------------------------------------
  void AnisoMeshSize::scale(double factor)
  {
    MAdMsgSgl::instance().error(__LINE__,__FILE__,"Not implemented");
  }

  // -------------------------------------------------------------------
  double AnisoMeshSize::getMeanLength() const
  {
    doubleMatrix V = doubleMatrix(3,3);
    doubleVector S = doubleVector(3);
    M.eig(V,S,false);
    double mean = 0.;
    for (int i=0; i<3; i++) mean += 1. / sqrt( S(i) );
    return ( MAdTHIRD * mean );
  }

  // -------------------------------------------------------------------
  double AnisoMeshSize::getMinLength() const
  {
    doubleMatrix V = doubleMatrix(3,3);
    doubleVector S = doubleVector(3);
    M.eig(V,S,false);
    int i = indexOfMax(S(0),S(1),S(2));
    return ( 1. / sqrt( S(i) ) );
  }

  // -------------------------------------------------------------------
  double AnisoMeshSize::getMaxLength() const
  {
    doubleMatrix V = doubleMatrix(3,3);
    doubleVector S = doubleVector(3);
    M.eig(V,S,false);
    int i = indexOfMin(S(0),S(1),S(2));
    return ( 1. / sqrt( S(i) ) );
  }

  // -------------------------------------------------------------------
  // get the square of the norm in the metric
  double AnisoMeshSize::normSq(const double vec[3]) const
  {
    return dot(vec,M,vec);
  }

  // -------------------------------------------------------------------
  // get the square of desired edge length along 'vec'
  double AnisoMeshSize::lengthSqInDir(const double vec[3]) const
  {
    double tmp[3];
    normalizeVec(vec,tmp);
    return ( 1. / dot(tmp,M,tmp) );
  }
    
  // -------------------------------------------------------------------
  // get the cosine of the angle with the direction of the minimal size
  double AnisoMeshSize::angleWithDir0(const double dir[3]) const
  {
    double tmp[3];
    normalizeVec(dir,tmp);
    double dir0[3];
    direction(0,dir0);
    normalizeVec(dir0,dir0);
    return dotProd(tmp,dir0);
  }

  // -------------------------------------------------------------------
  void AnisoMeshSize::print(string name) const
  {
    std::cout<<"Printing AnisoMeshSize \'"<<name<<"\' ("<<this<<")\n";
    M.print("Mesh size metric");
  }

  // -------------------------------------------------------------------

}
