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
// Authors: Gaetan Compere, Jean-Francois Remacle
// -------------------------------------------------------------------

#include "MAdMetric.h"

#include <algorithm>

namespace MAd {

  MAdMetric::MAdMetric(const double l1, // (h_1)^-2
                       const double l2,
                       const double l3,
                       const double t1[3],
                       const double t2[3],
                       const double t3[3])
    {
      // M = e^1 * diag * e^1^t
      // where the elements of diag are l_i = h_i^-2
      // and the rows of e are the UNIT and ORTHOGONAL directions

      double e[3][3];
      e[0][0] = t1[0]; e[0][1] = t1[1]; e[0][2] = t1[2];
      e[1][0] = t2[0]; e[1][1] = t2[1]; e[1][2] = t2[2];
      e[2][0] = t3[0]; e[2][1] = t3[1]; e[2][2] = t3[2];
      invertMat(e);
    
      double tmp[3][3];
      tmp[0][0] = l1 * e[0][0]; tmp[0][1] = l2 * e[0][1]; tmp[0][2] = l3 * e[0][2];
      tmp[1][0] = l1 * e[1][0]; tmp[1][1] = l2 * e[1][1]; tmp[1][2] = l3 * e[1][2];
      tmp[2][0] = l1 * e[2][0]; tmp[2][1] = l2 * e[2][1]; tmp[2][2] = l3 * e[2][2];
      
      transpose(e);

      _val[0] = tmp[0][0] * e[0][0] + tmp[0][1] * e[1][0] + tmp[0][2] * e[2][0];
      _val[1] = tmp[1][0] * e[0][0] + tmp[1][1] * e[1][0] + tmp[1][2] * e[2][0];
      _val[2] = tmp[1][0] * e[0][1] + tmp[1][1] * e[1][1] + tmp[1][2] * e[2][1];
      _val[3] = tmp[2][0] * e[0][0] + tmp[2][1] * e[1][0] + tmp[2][2] * e[2][0];
      _val[4] = tmp[2][0] * e[0][1] + tmp[2][1] * e[1][1] + tmp[2][2] * e[2][1];
      _val[5] = tmp[2][0] * e[0][2] + tmp[2][1] * e[1][2] + tmp[2][2] * e[2][2];
    }
//     MAdMetric(const double l1, // (h_1)^-2
//               const double l2,
//               const double l3,
//               const double t1[3],
//               const double t2[3],
//               const double t3[3]){
//       double t1b[3], t2b[3], t3b[3];
//       t1b[0] = t1[0] * l1; t2b[0] = t2[0] * l1; t3b[0] = t3[0] * l1;
//       t1b[1] = t1[1] * l2; t2b[1] = t2[1] * l2; t3b[1] = t3[1] * l2;
//       t1b[2] = t1[2] * l3; t2b[2] = t2[2] * l3; t3b[2] = t3[2] * l3;
//       _val[0] = dotProd (t1b,t1);
//       _val[1] = dotProd (t2b,t1);
//       _val[2] = dotProd (t2b,t2);
//       _val[3] = dotProd (t3b,t1);
//       _val[4] = dotProd (t3b,t2);    
//       _val[5] = dotProd (t3b,t3);
//     }

  void MAdMetric::print (const char *s) const
  {
    printf(" metric %s : %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E \n",s,
           (*this)(0,0),(*this)(1,1),(*this)(2,2),
           (*this)(0,1),(*this)(0,2),(*this)(1,2));
  }


  MAdMetric intersection (const MAdMetric &m1, const MAdMetric &m2)
  {
    MAdMetric im1 = m1.invert();
    doubleMatrix V(3,3);
    doubleVector S(3);
    im1 *= m2;
    im1.eig(V,S,true);
    double v0[3]; v0[0] = V(0,0); v0[1] = V(1,0); v0[2] = V(2,0);
    double v1[3]; v1[0] = V(0,1); v1[1] = V(1,1); v1[2] = V(2,1);
    double v2[3]; v2[0] = V(0,2); v2[1] = V(1,2); v2[2] = V(2,2);
    double l0 = std::max(dot(v0,m1,v0),dot(v0,m2,v0));
    double l1 = std::max(dot(v1,m1,v1),dot(v1,m2,v1));
    double l2 = std::max(dot(v2,m1,v2),dot(v2,m2,v2));
    MAdMetric iv(l0,l1,l2,v0,v1,v2);
    return iv;
  }

  // (1-t) * m1 + t * m2
  MAdMetric interpolation (const MAdMetric &m1, 
                           const MAdMetric &m2, 
                           const double t)
  {
    MAdMetric im1 = m1.invert();
    MAdMetric im2 = m2.invert();
    im1 *= (1.-t);
    im2 *= t;
    im1 += im2;
    return im1.invert();
  }

  MAdMetric interpolation (const MAdMetric &m1, 
                           const MAdMetric &m2, 
                           const MAdMetric &m3, 
                           const double u,
                           const double v)
  {
    MAdMetric im1 = m1.invert();
    MAdMetric im2 = m2.invert();
    MAdMetric im3 = m3.invert();
    im1 *= (1.-u-v);
    im2 *= u;
    im3 *= v;
    im1 += im2;
    im1 += im3;
    return im1.invert();
  }

  MAdMetric interpolation (const MAdMetric &m1, 
                           const MAdMetric &m2, 
                           const MAdMetric &m3,  
                           const MAdMetric &m4, 
                           const double u,
                           const double v,
                           const double w)
  {
    MAdMetric im1 = m1.invert();
    MAdMetric im2 = m2.invert();
    MAdMetric im3 = m3.invert();
    MAdMetric im4 = m4.invert();
    im1 *= (1.-u-v-w);
    im2 *= u;
    im3 *= v;
    im4 *= w;
    im1 += im2;
    im1 += im3;
    im1 += im4;
    return im1.invert();
  }

}
