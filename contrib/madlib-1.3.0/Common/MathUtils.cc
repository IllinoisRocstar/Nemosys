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

#include "MathUtils.h"
#include "MAdMessage.h"
#include "MAdDefines.h"

#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>

#ifdef _HAVE_GSL_
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#endif

namespace MAd {

  // -------------------------------------------------------------------
  // VECTORS
  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  void diffVec(const double v1[3], const double v0[3], double v01[3])
  {
    v01[0] = v1[0] - v0[0];
    v01[1] = v1[1] - v0[1];
    v01[2] = v1[2] - v0[2];
  }

  // -------------------------------------------------------------------
  double dotProd(const double v1[3], const double v2[3])
  {
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
  }

  // -------------------------------------------------------------------
  void crossProd(const double v1[3], const double v2[3], double cp[3])
  {
    cp[0] = v1[1]*v2[2] - v1[2]*v2[1];
    cp[1] = v1[2]*v2[0] - v1[0]*v2[2];
    cp[2] = v1[0]*v2[1] - v1[1]*v2[0];
  }

  // -------------------------------------------------------------------
  void normalizeVec(const double v[3], double nv[3])
  {
    double lSq = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    if ( lSq <= MAdTOLSQ ) {
      MAdMsgSgl::instance().warning(__LINE__,__FILE__,
                                    "Normalization called for a zero vector");
      nv[0] = v[0]; 
      nv[1] = v[1]; 
      nv[2] = v[2];
    }
    double invNorm = 1. / sqrt( lSq );
    nv[0] = v[0] * invNorm ; 
    nv[1] = v[1] * invNorm ; 
    nv[2] = v[2] * invNorm ;
  }

  // -------------------------------------------------------------------
  void printVec(const double vec[3], const char * name)
  {
    printf("\nPrinting vector %s\n",name);
    for (int i=0; i<3; i++) printf("  %g",vec[i]);
    printf("\n\n");
  }

  // -------------------------------------------------------------------
  bool isNanVec (const double vec[3])
  {
    if ( isnan(vec[0]) || isnan(vec[1]) || isnan(vec[2]) ) return true;
    return false;
  }

  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  // MATRICES
  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  void transpose(double mat[3][3])
  {
    double tmp;
    tmp = mat[0][1]; mat[0][1] = mat[1][0]; mat[1][0] = tmp;
    tmp = mat[0][2]; mat[0][2] = mat[2][0]; mat[2][0] = tmp;
    tmp = mat[1][2]; mat[1][2] = mat[2][1]; mat[2][1] = tmp;
  }

  // -------------------------------------------------------------------
  void transpMat(const double mat[3][3], double matT[3][3])
  {
    for(int i=0; i<3; i++) for(int j=0; j<3; j++)
      matT[i][j] = mat[j][i];
  }

  // -------------------------------------------------------------------
  double detMat(const double mat[3][3])
  {
    return (mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) -
            mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) +
            mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]));
  }

  // -------------------------------------------------------------------
  double traceMat(const double mat[3][3])
  {
    return mat[0][0] + mat[1][1] + mat[2][2];
  }

  // -------------------------------------------------------------------
  double traceMatMat(const double mat[3][3])
  {
    double a00 =  mat[0][0] * mat[0][0] + mat[1][0] * mat[0][1] + mat[2][0] * mat[0][2]; 
    double a11 =  mat[1][0] * mat[0][1] + mat[1][1] * mat[1][1] + mat[1][2] * mat[2][1]; 
    double a22 =  mat[2][0] * mat[0][2] + mat[2][1] * mat[1][2] + mat[2][2] * mat[2][2];
  
    return a00 + a11 + a22;
  }

  // -------------------------------------------------------------------
  void vecMat(const double vec[3], const double mat[3][3], double res[3])
  {
    res[0] = mat[0][0] * vec[0] + mat[1][0] * vec[1] + mat[2][0] * vec[2];
    res[1] = mat[0][1] * vec[0] + mat[1][1] * vec[1] + mat[2][1] * vec[2];
    res[2] = mat[0][2] * vec[0] + mat[1][2] * vec[1] + mat[2][2] * vec[2];
  }

  // -------------------------------------------------------------------
  void matVec(const double mat[3][3], const double vec[3], double res[3])
  {
    res[0] = mat[0][0] * vec[0] + mat[0][1] * vec[1] + mat[0][2] * vec[2];
    res[1] = mat[1][0] * vec[0] + mat[1][1] * vec[1] + mat[1][2] * vec[2];
    res[2] = mat[2][0] * vec[0] + mat[2][1] * vec[1] + mat[2][2] * vec[2];
  }

  // -------------------------------------------------------------------
  void matMat(const double mat1[3][3], const double mat2[3][3], double res[3][3])
  {
    res[0][0] = mat1[0][0] * mat2[0][0] + mat1[0][1] * mat2[1][0] + mat1[0][2] * mat2[2][0];
    res[0][1] = mat1[0][0] * mat2[0][1] + mat1[0][1] * mat2[1][1] + mat1[0][2] * mat2[2][1];
    res[0][2] = mat1[0][0] * mat2[0][2] + mat1[0][1] * mat2[1][2] + mat1[0][2] * mat2[2][2];

    res[1][0] = mat1[1][0] * mat2[0][0] + mat1[1][1] * mat2[1][0] + mat1[1][2] * mat2[2][0];
    res[1][1] = mat1[1][0] * mat2[0][1] + mat1[1][1] * mat2[1][1] + mat1[1][2] * mat2[2][1];
    res[1][2] = mat1[1][0] * mat2[0][2] + mat1[1][1] * mat2[1][2] + mat1[1][2] * mat2[2][2];

    res[2][0] = mat1[2][0] * mat2[0][0] + mat1[2][1] * mat2[1][0] + mat1[2][2] * mat2[2][0];
    res[2][1] = mat1[2][0] * mat2[0][1] + mat1[2][1] * mat2[1][1] + mat1[2][2] * mat2[2][1];
    res[2][2] = mat1[2][0] * mat2[0][2] + mat1[2][1] * mat2[1][2] + mat1[2][2] * mat2[2][2];
  }

  // -------------------------------------------------------------------
  double vecMatVec(const double mat[3][3], const double vec[3])
  {
    double tmp[3];
    matVec(mat,vec,tmp);
    return dotProd(vec,tmp);
  }

  // -------------------------------------------------------------------
  double inverseMat(const double mat[3][3], double inv[3][3])
  {
    double det = detMat(mat);
    if(det) {
      double idet = 1. / det;
      inv[0][0] =  (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) * idet;
      inv[1][0] = -(mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) * idet;
      inv[2][0] =  (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]) * idet;
      inv[0][1] = -(mat[0][1] * mat[2][2] - mat[0][2] * mat[2][1]) * idet;
      inv[1][1] =  (mat[0][0] * mat[2][2] - mat[0][2] * mat[2][0]) * idet;
      inv[2][1] = -(mat[0][0] * mat[2][1] - mat[0][1] * mat[2][0]) * idet;
      inv[0][2] =  (mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1]) * idet;
      inv[1][2] = -(mat[0][0] * mat[1][2] - mat[0][2] * mat[1][0]) * idet;
      inv[2][2] =  (mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]) * idet;
    }
    else{
      MAdMsgSgl::instance().warning(__LINE__,__FILE__,"Singular matrix");
      for(int i=0; i<3; i++) for(int j=0; j<3; j++) inv[i][j] = 0.;
    }
    return det;
  }

  // -------------------------------------------------------------------
  double inverseMat22 (const double mat[2][2], double inv[2][2])
  {
    double det = mat[0][0] * mat[1][1] -  mat[0][1] * mat[1][0];
    if (det) {
      double idet = 1. / det;
      inv[0][0] =  idet * mat[1][1];
      inv[0][1] = -idet * mat[0][1];
      inv[1][0] = -idet * mat[1][0];
      inv[1][1] =  idet * mat[0][0];
    }
    else{
      MAdMsgSgl::instance().warning(__LINE__,__FILE__,"Singular matrix");
      for(int i=0; i<2; i++) for(int j=0; j<2; j++) inv[i][j] = 0.;
    }
    return det;
  }

  // -------------------------------------------------------------------
  double invertMat(double mat[3][3])
  {
    double inv[3][3];
    double det = inverseMat(mat,inv);
    for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {
        mat[i][j] = inv[i][j];
      }
    }
    return det;
  }

  // -------------------------------------------------------------------
  bool system(const double A[3][3], const double b[3], double res[3], double * det)
  {
    *det = detMat(A);
    if(*det == 0.0) {
      res[0] = res[1] = res[2] = 0.0;
      return false;
    }

    res[0] = b[0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) -
      A[0][1] * (b[1] * A[2][2] - A[1][2] * b[2]) +
      A[0][2] * (b[1] * A[2][1] - A[1][1] * b[2]);
    res[1] = A[0][0] * (b[1] * A[2][2] - A[1][2] * b[2]) -
      b[0] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
      A[0][2] * (A[1][0] * b[2] - b[1] * A[2][0]);
    res[2] = A[0][0] * (A[1][1] * b[2] - b[1] * A[2][1]) -
      A[0][1] * (A[1][0] * b[2] - b[1] * A[2][0]) +
      b[0] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);

    double idet = 1. / (*det);
    for(int i=0; i<3; i++) res[i] *= idet;

    return true;
  }

  // -------------------------------------------------------------------
  void printMat(const double mat[3][3], const char * name)
  {
    printf("\nPrinting matrix %s\n",name);
    for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {
        printf("  %g",mat[i][j]);
      }
      printf("\n");
    }
    printf("\n");  
  }

  // -------------------------------------------------------------------
  void meanRow33(const double mat[3][3], double mean[3])
  {
    mean[0] = MAdTHIRD * ( mat[0][0] + mat[1][0] + mat[2][0] );
    mean[1] = MAdTHIRD * ( mat[0][1] + mat[1][1] + mat[2][1] );
    mean[2] = MAdTHIRD * ( mat[0][2] + mat[1][2] + mat[2][2] );
  }

  // -------------------------------------------------------------------
  void meanRow43(const double mat[4][3], double mean[3])
  {
    mean[0] = 0.25 * ( mat[0][0] + mat[1][0] + mat[2][0] + mat[3][0] );
    mean[1] = 0.25 * ( mat[0][1] + mat[1][1] + mat[2][1] + mat[3][1] );
    mean[2] = 0.25 * ( mat[0][2] + mat[1][2] + mat[2][2] + mat[3][2] );
  }

  // -------------------------------------------------------------------
  bool isNanMat (const double mat[3][3])
  {
    if ( isNanVec(mat[0]) || isNanVec(mat[1]) || isNanVec(mat[2]) ) return true;
    return false;
  }

  // -------------------------------------------------------------------
  // Miscellaneous
  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  void sort(double lst[3])
  {
    for (int i=0; i<3; i++) {
      int k = i;
      double mx = lst[i];
      for (int j=i+1; j<3; j++) if (lst[j] >= mx) { mx = lst[j]; k = j; }
      if (k != i) { lst[k] = lst[i]; lst[i] = mx; }
    }
  }

  // -------------------------------------------------------------------
  void cubicRoots(const double coef[4], double real[3], double imag[3]) // GCTODO: to be rewritten
  {
    double a = coef[3];
    double b = coef[2];
    double c = coef[1];
    double d = coef[0];

    if(!a || !d){
      printf("Error: Degenerate cubic: use a second degree solver!\n");
      return;
    }

    b /= a;
    c /= a;
    d /= a;
  
    double q = (3.0*c - (b*b))/9.0;
    double r = -(27.0*d) + b*(9.0*c - 2.0*(b*b));
    r /= 54.0;

    double discrim = q*q*q + r*r;
    imag[0] = 0.0; // The first root is always real.
    double term1 = (b/3.0);

    if (discrim > 0) { // one root is real, two are complex
      double s = r + sqrt(discrim);
      s = ((s < 0) ? -pow(-s, (1.0/3.0)) : pow(s, (1.0/3.0)));
      double t = r - sqrt(discrim);
      t = ((t < 0) ? -pow(-t, (1.0/3.0)) : pow(t, (1.0/3.0)));
      real[0] = -term1 + s + t;
      term1 += (s + t)/2.0;
      real[1] = real[2] = -term1;
      term1 = sqrt(3.0)*(-t + s)/2;
      imag[1] = term1;
      imag[2] = -term1;
      return;
    }

    // The remaining options are all real
    imag[1] = imag[2] = 0.0;
  
    double r13;
    if (discrim == 0){ // All roots real, at least two are equal.
      r13 = ((r < 0) ? -pow(-r,(1.0/3.0)) : pow(r,(1.0/3.0)));
      real[0] = -term1 + 2.0*r13;
      real[1] = real[2] = -(r13 + term1);
      return;
    }

    // Only option left is that all roots are real and unequal (to get
    // here, q < 0)
    q = -q;
    double dum1 = q*q*q;
    dum1 = acos(r/sqrt(dum1));
    r13 = 2.0*sqrt(q);
    real[0] = -term1 + r13*cos(dum1/3.0);
    real[1] = -term1 + r13*cos((dum1 + 2.0*M_PI)/3.0);
    real[2] = -term1 + r13*cos((dum1 + 4.0*M_PI)/3.0);
  }


  // -------------------------------------------------------------------
  // solve x^2 + b x + c = 0
  // x[2] is always set to be zero
  // long FindQuadraticRoots(const double b, const double c, double x[3]) // GCTODO: to be rewritten
  // {
  //   //    printf("Quadratic roots\n");
  //   x[2]=0.0;
  //   double delt=b*b-4.*c;
  //   if( delt >=0 ) {
  //     delt=sqrt(delt);
  //     x[0]=(-b+delt)/2.0;
  //     x[1]=(-b-delt)/2.0;
  //     return 3;
  //   }
    
  //   printf("Imaginary roots, impossible, delt=%f\n",delt);
  //   return 1;
  // }
  
  // -------------------------------------------------------------------
  // solve x^3 + a1 x^2 + a2 x + a3 = 0
  // long FindCubicRoots(const double coeff[4], double x[3]) // GCTODO: to be rewritten
  // {
  //   double a1 = coeff[2] / coeff[3];
  //   double a2 = coeff[1] / coeff[3];
  //   double a3 = coeff[0] / coeff[3];
    
  //   if( fabs(a3)<1.0e-8 ) 
  //     return FindQuadraticRoots(a1,a2,x);
    
  //   double Q = (a1 * a1 - 3 * a2) / 9.;
  //   double R = (2. * a1 * a1 * a1 - 9. * a1 * a2 + 27. * a3) / 54.;
  //   double Qcubed = Q * Q * Q;
  //   double d = Qcubed - R * R;
    
  //   //    printf ("d = %22.15e Q = %12.5E R = %12.5E Qcubed %12.5E\n",d,Q,R,Qcubed);

  //   /// three roots, 2 equal 
  //     if(Qcubed == 0.0 || fabs ( Qcubed - R * R ) < 1.e-8 * (fabs ( Qcubed) + fabs( R * R)) )
  //       {
  //         double theta;
  //         if (Qcubed <= 0.0)theta = acos(1.0);
  //         else if (R / sqrt(Qcubed) > 1.0)theta = acos(1.0); 
  //         else if (R / sqrt(Qcubed) < -1.0)theta = acos(-1.0); 
  //         else theta = acos(R / sqrt(Qcubed));
  //         double sqrtQ = sqrt(Q);
  //         //      printf("sqrtQ = %12.5E teta=%12.5E a1=%12.5E\n",sqrt(Q),theta,a1);
  //         x[0] = -2 * sqrtQ * cos( theta           / 3) - a1 / 3;
  //         x[1] = -2 * sqrtQ * cos((theta + 2 * M_PI) / 3) - a1 / 3;
  //         x[2] = -2 * sqrtQ * cos((theta + 4 * M_PI) / 3) - a1 / 3;
  //         return (3);
  //       }

  //     // Three real roots 
  //     if (d >= 0.0) {
  //       double theta = acos(R / sqrt(Qcubed));
  //       double sqrtQ = sqrt(Q);
  //       x[0] = -2 * sqrtQ * cos( theta           / 3) - a1 / 3;
  //       x[1] = -2 * sqrtQ * cos((theta + 2 * M_PI) / 3) - a1 / 3;
  //       x[2] = -2 * sqrtQ * cos((theta + 4 * M_PI) / 3) - a1 / 3;
  //       return (3);
  //     }
    
  //     // One real root 
  //     else {
  //       printf("IMPOSSIBLE !!!\n");

  //       double e = pow(sqrt(-d) + fabs(R), 1. / 3.);
  //       if (R > 0)
  //         e = -e;
  //       x[0] = (e + Q / e) - a1 / 3.;
  //       return (1);
  //     }
  // }

  // -------------------------------------------------------------------
  int indexOfMin(double v0, double v1, double v2)
  {
    double minV = v0;
    int index   = 0;
    if ( v1 < minV ) { minV = v1; index = 1; }
    if ( v2 < minV ) { minV = v2; index = 2; }
    return index;
  }

  // -------------------------------------------------------------------
  int indexOfMax(double v0, double v1, double v2)
  {
    double maxV = v0;
    int index   = 0;
    if ( v1 > maxV ) { maxV = v1; index = 1; }
    if ( v2 > maxV ) { maxV = v2; index = 2; }
    return index;
  }

  // -------------------------------------------------------------------
  double Tri_area(const double xyz[3][3])
  {
    double e0[3], e1[3];
    diffVec(xyz[1],xyz[0],e0);
    diffVec(xyz[2],xyz[0],e1);
    double nor[3];
    crossProd(e0,e1,nor);
    return ( 0.5 * dotProd(nor,nor) );
  }

  // -------------------------------------------------------------------
  double distToLineSq(const double p0[3], const double p1[3], 
                      const double xyz[3], double proj2pt[3], bool * onSegment)
  {
    double v01[3]; diffVec(p1,p0,v01);
    double v0x[3]; diffVec(xyz,p0,v0x);
    double t = dotProd(v01,v0x) / dotProd(v01,v01);
    
    if ( onSegment ) {
      *onSegment = false;
      if ( t >= 0. && t <= 1. ) *onSegment = true;
    }

    for (int i=0; i<3; i++) {
      proj2pt[i] = xyz[i] - ( (1.-t) * p0[i] + t * p1[i] );
    }

    return dotProd(proj2pt,proj2pt);
  }

  // -------------------------------------------------------------------
//   double distToLineSq(const double p0[3], const double p1[3], 
//                       const double xyz[3], bool * onSegment)
//   {
//     double v01[3]; diffVec(p1,p0,v01);
//     double vx0[3]; diffVec(p0,xyz,vx0);
//     double d01Sq = dotProd(v01,v01);
//     double dx0Sq = dotProd(vx0,vx0);
//     double prod  = dotProd(vx0,v01);

//     if ( onSegment ) {
//       *onSegment = false;
//       double vx1[3]; diffVec(p1,xyz,vx1);
//       if ( prod <= 0. && dotProd(vx1,v01) >= 0. ) *onSegment = true;
//     }

//     return std::max( ( dx0Sq * d01Sq - prod * prod ) / d01Sq, 0. );
//   }
  
  // -------------------------------------------------------------------
  //! Returns the coordinates (u,v) of the point in the parent element
  void Tri_linearParams(const double tri[3][3], const double xyz[3], 
                        double res[2])
  {
    double v0X[3], v01[3], v02[3];
    diffVec(xyz,tri[0],v0X);
    diffVec(tri[1],tri[0],v01);
    diffVec(tri[2],tri[0],v02);

    double n[3], nTmp[3];
    crossProd(v01,v02,n);
    double ASqInv = 1. / dotProd(n,n);

    crossProd(v0X,v02,nTmp);
    double A0X2Sq = dotProd(nTmp,nTmp);
    res[0] = sqrt( A0X2Sq * ASqInv );
    if ( dotProd(n,nTmp) < 0. ) res[0] = -res[0];

    crossProd(v01,v0X,nTmp);
    double A01XSq = dotProd(nTmp,nTmp);
    res[1] = sqrt( A01XSq * ASqInv );
    if ( dotProd(n,nTmp) < 0. ) res[1] = -res[1];

//     double n[3];
//     crossProd(v01,v02,n);
//     double A = sqrt ( dotProd(n,n) );
//     crossProd(v0X,v02,n);
//     double A1 = sqrt ( dotProd(n,n) );
//     crossProd(v01,v0X,n);
//     double A2 = sqrt ( dotProd(n,n) );
    
//     res[0] = A1/A;
//     res[1] = A2/A;
  }

  // -------------------------------------------------------------------
  //! Gets the projection of a point on a plane defined by 3 points
  //! and a bit representing the zone of the plane on wich the projection 
  //! falls given by:
  //!
  //!            010=2   | 011=3  /  001=1
  //!                    |       /
  //!        ------------+--e2--+-----------
  //!                  v0|     /v2
  //!                    | 7  /
  //!            110=6   e0  e1
  //!                    |  /
  //!                    | /    101=5
  //!                    |/
  //!                  v1+
  //!                   /|
  //!                  / |
  //!                   4
  //!
  //! and 8, 9 or 10 if it falls on v0, v1 or v2 (respectively)
  //!
  //! Computes the square area of the faces of the tetrahedron defined by the 
  //! triangle and the point, projected on the plane of the triangle, 
  //! with the triangle being the first, and then 
  //! the 3 faces being ordered the same way as the edges e0, e1, e2 in 
  //! the triangle.
  //!
  void pointToTriangle(const double tri[3][3], 
                       const double xyz[3], 
                       double proj[3],
                       int * bit,
                       double area[4])
  {
    pointToPlane(tri,xyz,proj);
    
    // find if the projection coincides with any of the 3 points 
    for(int i=0; i<3; i++) {
      if( distanceSq(proj,tri[i]) < MAdTOL ) {
        *bit=8+i;
//         area[0] = triArea(tri);
//         area[1] = area[2] = area[3] = 0.;
//#warning "implement this"
        return;
      }
    }
    
    // find normal to the triangle
    double v01[3], v02[3], norm[3];
    diffVec(tri[1],tri[0],v01);
    diffVec(tri[2],tri[0],v02);
    crossProd(v01,v02,norm);

    double ri[3], rj[3], rk[3];
    diffVec(proj,tri[0],ri);
    diffVec(proj,tri[1],rj);
    diffVec(proj,tri[2],rk);
    
    // determine on which side of the edges does the point lie.
    // First get normal vectors 
    double temp[3];
    double normi[3], normj[3], normk[3];
    crossProd(v01,ri,normi);
    diffVec(tri[2],tri[1],temp);
    crossProd(temp,rj,normj);
    diffVec(tri[0],tri[2],temp);
    crossProd(temp,rk,normk);
    
    double mag[3];
    mag[0] = dotProd(normi,norm);
    mag[1] = dotProd(normj,norm);
    mag[2] = dotProd(normk,norm);
    
    if (area) {
      area[0] = 0.5 * dotProd(norm,norm);
      area[1] = 0.5 * dotProd(normi,normi);
      area[2] = 0.5 * dotProd(normj,normj);
      area[3] = 0.5 * dotProd(normk,normk);
    }

    int filter[] = {1,2,4};
    *bit=0;
    for(int i=0; i<3; i++){
      if( mag[i] > 0. ) {
        *bit = *bit | filter[i];
      }
    }

  }

  // -------------------------------------------------------------------
  //! Gets the distance from a point to a triangle. Consider the 
  //! projection of the point to the plane of the triangle and a bit
  //! taking the following values according to the zone in which the 
  //! projection falls:
  //!
  //!            010=2   | 011=3|
  //!                    |      |     001=1
  //!        ------------+--e2--+v2
  //!                  v0|     / `
  //!                    | 7  /   `
  //!            110=6   e0  e1    `
  //!                    |  /       `
  //!                    | /  101=5
  //!                  v1|/
  //!        ------------+
  //!                     `
  //!             100=4    `
  //!                       `
  //!
  //! and 8, 9 or 10 if it falls on v0, v1 or v2 (respectively).
  //! The distance will be 
  //!    * the distance to v0, v1, v2 if the bit = 2,4,1 resp.,
  //!    * the distance to e0, e1, e2 if the bit = 6,5,3 resp.,
  //!    * the distance to the projection point if the bit >= 7
  //!
  //! Computes the vector from the closest point on the tri to 'xyz'.
  //!
  //! GC note: returning the square distance only allows to find the 
  //!          distance with a precision of the square of the machine precision.
  //!          It could be a problem when computing Laplacian of the distance for 
  //!          the curvature on highly anisotropic elements.
  double distToTriangleSq(const double tri[3][3], 
                          const double xyz[3],
                          double vec[3])
  {
    double distSq = MAdBIG;

    double proj[3];
    pointToPlane(tri,xyz,proj);

    double par[2];
    Tri_linearParams(tri,proj,par);
    
    // -- see if the closest point is inside the triangle ---
    if (  par[0] >= 0. && par[1] >= 0. && 
          par[0] + par[1] <= 1. )
      {
        diffVec(xyz,proj,vec);
        distSq = dotProd(vec,vec);
      }
    // --- else ---
    else
      {
        // --- see if the closest point is on an edge of the triangle ---
        bool onSeg=false, onSegTest;
        double testD;
        double testVec[3];
        for (int iE=0; iE<3; iE++) {
          testD = distToLineSq(tri[iE],tri[(iE+1)%3],xyz,testVec,&onSegTest);
          if ( onSegTest ) {
            if ( testD < distSq ) {
              distSq = testD;
              vec[0] = testVec[0]; vec[1] = testVec[1]; vec[2] = testVec[2];
            }
            onSeg = true;
          }
        }

        // --- otherwise it is a summit of the triangle
        if ( !onSeg ) {
          for (int iV=0; iV<3; iV++) {
            diffVec(xyz,tri[iV],testVec);
            testD = dotProd(testVec,testVec);
            if ( testD < distSq ) {
              distSq = testD;
              vec[0] = testVec[0]; vec[1] = testVec[1]; vec[2] = testVec[2];
            }
          }
        }
      }

    return distSq;
  }

  // -------------------------------------------------------------------
  //! Gets the projection of a point on a plane defined by 3 points
  void pointToPlane(const double plane[3][3], 
                    const double xyz[3], 
                    double proj[3])
  {
    double v01[3], v02[3];
    diffVec(plane[1],plane[0],v01);
    diffVec(plane[2],plane[0],v02);
    double normal[3];
    crossProd(v01,v02,normal);
    double areaSq = dotProd(normal,normal);
  
    double v0X[3];
    diffVec(xyz,plane[0],v0X);
    double ASqX = dotProd(v0X,normal);
  
    double ratio = ASqX / areaSq;
    for (int i=0; i<3; i++) proj[i] = xyz[i] - ratio * normal[i];
  }

  // -------------------------------------------------------------------
  double distanceSq(const double xyz0[3], const double xyz1[3])
  {
    return ( ( xyz1[0] - xyz0[0] ) * ( xyz1[0] - xyz0[0] ) +
             ( xyz1[1] - xyz0[1] ) * ( xyz1[1] - xyz0[1] ) +
             ( xyz1[2] - xyz0[2] ) * ( xyz1[2] - xyz0[2] )  );
  }

  // -------------------------------------------------------------------

} // End of namespace MAd
