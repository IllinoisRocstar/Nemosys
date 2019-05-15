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

#ifndef _H_MADMETRIC
#define _H_MADMETRIC

#include "MAdMatrix.h"
#include "MathUtils.h"

namespace MAd {

  // class for symmetric positive definite 3x3 matrix
  class MAdMetric {

  protected:
    // lower diagonal storage
    // 00 10 11 20 21 22 
    double _val[6];

  public:

    inline int getIndex(int i, int j) const{
      static int _index[3][3] = {{0,1,3},{1,2,4},{3,4,5}};
      return _index[i][j];
    }
    void getMat (doubleMatrix & mat) const{
      for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
          mat(i,j) = _val[getIndex(i,j)];			     
        }
      }
    }
    void getMat (double mat[3][3]) const{
      for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
          mat[i][j] = _val[getIndex(i,j)];			     
        }
      }
    }
    void setMat (const double mat[3][3]){
      for (int i=0;i<3;i++)
        for (int j=i;j<3;j++)
          _val[getIndex(i,j)] = mat[i][j];			     
    }
    MAdMetric(const MAdMetric& m){
      for (int i=0; i<6; i++) _val[i]=m._val[i];
    }
    // default constructor, identity 
    MAdMetric(const double v = 1.0){
      _val[0] = _val[2] = _val[5] = v;
      _val[1] = _val[3] = _val[4] = 0.0;
    }
    MAdMetric(const double l1, // (h_1)^-2
              const double l2,
              const double l3,
              const double t1[3],
              const double t2[3],
              const double t3[3]);
    inline double &operator()(int i, int j){ 
      return _val[getIndex(i,j)];
    }
    inline double operator()(int i, int j) const{ 
      return _val[getIndex(i,j)];
    }  
    MAdMetric invert () const {
      double m[3][3];
      getMat(m);
      invertMat(m);
      MAdMetric ithis;
      ithis.setMat(m);
      return ithis;
    }
    MAdMetric operator + (const MAdMetric &other) const {
      MAdMetric res(*this);
      for (int i=0;i<6;i++) res._val[i] += other._val[i];
      return res;
    }
    MAdMetric& operator += (const MAdMetric &other)  {
      for (int i=0;i<6;i++) _val[i] += other._val[i];
      return *this;
    }
    MAdMetric& operator *= (const double &other) {
      for (int i=0;i<6;i++) _val[i] *= other;
      return *this;
    }
    MAdMetric& operator *= (const MAdMetric &other) {
      double m1[3][3], m2[3][3], m3[3][3];
      getMat(m1);
      other.getMat(m2);
      matMat(m1,m2,m3);
      setMat(m3);
      return *this;
    }
    MAdMetric transform (double V[3][3]){
      double m[3][3], result[3][3], temp[3][3];
      getMat(m);
      transpose(V);
      matMat(V,m,temp);
      matMat(temp,V,result);
      MAdMetric a; a.setMat(result);
      return a;
    }
    // s: true if eigenvalues are sorted (from max to min of the absolute value of the REAL part)
    void eig (doubleMatrix &V, doubleVector &S, bool s=false) const {
      doubleMatrix me(3,3),right(3,3);
      doubleVector im(3);
      getMat(me);
      me.eig(V,S,im,right,s);
    }
    void print(const char *) const;
  };

  // scalar product with respect to the metric
  inline double dot(const double a[3], const MAdMetric &m, const double b[3])
  { return 
      b[0] * ( m(0,0) * a[0] + m(1,0) * a[1] + m(2,0) * a[2] ) + 
      b[1] * ( m(0,1) * a[0] + m(1,1) * a[1] + m(2,1) * a[2] ) + 
      b[2] * ( m(0,2) * a[0] + m(1,2) * a[1] + m(2,2) * a[2] ) ;}

  // compute the largest inscribed ellipsoid...
  MAdMetric intersection (const MAdMetric &m1, 
                          const MAdMetric &m2);
  MAdMetric interpolation (const MAdMetric &m1, 
                           const MAdMetric &m2, 
                           const double t);
  MAdMetric interpolation (const MAdMetric &m1, 
                           const MAdMetric &m2, 
                           const MAdMetric &m3, 
                           const double u,
                           const double v);
  MAdMetric interpolation (const MAdMetric &m1, 
                           const MAdMetric &m2, 
                           const MAdMetric &m3, 
                           const MAdMetric &m4, 
                           const double u,
                           const double v,
                           const double w);

}

#endif
