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

#ifndef _H_MATRIX_MAD
#define _H_MATRIX_MAD

#include "MAdMessage.h"

#include <assert.h>
#include <math.h>
#include <iostream>
#include <stdio.h>

#if defined(_HAVE_GSL_)
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#endif

// -------------------------------------------------------------------
#if defined(HAVE_BLAS)
extern "C" {
  void dgemm_(const char *transa, const char *transb, int *m, int *n, int *k, 
              double *alpha, double *a, int *lda, double *b, int *ldb, 
              double *beta, double *c, int *ldc);
  void dgemv_(const char *trans, int *m, int *n, double *alpha, double *a, 
              int *lda, double *x, int *incx, double *beta, double *y, int *incy);
}
#endif

// -------------------------------------------------------------------
#if defined(HAVE_LAPACK)
extern "C" {
  void dgesv_(int *N, int *nrhs, double *A, int *lda, int *ipiv, 
              double *b, int *ldb, int *info);
  void dgetrf_(int *M, int *N, double *A, int *lda, int *ipiv, int *info);
  void dgesvd_(const char* jobu, const char *jobvt, int *M, int *N,
               double *A, int *lda, double *S, double* U, int *ldu,
               double *VT, int *ldvt, double *work, int *lwork, int *info);
  void dgeev_(const char *jobvl, const char *jobvr, 
	      int *n, double *a, int *lda, 
	      double *wr, double *wi, 
	      double *vl, int *ldvl, 
	      double *vr, int *ldvr, 
	      double *work, int *lwork,
	      int *info); 
}
#endif

// -------------------------------------------------------------------

namespace MAd {

  // -------------------------------------------------------------------
  // Basic vector / matrix
  // -------------------------------------------------------------------
  template <class SCALAR>
  class MAd_Vector
  {
  private:
    int r;
  public:
    inline int size() const { return r; }
    SCALAR *data;
    ~MAd_Vector() { if(data) delete [] data; }
    MAd_Vector() : r(0)
    {
      data = 0;
    }
    MAd_Vector(int R) : r(R)
    {
      data = new SCALAR[r];
      scale(0);
    }
    MAd_Vector(const MAd_Vector<SCALAR> &other) : r(other.r)
    {
      data = new SCALAR[r];
      for (int i = 0; i < r; ++i) data[i] = other.data[i];
    }
    inline MAd_Vector<SCALAR> operator= (const MAd_Vector<SCALAR> &other)
    {
      if ( this != &other ) {
        r = other.r;
        if (data) delete [] data;
        data = new SCALAR[r];
        for (int i = 0; i < r; ++i) data[i] = other.data[i];
      }
      return *this;
    }
    inline SCALAR operator () (int i) const
    {
      return data[i];
    }
    inline SCALAR & operator () (int i)
    {
      return data[i];
    }
    inline double norm()
    {
      double n = 0.;
      for(int i = 0; i < r; ++i) n += data[i] * data[i];
      return sqrt(n);
    }
    inline void scale(const SCALAR s)
    {
      for (int i = 0; i < r; ++i) data[i] *= s;
    }
    inline void set_all(const double &m) 
    {
      for (int i = 0; i < r; ++i) data[i] = m;
    }
    inline void add(const MAd_Vector &v) 
    {
      for (int i = 0; i < r; ++i) data[i] += v(i);    
    }
    void print(std::string name) const
    {
      printf("Vector %s:\n  ",name.c_str());
      for (int i = 0; i < r; ++i) printf("%12.5E ",data[i]);
      printf("\n");
    }
  };

  // -------------------------------------------------------------------
  template <class SCALAR>
  class MAd_Matrix
  {
  private:
    int _r, _c; // r = nb rows, c = nb columns
    SCALAR *_data;

  public:

    MAd_Matrix(int R,int C) : _r(R), _c(C)
    {
      _data = new SCALAR[_r * _c];
      scale(0.);
    }
    MAd_Matrix(const MAd_Matrix<SCALAR> &other) : _r(other._r), _c(other._c)
    {
      _data = new double[_r * _c];
      memcpy(other);
    }
    MAd_Matrix() : _r(0), _c(0), _data(0) {}
    ~MAd_Matrix() { delete [] _data; }

  public:

    void print()const {
			printf("matrix(%d,%d)",_r,_c);
			for(int i=0;i<_r;i++){
				for(int j=0;j<_c;j++){
					printf("%0.3e\t",_data[i + _r * j]);
				}
				printf("\n");
			}
			printf("--------\n");
		}

    inline int size1() const { return _r; }
    inline int size2() const { return _c; }

    inline SCALAR operator () (int i, int j) const
    {
      return _data[i + _r * j];
    }
    inline SCALAR & operator () (int i, int j)
    {
      return _data[i + _r * j];
    }

    MAd_Matrix<SCALAR> & operator = (const MAd_Matrix<SCALAR> &other)
    {
      if(this != &other){
        _r = other._r; 
        _c = other._c;
        _data = new SCALAR[_r * _c];
        memcpy(other);
      }
      return *this;
    }
    void memcpy(const MAd_Matrix &other)
    {
      for (int i = 0; i < _r * _c; ++i) _data[i] = other._data[i];
    }

    void copy(const MAd_Matrix<SCALAR> &a, int i0, int ni, int j0, int nj, 
              int desti0, int destj0)
    {
      for(int i = i0, desti = desti0; i < i0 + ni; i++, desti++)
        for(int j = j0, destj = destj0; j < j0 + nj; j++, destj++)
          (*this)(desti, destj) = a(i, j);
    }

    // c = c + data * b
    inline void mult(const MAd_Matrix<SCALAR> &b, MAd_Matrix<SCALAR> &c)
    {
      c.scale(0.);
      for(int i = 0; i < _r; i++)
        for(int j = 0; j < b.size2(); j++)
          for(int k = 0; k < _c; k++)
            c._data[i + _r * j] += (*this)(i, k) * b(k, j);
    }

    // y = y + data * x
    inline void mult(const MAd_Vector<SCALAR> &x, MAd_Vector<SCALAR> &y)
    {
      y.scale(0.);
      for(int i = 0; i < _r; i++)
        for(int j = 0; j < _c; j++)
          y._data[i] += (*this)(i, j) * x(j);
    }

    // data = alpha * ( a * b ) + beta * data
    inline void blas_dgemm(MAd_Matrix<SCALAR> &a, MAd_Matrix<SCALAR> &b, 
                           SCALAR alpha=1., SCALAR beta=1.)
    {
      MAd_Matrix<SCALAR> temp(a.size1(), b.size2()); // temp = 0;
      a.mult(b, temp); // temp = a * b
      temp.scale(alpha); // temp = alpha * ( a * b )
      scale(beta);
      add(temp);
    }

    inline void set_all(const double &m) 
    {
      for (int i = 0; i < _r * _c; ++i) _data[i] = m;
    }
    inline void setValues(const SCALAR *M[]) 
    {
      for (int i = 0; i < _r ; ++i) 
        for (int j = 0; j < _c ; ++j) _data[i + _r * j] = M[i][j];
    }
    inline void scale(const double s)
    {
      if(s == 0.)
        for(int i = 0; i < _r * _c; ++i) _data[i] = 0.;
      else
        for(int i = 0; i < _r * _c; ++i) _data[i] *= s;
    }
    inline void add(const double &a) 
    {
      for(int i = 0; i < _r * _c; ++i) _data[i] += a;
    }
    inline void add(const MAd_Matrix<SCALAR> &m) 
    {
      for(int i = 0; i < size1(); i++)
        for(int j = 0; j < size2(); j++)
          (*this)(i, j) += m(i, j);
    }

    inline MAd_Matrix<SCALAR> transpose()
    {
      MAd_Matrix<SCALAR> T(size2(), size1());
      for(int i = 0; i < size1(); i++)
        for(int j = 0; j < size2(); j++)
          T(j, i) = (*this)(i, j);
      return T;
    }
    inline bool lu_solve(const MAd_Vector<SCALAR> &rhs, MAd_Vector<SCALAR> &result)
    {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,"LU factorization requires LAPACK");
      return false;
    }
    MAd_Matrix<SCALAR> cofactor(int i, int j) const 
    {
      int ni = size1();
      int nj = size2();
      MAd_Matrix<SCALAR> cof(ni - 1, nj - 1);
      for(int I = 0; I < ni; I++){
        for(int J = 0; J < nj; J++){
          if(J != j && I != i)
            cof(I < i ? I : I - 1, J < j ? J : J - 1) = (*this)(I, J);
        }
      }
      return cof;
    }
    SCALAR determinant() const
    {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,"Determinant computation requires LAPACK");
      return 0.;
    }
    bool svd(MAd_Matrix<SCALAR> &V, MAd_Vector<SCALAR> &S)
    {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,"Singular value decomposition requires LAPACK");
      return false;
    } 
    bool eig(MAd_Matrix<double> &VL, // left eigenvectors 
             MAd_Vector<double> &DR, // Real part of eigenvalues
             MAd_Vector<double> &DI, // Im part of eigenvalues
             MAd_Matrix<double> &VR,
             bool sortRealPart=false )
    {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,"Eigen vectors computation requires LAPACK");
      return false;
    } 
    void print(std::string name) const
    {
      printf("Matrix %s (%d, %d):\n  ",name.c_str(),_r,_c);
      for(int i = 0; i < _r ; ++i) {
        for(int j = 0; j < _c ; ++j) printf("%12.5E ",(*this)(i, j));
        printf("\n");
      }
    }
  };


// -------------------------------------------------------------------
// With BLAS / LAPACK
// -------------------------------------------------------------------

  // -------------------------------------------------------------------
  template <class SCALAR>
  class MAd_BLASLAPACK_Vector
  {
  private:
    int r;
  public:
    inline int size() const { return r; }
    SCALAR *data;
    ~MAd_BLASLAPACK_Vector() { delete [] data; }
    MAd_BLASLAPACK_Vector() : r(0)
    {
      data = 0;
    }
    MAd_BLASLAPACK_Vector(int R) : r(R)
    {
      data = new SCALAR[r];
      scale(0);
    }
    MAd_BLASLAPACK_Vector(const MAd_BLASLAPACK_Vector<SCALAR> &other) : r(other.r)
    {
      data = new double[r];
      for (int i = 0; i < r; ++i) data[i] = other.data[i];
    }
    inline SCALAR operator () (int i) const
    {
      return data[i];
    }
    inline SCALAR & operator () (int i)
    {
      return data[i];
    }
    inline double norm()
    {
      double n = 0.;
      for(int i = 0; i < r; ++i) n += data[i] * data[i];
      return sqrt(n);
    }
    inline void scale(const SCALAR s)
    {
      for (int i = 0; i < r; ++i) data[i] *= s;
    }
    inline void set_all(const double &m) 
    {
      for (int i = 0; i < r; ++i) data[i] = m;
    }
    inline void add(const MAd_BLASLAPACK_Vector &v) 
    {
      for (int i = 0; i < r; ++i) data[i] += v(i);    
    }
    void print(std::string name) const
    {
      printf("Vector %s:\n  ",name.c_str());
      for (int i = 0; i < r; ++i) printf("%12.5E ",data[i]);
      printf("\n");
    }
  };

  // -------------------------------------------------------------------
  template <class SCALAR>
  class MAd_BLASLAPACK_Matrix
  {
  private:
    int _r, _c; // r = nb rows, c = nb columns
    SCALAR *_data;

  public:

    MAd_BLASLAPACK_Matrix(int R,int C) : _r(R), _c(C)
    {
      _data = new SCALAR[_r * _c];
      scale(0.);
    }
    MAd_BLASLAPACK_Matrix(const MAd_BLASLAPACK_Matrix<SCALAR> &other) : _r(other._r), _c(other._c)
    {
      _data = new double[_r * _c];
      memcpy(other);
    }
    MAd_BLASLAPACK_Matrix() : _r(0), _c(0), _data(0) {}
    ~MAd_BLASLAPACK_Matrix() { delete [] _data; }

  public:

    inline int size1() const { return _r; }
    inline int size2() const { return _c; }

    inline SCALAR operator () (int i, int j) const
    {
      return _data[i + _r * j];
    }
    inline SCALAR & operator () (int i, int j)
    {
      return _data[i + _r * j];
    }

    MAd_BLASLAPACK_Matrix<SCALAR> & operator = (const MAd_BLASLAPACK_Matrix<SCALAR> &other)
    {
      if(this != &other){
        _r = other._r; 
        _c = other._c;
        _data = new SCALAR[_r * _c];
        memcpy(other);
      }
      return *this;
    }
    void memcpy(const MAd_BLASLAPACK_Matrix &other)
    {
      for (int i = 0; i < _r * _c; ++i) _data[i] = other._data[i];
    }

    void copy(const MAd_BLASLAPACK_Matrix<SCALAR> &a, int i0, int ni, int j0, int nj, 
              int desti0, int destj0)
    {
      for(int i = i0, desti = desti0; i < i0 + ni; i++, desti++)
        for(int j = j0, destj = destj0; j < j0 + nj; j++, destj++)
          (*this)(desti, destj) = a(i, j);
    }

    // c = c + data * b
    inline void mult(const MAd_BLASLAPACK_Matrix<SCALAR> &b, MAd_BLASLAPACK_Matrix<SCALAR> &c)
    {
#if defined(HAVE_BLAS)
      int M = c.size1(), N = c.size2(), K = _c;
      int LDA = _r, LDB = b.size1(), LDC = c.size1();
      double alpha = 1., beta = 0.;
      dgemm_("N", "N", &M, &N, &K, &alpha, _data, &LDA, b._data, &LDB, 
             &beta, c._data, &LDC);
#else
      c.scale(0.);
      for(int i = 0; i < _r; i++)
        for(int j = 0; j < b.size2(); j++)
          for(int k = 0; k < _c; k++)
            c._data[i + _r * j] += (*this)(i, k) * b(k, j);
#endif
    }

    // y = y + data * x
    inline void mult(const MAd_BLASLAPACK_Vector<SCALAR> &x, MAd_BLASLAPACK_Vector<SCALAR> &y)
    {
#if defined(HAVE_BLAS)
      int M = _r, N = _c, LDA = _r, INCX = 1, INCY = 1;
      double alpha = 1., beta = 0.;
      dgemv_("N", &M, &N, &alpha, _data, &LDA, x._data, &INCX,
             &beta, y._data, &INCY);
#else
      y.scale(0.);
      for(int i = 0; i < _r; i++)
        for(int j = 0; j < _c; j++)
          y._data[i] += (*this)(i, j) * x(j);
#endif
    }

    // data = alpha * ( a * b ) + beta * data
    inline void blas_dgemm(MAd_BLASLAPACK_Matrix<SCALAR> &a, MAd_BLASLAPACK_Matrix<SCALAR> &b, 
                           SCALAR alpha=1., SCALAR beta=1.)
    {
#if defined(HAVE_BLAS)
      int M = size1(), N = size2(), K = a.size2();
      int LDA = a.size1(), LDB = b.size1(), LDC = size1();
      dgemm_("N", "N", &M, &N, &K, &alpha, a._data, &LDA, b._data, &LDB, 
             &beta, _data, &LDC);
#else
      MAd_BLASLAPACK_Matrix<SCALAR> temp(a.size1(), b.size2()); // temp = 0;
      a.mult(b, temp); // temp = a * b
      temp.scale(alpha); // temp = alpha * ( a * b )
      scale(beta);
      add(temp);
#endif
    }

    inline void set_all(const double &m) 
    {
      for (int i = 0; i < _r * _c; ++i) _data[i] = m;
    }
    inline void setValues(const SCALAR *M[]) 
    {
      for (int i = 0; i < _r ; ++i) 
        for (int j = 0; j < _c ; ++j) _data[i + _r * j] = M[i][j];
    }
    inline void scale(const double s)
    {
      if(s == 0.)
        for(int i = 0; i < _r * _c; ++i) _data[i] = 0.;
      else
        for(int i = 0; i < _r * _c; ++i) _data[i] *= s;
    }
    inline void add(const double &a) 
    {
      for(int i = 0; i < _r * _c; ++i) _data[i] += a;
    }
    inline void add(const MAd_BLASLAPACK_Matrix<SCALAR> &m) 
    {
      for(int i = 0; i < size1(); i++)
        for(int j = 0; j < size2(); j++)
          (*this)(i, j) += m(i, j);
    }

    inline MAd_BLASLAPACK_Matrix<SCALAR> transpose()
    {
      MAd_BLASLAPACK_Matrix<SCALAR> T(size2(), size1());
      for(int i = 0; i < size1(); i++)
        for(int j = 0; j < size2(); j++)
          T(j, i) = (*this)(i, j);
      return T;
    }
    inline bool lu_solve(const MAd_BLASLAPACK_Vector<SCALAR> &rhs, MAd_BLASLAPACK_Vector<SCALAR> &result)
    {
#if defined(HAVE_LAPACK)
      int N = size1(), nrhs = 1, lda = N, ldb = N, info;
      int *ipiv = new int[N];
      for(int i = 0; i < N; i++) result(i) = rhs(i);
      dgesv_(&N, &nrhs, _data, &lda, ipiv, result.data, &ldb, &info);
      delete [] ipiv;
      if(info == 0) return true;
      MAdMsgSgl::instance().error(__LINE__,__FILE__,"Problem in LAPACK LU (info=%d)", info);
#else
      MAdMsgSgl::instance().error(__LINE__,__FILE__,"LU factorization requires LAPACK");
#endif
      return false;
    }
    MAd_BLASLAPACK_Matrix<SCALAR> cofactor(int i, int j) const 
    {
      int ni = size1();
      int nj = size2();
      MAd_BLASLAPACK_Matrix<SCALAR> cof(ni - 1, nj - 1);
      for(int I = 0; I < ni; I++){
        for(int J = 0; J < nj; J++){
          if(J != j && I != i)
            cof(I < i ? I : I - 1, J < j ? J : J - 1) = (*this)(I, J);
        }
      }
      return cof;
    }
    SCALAR determinant() const
    {
#if defined(HAVE_LAPACK)
      MAd_BLASLAPACK_Matrix<SCALAR> tmp(*this);
      int M = size1(), N = size2(), lda = size1(), info;
      int *ipiv = new int[std::min(M, N)];
      dgetrf_(&M, &N, tmp._data, &lda, ipiv, &info);
      SCALAR det = 1.;
      for(int i = 0; i < size1(); i++){
        det *= tmp(i, i);
        if(ipiv[i] != i + 1) det = -det;
      }
      return det;
#else
      MAdMsgSgl::instance().error(__LINE__,__FILE__,"Determinant computation requires LAPACK");
      return 0.;
#endif
    }
    bool svd(MAd_BLASLAPACK_Matrix<SCALAR> &V, MAd_BLASLAPACK_Vector<SCALAR> &S)
    {
#if defined(HAVE_LAPACK)
      MAd_BLASLAPACK_Matrix<SCALAR> VT(V.size2(), V.size1());
      int M = size1(), N = size2(), LDA = size1(), LDVT = VT.size1(), info;
      int LWORK = std::max(3 * std::min(M, N) + std::max(M, N), 5 * std::min(M, N));
      MAd_BLASLAPACK_Vector<SCALAR> WORK(LWORK);
      dgesvd_("O", "A", &M, &N, _data, &LDA, S._data, _data, &LDA,
              VT._data, &LDVT, WORK._data, &LWORK, &info);
      V = VT.transpose();
      if(info == 0) return true;
      MAdMsgSgl::instance().error(__LINE__,__FILE__,"Problem in LAPACK SVD (info=%d)", info);
#else
      MAdMsgSgl::instance().error(__LINE__,__FILE__,"Singular value decomposition requires LAPACK");
#endif
      return false;
    }
    bool eig(MAd_BLASLAPACK_Matrix<double> &VL, // left eigenvectors 
             MAd_BLASLAPACK_Vector<double> &DR, // Real part of eigenvalues
             MAd_BLASLAPACK_Vector<double> &DI, // Im part of eigenvalues
             MAd_BLASLAPACK_Matrix<double> &VR,
             bool sortRealPart=false ) // if true: sorted from max '|DR|' to min '|DR|'
#if defined(HAVE_LAPACK)
      ;
#else
    {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,"Eigen vectors computation requires LAPACK");
      return false;
    }
#endif
    void print(std::string name) const
    {
      printf("Matrix %s (%d, %d):\n  ",name.c_str(),_r,_c);
      for(int i = 0; i < _r ; ++i) {
        for(int j = 0; j < _c ; ++j) printf("%12.5E ",(*this)(i, j));
        printf("\n");
      }
    }
  };

// -------------------------------------------------------------------
// With GSL
// -------------------------------------------------------------------
#ifdef _HAVE_GSL_

  class MAd_GSL_Vector
  {
  private:
    int _r;
    gsl_vector *_data;
    friend class MAd_GSL_Matrix;

  public:

    MAd_GSL_Vector(int r) : _r(r)
    {
      _data = gsl_vector_calloc(_r);
    }
    MAd_GSL_Vector(const MAd_GSL_Vector &other) : _r(other._r)
    {
      _data = gsl_vector_calloc(_r);
      gsl_vector_memcpy(_data, other._data);
    }
    ~MAd_GSL_Vector() { gsl_vector_free(_data); }

    inline int size() const { return _r; }

    inline double operator () (int i) const
    {
      return gsl_vector_get(_data, i);
    }
    inline double & operator () (int i)
    {
      return *gsl_vector_ptr(_data, i);
    }
    inline double norm()
    {
      return gsl_blas_dnrm2(_data);
    }
    inline void scale(const double s)
    {
      if(s == 0.) gsl_vector_set_zero(_data);
      else gsl_vector_scale(_data, s);
    }
    inline void set_all(const double &m) 
    {
      gsl_vector_set_all(_data, m);
    }
    inline void add(const MAd_GSL_Vector &v) 
    {
      gsl_vector_add (_data, v._data);
    }
    void print(std::string name) const
    {
      printf("Vector %s:\n  ",name.c_str());
      for (int i = 0; i < _r; ++i) printf("%12.5E ",gsl_vector_get(_data, i));
      printf("\n");
    }
  };

  // -------------------------------------------------------------------
  class MAd_GSL_Matrix
  {
  private:
    gsl_matrix *_data;
  
  public:
    MAd_GSL_Matrix(int r, int  c) { _data = gsl_matrix_calloc(r, c); }
    MAd_GSL_Matrix(const MAd_GSL_Matrix &other) : _data(0)
    {
      if(_data) gsl_matrix_free(_data);
      _data = gsl_matrix_calloc(other._data->size1, other._data->size2);
      gsl_matrix_memcpy(_data, other._data);
    }
    MAd_GSL_Matrix() : _data(0) {}
    ~MAd_GSL_Matrix() { if(_data && _data->owner == 1) gsl_matrix_free(_data); }

  public:
    inline int size1() const { return _data->size1; }
    inline int size2() const { return _data->size2; }
    MAd_GSL_Matrix & operator = (const MAd_GSL_Matrix &other)
    {
      if(&other != this){
        if(_data) gsl_matrix_free(_data);
        _data = gsl_matrix_calloc(other._data->size1, other._data->size2);
        gsl_matrix_memcpy(_data, other._data);
      }
      return *this;
    }
    void memcpy(const MAd_GSL_Matrix &other)
    {
      gsl_matrix_memcpy(_data, other._data);
    }
    inline double operator () (int i, int j) const
    {
      return gsl_matrix_get(_data, i, j);
    }
    inline double & operator () (int i, int j)
    {
      return *gsl_matrix_ptr(_data, i, j);
    }
    void copy(const MAd_GSL_Matrix &a, int i0, int ni, int j0, int nj, 
              int desti0, int destj0)
    {
      for(int i = i0, desti = desti0; i < i0 + ni; i++, desti++)
        for(int j = j0, destj = destj0; j < j0 + nj; j++, destj++)
          (*this)(desti, destj) = a(i, j);
    }
    inline void mult(const MAd_GSL_Matrix &b, MAd_GSL_Matrix &c)
    {
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., _data, b._data, 0., c._data);
    }
    inline void mult(const MAd_GSL_Vector &x, MAd_GSL_Vector &y)
    {
      gsl_blas_dgemv(CblasNoTrans, 1., _data, x._data, 0., y._data);
    }
    inline void blas_dgemm(MAd_GSL_Matrix &a, MAd_GSL_Matrix &b, 
                           double alpha=1., double beta=1.)
    {      
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, alpha, a._data, b._data, beta, _data);
    }
    inline void set_all(const double &m) 
    {
      gsl_matrix_set_all(_data, m);
    }
    inline void setValues(const double *M[]) 
    {
      for (int i = 0; i < size1() ; ++i) 
        for (int j = 0; j < size2() ; ++j) *gsl_matrix_ptr(_data, i, j) = M[i][j];
    }
    inline void scale(const double s) 
    {
      if(s == 0.) gsl_matrix_set_zero(_data);
      else gsl_matrix_scale(_data, s);
    }
    inline void add(const double &a) 
    {
      gsl_matrix_add_constant(_data, a);
    }
    inline void add(const MAd_GSL_Matrix &m) 
    {
      gsl_matrix_add(_data, m._data);
    }
    inline MAd_GSL_Matrix transpose()
    {
      MAd_GSL_Matrix T(size2(), size1());
      for(int i = 0; i < size1(); i++)
        for(int j = 0; j < size2(); j++)
          T(j, i) = (*this)(i, j);
      return T;
    }
    inline bool lu_solve(const MAd_GSL_Vector &rhs, MAd_GSL_Vector &result)
    {
      int s;
      gsl_permutation *p = gsl_permutation_alloc(size1());
      gsl_linalg_LU_decomp(_data, p, &s);
      gsl_linalg_LU_solve(_data, p, rhs._data, result._data);
      gsl_permutation_free(p);
      return true;
    }
    MAd_GSL_Matrix cofactor(int i, int j) const 
    {
      int ni = size1();
      int nj = size2();
      MAd_GSL_Matrix cof(ni - 1, nj - 1);
      for(int I = 0; I < ni; I++){
        for(int J = 0; J < nj; J++){
          if(J != j && I != i)
            cof(I < i ? I : I - 1, J < j ? J : J - 1) = (*this)(I, J);
        }
      }
      return cof;
    }
    double determinant() const 
    {
      MAd_GSL_Matrix tmp = *this;
      gsl_permutation *p = gsl_permutation_alloc(size1());
      int s;
      gsl_linalg_LU_decomp(tmp._data, p, &s);
      gsl_permutation_free(p);
      return gsl_linalg_LU_det(tmp._data, s);
    } 
    bool svd(MAd_GSL_Matrix &V, MAd_GSL_Vector &S)
    {
      MAd_GSL_Vector tmp(S.size());
      gsl_linalg_SV_decomp(_data, V._data, S._data, tmp._data);
      return true;
    }
    inline void invert ()
    {
      int s;
      gsl_permutation *p = gsl_permutation_alloc (size1());
      gsl_linalg_LU_decomp(_data, p, &s);
      gsl_matrix *data_inv = gsl_matrix_calloc(size1(), size2());
      gsl_linalg_LU_invert(_data, p, data_inv) ;
      gsl_matrix_memcpy(_data, data_inv);
      gsl_matrix_free(data_inv);
      gsl_permutation_free(p);
    }
    inline bool invertSecure(double &det)
    {
      int s;
      gsl_permutation *p = gsl_permutation_alloc (size1());
      gsl_linalg_LU_decomp(_data, p, &s);
      det = gsl_linalg_LU_det(_data, s);
      gsl_matrix *data_inv = gsl_matrix_calloc(size1(), size2());
      gsl_linalg_LU_invert(_data, p, data_inv);
      gsl_matrix_memcpy(_data, data_inv);
      gsl_matrix_free(data_inv);
      gsl_permutation_free(p);
      return (det != 0.);
    }
    bool eig(MAd_GSL_Matrix &VL, // left eigenvectors 
             MAd_GSL_Vector &DR, // Real part of eigenvalues
             MAd_GSL_Vector &DI, // Im part of eigenvalues
             MAd_GSL_Matrix &VR,
             bool sortRealPart=false )
    {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,"Eigen vectors computation requires LAPACK");
      return false;
    }
    void print(std::string name) const
    {
      printf("Matrix %s (%d, %d):\n  ",name.c_str(),_data->size1,_data->size2);
      for(int i = 0; i < _data->size1 ; ++i) {
        for(int j = 0; j < _data->size2 ; ++j) printf("%12.5E ",(*this)(i, j));
        printf("\n");
      }
    }
  };

#endif

  // -------------------------------------------------------------------
  // -------------------------------------------------------------------

#if defined(HAVE_LAPACK)
  typedef MAd_BLASLAPACK_Matrix<double> doubleMatrix;
  typedef MAd_BLASLAPACK_Vector<double> doubleVector;
#else
 #if defined(_HAVE_GSL_)
   typedef MAd_GSL_Matrix doubleMatrix;
   typedef MAd_GSL_Vector doubleVector;
 #else
   typedef MAd_Matrix<double> doubleMatrix;
   typedef MAd_Vector<double> doubleVector;
 #endif
#endif

  // -------------------------------------------------------------------

  // Should be used for operations on small matrices
  typedef MAd_Matrix<double> smallMatrix;
  typedef MAd_Vector<double> smallVector;

  // -------------------------------------------------------------------

}

#endif
