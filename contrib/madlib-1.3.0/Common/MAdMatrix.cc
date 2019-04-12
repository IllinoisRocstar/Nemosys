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

#include "MAdMatrix.h"
#include <stdio.h>

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

namespace MAd {

#if defined(HAVE_LAPACK)

  template<> 
  bool MAd_BLASLAPACK_Matrix<double>::eig(MAd_BLASLAPACK_Matrix<double> &VL, // left eigenvectors 
                                          MAd_BLASLAPACK_Vector<double> &DR, // Real part of eigenvalues
                                          MAd_BLASLAPACK_Vector<double> &DI, // Im part of eigenvalues
                                          MAd_BLASLAPACK_Matrix<double> &VR,
                                          bool sortRealPart ) // if true: sorted from max '|DR|' to min '|DR|'
  {
    int N = size1(), info;
    int LWORK = 10*N;
    double * work = new double[LWORK];
      
    dgeev_("V","V",
           &N,_data,
           &N,DR.data,DI.data,
           VL._data,&N,
           VR._data,&N,
           work,&LWORK,&info);
      
    delete [] work;
      
    if(info == 0)
      {
        if (sortRealPart) {
          double tmp[8];
          // do permutations
          for (int i=0; i<(size1()-1); i++) {
            int maxR = i;
            for (int j=i+1; j<size1(); j++) if ( fabs(DR(j)) > fabs(DR(maxR)) ) maxR = j;
            if ( maxR != i )
              {
                tmp[0] = DR(i); tmp[1] = DI(i);
                tmp[2] = VL(0,i); tmp[3] = VL(1,i); tmp[4] = VL(2,i);
                tmp[5] = VR(0,i); tmp[6] = VR(1,i); tmp[7] = VR(2,i);
                  
                DR(i) = DR(maxR); DI(i) = DI(maxR);
                VL(0,i) = VL(0,maxR); VL(1,i) = VL(1,maxR); VL(2,i) = VL(2,maxR);
                VR(0,i) = VR(0,maxR); VR(1,i) = VR(1,maxR); VR(2,i) = VR(2,maxR);
                  
                DR(maxR) = tmp[0]; DI(maxR) = tmp[1];
                VL(0,maxR) = tmp[2]; VL(1,maxR) = tmp[3]; VL(2,maxR) = tmp[4];
                VR(0,maxR) = tmp[5]; VR(1,maxR) = tmp[6]; VR(2,maxR) = tmp[7];
              }
          }
        }
        return true;
      }
    if(info > 0)
      MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                  "QR Algorithm failed to compute all the eigenvalues (info = %d)", info);
    else
      MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                  "Wrong %d-th argument in eig (info = %d)", -info);

    return false;
  }

#endif
}
