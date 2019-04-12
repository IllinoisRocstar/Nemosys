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

#ifndef _H_MATHUTILS
#define _H_MATHUTILS

namespace MAd {

  // -------------------------------------------------------------------
  // --- Vectors ---

  void   diffVec      (const double [3], const double [3], double [3]);
  double dotProd      (const double [3], const double [3]);
  void   crossProd    (const double [3], const double [3], double [3]);
  void   normalizeVec (const double [3], double[3]);
  void   printVec     (const double [3], const char * name="");
  bool   isNanVec     (const double [3]);

  // -------------------------------------------------------------------
  // --- Matrices ---

  // 3 x 3
  void   transpose    (double [3][3]);
  void   transpMat    (const double [3][3], double [3][3]);
  double detMat       (const double [3][3]);
  double traceMat     (const double [3][3]);
  double traceMatMat  (const double [3][3]);
  void   vecMat       (const double [3],    const double [3][3], double [3]);
  void   matVec       (const double [3][3], const double [3],    double [3]);
  double vecMatVec    (const double [3][3], const double [3]);
  void   matMat       (const double [3][3], const double [3][3], double [3][3]);
  double inverseMat   (const double [3][3], double [3][3]);
  double invertMat    (double [3][3]);
  bool   system       (const double [3][3], const double [3], double [3], double *);
  void   printMat     (const double [3][3], const char * name="");
  void   meanRow33    (const double [3][3], double [3]);
  void   meanRow43    (const double [4][3], double [3]);
  bool   isNanMat     (const double [3][3]);

  // 2 x 2
  double inverseMat22 (const double [2][2], double [2][2]);

  // -------------------------------------------------------------------

  void sort(double [3]);
  void cubicRoots(const double [4], double [3], double [3]);
  int indexOfMin(double, double, double);
  int indexOfMax(double, double, double);
  double Tri_area(const double [3][3]);
  void Tri_linearParams(const double tri[3][3], const double xyz[3], double res[2]);
  double distToLineSq(const double p0[3], const double p1[3], 
                      const double xyz[3], double proj2pt[3], bool *onSegment=0);
  double distToTriangleSq(const double tri[3][3], const double xyz[3], double vec[3]);
  void pointToTriangle(const double tri[3][3], const double xyz[3], 
                       double proj[3], int * bit, double area[4]=0);
  void pointToPlane(const double [3][3], const double [3], double[3]);
  double distanceSq(const double[3], const double[3]);

  // -------------------------------------------------------------------

}

#endif
