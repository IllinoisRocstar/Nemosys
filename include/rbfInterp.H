/* Radial Basis Function (RBF) interpolation class declatation
   for n-dimensional sapce. Based on the rbf_interp_nd implem-
   ntation provided in:
   https://people.sc.fsu.edu/~jburkardt/cpp_src/rbf_interp_nd/rbf_interp_nd.html
*/
#ifndef _RBFINTERP_H_
#define _RBFINTERP_H_

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstring>

using namespace std;

#include "rbf_interp_nd.hpp"
#include "r8lib.hpp"

// types
enum rbf_interp_t {MULQUAD, INVMULQUAD, THNSPLINE, GAUSS};


class RBFInterpolant {

// public members
public:
   RBFInterpolant(int dim, int nPnt, rbf_interp_t type, double r0):
   dim(dim), nPnt(nPnt), type(type), wCalced(false), r0(r0)
   {};
   ~RBFInterpolant()
    {
      if (w)
        delete [] w;
    };

    void setPointCoords(double* inPntCrds);
    void setPointData(double* inPntData);
    double* interpolate(int ni, double* xi);

// private members
private:
    void calcWeights();
   
// private members   
private:
  int dim;              // space dimensions
  int nPnt;             // number of data points
  rbf_interp_t type;    // type of RBF interpolation
  double r0;            // the r0 parameter for interpolation
  double* pntCrds;      // array of point coordinates
  double* pntData;      // array of point data
  double* w;            // weights for interpolation
  bool wCalced;         // weight calculation switch
};

#endif
