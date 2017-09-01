/* 
  Basic interpolation class declatation for 3-nDimensional sapce.
*/
#ifndef _BASEINTERP_H_
#define _BASEINTERP_H_

// standard
#include<vector>

// third party
#include "ANN.h"

// types

class basicInterpolant {

// public members
public:
   basicInterpolant(int nDim, int nPnt, int nNib, std::vector<double>& pntCrds):
   nDim(nDim), nPnt(nPnt), nNib(nNib), wCalced(false), treeExist(false),
   w(NULL), pntNibIdx(NULL), kdTree(NULL)
   {
     buildPointKDTree(pntCrds);
   };
   ~basicInterpolant()
    {
      if (w) delete [] w;
      if (pntNibIdx) delete [] pntNibIdx;
      if (kdTree) delete kdTree; 
    };
   
    void setPointData(std::vector<double>& inPntData);
    void interpolate(int ni, std::vector<double>& xi, 
                     std::vector<double>& pntData, 
                     std::vector<double>& newPntData,
                     int verb = 0);

    void interpolate(int ni, std::vector<double>& xi, 
                     std::vector<double>& pntData, 
                     std::vector<double>& newPntData,
                     double tol, int verb = 0);
    void clearCache() 
    {wCalced = false;};

// private members
private:
    void calcWeights();
    void buildPointKDTree(std::vector<double>& pntCrds);
   
// private members   
private:
  int nDim;              // space dimensions
  int nPnt;             // number of data points
  int nNib;             // number of neighbours used for the interpolation
  double* w;            // weights for interpolation
  bool wCalced;         // weight calculation switch
  int* pntNibIdx;
  // search support data structures
  bool treeExist;
  ANNkd_tree* kdTree;
};

#endif
