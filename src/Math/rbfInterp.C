/* Implementation of RBF interpolation class */

#include "rbfInterp.H"

/* 
   sets data point coordinates 
   input:
      inPntCrds[dim*nPnt] : coordinates of the data points
*/
void RBFInterpolant::setPointCoords(double* inPntCrds)
{
  pntCrds = inPntCrds;
}

/* 
   sets data point values 
   input:
      inPntdata[nPnt] :  values of the data points
*/
void RBFInterpolant::setPointData(double* inPntData)
{
  // evertime called weights should be upadted
  wCalced = false;
  pntData = inPntData;
}

/* 
  Calculate the weights for interpolation 
*/
void RBFInterpolant::calcWeights()
{
  // return if calculated
  if (wCalced)
    return;
  // calculate weights
  switch(type) {
  case MULQUAD:
     w = rbf_weight ( dim, nPnt, pntCrds, r0, phi1, pntData );
     break;
  case INVMULQUAD:
     w = rbf_weight ( dim, nPnt, pntCrds, r0, phi2, pntData );
     break;
  case THNSPLINE:
     w = rbf_weight ( dim, nPnt, pntCrds, r0, phi3, pntData );
     break;
  case GAUSS:
     w = rbf_weight ( dim, nPnt, pntCrds, r0, phi4, pntData );
     break;
  default:
     std::cerr << "Unknown interpolation type is not supported!\n";
  }
  // print weights
  r8vec_print ( nPnt, w, "  Weight vector:" );

  wCalced = true;
  return;
}
/*
   interpolates values for given point coordinates
   input:
      ni : number of interpolation points
      xi[nDim*ni] : point coordinates
      fi[nPnt] : interpolation values
*/
double* RBFInterpolant::interpolate(int ni, double* xi)
{
  // calculate weights if needed
  if (!wCalced)
    calcWeights();
  // performa interpolation 
  switch(type) {
  case MULQUAD:
     return( rbf_interp_nd ( dim, nPnt, pntCrds, r0, phi1, w, ni, xi ) );
     break;
  case INVMULQUAD:
     return( rbf_interp_nd ( dim, nPnt, pntCrds, r0, phi2, w, ni, xi ) );
     break;
  case THNSPLINE:
     return( rbf_interp_nd ( dim, nPnt, pntCrds, r0, phi3, w, ni, xi ) );
     break;
  case GAUSS:
     return( rbf_interp_nd ( dim, nPnt, pntCrds, r0, phi4, w, ni, xi ) );
     break;
  default:
     std::cerr << "Unknown interpolation type is not supported!\n";
     return(NULL);
  }
}
