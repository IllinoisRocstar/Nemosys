/* Implementation of RBF interpolation class */

#include "baseInterp.H"

/*
   interpolates values for given point coordinates
   input:
      ni : number of interpolation points
      xi[nDim*ni] : point coordinates

   output:
      returns fi[ni] interpolation values
*/
void basicInterpolant::interpolate(int ni, std::vector<double>& xi, 
              std::vector<double>& pntData, std::vector<double>& newPntData)
{
  // calculating neighbouring indices and weights 
  if (!wCalced) {
    pntNibIdx = new int[ni*nNib];
    w = new double[ni*nNib];
    for (int iPnt=0; iPnt<ni; iPnt++)
    {
      ANNpoint     qryPnt;
      ANNidxArray  nnIdx;
      ANNdistArray dists;
      qryPnt = annAllocPt(nDim);
      qryPnt[0] = xi[iPnt*nDim];
      qryPnt[1] = xi[iPnt*nDim+1];
      qryPnt[2] = xi[iPnt*nDim+2];
      nnIdx  = new ANNidx[nNib];
      dists  = new ANNdist[nNib];
      kdTree->annkSearch(qryPnt, nNib, nnIdx, dists);
      double totW = 0.0;
      // correcting for zero distance points
      for (int iNib=0; iNib<nNib; iNib++)
      {
        if (dists[iNib] == 0) dists[iNib] = 1.0e-14;
        totW +=1.0/dists[iNib];
      }
      for (int iNib=0; iNib<nNib; iNib++) 
      {
        pntNibIdx[iPnt*nNib+iNib] = nnIdx[iNib];
        w[iPnt*nNib+iNib] = (1.0/dists[iNib])/totW;
      }
     }
     wCalced = true;
  }  
  // performing the interpolation
  for (int iPnt=0; iPnt<ni; iPnt++)
  {
    newPntData.push_back(0.0);
    for (int iNib=0; iNib<nNib; iNib++)
      newPntData[iPnt] += pntData[pntNibIdx[iPnt*nNib+iNib] ] 
                   *w[iPnt*nNib+iNib];
  }
}

/* Builds kd-Tree */
void basicInterpolant::buildPointKDTree(std::vector<double>& pntCrds)
{
  ANNpointArray pntCrd;
  pntCrd = annAllocPts(nPnt, nDim);
  // clearing old instance
  if (kdTree)
    delete kdTree;
  // filling up vertex coordinate array for the current mesh
  for (int iPnt=0; iPnt<nPnt; iPnt++)
  {
     pntCrd[iPnt][0] = pntCrds[iPnt*nDim];
     pntCrd[iPnt][1] = pntCrds[iPnt*nDim+1];
     pntCrd[iPnt][2] = pntCrds[iPnt*nDim+2];
  }
  // building kdTree
  kdTree = new ANNkd_tree(pntCrd, nPnt, nDim);
  treeExist = true;
}
