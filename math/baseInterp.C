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
              std::vector<double>& pntData, std::vector<double>& newPntData,
              int verb)
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
    
      // using tolerance to exclude points with dist>tol
      // by setting them to -1
      //for (int i = 0; i < nNib; ++i) {
      //  if (dists[i] > tol)

      double totW = 0.0;
      // searching for zero-distance points
      int iNibZeroDist = -1;
      for (int iNib=0; iNib<nNib; iNib++)
        if (dists[iNib] == 0)
        {
          iNibZeroDist = iNib;
          break; // for loop
        }
      //std::cout << "iNibZeroDist = " << iNibZeroDist << std::endl;
      if (iNibZeroDist != -1) {
         // query point is repeating
        for (int iNib=0; iNib<nNib; iNib++) {
          pntNibIdx[iPnt*nNib+iNib] = nnIdx[iNib];
          w[iPnt*nNib+iNib] = 0.0;
        }
         w[iPnt*nNib + iNibZeroDist] = 1.0;
      } 
      else {
        // query point is not repeating
        for (int iNib=0; iNib<nNib; iNib++) {
          totW +=1.0/dists[iNib];
        }
        for (int iNib=0; iNib<nNib; iNib++) {
          pntNibIdx[iPnt*nNib+iNib] = nnIdx[iNib];
          w[iPnt*nNib+iNib] = (1.0/dists[iNib])/totW;
        }
      }
    }
     wCalced = true;
  }  
  // performing the interpolation
  for (int iPnt=0; iPnt<ni; iPnt++)
  {
    //newPntData.push_back(0.0);
    newPntData.resize(ni,0.0);
    for (int iNib=0; iNib<nNib; iNib++)
    {
      newPntData[iPnt] += pntData[pntNibIdx[iPnt*nNib+iNib] ] 
                   *w[iPnt*nNib+iNib];
      if (verb>0)
        std::cout << "Nib Indx = "
                  << pntNibIdx[iPnt*nNib+iNib]
                  << " weight = "
                  << w[iPnt*nNib+iNib]
                  << std::endl;
    }
  }
}

// overload to add distance tol to be passed by user
// checks if dists[i] > tol
void basicInterpolant::interpolate(int ni, std::vector<double>& xi, 
              std::vector<double>& pntData, std::vector<double>& newPntData,
              double tol, int verb)

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
    
      // using tolerance to exclude points with dist>tol
      // by setting them to -1
      for (int i = 0; i < nNib; ++i) {
        if (dists[i] > tol)
          dists[i] = -1;
      }
      double totW = 0.0;
      // searching for zero-distance points
      int iNibZeroDist = -1;
      for (int iNib=0; iNib<nNib; iNib++) {
        if (dists[iNib] == 0)
        {
          iNibZeroDist = iNib;
          break; // for loop
        }
      }
      //std::cout << "iNibZeroDist = " << iNibZeroDist << std::endl;
      if (iNibZeroDist != -1) {
         // query point is repeating
        for (int iNib=0; iNib<nNib; iNib++) {
          pntNibIdx[iPnt*nNib+iNib] = nnIdx[iNib];
          w[iPnt*nNib+iNib] = 0.0;
        }
         w[iPnt*nNib + iNibZeroDist] = 1.0;
      } 
      else {
        // query point is not repeating
        for (int iNib=0; iNib<nNib; iNib++) {
          // checking if distance is within tol
          if (dists[iNib] != -1)
            totW +=1.0/dists[iNib];
        }
        for (int iNib=0; iNib<nNib; iNib++) {
          pntNibIdx[iPnt*nNib+iNib] = nnIdx[iNib];
            if (dists[iNib] != -1) 
              w[iPnt*nNib+iNib] = (1.0/dists[iNib])/totW;
            // if dist to query outside tolerance, weight=0
            else
              w[iPnt*nNib+iNib] = 0.0;
        }
      }
    }
     wCalced = true;
  }  
  // performing the interpolation
  for (int iPnt=0; iPnt<ni; iPnt++)
  {
    //newPntData.push_back(0.0);
    newPntData.resize(ni,0.0);
    for (int iNib=0; iNib<nNib; iNib++)
    {
      newPntData[iPnt] += pntData[pntNibIdx[iPnt*nNib+iNib] ] 
                   *w[iPnt*nNib+iNib];
      if (verb>0)
        std::cout << "Nib Indx = "
                  << pntNibIdx[iPnt*nNib+iNib]
                  << " weight = "
                  << w[iPnt*nNib+iNib]
                  << std::endl;
    }
  }
}


// overload to add distance tol to be passed by user and exception handling
// checks if dists[i] > tol
// checks if N_tol(x) has all 0 data and throws exception 
// for plane points outside of inclusion
void basicInterpolant::interpolate(int ni, std::vector<double>& xi,
                                   std::vector<sphere>& spheres, 
                                   std::vector<double>& maskData,
                                   std::vector<double>& pntData, std::vector<double>& newPntData,
                                   double tol, int verb)

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

      // using tolerance to exclude points with dist>tol
      // by setting them to -1
      for (int iNib = 0; iNib < nNib; ++iNib) {
       
         if (dists[iNib] > tol)
          dists[iNib] = -1;
      }
      
      double totW = 0.0;
      // searching for zero-distance points
      int iNibZeroDist = -1;
      for (int iNib=0; iNib<nNib; iNib++) {
        if (dists[iNib] == 0)
        {
          iNibZeroDist = iNib;
          break; // for loop
        }
      }
      //std::cout << "iNibZeroDist = " << iNibZeroDist << std::endl;
      if (iNibZeroDist != -1) {
         // query point is repeating
        for (int iNib=0; iNib<nNib; iNib++) {
          pntNibIdx[iPnt*nNib+iNib] = nnIdx[iNib];
          w[iPnt*nNib+iNib] = 0.0;
        }
         w[iPnt*nNib + iNibZeroDist] = 1.0;
      } 
      else {
        // query point is not repeating
        for (int iNib=0; iNib<nNib; iNib++) {
          // checking if distance is within tol and weather point is in inclusion
          pntNibIdx[iPnt*nNib+iNib] = nnIdx[iNib];
          if (dists[iNib] != -1 || maskData[pntNibIdx[iPnt*nNib+iNib]] == 0.0)
            totW +=1.0/dists[iNib];
        }
        for (int iNib=0; iNib<nNib; iNib++) {
          //pntNibIdx[iPnt*nNib+iNib] = nnIdx[iNib];
            if (dists[iNib] != -1 || maskData[pntNibIdx[iPnt*nNib+iNib]] == 0.0) 
              w[iPnt*nNib+iNib] = (1.0/dists[iNib])/totW;
            // if dist to query outside tolerance, weight=0
            else
              w[iPnt*nNib+iNib] = 0.0;
        }
      }
    }
     wCalced = true;
  }


          
  // performing the interpolation
  for (int iPnt=0; iPnt<ni; iPnt++)
  {
    // if point on plane outside of inclusion is surrounded by neighbors inside 
    // inclusions, exception must be thrown
    std::vector<double> point;
    point.push_back(xi[iPnt*nDim]);
    point.push_back(xi[iPnt*nDim+1]);
    point.push_back(xi[iPnt*nDim+2]);
    bool in_sphere = false;
    for (int i = 0; i < spheres.size(); ++i) {
      if (in_sphere=spheres[i].in_sphere(point))
        break;
    }
    if (!in_sphere) {
      bool all_inclusions=true;
      for (int iNib = 0; iNib<nNib; ++iNib) {
        if (maskData[pntNibIdx[iPnt*nNib+iNib]] == 0.0) {
          all_inclusions=false;
          break;
        }
      }
      if (all_inclusions) {
        std::cerr << "All Neighbors of non-inclusion point " << iPnt << " are in inclusions!" << std::endl
                  << "Check point at: " << point[0] << " " << point[1] << " " << point[2] << std::endl
                  << "Refine RocLB mesh or use coarser planar mesh" << std::endl;
        exit(5);
      }
    } // finished check
    
    newPntData.resize(ni,0.0);
    for (int iNib=0; iNib<nNib; iNib++)
    {
      newPntData[iPnt] += pntData[pntNibIdx[iPnt*nNib+iNib] ] 
                   *w[iPnt*nNib+iNib];
      if (verb>0)
        std::cout << "Nib Indx = "
                  << pntNibIdx[iPnt*nNib+iNib]
                  << " weight = "
                  << w[iPnt*nNib+iNib]
                  << std::endl;
    }
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
