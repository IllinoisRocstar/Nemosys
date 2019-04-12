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

#ifndef _H_DISTANCETOPOINTS
#define _H_DISTANCETOPOINTS

#include <stdio.h>

#ifdef _HAVE_ANN_
#include "ANN.h"
#endif

namespace MAd {

  // -------------------------------------------------------------------
  class MAd_searchTool {

  public:
    MAd_searchTool() {
      printf("Error: no search tool implemented: flag _HAVE_ANN_ not set\n");
      throw;
    }
    void reset() { throw; }
    void allocatePoints(int n) { throw; }
    void addPoint(int index, double xyz[3]) { throw; }
    void allocateTree(int n) { throw; }
    double computeDistanceSq(const double xyz[3], int * id=NULL) const { throw; }
  };

  // -------------------------------------------------------------------
#ifdef _HAVE_ANN_
  class ANN_searchTool {

  public:

    ANN_searchTool()
    {
      kdTree = NULL;
      points = NULL;
    }

    ~ANN_searchTool()
    {
      if (kdTree) { delete kdTree; kdTree = NULL; }
      if (points) { annDeallocPts(points); points = NULL; }
      annClose();
    }

    void reset()
    {
      if (kdTree) { delete kdTree; kdTree = NULL; }
      if (points) { annDeallocPts(points); points = NULL; }
    }

    void allocatePoints(int n)
    {
      if (points) annDeallocPts(points);
      points = annAllocPts(n,3);
    }

    void addPoint(int index, double xyz[3])
    {
      for (int i=0; i<3; i++) points[index][i] = xyz[i];
    }

    void allocateTree(int n)
    {
      if (kdTree) delete kdTree;
      kdTree = new ANNkd_tree(points, n, 3);
    }

    double computeDistanceSq(const double xyz[3], int * id=NULL) const
    {
      if (!kdTree) {
        printf("Error: no research tree built\n");
        throw;
      }

      double xyz_copy[3]; // because the first argument of annkSearch is not const
      for (int i=0; i<3; i++) xyz_copy[i] = xyz[i];
      int maxpts = 1;
      // the use of allocation is done due to ms visual c++ compiler 
      // that do not support c99 standard (it uses the c95 et c++98 standards)
      ANNidx*  index= new ANNidx[maxpts];
      ANNdist* distSq= new ANNdist[maxpts];
      kdTree->annkSearch(xyz_copy, maxpts, index, distSq);
      double theDistSq = distSq[0];
      if ( id ) *id = index[0];
      delete [] index;
      delete [] distSq;
      return theDistSq;
    }

    // Compute gradient of the square distance to the cloud of points 
    // by non-centered differences. Central differences would lead to 
    // zero gradients for points located on the walls (distance is not signed).
    double gradDistanceSq(const double xyz[3], double grad[3]) const
    {
      double dSq = computeDistanceSq(xyz);
      double eps = 1.e-3;
      double dEps[2];

      int maxpts = 1;
      ANNidx*  index= new ANNidx[maxpts];
      ANNdist* distSq= new ANNdist[maxpts];
      double tmp[3];
      tmp[0] = xyz[0]; tmp[1] = xyz[1]; tmp[2] = xyz[2];
      for (int i=0; i<3; i++) {
        tmp[i] += eps;
        kdTree->annkSearch(tmp, maxpts, index, distSq);
        dEps[0] = distSq[0];
        tmp[i] -= 2*eps;
        kdTree->annkSearch(tmp, maxpts, index, distSq);
        dEps[1] = distSq[1];
        tmp[i] += eps;
        grad[i] = 1000. * ( dEps[1] - dEps[0] );
      }
      
      return dSq;
    }

    private:

      ANNkd_tree * kdTree;
      ANNpoint * points;
  };

  typedef ANN_searchTool SearchTool;

  // -------------------------------------------------------------------
#else

  typedef MAd_searchTool SearchTool;

  // -------------------------------------------------------------------
#endif

  // -------------------------------------------------------------------

}

#endif
