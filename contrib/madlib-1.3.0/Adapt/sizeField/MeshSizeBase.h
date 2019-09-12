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

#ifndef _H_MESHSIZEBASE
#define _H_MESHSIZEBASE

#include "MAdMetric.h"

#include <string>

namespace MAd {

  // -------------------------------------------------------------------
  enum MeshSizeType {
    ISOTROPIC,
    ANISOTROPIC
  };

  // -------------------------------------------------------------------
  typedef class MeshSizeBase * pMSize;

  // -------------------------------------------------------------------
  pMSize MS_copy(const pMSize pMS);
  pMSize MS_intersect(const pMSize pMS0, const pMSize pMS1);
  // interpolate on an edge
  pMSize MS_interpolate(const pMSize, const pMSize, double);
  // interpolate on a triangle
  pMSize MS_interpolate(const pMSize, const pMSize, const pMSize, 
                        double, double);
  // interpolate on a tetrahedron
  pMSize MS_interpolate(const pMSize, const pMSize, const pMSize, const pMSize, 
                        double, double, double);

  // -------------------------------------------------------------------
  // This class stores a mesh size (prescribed edge length)
  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  class MeshSizeBase
  {
  public:
  
    MeshSizeBase() {};
    virtual ~MeshSizeBase() {};

  public:

    virtual MeshSizeType getType() const = 0;
    virtual MAdMetric getMetric() const = 0;

    // this becomes the minimum size between two sizes
    virtual void   intersect(const pMSize pMS0, const pMSize pMS1) = 0;

    // this becomes the interpolation between two sizes
    virtual void   interpolate(const pMSize, const pMSize, double) = 0;
    // this becomes the interpolation between three sizes
    virtual void   interpolate(const pMSize, const pMSize, const pMSize, 
                               double, double) = 0;
    // this becomes the interpolation between four sizes
    virtual void   interpolate(const pMSize, const pMSize, const pMSize, const pMSize, 
                               double, double, double) = 0;

    // get a principle direction (eigenvectors of the metric)
    virtual double direction(int i, double dir[3]) const = 0;

    // get the desired size in a principle direction (eigenvalues)
    virtual double size(int i=0) const = 0;

    // get the desired size in the 3 directions (eigenvalues)
    virtual void   sizes(double _h[3]) const = 0;

    // get the square norm of the vector in the metric
    virtual double normSq(const double vec[3]) const = 0;

    // get the desired size in the given direction
    virtual double lengthSqInDir(const double dir[3]) const = 0;

    // get the mean, min, max sizes
    virtual double getMeanLength() const = 0;
    virtual double getMinLength()  const = 0;
    virtual double getMaxLength()  const = 0;

    // scale the size
    virtual void   scale(double) = 0;

    // print info
    virtual void   print(std::string name="") const = 0;

  };

  // -------------------------------------------------------------------

}

#endif
