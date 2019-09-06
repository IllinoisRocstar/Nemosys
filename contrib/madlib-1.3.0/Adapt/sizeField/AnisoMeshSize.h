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

#ifndef _H_ANISOMESHSIZE
#define _H_ANISOMESHSIZE

#include "MeshSizeBase.h"

namespace MAd {

  // -------------------------------------------------------------------
  // This class stores an anisotropic mesh size (prescribed edge length)
  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  class AnisoMeshSize : public MeshSizeBase {

  public:
  
    AnisoMeshSize(double e[3][3], double _h[3]);
    AnisoMeshSize(double e[3], double hDir, double hTg);
    AnisoMeshSize(double h=1.);
    AnisoMeshSize(const AnisoMeshSize &);
    ~AnisoMeshSize() {};

  public:

    MeshSizeType getType() const;
    MAdMetric getMetric() const { return M; }

    void   intersect(const pMSize pMS0, const pMSize pMS1);

    void   interpolate(const pMSize, const pMSize, double);
    void   interpolate(const pMSize, const pMSize, const pMSize, 
                       double, double);
    void   interpolate(const pMSize, const pMSize, const pMSize, const pMSize, 
                       double, double, double);

    double direction(int i, double dir[3]) const;
    double size(int i=0) const;
    void   sizes(double _h[3]) const;

    // get the square of the norm in the metric
    double normSq(const double vec[3]) const;

    double lengthSqInDir(const double dir[3]) const ;
    
    // get the cosine of the angle with the direction of the minimal size
    double angleWithDir0(const double dir[3]) const;

    double getMeanLength() const;
    double getMinLength() const;
    double getMaxLength() const;

    void   scale(int, double);
    void   scale(double);

    void   print(std::string name="") const;
  
  private:

    MAdMetric M;
  };

  // -------------------------------------------------------------------

}

#endif
