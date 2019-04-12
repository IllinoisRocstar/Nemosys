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

#ifndef _H_ISOMESHSIZE
#define _H_ISOMESHSIZE

#include "MeshSizeBase.h"

namespace MAd {

  // -------------------------------------------------------------------
  // This class stores an isotropic mesh size (prescribed edge length)
  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  class IsoMeshSize : public MeshSizeBase {

  public:
  
    IsoMeshSize(double _h=-1.);
    IsoMeshSize(const IsoMeshSize &);
    ~IsoMeshSize() {};

  public:

    MeshSizeType getType() const;
    MAdMetric getMetric() const { return MAdMetric(1./(h*h)); }

    void   intersect(const pMSize pMS0, const pMSize pMS1);
    void   intersect(const pMSize pMS);

    void   interpolate(const pMSize, const pMSize, double);
    void   interpolate(const pMSize, const pMSize, const pMSize, 
                       double, double);
    void   interpolate(const pMSize, const pMSize, const pMSize, const pMSize, 
                       double, double, double);

    double direction(int i, double dir[3]) const;
    double size(int i=0) const;
    void   sizes(double _h[3]) const;
    double normSq(const double vec[3]) const;
    double lengthSqInDir(const double dir[3]) const ;

    double getMeanLength() const;
    double getMinLength() const;
    double getMaxLength() const;

    void   setSize(double);

    void   scale(int, double);
    void   scale(double);

    void   print(std::string name="") const;
  
  private:

    double h;
  };

  // -------------------------------------------------------------------

}

#endif
