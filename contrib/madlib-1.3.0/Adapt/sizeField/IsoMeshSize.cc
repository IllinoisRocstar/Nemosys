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

#include "IsoMeshSize.h"
#include "MathUtils.h"
#include "MAdMessage.h"

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
using std::string;

namespace MAd {

  // -------------------------------------------------------------------
  IsoMeshSize::IsoMeshSize(double _h):
    MeshSizeBase()
  {
    h = _h;
  }

  // -------------------------------------------------------------------
  IsoMeshSize::IsoMeshSize(const IsoMeshSize &pm):
    MeshSizeBase()
  {
    h = pm.size();
  }

  // -------------------------------------------------------------------
  MeshSizeType IsoMeshSize::getType() const 
  {
    return ISOTROPIC;
  }

  // -------------------------------------------------------------------
  void IsoMeshSize::intersect(const pMSize pMS0, const pMSize pMS1)
  {
    if ( pMS0->getType() != ISOTROPIC || 
         pMS1->getType() != ISOTROPIC    ) {
      printf("Error: intersecting anisotropic sizes in isotropic mesh size\n");
      exit(1);
    }

    double h0 = ((IsoMeshSize*)pMS0)->size();
    double h1 = ((IsoMeshSize*)pMS1)->size();
    h = std::min(h0,h1);
  }

  // -------------------------------------------------------------------
  void IsoMeshSize::interpolate(const pMSize pMS0, const pMSize pMS1, 
                                double param)
  {
    if ( pMS0->getType() != ISOTROPIC || 
         pMS1->getType() != ISOTROPIC    ) {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,"Anisotropic size(s)");
    }

    double h0 = ((IsoMeshSize*)pMS0)->size();
    double h1 = ((IsoMeshSize*)pMS1)->size();
    h = h0 + param * ( h1 - h0 );
  }

  // -------------------------------------------------------------------
  void IsoMeshSize::interpolate(const pMSize pMS0, const pMSize pMS1, 
                                const pMSize pMS2, double u, double v)
  {
    if ( pMS0->getType() != ISOTROPIC || 
         pMS1->getType() != ISOTROPIC || 
         pMS2->getType() != ISOTROPIC    ) {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,"Anisotropic size(s)");
    }

    double h0 = ((IsoMeshSize*)pMS0)->size();
    double h1 = ((IsoMeshSize*)pMS1)->size();
    double h2 = ((IsoMeshSize*)pMS2)->size();
    h = (1.-u-v) * h0 + u * h1 + v * h2;
  }

  // -------------------------------------------------------------------
  void IsoMeshSize::interpolate(const pMSize pMS0, const pMSize pMS1, 
                                const pMSize pMS2, const pMSize pMS3,
                                double u, double v, double w)
  {
    if ( pMS0->getType() != ISOTROPIC || 
         pMS1->getType() != ISOTROPIC || 
         pMS2->getType() != ISOTROPIC || 
         pMS3->getType() != ISOTROPIC    ) {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,"Anisotropic size(s)");
    }

    double h0 = ((IsoMeshSize*)pMS0)->size();
    double h1 = ((IsoMeshSize*)pMS1)->size();
    double h2 = ((IsoMeshSize*)pMS2)->size();
    double h3 = ((IsoMeshSize*)pMS3)->size();
    h = (1.-u-v-w) * h0 + u * h1 + v * h2 + w * h3;
  }

  // -------------------------------------------------------------------
  double IsoMeshSize::direction(int i, double dir[3]) const 
  {
    for (int c=0; c<3; c++) dir[c] = 0.;
    dir[i] = 1.;
    return h;
  }

  // -------------------------------------------------------------------
  double IsoMeshSize::size(int i) const {
    return h;
  }

  // -------------------------------------------------------------------
  void IsoMeshSize::sizes(double _h[3]) const {
    _h[0] = _h[1] = _h[2] = h;
  }

  // -------------------------------------------------------------------
  double IsoMeshSize::getMeanLength() const 
  {
    return h;
  }

  // -------------------------------------------------------------------
  double IsoMeshSize::getMinLength() const 
  {
    return h;
  }

  // -------------------------------------------------------------------
  double IsoMeshSize::getMaxLength() const 
  {
    return h;
  }

  // -------------------------------------------------------------------
  void IsoMeshSize::setSize(double _h)
  {
    h = _h;
  }

  // -------------------------------------------------------------------
  void IsoMeshSize::scale(double factor)
  {
    h *= factor;
  }

  // -------------------------------------------------------------------
  double IsoMeshSize::normSq(const double vec[3]) const
  {
    return ( dotProd(vec,vec) / (h*h) );
  }

  // -------------------------------------------------------------------
  double IsoMeshSize::lengthSqInDir(const double vec[3]) const
  {
    return h*h;
  }

  // -------------------------------------------------------------------
  void IsoMeshSize::print(string name) const
  {
    std::cout<<"Printing IsoMeshSize \'"<<name<<"\' ("<<this<<")\n";
    std::cout<<"  h:\t"<<h<<"\n";
  }

  // -------------------------------------------------------------------

}
