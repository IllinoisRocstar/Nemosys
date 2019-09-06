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

#include "MeshSizeBase.h"
#include "IsoMeshSize.h"
#include "AnisoMeshSize.h"

namespace MAd {

  // -------------------------------------------------------------------
  pMSize MS_copy(const pMSize pMS)
  {
    MeshSizeType type = pMS->getType();
    switch (type) {
    case ISOTROPIC:
      return new IsoMeshSize((IsoMeshSize&)(*pMS));
      break;
    case ANISOTROPIC:
      return new AnisoMeshSize((AnisoMeshSize&)(*pMS));
      break;
    default:
      throw;
    }
  }

  // -------------------------------------------------------------------
  pMSize MS_intersect(const pMSize pMS0, const pMSize pMS1)
  {
    pMSize newS;
    if ( pMS0->getType() == ISOTROPIC && 
         pMS1->getType() == ISOTROPIC    ) {
      newS = new IsoMeshSize();
    }
    else {
      newS = new AnisoMeshSize();
    }
    newS->intersect(pMS0,pMS1);

    return newS;
  }

  // -------------------------------------------------------------------
  pMSize MS_interpolate(const pMSize pMS0, const pMSize pMS1, double param)
  {
    pMSize newS;
    if ( pMS0->getType() == ISOTROPIC && 
         pMS1->getType() == ISOTROPIC    ) {
      newS = new IsoMeshSize();
    }
    else {
      newS = new AnisoMeshSize();
    }
    newS->interpolate(pMS0,pMS1,param);
    return newS;
  }

  // -------------------------------------------------------------------
  pMSize MS_interpolate(const pMSize pMS0, const pMSize pMS1, 
                        const pMSize pMS2, double u, double v)
  {
    pMSize newS;
    if ( pMS0->getType() == ISOTROPIC &&
         pMS1->getType() == ISOTROPIC &&
         pMS2->getType() == ISOTROPIC    ) {
      newS = new IsoMeshSize();
    }
    else {
      newS = new AnisoMeshSize();
    }
    newS->interpolate(pMS0,pMS1,pMS2,u,v);
    return newS;
  }

  // -------------------------------------------------------------------
  pMSize MS_interpolate(const pMSize pMS0, const pMSize pMS1, 
                        const pMSize pMS2, const pMSize pMS3,
                        double u, double v, double w)
  {
    pMSize newS;
    if ( pMS0->getType() == ISOTROPIC &&
         pMS1->getType() == ISOTROPIC &&
         pMS2->getType() == ISOTROPIC &&
         pMS3->getType() == ISOTROPIC    ) {
      newS = new IsoMeshSize();
    }
    else {
      newS = new AnisoMeshSize();
    }
    newS->interpolate(pMS0,pMS1,pMS2,pMS3,u,v,w);
    return newS;
  }

  // -------------------------------------------------------------------

}
