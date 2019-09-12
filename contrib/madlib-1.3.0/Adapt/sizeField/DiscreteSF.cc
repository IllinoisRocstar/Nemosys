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

#include "DiscreteSF.h"
#include "IsoMeshSize.h"
#include "AnisoMeshSize.h"

using namespace MAd;

// -------------------------------------------------------------------

namespace MAd {

  // -------------------------------------------------------------------
  DiscreteSF::DiscreteSF(pMesh m, std::string name): SizeFieldBase(name)
  {
    mesh = m;
    dim = M_dim(mesh);

    pMSizeFieldId = MD_newMeshDataId("");
  }
 
  // -------------------------------------------------------------------
  DiscreteSF::~DiscreteSF()
  {
    MD_deleteMeshDataId(pMSizeFieldId);
  }

  // -------------------------------------------------------------------
  void DiscreteSF::deleteSize(pEntity entity)
  {
    void * temp;
    if( EN_getDataPtr(entity,pMSizeFieldId,&temp) ) {
      EN_deleteData(entity,pMSizeFieldId);
      delete (pMSize)temp; 
    }
  }

  // -------------------------------------------------------------------
  void DiscreteSF::setSize(pEntity ent, 
                           double dirs[3][3],
                           double h[3])
  {
    pMSize pS = new AnisoMeshSize(dirs,h);
    setSize(ent,pS);
  }


  // -------------------------------------------------------------------
  void DiscreteSF::setSize(pEntity ent, double h)
  {
    pMSize pS = new IsoMeshSize(h);
    setSize(ent,pS);
  }

  // -------------------------------------------------------------------
  void DiscreteSF::setSize(pEntity pEnt, pMSize pS)
  {
    void * temp;
    if( EN_getDataPtr(pEnt,pMSizeFieldId,&temp) ) {
      delete (pMSize)temp;
      EN_modifyDataPtr(pEnt,pMSizeFieldId,pS);
    }
    else {
      EN_attachDataPtr(pEnt,pMSizeFieldId,pS);
    }
  }

  // -------------------------------------------------------------------
  void DiscreteSF::setSize(int id, double h)
  {
    pVertex pv = M_vertex(mesh,id);
    pMSize pS = new IsoMeshSize(h);
    setSize((pEntity)pv,pS);
  }

  // -------------------------------------------------------------------

}
