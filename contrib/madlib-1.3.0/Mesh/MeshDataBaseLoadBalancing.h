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
// Authors: Cecile Dobrzynski, Jean-Francois Remacle, Gaetan Compere
// -------------------------------------------------------------------

#ifndef H_MESHDATABASELOADBALANCE
#define H_MESHDATABASELOADBALANCE

#ifdef PARALLEL

namespace MAd {

#ifdef DEBUG
  void checkRemotePointer(pMesh mesh, pMeshDataId tagData );
  void checkRemotePointer2(pMesh mesh, pMeshDataId tagData );
  void checkRemotePointerChange(pMesh mesh, pMeshDataId tagData,pMeshDataId tagNew, pMeshDataId tagChange );
  void checkNew(pMesh mesh, pMeshDataId tagData,pMeshDataId tagChange);
#endif
  void MarkedEltVertex(pMesh mesh,pMeshDataId tagElt);
  void CommOldInterfaces(pMesh mesh);

  void UpdateInterfaces(pMesh mesh, MDB_DataExchanger &de);
  void SendVertex(pMesh mesh, MDB_DataExchanger &de);
  void SendElt(pMesh mesh,pMeshDataId tagElt, MDB_DataExchanger &de);
  void DeleteEntitiesAndData(pMesh mesh,pMeshDataId tagElt, MDB_DataExchanger &de );
  void MarkEltSubEntities(pMesh mesh, pMeshDataId tagElt);
  void MigrateEntitiesAndData(pMesh mesh, pMeshDataId tagDest, MDB_DataExchanger &de);

}

#endif
#endif
