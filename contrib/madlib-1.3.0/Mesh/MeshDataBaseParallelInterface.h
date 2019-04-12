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
// Authors: J.-F. Remacle, C. Dobrzynski, K. Hillewaert, G. Compere
// -------------------------------------------------------------------

#ifndef H_MESHDATABASEPARALLELINTEFACE
#define H_MESHDATABASEPARALLELINTEFACE

#include "MeshDataBaseInterface.h"
#include "MeshDataBaseComm.h"
#include "MeshDataBaseCommPeriodic.h"

#include <vector>

namespace MAd {

  class MDB_DataExchanger;
  class MDB_DataExchangerPeriodic;

  // -------------------------------------------------------------------

#ifdef PARALLEL
  void M_writeParallel(pMesh, const char*, int version=2);
#endif 

  // -------------------------------------------------------------------

  /*! \brief Return true if on a parallel or periodic interface \ingroup internal */ 
  bool EN_isInterface(pEntity pv);

  /*! \brief Return true if on a parallel or periodic interface \ingroup parallel */ 
  bool V_isInterface(pVertex pv);

  /*! \brief Return true if pv is a (periodic) copy of vertex with tag id \ingroup parallel */
  bool V_corresponds(pVertex pv,int id);

  /*! \brief Return size of the list and list is an array containing 
    proc numbers where pv exist except calling proc \ingroup parallel */ 
  int  V_listInterface(pVertex pv, std::vector<int>* list);

  /*! \brief  Return true if edge is //possibly// on a parallel or periodic interface \ingroup internal */ 
  bool E_isPotentialInterface(pMesh, pEdge);

  /*! \brief  Return true if on a parallel or periodic interface \ingroup parallel */ 
  bool E_isInterface(pEdge);

  /*! \brief Return true if the edge corresponds either directly or periodically \ingroup parallel */
  bool E_corresponds(pEdge,int,int);

  /*! \brief  Return true if on a parallel interface, distProc is 
    (one of) the proc sharing the edge and distVt is a list 
    of the pointers to the vertices on the distProc mesh \ingroup parallel */ 
  bool E_isInterface(pMesh, pEdge, int * distProc,
                     std::vector<pVertex>* distVt);
  
  /*! \brief  Return true if face is potentially on a parallel interface \ingroup internal */ 
  bool F_isPotentialInterface(pMesh, pFace);

  /*! \brief  Return true if face is on a parallel interface \ingroup parallel */ 
  bool F_isInterface(pFace);

  /*! \brief  Return true if on a parallel interface, distProc is the 
    proc sharing the face and distVt is a list of the pointers 
    to the vertices on the distant mesh \ingroup parallel */ 
  bool F_isInterface(pMesh, pFace, int * distProc,
                     std::vector<pVertex>* distVt);
  /*! \brief  Return true if on a parallel interface (always false) \ingroup parallel */ 
  bool R_isInterface(pRegion pr);

  // -------------------------------------------------------------------

  /*! \brief Fill the attached pointer containing the distant proc
    numbers and the distant node pointers \ingroup parallel */ 
  void V_createInfoInterface(pMesh mesh, pMeshDataId tagVertex);

  /*! \brief Fill the attached pointer containing the distant proc
    numbers and the distant edge pointers \ingroup parallel */ 
  void E_createInfoInterface(pMesh mesh, pMeshDataId tagEdge);

  /*! \brief Fill the attached pointer containing the distant proc
    numbers and the distant face pointers \ingroup parallel */ 
  void F_createInfoInterface(pMesh mesh, pMeshDataId tagFace);

//   void UpdateIDGlobal(pMesh mesh, int IdGlobal);

  // -------------------------------------------------------------------
#ifdef PARALLEL

  // interface elt migration
  void Balance(pMesh mesh,MDB_DataExchanger &de);
  void Balance2(pMesh mesh,MDB_DataExchanger &de);
  void BalanceRandom(pMesh mesh,MDB_DataExchanger &de);
  int BalanceManifold(pMesh mesh,MDB_DataExchanger &de);

#ifdef _HAVE_PARMETIS_
  // load balancing with metis
  void BalanceMetis(pMesh mesh,MDB_DataExchanger &de);
  void BalanceMetis2(pMesh mesh,MDB_DataExchanger &de);
#endif

  // migration of elt tagged tagElt
  void loadBalancing(pMesh mesh,pMeshDataId tagElt, MDB_DataExchanger &de);  
  void loadBalancing2(pMesh mesh,pMeshDataId tagElt, MDB_DataExchanger &de);
#endif

  //Move periodic interfaces
  void BalancePeriodic(pMesh mesh,int dim,MDB_DataExchanger &de,
                       MDB_DataExchangerPeriodic &deperiodic,std::vector<std::vector<int> >& transfo);

// -------------------------------------------------------------------

//migration of periodic elt
void PeriodicInterfaceMigration(MAd::pMesh mesh,MAd::pMeshDataId tagElt,MAd::pMeshDataId tagMove,
                                MAd::pMeshDataId tagTransfo, MDB_DataExchanger &de,
                                MDB_DataExchangerPeriodic &deperiodic); 


void GroupPeriodicTetra(MAd::pMesh mesh, MDB_DataExchanger &de,
                        MDB_DataExchangerPeriodic &deperiodic);


// -------------------------------------------------------------------

}

#endif
