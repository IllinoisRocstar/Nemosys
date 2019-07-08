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

#include "madlib_export.h"

#include "MeshDataBaseInterface.h"
#include "MeshDataBaseComm.h"
#include "MeshDataBaseCommPeriodic.h"

#include <vector>

namespace MAd {

  class MDB_DataExchanger;
  class MDB_DataExchangerPeriodic;

  // -------------------------------------------------------------------

#ifdef PARALLEL
  void MADLIB_EXPORT M_writeParallel(pMesh, const char*, int version=2);
#endif 

  // -------------------------------------------------------------------

  /*! \brief Return true if on a parallel or periodic interface \ingroup internal */ 
  bool MADLIB_EXPORT EN_isInterface(pEntity pv);

  /*! \brief Return true if on a parallel or periodic interface \ingroup parallel */ 
  bool MADLIB_EXPORT V_isInterface(pVertex pv);

  /*! \brief Return true if pv is a (periodic) copy of vertex with tag id \ingroup parallel */
  bool MADLIB_EXPORT V_corresponds(pVertex pv,int id);

  /*! \brief Return size of the list and list is an array containing 
    proc numbers where pv exist except calling proc \ingroup parallel */ 
  int  MADLIB_EXPORT V_listInterface(pVertex pv, std::vector<int>* list);

  /*! \brief  Return true if edge is //possibly// on a parallel or periodic interface \ingroup internal */ 
  bool MADLIB_EXPORT E_isPotentialInterface(pMesh, pEdge);

  /*! \brief  Return true if on a parallel or periodic interface \ingroup parallel */ 
  bool MADLIB_EXPORT E_isInterface(pEdge);

  /*! \brief Return true if the edge corresponds either directly or periodically \ingroup parallel */
  bool MADLIB_EXPORT E_corresponds(pEdge,int,int);

  /*! \brief  Return true if on a parallel interface, distProc is 
    (one of) the proc sharing the edge and distVt is a list 
    of the pointers to the vertices on the distProc mesh \ingroup parallel */ 
  bool MADLIB_EXPORT E_isInterface(pMesh, pEdge, int * distProc,
                                   std::vector<pVertex>* distVt);
  
  /*! \brief  Return true if face is potentially on a parallel interface \ingroup internal */ 
  bool MADLIB_EXPORT F_isPotentialInterface(pMesh, pFace);

  /*! \brief  Return true if face is on a parallel interface \ingroup parallel */ 
  bool MADLIB_EXPORT F_isInterface(pFace);

  /*! \brief  Return true if on a parallel interface, distProc is the 
    proc sharing the face and distVt is a list of the pointers 
    to the vertices on the distant mesh \ingroup parallel */ 
  bool MADLIB_EXPORT F_isInterface(pMesh, pFace, int * distProc,
                                   std::vector<pVertex>* distVt);
  /*! \brief  Return true if on a parallel interface (always false) \ingroup parallel */ 
  bool MADLIB_EXPORT R_isInterface(pRegion pr);

  // -------------------------------------------------------------------

  /*! \brief Fill the attached pointer containing the distant proc
    numbers and the distant node pointers \ingroup parallel */ 
  void MADLIB_EXPORT V_createInfoInterface(pMesh mesh, pMeshDataId tagVertex);

  /*! \brief Fill the attached pointer containing the distant proc
    numbers and the distant edge pointers \ingroup parallel */ 
  void MADLIB_EXPORT E_createInfoInterface(pMesh mesh, pMeshDataId tagEdge);

  /*! \brief Fill the attached pointer containing the distant proc
    numbers and the distant face pointers \ingroup parallel */ 
  void MADLIB_EXPORT F_createInfoInterface(pMesh mesh, pMeshDataId tagFace);

//   void MADLIB_EXPORT UpdateIDGlobal(pMesh mesh, int IdGlobal);

  // -------------------------------------------------------------------
#ifdef PARALLEL

  // interface elt migration
  void MADLIB_EXPORT Balance(pMesh mesh,MDB_DataExchanger &de);
  void MADLIB_EXPORT Balance2(pMesh mesh,MDB_DataExchanger &de);
  void MADLIB_EXPORT BalanceRandom(pMesh mesh,MDB_DataExchanger &de);
  int  MADLIB_EXPORT BalanceManifold(pMesh mesh,MDB_DataExchanger &de);

#ifdef _HAVE_PARMETIS_
  // load balancing with metis
  void MADLIB_EXPORT BalanceMetis(pMesh mesh,MDB_DataExchanger &de);
  void MADLIB_EXPORT BalanceMetis2(pMesh mesh,MDB_DataExchanger &de);
#endif

  // migration of elt tagged tagElt
  void MADLIB_EXPORT loadBalancing(pMesh mesh,pMeshDataId tagElt, MDB_DataExchanger &de);
  void MADLIB_EXPORT loadBalancing2(pMesh mesh,pMeshDataId tagElt, MDB_DataExchanger &de);
#endif

  //Move periodic interfaces
  void MADLIB_EXPORT BalancePeriodic(pMesh mesh,int dim,MDB_DataExchanger &de,
                                     MDB_DataExchangerPeriodic &deperiodic,
                                     std::vector<std::vector<int> >& transfo);

// -------------------------------------------------------------------

//migration of periodic elt
void MADLIB_EXPORT PeriodicInterfaceMigration(MAd::pMesh mesh,
                                              MAd::pMeshDataId tagElt,
                                              MAd::pMeshDataId tagMove,
                                              MAd::pMeshDataId tagTransfo,
                                              MDB_DataExchanger &de,
                                              MDB_DataExchangerPeriodic &deperiodic); 


void MADLIB_EXPORT GroupPeriodicTetra(MAd::pMesh mesh, MDB_DataExchanger &de,
                                      MDB_DataExchangerPeriodic &deperiodic);


// -------------------------------------------------------------------

}

#endif
