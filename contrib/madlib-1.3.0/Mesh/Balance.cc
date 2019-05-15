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

#include "MeshDataBase.h"
#include "MeshDataBaseInterface.h"
#include "MeshDataBaseComm.h"
#include "MeshDataBaseParallelInterface.h"
#include "Mark.h"
#include "MAdMessage.h"

#ifdef PARALLEL
#include "mpi.h"
#include "autopack.h"
#ifdef _HAVE_PARMETIS_
#include "metisAdaptiveRepart.h"
#endif
#endif

#include <stdio.h>

using namespace MAd;

namespace MAd {

#ifdef PARALLEL
  // -------------------------------------------------------------------
  void Balance(pMesh mesh,MDB_DataExchanger &de) {

    int dim = (mesh->tets.empty()) ? 2 : 3;
    pMeshDataId tagElt= MD_newMeshDataId("EltDestination"); //dest = (int - 1)
  
    // --- Mark the elements to be moved and their destination ---
    if(dim==2) MarkTriangles (mesh, tagElt);
    else       MarkTets      (mesh, tagElt);
    //  if(dim==2) MarkTrianglesSmooth (mesh, tagElt);
    //  else       MarkTetsSmooth      (mesh, tagElt);

    // --- Move marked elements ---
    loadBalancing(mesh, tagElt, de);
  
    // --- Check that the mesh is not empty (debug) ---
    if ( dim==2 ) assert( !( mesh->triangles.empty() ) );
    if ( dim==3 ) assert( !( mesh->tets.empty() ) );

    // ----------------------------------------------
    // ------ Tagging inter-partition nodes
    // ----------------------------------------------

    pMeshDataId tag = MD_lookupMeshDataId("RemotePoint");
  
    V_createInfoInterface(mesh,tag);
    E_createInfoInterface(mesh,tag);
    F_createInfoInterface(mesh,tag);

    return;
  }
  // -------------------------------------------------------------------
  void Balance2(pMesh mesh,MDB_DataExchanger &de) {

    int dim = (mesh->tets.empty()) ? 2 : 3;
    pMeshDataId tagElt= MD_newMeshDataId("EltDestination"); //dest = (int - 1)

    // --- Mark the elements to be moved and their destination ---
    if(dim==2) MarkTriangles (mesh, tagElt);
    else       MarkTets      (mesh, tagElt);

    // --- Move marked elements ---
    loadBalancing2(mesh, tagElt, de);

    // --- Check that the mesh is not empty (debug) ---
    if ( dim==2 ) assert( !( mesh->triangles.empty() ) );
    if ( dim==3 ) assert( !( mesh->tets.empty() ) );
    // ----------------------------------------------
    // ------ Tagging inter-partition nodes
    // ----------------------------------------------

    pMeshDataId tag = MD_lookupMeshDataId("RemotePoint");
 
    V_createInfoInterface(mesh,tag);
    E_createInfoInterface(mesh,tag);
    F_createInfoInterface(mesh,tag);

    return;
  }

  // -------------------------------------------------------------------
  int BalanceManifold(pMesh mesh,MDB_DataExchanger &de) {
  
    int dim = (mesh->tets.empty()) ? 2 : 3;  
    pMeshDataId tagElt= MD_newMeshDataId("EltDestination"); //dest = (int - 1)
  
    int nmanifold = MarkTetsManifold(mesh,tagElt);
    
    loadBalancing(mesh,tagElt,de);

    if ( dim==3 ) assert( !(mesh->tets.empty()) );
 
    return(nmanifold);
  }

  // -------------------------------------------------------------------
  void BalanceRandom(pMesh mesh, MDB_DataExchanger &de) {

    int dim = (mesh->tets.empty()) ? 2 : 3;
    pMeshDataId tagElt = MD_newMeshDataId("EltDestination"); //dest = (int - 1)

    if(dim==2) MarkTrianglesRandom(mesh,tagElt);
    else       MarkTetsRandom(mesh,tagElt);

    loadBalancing(mesh,tagElt,de);

    if ( dim==3 ) assert( !(mesh->tets.empty()) );
  }

  // -------------------------------------------------------------------
#ifdef _HAVE_PARMETIS_
  void BalanceMetis(pMesh mesh,MDB_DataExchanger &de) {
  
    int dim = (mesh->tets.empty()) ? 2 : 3;
  
    //std::cout<<"metisAdaptiveRepart"<<mesh->nbPoints<<" "<<mesh->nbTets<<std::endl;
    pMeshDataId tagElt    = MD_newMeshDataId("EltDestination"); //dest = (int - 1)
    std::cout <<"USING PARMETIS"<<std::endl;
    metisAdaptiveRepart(mesh,tagElt);

    loadBalancing(mesh,tagElt,de);

    //std::cout<<"end metisAdaptiveRepart"<<mesh->nbPoints<<" "<<mesh->nbTets<<std::endl;
  
    if ( dim==3 ) assert( !(mesh->tets.empty()) );
  
    return;
  }
#endif

  // -------------------------------------------------------------------
#ifdef _HAVE_PARMETIS_
  void BalanceMetis2(pMesh mesh, MDB_DataExchanger &de) {
  
    int dim = M_dim(mesh);
  
    pMeshDataId tagElt = MD_newMeshDataId("EltDestination"); //dest = (int - 1)
    metisAdaptiveRepart(mesh,tagElt);
    loadBalancing2(mesh,tagElt,de);

    if ( M_dim(mesh) != dim ) {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                  "The dimension of the mesh has been reduced from %d to %d during a balancing operation",
                                  dim, M_dim(mesh));
    }

    pMeshDataId tag = MD_lookupMeshDataId("RemotePoint");
    V_createInfoInterface(mesh,tag);
    E_createInfoInterface(mesh,tag);
    F_createInfoInterface(mesh,tag);
  }
#endif

  // -------------------------------------------------------------------
#endif

  // -------------------------------------------------------------------
  void BalancePeriodic(pMesh mesh,int dim,MDB_DataExchanger &de,
                       MDB_DataExchangerPeriodic &deperiodic,std::vector<std::vector<int> >& transfo) {
    pMeshDataId tagMove = MD_newMeshDataId("TagMovePeriodic");
    pMeshDataId tagElt  = MD_newMeshDataId("EltDestination"); //dest = (int - 1)
    /*for 3-periodic cases*/
    pMeshDataId tagTransfo  = MD_newMeshDataId("Transformation"); 
  
    if(dim==2) MarkPeriodicTriangles(mesh,transfo,tagElt,tagMove,tagTransfo);
    else       MarkPeriodicTets(mesh,transfo,tagElt,tagMove,tagTransfo);  
  
    PeriodicInterfaceMigration(mesh,tagElt,tagMove,tagTransfo,de,deperiodic);
  
    MD_deleteMeshDataId(tagMove);
    MD_deleteMeshDataId(tagElt);
    MD_deleteMeshDataId(tagTransfo);
  
    return;
  }

  // -------------------------------------------------------------------

}
