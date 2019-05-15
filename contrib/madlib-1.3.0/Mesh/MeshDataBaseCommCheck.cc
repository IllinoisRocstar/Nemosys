// -------------------------------------------------------------------
// MAdLib - Copyright (C) 2008-2009 Universite catholique de Louvain
//
// See the Copyright.txt and License.txt files for license information. 
// You should have received a copy of these files along with MAdLib. 
// If not, see <http://www.madlib.be/license/>
//
// Please report all bugs and problems to <contrib@madlib.be>
//
// Authors: Koen Hillewaert
// -------------------------------------------------------------------

#include "MeshDataBaseCommCheck.h"
#include "MeshDataBaseInterface.h"
#include "MeshDataBaseParallelInterface.h"
#include "MAdMessage.h"
#include "MeshDataBase.h"

#include <sstream>
using std::ostringstream;
using std::ostream;
using std::cout;
using std::endl;
#include <iostream>
#include <stdlib.h>

#ifdef PARALLEL
#include "mpi.h"
#endif

namespace MAd {

  // -----------------------------------------------------------------------------

  MDB_CommCheck::MDB_CommCheck():
    MDB_DataExchanger(1354),numRecvs(0),numSends(0),numErrors(0) {}

  // -----------------------------------------------------------------------------

  void MDB_CommCheck::printStatus(ostream& out) const {

    out << "Number of sends " << numSends
        << ", number of receives " << numRecvs
        << ", number of errors " << numErrors << std::endl;
  }


  // -----------------------------------------------------------------------------

  void* MDB_CommCheck::sendData(pEntity pe,int iProc,int& size) {  
  
    int type = EN_type(pe);

    void* buf = NULL;
    size = sizeof(pEntity);
  
    switch (type) {
    
    case 0:
      {
        size += 2*sizeof(int);
        buf = malloc(size);

        int* ibuf = reinterpret_cast<int*>(buf);

        *(ibuf++) = 1;
        *(ibuf++) = EN_id(pe);
        
        pEntity* ebuf = reinterpret_cast<pEntity*>(ibuf);
        *ebuf = pe;
        
        return buf;
      
        break;
      }
    case 1:
      {
        pEdge pl = (pEdge) pe;
        size += 3 * sizeof(int);
        buf = malloc(size);

        int *ibuf = reinterpret_cast<int*>(buf);

        *(ibuf++) = 2;
        *(ibuf++) = EN_id((pEntity) E_vertex(pl,0));
        *(ibuf++) = EN_id((pEntity) E_vertex(pl,1));      
        
        pEntity* ebuf = reinterpret_cast<pEntity*>(ibuf);
        *ebuf = pl;
        
        break;
      }
    case 2:
      {
        pFace pf = (pFace) pe;
        size += (1 + F_numVertices(pf)) * sizeof(int);
        buf  = malloc(size);
        
        int *ibuf = reinterpret_cast<int*>(buf);
        
        *(ibuf++) = F_numVertices(pf);
        for (int i=0;i<F_numVertices(pf);i++) *(ibuf++) = EN_id((pEntity) F_vertex(pf,i));
        
        pEntity* ebuf = reinterpret_cast<pEntity*>(ibuf);
        *ebuf = pf;
        
        break;
      }
    }
    return buf;
  }

  // -----------------------------------------------------------------------------

  void MDB_CommCheck::receiveData(pEntity pe,int from,void* buf) {
  
    int myrank = 0;
#ifdef PARALLEL
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
#endif

    int type = EN_type(pe);

    int* ibuf = reinterpret_cast<int*>(buf);
    int  np   = *(ibuf++);
    

    switch (type) {

    case 0: 
      {
        int id = *(ibuf++);
        
        if (!V_corresponds((pVertex) pe,id)) {
          
          MAdMsgSgl::instance().error(__LINE__,__FILE__,
                   "Non-corresponding communication vertex from proc %d : %d vs %d",
                   from,EN_id(pe),id);
        }
        
        pEntity* pbuf = reinterpret_cast<pEntity*>(ibuf);
        pEntity  corr = *pbuf;
        
        std::pair<int,pEntity> comm = std::make_pair(from,corr);
        
        if (communications[pe].find(comm) != communications[pe].end()) {
          MAdMsgSgl::instance().error(__LINE__,__FILE__,
                   "Multiple identical communications for vertex %d on %d from %d",
                   EN_id(pe),myrank,from);
          
        }
        communications[pe].insert(comm);
        break;
      }
    case 1:
      {
        pEdge pl = (pEdge) pe;
        bool check = true;
      
        int idSend[2];
      
        for (int i=0;i<2;i++){
          idSend[i] = *(ibuf++);
          check = check && V_corresponds(E_vertex(pl,i),idSend[i]);
        }
        
        if (!check) {        
          
          MAdMsgSgl::instance().error(__LINE__,__FILE__,
                   "Non-corresponding edge from proc %d on proc %d : (%d,%d) vs (%d,%d)",
                   from,myrank,
                   idSend[0],
                   idSend[1],
                   EN_id(E_vertex(pl,0)),
                   EN_id(E_vertex(pl,1)));
        }

        
        pEntity* pbuf = reinterpret_cast<pEntity*>(ibuf);
        pEntity  corr = *pbuf;
        
        std::pair<int,pEntity> comm = std::make_pair(from,corr);

        if (communications[pe].find(comm) != communications[pe].end()) {
          MAdMsgSgl::instance().error(__LINE__,__FILE__,
                   "Multiple identical communications for edge %d-%d on %d from %d",
                   EN_id((pEntity) E_vertex(pl,0)),
                   EN_id((pEntity) E_vertex(pl,1))
                   ,myrank,from);
        }
        communications[pe].insert(comm);
        break;
      }
    case 2:
      {
        pFace pf = (pFace) pe;
        bool check = (np == F_numVertices(pf));
        
        int idSend[4];
      
        for (int i=0;i<np ;i++){
          idSend[i] = *(ibuf++);
          check = check && V_corresponds(F_vertex(pf,i),idSend[i]);
        }
        
        if (!check) {
          
          ostringstream recv;
          for (int i=0;i<np;i++) recv << " " << idSend[i];
          ostringstream send;
          for (int i=0;i<F_numVertices(pf);i++) send << " " << EN_id((pEntity) F_vertex(pf,i));
        
          MAdMsgSgl::instance().error(__LINE__,__FILE__,
                   "Non-corresponding face from proc %d - (%s) vs (%s)",
                   recv.str().c_str(),send.str().c_str());
        }

        
        pEntity* pbuf = reinterpret_cast<pEntity*>(ibuf);
        pEntity  corr = *pbuf;
        
        std::pair<int,pEntity> comm = std::make_pair(from,corr);

        if (communications[pe].find(comm) != communications[pe].end()) {
          MAdMsgSgl::instance().error(__LINE__,__FILE__,
                   "Multiple identical communications for face %d-%d-%d on %d from %d",
                   EN_id((pEntity) F_vertex(pf,0)),
                   EN_id((pEntity) F_vertex(pf,1)),
                   EN_id((pEntity) F_vertex(pf,2)),
                   myrank,from);
        }
        communications[pe].insert(comm);
        break;
      } 
    }
  }
}
