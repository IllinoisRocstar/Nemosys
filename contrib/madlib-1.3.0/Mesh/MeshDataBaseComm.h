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

#ifndef _MeshDataBaseCOMM_H
#define _MeshDataBaseCOMM_H

#include "MeshDataBaseInterface.h"

/*
  MDB_DataExchanger is a class that allow user to exchange
  data's between partitions and/or periodic boundaries. User
  has to provide a tag for non synchronous communications. 
  Basically, choose a tag that is higher that 200 for being
  sure that it does not interact with internal messages.
  User has to allocate buffers using malloc.

  User will receive its datas.
*/

namespace MAd {

  class MDB_DataExchanger
  {
  private:
    int tagComm;
  public :
    // get a tag for the communication
    int tag() const { return tagComm; }
    // user allocates sends a message of _size size related to mesh entity pe to proc iProc
    virtual void * sendData (pEntity pe,       // in
                             int iProcDest ,   // in
                             int &_size ) = 0; // out
    // mesh entity pe recieves data *buf form proc iProc.
    // The user shall NOT delete the message !!
    virtual void receiveData (pEntity pe,       //in
                              int iProcSender , //in
                              void *buf ) = 0;  //in
    // In case the user has to delete a data when 'pe' is deleted
    virtual void deleteExternalData (pEntity pe) const {}
    MDB_DataExchanger(int _tag): tagComm(_tag) {}
    virtual ~MDB_DataExchanger() {}
  };
  void exchangeDataOnVertices (pMesh m, MDB_DataExchanger &de );
  void exchangeDataOnEdges    (pMesh m, MDB_DataExchanger &de );
  void exchangeDataOnFaces    (pMesh m, MDB_DataExchanger &de );
  void exchangeDataOnTriangles(pMesh m, MDB_DataExchanger &de );
  void exchangeDataOnQuads    (pMesh m, MDB_DataExchanger &de );

}

#endif
