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
// Authors: Koen Hillewaert
// -------------------------------------------------------------------

#ifndef MESHDATABASEPARALLELCHECK_H
#define MESHDATABASEPARALLELCHECK_H

#include "MeshDataBaseComm.h"

#include <iostream>
#include <map>
#include <set>


namespace MAd {

  /*! \brief Class for checking orientation of corresponding entities \ingroup parallel */

  class MDB_CommCheck : public MDB_DataExchanger {

  public:

    MDB_CommCheck();

    void printStatus(std::ostream&) const;
  
    virtual void* sendData(pEntity,int,int&);
    virtual void  receiveData(pEntity,int,void*);

  private:

    int numRecvs;
    int numSends;
    int numErrors;

    std::map<pEntity,std::set<std::pair<int,pEntity> > > communications;
    
  
  };

}

#endif
