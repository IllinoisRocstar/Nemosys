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

#ifndef _MeshDataBaseCommPerio_H
#define _MeshDataBaseCommPerio_H

/*
  MDB_DataExchangerPeriodic is a class that allow user to define function for
  periodic boundaries.
*/

namespace MAd {

  class MDB_DataExchangerPeriodic
  {
  public :
    virtual int nbRefPeriodic() const = 0;
    virtual void   fperiodic (const int inv,
                              const double X, 
                              const double Y,
                              const double Z,
                              int numref,      // number of the transformation
                              double *Xnew,
                              double *Ynew,
                              double *Znew) = 0;
    virtual ~MDB_DataExchangerPeriodic() {};			    
  };

}

#endif
