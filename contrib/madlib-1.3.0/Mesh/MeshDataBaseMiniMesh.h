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
// Authors: Jean-Francois Remacle, Gaetan Compere
// -------------------------------------------------------------------

#ifndef _MESHDATABASE_MINIMESH_
#define _MESHDATABASE_MINIMESH_

#include <stdio.h>

namespace MAd {

  class MDB_MiniMesh
  {
  protected:
    int nbVertex;
    double * coords;
    int nbInfos;
    int *infos;
    int *ids;
    void _flush ( FILE *F ) const;
    void _load  ( FILE *F );
    size_t size_buf () const;
  public:
    virtual ~MDB_MiniMesh ();
    MDB_MiniMesh () : nbVertex (0) , coords(0), nbInfos (0), infos(0), ids(0) {}
  };

}

#endif
