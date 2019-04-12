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

#include "MeshDataBaseMiniMesh.h"

namespace MAd {

  size_t MDB_MiniMesh::size_buf () const
  {
    return sizeof(int) * ( 2 + nbInfos ) + sizeof(double) * 3 * nbVertex + sizeof(int) * nbVertex;
  }

  MDB_MiniMesh::~MDB_MiniMesh()
  {
    if (nbVertex)
      {
        if (coords) delete [] coords;
        delete [] infos;
        delete [] ids;
      }
  }

  void MDB_MiniMesh::_load ( FILE *f )
  {
    if (nbVertex != 0)
      {
        if (coords) { delete [] coords; coords=NULL; }
        delete [] infos;
        delete [] ids;
      }
    fread ( &nbVertex , sizeof(int)   , 1, f );
    fread ( &nbInfos  , sizeof(int)   , 1, f );
    coords = new double [3 * nbVertex];
    infos  = new int    [nbInfos     ];
    ids    = new int    [nbVertex    ];
    fread ( coords    , sizeof(double), nbVertex * 3, f );
    fread ( infos     , sizeof(int)   , nbInfos, f );
    fread ( ids       , sizeof(int)   , nbVertex, f );
  }

  void MDB_MiniMesh::_flush ( FILE *f ) const
  {
    if (nbVertex == 0) throw;
    fwrite ( &nbVertex , sizeof(int)   , 1, f );
    fwrite ( &nbInfos  , sizeof(int)   , 1, f );
    fwrite ( coords    , sizeof(double), nbVertex * 3, f );
    fwrite ( infos     , sizeof(int)   , nbInfos, f );
    fwrite ( ids       , sizeof(int)   , nbVertex, f );
  }

}
