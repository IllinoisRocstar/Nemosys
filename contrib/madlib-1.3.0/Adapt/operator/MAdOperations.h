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
// Authors: Gaetan Compere, Jean-Francois Remacle
// -------------------------------------------------------------------

#ifndef _H_MADOPERATIONS
#define _H_MADOPERATIONS

namespace MAd {

  // -------------------------------------------------------------------
  enum operationType {
    MAd_UNKNOWNOPERATION = 0,
    MAd_VERTEXMOVE       = 1,
    MAd_ESPLIT           = 2,
    MAd_ECOLLAPSE        = 3,
    MAd_FCOLLAPSE        = 4,
    MAd_DESPLTCLPS       = 5,
    MAd_ESWAP            = 6,
    MAd_FSWAP            = 7,
    MAd_RREMOVE          = 8
  };

  // -------------------------------------------------------------------

}

#endif
