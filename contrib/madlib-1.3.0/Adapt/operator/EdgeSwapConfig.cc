// -------------------------------------------------------------------
// MAdLib - Copyright (C) 2008-2009 Universite catholique de Louvain
//
// See the Copyright.txt and License.txt files for license information. 
// You should have received a copy of these files along with MAdLib. 
// If not, see <http://www.madlib.be/license/>
//
// Please report all bugs and problems to <contrib@madlib.be>
//
// Authors: Arnaud Francois, Gaetan Compere, Jean-Francois Remacle
// -------------------------------------------------------------------

#include "EdgeSwapConfig.h"
#include <iostream>

using namespace MAd;

// -------------------------------------------------------------------
// Edge swap template with 3 faces connected to the edge
const int EdgeSwap3::triangles[1][3] = { { 0,1,2 } };
const int EdgeSwap3::triangulations[1][1] = {{0}};

// Edge swap template with 4 faces connected to the edge
const int EdgeSwap4::triangles[4][3] = { {0,1,2}, {0,2,3}, {0,1,3}, {1,2,3} };
const int EdgeSwap4::triangulations[2][2] = { {0, 1}, {2, 3} };

// Edge swap template with 5 faces connected to the edge
const int EdgeSwap5::triangles[10][3]= { {0,1,2}, {0,2,3}, {0,3,4}, {0,1,4}, {1,3,4},
                                         {1,2,3}, {2,3,4}, {0,2,4}, {0,1,3}, {1,2,4} };
const int EdgeSwap5::triangulations[5][3] = { {0,1,2}, {3,4,5}, {0,6,7}, {2,5,8}, {3,6,9} } ;

// Edge swap template with 6 faces connected to the edge
const int EdgeSwap6::triangles[20][3]= { {0,1,2}, {0,2,3}, {0,3,4}, {0,4,5}, 
      {0,2,5}, {2,4,5}, {2,3,4}, {0,3,5}, {3,4,5}, {0,2,4}, {2,3,5}, {1,2,3},
      {0,1,3}, {0,1,5}, {1,4,5}, {1,3,4}, {0,1,4}, {1,3,5}, {1,2,4}, {1,2,5} };
const int EdgeSwap6::triangulations[14][4] = { {0,1,2,3}, {0,4,5,6}, {0,1,7,8},
      {0,3,6,9}, {0,4,8,10}, {2,3,11,12}, {11,13,14,15}, {7,8,11,12}, {3,11,15,16},
      {8,11,13,17}, {6,13,14,18}, {3,6,16,18}, {5,6,13,19}, {8,10,13,19} };

// Edge swap template with 7 faces connected to the edge
const int EdgeSwap7::triangles[35][3] = { {0,1,2}, {0,2,3}, {0,3,4}, {0,4,5}, 
      {0,5,6}, {0,3,6}, {3,5,6}, {3,4,5}, {0,4,6}, {4,5,6}, {0,3,5}, {3,4,6},
      {0,2,4}, {2,3,4}, {0,2,6}, {2,5,6}, {2,4,5}, {0,2,5}, {2,4,6}, {2,3,5},
      {2,3,6}, {0,1,3}, {1,2,3}, {0,1,4}, {1,3,4}, {0,1,6}, {1,5,6}, {1,4,5},
      {0,1,5}, {1,4,6}, {1,3,5}, {1,3,6}, {1,2,4}, {1,2,5}, {1,2,6} };
const int EdgeSwap7::triangulations[42][5] = { {0,1,2,3,4}, {0,1,5,6,7}, 
      {0,1,2,8,9}, {0,1,4,7,10}, {0,1,5,9,11}, {0,3,4,12,13}, {0,13,14,15,16},
      {0,8,9,12,13}, {0,4,13,16,17}, {0,9,13,14,18}, {0,7,14,15,19},{0,4,7,17,19},
      {0,6,7,14,20}, {0,9,11,14,20}, {2,3,4,21,22}, {5,6,7,21,22}, {2,8,9,21,22},
      {4,7,10,21,22}, {5,9,11,21,22}, {3,4,22,23,24}, {22,24,25,26,27}, 
      {8,9,22,23,24}, {4,22,24,27,28}, {9,22,24,25,29}, {7,22,25,26,30},
      {4,7,22,28,30}, {6,7,22,25,31}, {9,11,22,25,31}, {3,4,13,23,32},
      {13,25,26,27,32}, {8,9,13,23,32}, {4,13,27,28,32}, {9,13,25,29,32},
      {13,16,25,26,33}, {4,13,16,28,33}, {13,15,16,25,34}, {9,13,18,25,34},
      {7,19,25,26,33}, {4,7,19,28,33}, {7,15,19,25,34}, {6,7,20,25,34}, 
      {9,11,20,25,34} };

// -------------------------------------------------------------------
void EdgeSwapConfiguration::set( int i )
{
  switch( i )
  {
    case 0: c = &cNull; break;
    case 3: c = &c3; break;
    case 4: c = &c4; break;
    case 5: c = &c5; break;
    case 6: c = &c6; break;
    case 7: c = &c7; break;
    default:
      std::cerr << "Error: Swap configuration not implemented for n="<< i << std::endl;
      throw;
  }
}

// -------------------------------------------------------------------
