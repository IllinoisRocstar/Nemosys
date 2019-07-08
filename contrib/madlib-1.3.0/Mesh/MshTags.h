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

#ifndef MSHTAGS__H
#define MSHTAGS__H

#include <string>
#include <sstream>

// Element types in .msh file format
#define MSH_LIN_2  1
#define MSH_TRI_3  2
#define MSH_QUA_4  3
#define MSH_TET_4  4
#define MSH_HEX_8  5
#define MSH_PRI_6  6
#define MSH_PYR_5  7
#define MSH_LIN_3  8
#define MSH_TRI_6  9
#define MSH_QUA_9  10
#define MSH_TET_10 11
#define MSH_HEX_27 12
#define MSH_PRI_18 13
#define MSH_PYR_14 14
#define MSH_PNT    15
#define MSH_QUA_8  16
#define MSH_HEX_20 17
#define MSH_PRI_15 18
#define MSH_PYR_13 19
#define MSH_TRI_9  20
#define MSH_TRI_10 21
#define MSH_TRI_12 22
#define MSH_TRI_15 23
#define MSH_TRI_15I 24
#define MSH_TRI_21 25
#define MSH_LIN_4  26
#define MSH_LIN_5  27
#define MSH_LIN_6  28
#define MSH_TET_20 29
#define MSH_TET_35 30
#define MSH_TET_56 31
// AEG: This is changed in the official Gmsh declaration.
//#define MSH_TET_34 32
#define MSH_TET_34 79
// /AEG
#define MSH_UNDEFINED_ELEM 33
#define MSH_MAX_ELEMENT_NODES 56

namespace MAd {

  inline
  int getNumVerticesForElementTypeMSH(int type)
  {
    switch (type) {
    case MSH_PNT    : return 1;
    case MSH_LIN_2  : return 2;
    case MSH_LIN_3  : return 2 + 1;
    case MSH_LIN_4  : return 2 + 2;
    case MSH_LIN_5  : return 2 + 3;
    case MSH_LIN_6  : return 2 + 4;
    case MSH_TRI_3  : return 3;
    case MSH_TRI_6  : return 3 + 3;
    case MSH_TRI_9  : return 3 + 6;
    case MSH_TRI_10  : return 3 + 6 + 1;
    case MSH_TRI_12  : return 3 + 9;
    case MSH_TRI_15  : return 3 + 9 + 3;
    case MSH_TRI_15I : return 3 + 12;
    case MSH_TRI_21 : return 3 + 12 + 6;
    case MSH_QUA_4  : return 4;
    case MSH_QUA_8  : return 4 + 4;
    case MSH_QUA_9  : return 4 + 4 + 1;
    case MSH_TET_4  : return 4;
    case MSH_TET_10 : return 4 + 6;
    case MSH_HEX_8  : return 8;
    case MSH_HEX_20 : return 8 + 12;
    case MSH_HEX_27 : return 8 + 12 + 6 + 1;
    case MSH_PRI_6  : return 6;
    case MSH_PRI_15 : return 6 + 9;
    case MSH_PRI_18 : return 6 + 9 + 3;
    case MSH_PYR_5  : return 5;
    case MSH_PYR_13 : return 5 + 8;
    case MSH_PYR_14 : return 5 + 8 + 1;
    case MSH_TET_20 : return 20;
    case MSH_TET_35 : return 35;
    case MSH_TET_34 : return 34;
    case MSH_TET_56 : return 56;
    default: 
      printf("Error (MshTags.h): Unknown type of element %d\n", type);
    }
    return 0;
  }

  inline
  int getDimForElementTypeMSH(int type)
  {
    switch (type) {
    case MSH_PNT    : return 0;
    case MSH_LIN_2  : 
    case MSH_LIN_3  : 
    case MSH_LIN_4  : 
    case MSH_LIN_5  : 
    case MSH_LIN_6  : return 1;
    case MSH_TRI_3  : 
    case MSH_TRI_6  : 
    case MSH_TRI_9  : 
    case MSH_TRI_10  : 
    case MSH_TRI_12  : 
    case MSH_TRI_15  : 
    case MSH_TRI_15I : 
    case MSH_TRI_21 : 
    case MSH_QUA_4  : 
    case MSH_QUA_8  : 
    case MSH_QUA_9  : return 2;
    case MSH_TET_4  : 
    case MSH_TET_10 : 
    case MSH_HEX_8  : 
    case MSH_HEX_20 : 
    case MSH_HEX_27 : 
    case MSH_PRI_6  : 
    case MSH_PRI_15 : 
    case MSH_PRI_18 : 
    case MSH_PYR_5  : 
    case MSH_PYR_13 : 
    case MSH_PYR_14 : 
    case MSH_TET_20 : 
    case MSH_TET_35 : 
    case MSH_TET_34 : 
    case MSH_TET_56 : return 3;
    default: 
      printf("Error (MshTags.h): Unknown type of element %d\n", type);
    }
    return -1;
  }

  inline
  std::string getElementName(int type) {

    std::string name;
  
    switch (type) {
    case MSH_PNT     : name = "Point"; break;
    case MSH_LIN_2   : name = "Linear edge"; break;
    case MSH_LIN_3   : name = "Quadratic edge"; break;
    case MSH_LIN_4   : name = "Cubic edge"; break;
    case MSH_LIN_5   : name = "Quartic edge"; break;
    case MSH_LIN_6   : name = "Pentic edge"; break;
    case MSH_TRI_3   : name = "Linear triangle"; break;
    case MSH_TRI_6   : name = "Quadratic triangle"; break;
    case MSH_TRI_9   : name = "Cubic serendipity triangle"; break;
    case MSH_TRI_10  : name = "Cubic triangle"; break;
    case MSH_TRI_12  : name = "Quartic serendipity triangle"; break;
    case MSH_TRI_15  : name = "Quartic triangle"; break;
    case MSH_TRI_15I : name = "Pentic serendipity triangle"; break;
    case MSH_TRI_21  : name = "Pentic triangle"; break;
    case MSH_QUA_4   : name = "Bilinear Quadrangle"; break;
    case MSH_QUA_8   : name = "Quadratic serendipity quadrangle"; break;
    case MSH_QUA_9   : name = "Quadratic quadrangle"; break;
    case MSH_TET_4   : name = "Linear tetrahedron"; break;
    case MSH_TET_10  : name = "Quadratic tetrahedron"; break;
    case MSH_HEX_8   : name = "Trilinear hexahedron"; break;
    case MSH_HEX_20  : name = "Quadratic edge serendipity hexahedron"; break;
    case MSH_HEX_27  : name = "Quadratic serendipity hexahedron"; break;
    case MSH_PRI_6   : name = "Linear prism"; break;
    case MSH_PRI_15  : name = "Quadratic edge serendipity prism"; break;
    case MSH_PRI_18  : name = "Quadratic serendipity prism"; break;
    case MSH_PYR_5   : name = "Linear pyramid"; break;
    case MSH_PYR_13  : name = "Quadratic edge serendipity pyramid"; break;
    case MSH_PYR_14  : name = "Quadratic serendipty pyramid"; break;
    case MSH_TET_20  : name = "Cubic tetrahedron";break;
    case MSH_TET_35  : name = "Quartic tetrahedron";break;
    case MSH_TET_34  : name = "Quartic serendipity tetrahedron";break;
    case MSH_TET_56  : name = "Pentic tetrahedron";break;
    default:
      std::stringstream ss; break;
      ss <<  "Unknown type of element (tag " << type << ")"; break;
      name = ss.str(); 
    }
    return name;
  }

}

#endif
