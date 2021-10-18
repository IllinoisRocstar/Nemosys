#ifndef NEMOSYS_GMSHTYPES_H_
#define NEMOSYS_GMSHTYPES_H_

/**
 * The Gmsh types are taken from gmsh/GmshDefines.h header from Gmsh 4.4.1
 *
 * Useful regex when adding in bulk: (?<=MSH_\w*)\W*(?=[0-9]*$)
 *
 * Note: enums are UNSCOPED. This is to enable compatibility when using gmsh
 * functions as they take the int as the argument.
 */
enum GMSH_ALGO_TYPES {
  // Surface and Volume meshing algorithms
  ALGO_2D_MESHADAPT = 1,
  ALGO_2D_AUTO = 2,
  ALGO_2D_DELAUNAY = 5,
  ALGO_2D_FRONTAL = 6,  // gmsh default
  ALGO_2D_BAMG = 7,
  ALGO_2D_FRONTAL_QUAD = 8,
  ALGO_2D_PACK_PRLGRMS = 9,
  ALGO_2D_PACK_PRLGRMS_CSTR = 10,  // not in documentation
};
enum GMSH_ALGO_3D {
  ALGO_3D_DELAUNAY = 1,  // gmsh default
  ALGO_3D_FRONTAL = 4,
  ALGO_3D_MMG3D = 7,
  ALGO_3D_RTREE = 9,
  ALGO_3D_HXT = 10,
};
// Not in GmshDefines.h
enum GMSH_ALGO_RECOMBINATION {
  ALGO_RECOMBINE_SIMPLE = 0,
  ALGO_RECOMBINE_BLOSSOM = 1,  // gmsh default
  ALGO_RECOMBINE_SIMPLE_FULLQUAD = 2,
  ALGO_RECOMBINE_BLOSSOM_FULLQUAD = 3,
};
enum GMSH_MSH_TYPES {
  // Element types in .msh file format (numbers should not be changed)
  MSH_LIN_2 = 1,
  MSH_TRI_3 = 2,
  MSH_QUA_4 = 3,
  MSH_TET_4 = 4,
  MSH_HEX_8 = 5,
  MSH_PRI_6 = 6,
  MSH_PYR_5 = 7,
  MSH_LIN_3 = 8,
  MSH_TRI_6 = 9,
  MSH_QUA_9 = 10,
  MSH_TET_10 = 11,
  MSH_HEX_27 = 12,
  MSH_PRI_18 = 13,
  MSH_PYR_14 = 14,
  MSH_PNT = 15,
  MSH_QUA_8 = 16,
  MSH_HEX_20 = 17,
  MSH_PRI_15 = 18,
  MSH_PYR_13 = 19,
  MSH_TRI_9 = 20,
  MSH_TRI_10 = 21,
  MSH_TRI_12 = 22,
  MSH_TRI_15 = 23,
  MSH_TRI_15I = 24,
  MSH_TRI_21 = 25,
  MSH_LIN_4 = 26,
  MSH_LIN_5 = 27,
  MSH_LIN_6 = 28,
  MSH_TET_20 = 29,
  MSH_TET_35 = 30,
  MSH_TET_56 = 31,
  MSH_TET_22 = 32,
  MSH_TET_28 = 33,
  MSH_POLYG_ = 34,
  MSH_POLYH_ = 35,
  MSH_QUA_16 = 36,
  MSH_QUA_25 = 37,
  MSH_QUA_36 = 38,
  MSH_QUA_12 = 39,
  MSH_QUA_16I = 40,
  MSH_QUA_20 = 41,
  MSH_TRI_28 = 42,
  MSH_TRI_36 = 43,
  MSH_TRI_45 = 44,
  MSH_TRI_55 = 45,
  MSH_TRI_66 = 46,
  MSH_QUA_49 = 47,
  MSH_QUA_64 = 48,
  MSH_QUA_81 = 49,
  MSH_QUA_100 = 50,
  MSH_QUA_121 = 51,
  MSH_TRI_18 = 52,
  MSH_TRI_21I = 53,
  MSH_TRI_24 = 54,
  MSH_TRI_27 = 55,
  MSH_TRI_30 = 56,
  MSH_QUA_24 = 57,
  MSH_QUA_28 = 58,
  MSH_QUA_32 = 59,
  MSH_QUA_36I = 60,
  MSH_QUA_40 = 61,
  MSH_LIN_7 = 62,
  MSH_LIN_8 = 63,
  MSH_LIN_9 = 64,
  MSH_LIN_10 = 65,
  MSH_LIN_11 = 66,
  MSH_LIN_B = 67,
  MSH_TRI_B = 68,
  MSH_POLYG_B = 69,
  MSH_LIN_C = 70,
  // TETS COMPLETE (6->10)
  MSH_TET_84 = 71,
  MSH_TET_120 = 72,
  MSH_TET_165 = 73,
  MSH_TET_220 = 74,
  MSH_TET_286 = 75,
  // TETS INCOMPLETE (6->10)
  MSH_TET_34 = 79,
  MSH_TET_40 = 80,
  MSH_TET_46 = 81,
  MSH_TET_52 = 82,
  MSH_TET_58 = 83,
  //
  MSH_LIN_1 = 84,
  MSH_TRI_1 = 85,
  MSH_QUA_1 = 86,
  MSH_TET_1 = 87,
  MSH_HEX_1 = 88,
  MSH_PRI_1 = 89,
  MSH_PRI_40 = 90,
  MSH_PRI_75 = 91,
  // HEXES COMPLETE (3->9)
  MSH_HEX_64 = 92,
  MSH_HEX_125 = 93,
  MSH_HEX_216 = 94,
  MSH_HEX_343 = 95,
  MSH_HEX_512 = 96,
  MSH_HEX_729 = 97,
  MSH_HEX_1000 = 98,
  // HEXES INCOMPLETE (3->9)
  MSH_HEX_32 = 99,
  MSH_HEX_44 = 100,
  MSH_HEX_56 = 101,
  MSH_HEX_68 = 102,
  MSH_HEX_80 = 103,
  MSH_HEX_92 = 104,
  MSH_HEX_104 = 105,
  // PRISMS COMPLETE (5->9)
  MSH_PRI_126 = 106,
  MSH_PRI_196 = 107,
  MSH_PRI_288 = 108,
  MSH_PRI_405 = 109,
  MSH_PRI_550 = 110,
  // PRISMS INCOMPLETE (3->9)
  MSH_PRI_24 = 111,
  MSH_PRI_33 = 112,
  MSH_PRI_42 = 113,
  MSH_PRI_51 = 114,
  MSH_PRI_60 = 115,
  MSH_PRI_69 = 116,
  MSH_PRI_78 = 117,
  // PYRAMIDS COMPLETE (3->9)
  MSH_PYR_30 = 118,
  MSH_PYR_55 = 119,
  MSH_PYR_91 = 120,
  MSH_PYR_140 = 121,
  MSH_PYR_204 = 122,
  MSH_PYR_285 = 123,
  MSH_PYR_385 = 124,
  // PYRAMIDS INCOMPLETE (3->9)
  MSH_PYR_21 = 125,
  MSH_PYR_29 = 126,
  MSH_PYR_37 = 127,
  MSH_PYR_45 = 128,
  MSH_PYR_53 = 129,
  MSH_PYR_61 = 130,
  MSH_PYR_69 = 131,
  // Additional types
  MSH_PYR_1 = 132,
  MSH_PNT_SUB = 133,
  MSH_LIN_SUB = 134,
  MSH_TRI_SUB = 135,
  MSH_TET_SUB = 136,
  MSH_TET_16 = 137,
  MSH_TRI_MINI = 138,
  MSH_TET_MINI = 139,
  MSH_TRIH_4 = 140
};

#endif  // NEMOSYS_GMSHTYPES_H_
