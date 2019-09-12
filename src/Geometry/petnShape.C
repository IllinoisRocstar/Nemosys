#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <rocPackShape.H>
#include <petnShape.H>


namespace NEM {

namespace GEO {


petnShape::petnShape()
{
  // Nothing
}

petnShape::~petnShape()
{
  // Nothing
}


std::vector<std::vector<double>> petnShape::getVertices()
{
  // Adding vertex data
  std::vector<std::vector<double>> verts;
  verts.resize(18);

  verts[0].resize(3);
  verts[0][0] = 0.00;
  verts[0][1] = -0.45;
  verts[0][2] = 0.00;

  verts[1].resize(3);
  verts[1][0] = 0.20;
  verts[1][1] = -0.25;
  verts[1][2] = -0.20;

  verts[2].resize(3);
  verts[2][0] = 0.20;
  verts[2][1] = -0.35;
  verts[2][2] = 0.00;

  verts[3].resize(3);
  verts[3][0] = 0.00;
  verts[3][1] = -0.35;
  verts[3][2] = -0.20;

  verts[4].resize(3);
  verts[4][0] = 0.20;
  verts[4][1] = -0.25;
  verts[4][2] = 0.20;

  verts[5].resize(3);
  verts[5][0] = 0.00;
  verts[5][1] = -0.35;
  verts[5][2] = 0.20;

  verts[6].resize(3);
  verts[6][0] = -0.20;
  verts[6][1] = -0.25;
  verts[6][2] = 0.20;

  verts[7].resize(3);
  verts[7][0] = -0.20;
  verts[7][1] = -0.35;
  verts[7][2] = 0.00;

  verts[8].resize(3);
  verts[8][0] = -0.20;
  verts[8][1] = -0.25;
  verts[8][2] = -0.20;

  verts[9].resize(3);
  verts[9][0] = 0.00;
  verts[9][1] = 0.45;
  verts[9][2] = 0.00;

  verts[10].resize(3);
  verts[10][0] = 0.20;
  verts[10][1] = 0.25;
  verts[10][2] = -0.20;

  verts[11].resize(3);
  verts[11][0] = 0.20;
  verts[11][1] = 0.35;
  verts[11][2] = 0.00;

  verts[12].resize(3);
  verts[12][0] = 0.00;
  verts[12][1] = 0.35;
  verts[12][2] = -0.20;

  verts[13].resize(3);
  verts[13][0] = 0.20;
  verts[13][1] = 0.25;
  verts[13][2] = 0.20;

  verts[14].resize(3);
  verts[14][0] = 0.00;
  verts[14][1] = 0.35;
  verts[14][2] = 0.20;

  verts[15].resize(3);
  verts[15][0] = -0.20;
  verts[15][1] = 0.25;
  verts[15][2] = 0.20;

  verts[16].resize(3);
  verts[16][0] = -0.20;
  verts[16][1] = 0.35;
  verts[16][2] = 0.00;

  verts[17].resize(3);
  verts[17][0] = -0.20;
  verts[17][1] = 0.25;
  verts[17][2] = -0.20;

  return verts;
}

std::vector<std::vector<int>> petnShape::getFaces()
{
  // Adding face data
  std::vector<std::vector<int>> faces;
  faces.resize(12);

  faces[0].resize(6);
  faces[0][0] = 10;
  faces[0][1] = 11;
  faces[0][2] = 13;
  faces[0][3] = 4;
  faces[0][4] = 2;
  faces[0][5] = 1;

  faces[1].resize(6);
  faces[1][0] = 13;
  faces[1][1] = 14;
  faces[1][2] = 15;
  faces[1][3] = 6;
  faces[1][4] = 5;
  faces[1][5] = 4;

  faces[2].resize(6);
  faces[2][0] = 15;
  faces[2][1] = 16;
  faces[2][2] = 17;
  faces[2][3] = 8;
  faces[2][4] = 7;
  faces[2][5] = 6;

  faces[3].resize(4);
  faces[3][0] = 14;
  faces[3][1] = 13;
  faces[3][2] = 11;
  faces[3][3] = 9;

  faces[4].resize(4);
  faces[4][0] = 16;
  faces[4][1] = 15;
  faces[4][2] = 14;
  faces[4][3] = 9;

  faces[5].resize(4);
  faces[5][0] = 7;
  faces[5][1] = 8;
  faces[5][2] = 3;
  faces[5][3] = 0;

  faces[6].resize(4);
  faces[6][0] = 3;
  faces[6][1] = 1;
  faces[6][2] = 2;
  faces[6][3] = 0;

  faces[7].resize(4);
  faces[7][0] = 11;
  faces[7][1] = 10;
  faces[7][2] = 12;
  faces[7][3] = 9;

  faces[8].resize(4);
  faces[8][0] = 5;
  faces[8][1] = 6;
  faces[8][2] = 7;
  faces[8][3] = 0;

  faces[9].resize(6);
  faces[9][0] = 3;
  faces[9][1] = 8;
  faces[9][2] = 17;
  faces[9][3] = 12;
  faces[9][4] = 10;
  faces[9][5] = 1;

  faces[10].resize(4);
  faces[10][0] = 9;
  faces[10][1] = 12;
  faces[10][2] = 17;
  faces[10][3] = 16;

  faces[11].resize(4);
  faces[11][0] = 5;
  faces[11][1] = 0;
  faces[11][2] = 2;
  faces[11][3] = 4;

  return faces;
}

}

}
