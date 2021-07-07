#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "Geometry/rocPackShape.H"
#include "Geometry/hmxShape.H"


namespace NEM {

namespace GEO {


hmxShape::hmxShape()
{
  // Noething
}

hmxShape::~hmxShape()
{
  // Noething
}

std::vector<std::vector<double>> hmxShape::getVertices()
{
  // Adding vertex data
  std::vector<std::vector<double>> verts;
  verts.resize(18);

  verts[0].resize(3);
  verts[0][0] = -1.3973946670992046;
  verts[0][1] = 0.8294297402753090;
  verts[0][2] = 0.7394707591127669;

  verts[1].resize(3);
  verts[1][0] = -0.6657703416741856;
  verts[1][1] = -1.3009854568631263;
  verts[1][2] = 0.4264290990409822;

  verts[2].resize(3);
  verts[2][0] = 1.3973947035090559;
  verts[2][1] = -0.8294297602140067;
  verts[2][2] = -0.7394707784531516;

  verts[3].resize(3);
  verts[3][0] = 0.6657703744131530;
  verts[3][1] = 1.3009854476136657;
  verts[3][2] = -0.4264291168106977;

  verts[4].resize(3);
  verts[4][0] = 0.3835703684025527;
  verts[4][1] = -1.5974467033338673;
  verts[4][2] = -0.1431211054793107;

  verts[5].resize(3);
  verts[5][0] = -0.3835703227123897;
  verts[5][1] = 1.5974466904254172;
  verts[5][2] = 0.1431210806800808;

  verts[6].resize(3);
  verts[6][0] = -0.6302224609427972;
  verts[6][1] = -0.9781790576670774;
  verts[6][2] = -1.1521344230410329;

  verts[7].resize(3);
  verts[7][0] = -1.1453321714573643;
  verts[7][1] = 0.5217677283235427;
  verts[7][2] = -0.9317333328983531;

  verts[8].resize(3);
  verts[8][0] = -0.5198978031779152;
  verts[8][1] = -1.4790343116441396;
  verts[8][2] = -0.5407230456136918;

  verts[9].resize(3);
  verts[9][0] = -0.1637657754207528;
  verts[9][1] = -0.6248174097886718;
  verts[9][2] = -1.4265126187846873;

  verts[10].resize(3);
  verts[10][0] = -0.6788754859353187;
  verts[10][1] = 0.8751293762019474;
  verts[10][2] = -1.2061115286420081;

  verts[11].resize(3);
  verts[11][0] = -0.2376977860165109;
  verts[11][1] = 1.4193978378419216;
  verts[11][2] = -0.8240310520377883;

  verts[12].resize(3);
  verts[12][0] = 0.6788754551281021;
  verts[12][1] = -0.8751294321178871;
  verts[12][2] = 1.2061114962377961;

  verts[13].resize(3);
  verts[13][0] = 0.1637657886385562;
  verts[13][1] = 0.6248173573382960;
  verts[13][2] = 1.4265126656883920;

  verts[14].resize(3);
  verts[14][0] = 0.2376977724807246;
  verts[14][1] = -1.4193978708317796;
  verts[14][2] = 0.8240309982418957;

  verts[15].resize(3);
  verts[15][0] = 1.1453321406501467;
  verts[15][1] = -0.5217677842394817;
  verts[15][2] = 0.9317333004941414;

  verts[16].resize(3);
  verts[16][0] = 0.6302224741606013;
  verts[16][1] = 0.9781790052167011;
  verts[16][2] = 1.1521344699447380;

  verts[17].resize(3);
  verts[17][0] = 0.5198978428216525;
  verts[17][1] = 1.4790342828404699;
  verts[17][2] = 0.5407230876168944;

  return verts;

}

std::vector<std::vector<int>> hmxShape::getFaces()
{
  // Adding face data
  std::vector<std::vector<int>> faces;
  faces.resize(12);

  faces[0].resize(5);
  faces[0][0] = 0;
  faces[0][1] = 7;
  faces[0][2] = 6;
  faces[0][3] = 8;
  faces[0][4] = 1;

  faces[1].resize(5);
  faces[1][0] = 2;
  faces[1][1] = 9;
  faces[1][2] = 10;
  faces[1][3] = 11;
  faces[1][4] = 3;

  faces[2].resize(5);
  faces[2][0] = 0;
  faces[2][1] = 1;
  faces[2][2] = 14;
  faces[2][3] = 12;
  faces[2][4] = 13;

  faces[3].resize(5);
  faces[3][0] = 2;
  faces[3][1] = 3;
  faces[3][2] = 17;
  faces[3][3] = 16;
  faces[3][4] = 15;

  faces[4].resize(4);
  faces[4][0] = 6;
  faces[4][1] = 7;
  faces[4][2] = 10;
  faces[4][3] = 9;

  faces[5].resize(4);
  faces[5][0] = 12;
  faces[5][1] = 15;
  faces[5][2] = 16;
  faces[5][3] = 13;

  faces[6].resize(5);
  faces[6][0] = 2;
  faces[6][1] = 4;
  faces[6][2] = 8;
  faces[6][3] = 6;
  faces[6][4] = 9;

  faces[7].resize(5);
  faces[7][0] = 0;
  faces[7][1] = 5;
  faces[7][2] = 11;
  faces[7][3] = 10;
  faces[7][4] = 7;

  faces[8].resize(5);
  faces[8][0] = 2;
  faces[8][1] = 15;
  faces[8][2] = 12;
  faces[8][3] = 14;
  faces[8][4] = 4;

  faces[9].resize(5);
  faces[9][0] = 0;
  faces[9][1] = 13;
  faces[9][2] = 16;
  faces[9][3] = 17;
  faces[9][4] = 5;

  faces[10].resize(4);
  faces[10][0] = 1;
  faces[10][1] = 8;
  faces[10][2] = 4;
  faces[10][3] = 14;

  faces[11].resize(4);
  faces[11][0] = 3;
  faces[11][1] = 11;
  faces[11][2] = 5;
  faces[11][3] = 17;

  return faces;
}


}

}
