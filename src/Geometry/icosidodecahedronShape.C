#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <rocPackShape.H>
#include <icosidodecahedronShape.H>

namespace NEM {

namespace GEO {


icosidodecahedronShape::icosidodecahedronShape()
{
  // Nothing
}

icosidodecahedronShape::~icosidodecahedronShape()
{
  // Nothing
}

std::vector<std::vector<double>> icosidodecahedronShape::getVertices()
{
  // Adding vertex data
  std::vector<std::vector<double>> verts;
  verts.resize(30);

  for (int i=0; i<verts.size(); i++)
    verts[i].resize(3);

  verts[0][0] = 0.35682146441852025;
  verts[0][1] = -1.5265543675337354e-17;
  verts[0][2] = 0.9341725978266664;

  verts[1][0] = 0.7557598656204081;
  verts[1][1] = 0.3090165361727076;
  verts[1][2] = 0.5773491334111481;

  verts[2][0] = 0.4670852989148342;
  verts[2][1] = 0.8090157856843972;
  verts[2][2] = 0.35682146441852025;

  verts[3][0] = -0.17841073220926015;
  verts[3][1] = 0.30899953619822423;
  verts[3][2] = 0.9341725978266664;

  verts[4][0] = -0.6454960311240943;
  verts[4][1] = 0.49998224953720627;
  verts[4][2] = 0.5773491334111481;

  verts[5][0] = -0.9341705978296684;
  verts[5][1] = -1.5265543675337354e-17;
  verts[5][2] = 0.35682146441852025;

  verts[6][0] = -0.17841073220926015;
  verts[6][1] = -0.30899953619822423;
  verts[6][2] = 0.9341725978266664;

  verts[7][0] = -0.1102638344963139;
  verts[7][1] = -0.8089987857099138;
  verts[7][2] = 0.5773491334111481;

  verts[8][0] = 0.7557598656204081;
  verts[8][1] = -0.30899953619822423;
  verts[8][2] = 0.5773491334111481;

  verts[9][0] = 0.8660237001167221;
  verts[9][1] = -0.49998224953720627;
  verts[9][2] = -2.220442716412706e-17;

  verts[10][0] = -0.1102638344963139;
  verts[10][1] = 0.8089987857099138;
  verts[10][2] = 0.5773491334111481;

  verts[11][0] = -1.480295144275137e-17;
  verts[11][1] = 0.9999984990233793;
  verts[11][2] = -2.220442716412706e-17;

  verts[12][0] = -0.6454960311240943;
  verts[12][1] = -0.49998224953720627;
  verts[12][2] = 0.5773491334111481;

  verts[13][0] = -0.8660237001167221;
  verts[13][1] = -0.49999924951168967;
  verts[13][2] = -2.220442716412706e-17;

  verts[14][0] = 0.8660237001167221;
  verts[14][1] = 0.49999924951168967;
  verts[14][2] = -2.220442716412706e-17;

  verts[15][0] = 0.6454960311240943;
  verts[15][1] = 0.49998224953720627;
  verts[15][2] = -0.5773491334111481;

  verts[16][0] = 0.4670852989148342;
  verts[16][1] = -0.8089987857099138;
  verts[16][2] = 0.35682146441852025;

  verts[17][0] = -1.480295144275137e-17;
  verts[17][1] = -0.9999984990233793;
  verts[17][2] = -2.220442716412706e-17;

  verts[18][0] = 0.11026383449631387;
  verts[18][1] = -0.8089987857099138;
  verts[18][2] = -0.5773491334111481;

  verts[19][0] = -0.8660237001167221;
  verts[19][1] = 0.49998224953720627;
  verts[19][2] = -2.220442716412706e-17;

  verts[20][0] = -0.7557598656204081;
  verts[20][1] = 0.30899953619822423;
  verts[20][2] = -0.5773491334111481;

  verts[21][0] = 0.9341705978296684;
  verts[21][1] = -1.5265543675337354e-17;
  verts[21][2] = -0.35682146441852025;

  verts[22][0] = 0.6454960311240943;
  verts[22][1] = -0.49998224953720627;
  verts[22][2] = -0.5773491334111481;

  verts[23][0] = 0.1784107322092601;
  verts[23][1] = -0.30899953619822423;
  verts[23][2] = -0.9341725978266664;

  verts[24][0] = -0.4670852989148342;
  verts[24][1] = 0.8089987857099138;
  verts[24][2] = -0.35682146441852025;

  verts[25][0] = 0.11026383449631387;
  verts[25][1] = 0.8089987857099138;
  verts[25][2] = -0.5773491334111481;

  verts[26][0] = 0.1784107322092601;
  verts[26][1] = 0.30899953619822423;
  verts[26][2] = -0.9341725978266664;

  verts[27][0] = -0.4670852989148342;
  verts[27][1] = -0.8090157856843972;
  verts[27][2] = -0.35682146441852025;

  verts[28][0] = -0.7557598656204081;
  verts[28][1] = -0.3090165361727076;
  verts[28][2] = -0.5773491334111481;

  verts[29][0] = -0.35682146441852025;
  verts[29][1] = -1.5265543675337354e-17;
  verts[29][2] = -0.9341725978266664;

  return verts;
}

std::vector<std::vector<int>> icosidodecahedronShape::getFaces()
{
  // Adding face data
  std::vector<std::vector<int>> faces;
  faces.resize(32);

  faces[0].resize(3);
  faces[0][0] = 23;
  faces[0][1] = 29;
  faces[0][2] = 26;

  faces[1].resize(3);
  faces[1][0] = 15;
  faces[1][1] = 26;
  faces[1][2] = 25;

  faces[2].resize(3);
  faces[2][0] = 20;
  faces[2][1] = 29;
  faces[2][2] = 28;

  faces[3].resize(3);
  faces[3][0] = 18;
  faces[3][1] = 23;
  faces[3][2] = 22;

  faces[4].resize(5);
  faces[4][0] = 18;
  faces[4][1] = 27;
  faces[4][2] = 28;
  faces[4][3] = 29;
  faces[4][4] = 23;

  faces[5].resize(5);
  faces[5][0] = 20;
  faces[5][1] = 24;
  faces[5][2] = 25;
  faces[5][3] = 26;
  faces[5][4] = 29;

  faces[6].resize(5);
  faces[6][0] = 15;
  faces[6][1] = 21;
  faces[6][2] = 22;
  faces[6][3] = 23;
  faces[6][4] = 26;

  faces[7].resize(3);
  faces[7][0] = 27;
  faces[7][1] = 18;
  faces[7][2] = 17;

  faces[8].resize(3);
  faces[8][0] = 13;
  faces[8][1] = 28;
  faces[8][2] = 27;

  faces[9].resize(3);
  faces[9][0] = 9;
  faces[9][1] = 22;
  faces[9][2] = 21;

  faces[10].resize(3);
  faces[10][0] = 21;
  faces[10][1] = 15;
  faces[10][2] = 14;

  faces[11].resize(3);
  faces[11][0] = 11;
  faces[11][1] = 25;
  faces[11][2] = 24;

  faces[12].resize(3);
  faces[12][0] = 24;
  faces[12][1] = 20;
  faces[12][2] = 19;

  faces[13].resize(5);
  faces[13][0] = 9;
  faces[13][1] = 16;
  faces[13][2] = 17;
  faces[13][3] = 18;
  faces[13][4] = 22;

  faces[14].resize(5);
  faces[14][0] = 5;
  faces[14][1] = 19;
  faces[14][2] = 20;
  faces[14][3] = 28;
  faces[14][4] = 13;

  faces[15].resize(5);
  faces[15][0] = 2;
  faces[15][1] = 14;
  faces[15][2] = 15;
  faces[15][3] = 25;
  faces[15][4] = 11;

  faces[16].resize(5);
  faces[16][0] = 7;
  faces[16][1] = 12;
  faces[16][2] = 13;
  faces[16][3] = 27;
  faces[16][4] = 17;

  faces[17].resize(5);
  faces[17][0] = 4;
  faces[17][1] = 10;
  faces[17][2] = 11;
  faces[17][3] = 24;
  faces[17][4] = 19;

  faces[18].resize(5);
  faces[18][0] = 1;
  faces[18][1] = 8;
  faces[18][2] = 9;
  faces[18][3] = 21;
  faces[18][4] = 14;

  faces[19].resize(3);
  faces[19][0] = 7;
  faces[19][1] = 17;
  faces[19][2] = 16;

  faces[20].resize(3);
  faces[20][0] = 5;
  faces[20][1] = 13;
  faces[20][2] = 12;

  faces[21].resize(3);
  faces[21][0] = 4;
  faces[21][1] = 19;
  faces[21][2] = 5;

  faces[22].resize(3);
  faces[22][0] = 16;
  faces[22][1] = 9;
  faces[22][2] = 8;

  faces[23].resize(3);
  faces[23][0] = 2;
  faces[23][1] = 11;
  faces[23][2] = 10;

  faces[24].resize(3);
  faces[24][0] = 1;
  faces[24][1] = 14;
  faces[24][2] = 2;

  faces[25].resize(5);
  faces[25][0] = 3;
  faces[25][1] = 0;
  faces[25][2] = 1;
  faces[25][3] = 2;
  faces[25][4] = 10;

  faces[26].resize(5);
  faces[26][0] = 0;
  faces[26][1] = 6;
  faces[26][2] = 7;
  faces[26][3] = 16;
  faces[26][4] = 8;

  faces[27].resize(5);
  faces[27][0] = 6;
  faces[27][1] = 3;
  faces[27][2] = 4;
  faces[27][3] = 5;
  faces[27][4] = 12;

  faces[28].resize(3);
  faces[28][0] = 10;
  faces[28][1] = 4;
  faces[28][2] = 3;

  faces[29].resize(3);
  faces[29][0] = 0;
  faces[29][1] = 8;
  faces[29][2] = 1;

  faces[30].resize(3);
  faces[30][0] = 12;
  faces[30][1] = 7;
  faces[30][2] = 6;

  faces[31].resize(3);
  faces[31][0] = 3;
  faces[31][1] = 6;
  faces[31][2] = 0;

  return faces;
}

}
}