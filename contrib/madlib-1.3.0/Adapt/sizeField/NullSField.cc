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

#include "NullSField.h"
#include "MathUtils.h"

namespace MAd {

  // -------------------------------------------------------------------
  // Length squared computation
  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  double NullSField::SF_VV_lengthSq(const pVertex v0, const pVertex v1) const
  {
    double xyz[2][3];
    V_coord(v0,xyz[0]);
    V_coord(v1,xyz[1]);

    return SF_XYZ_lengthSq(xyz[0],xyz[1],NULL,NULL);
  }

  // -------------------------------------------------------------------
  double NullSField::SF_XYZ_lengthSq(const double xyz0[3], const double xyz1[3], 
                                     const pMSize, const pMSize) const
  {
    double vec[3];
    diffVec(xyz0,xyz1,vec);

    return dotProd(vec,vec); 
  }

  // -------------------------------------------------------------------
  // Area squared computation
  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  double NullSField::SF_F_areaSq(const pFace face) const
  { 
    double xyz[3][3];
    F_coordP1(face,xyz);

    return SF_XYZ_areaSq(xyz,NULL,NULL);
  }

  // -------------------------------------------------------------------
  double NullSField::SF_XYZ_areaSq(const double xyz[3][3], const pMSize, 
                                   const double norDir[3]) const
  {
    return XYZ_F_areaSq(xyz,norDir);
  }

  // -------------------------------------------------------------------
  // Volume computation
  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  double NullSField::SF_R_volume(const pRegion r) const
  {
    return R_volume(r);
  }

  // -------------------------------------------------------------------
  double NullSField::SF_XYZ_volume(const double xyz[4][3], 
                                   const pMSize) const
  {
    return R_XYZ_volume(xyz);
  }

  // -------------------------------------------------------------------
  // Center of edge computation
  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  double NullSField::SF_E_center(const pEdge edge, 
                                 double center[3], double * reducSq, 
                                 pMSize *) const
  {
    return SF_VV_center(E_vertex(edge,0),E_vertex(edge,1),center,reducSq,NULL);
  }

  // -------------------------------------------------------------------
  double NullSField::SF_VV_center(const pVertex v0, const pVertex v1, 
                                  double center[3], double * reducSq, 
                                  pMSize *) const
  {
    double xyz[2][3];
    V_coord(v0,xyz[0]);
    V_coord(v1,xyz[1]);

    for(int i=0; i<3; i++) center[i] = 0.5*(xyz[0][i]+xyz[1][i]);
    *reducSq = 0.5;
    return 0.5;
  }

  // -------------------------------------------------------------------

}
