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

#include "SizeFieldBase.h"
#include "MAdOutput.h"
#include "MathUtils.h"
#include "MeshSizeBase.h"

using std::string;

namespace MAd {

  // -------------------------------------------------------------------
  // Length squared computation
  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  double SizeFieldBase::SF_E_lengthSq(const pEdge edge) const
  {
    return SF_VV_lengthSq(E_vertex(edge,0),E_vertex(edge,1));
  }

  // -------------------------------------------------------------------
  // Center of edge computation
  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  double SizeFieldBase::SF_XYZ_center(const double xyz[2][3], const pMSize pS[2],
                                      double center[3], double * reducSq, 
                                      pMSize * cSize) const
  {
    if( !pS[0] || !pS[1] ) {
      printf("Error in SizeFieldBase::SF_XYZ_center: one of the sizes is NULL\n");
      throw;
    }

    double vec[3];
    diffVec(xyz[1],xyz[0],vec);

    double lenSq[2];
    lenSq[0] = pS[0]->lengthSqInDir(vec);
    lenSq[1] = pS[1]->lengthSqInDir(vec);
    double lenRatSqrt = pow(lenSq[1]/lenSq[0],0.25);

    double cParam = 1. / ( 1. + lenRatSqrt );
    *reducSq = 1. / ( 1./lenRatSqrt + lenRatSqrt + 2 );

//     double len[2];
//     len[0] = sqrt(pS[0]->lengthSqInDir(vec));
//     len[1] = sqrt(pS[1]->lengthSqInDir(vec));

//     double cParam = len[0] / ( len[0] + len[1] );

    for (int i=0; i<3; i++) {
      center[i] = cParam * xyz[1][i] + (1-cParam) * xyz[0][i];
    }

    *cSize = MS_interpolate(pS[0],pS[1],cParam);

    return cParam;
  }

  // -------------------------------------------------------------------
  void SizeFieldBase::printPosIsotropic(const pMesh mesh, const string name)
  {
    MAdGmshOutput(mesh, this, name.c_str(), OD_SIZEFIELD_MEAN);
  }

  // -------------------------------------------------------------------
  void SizeFieldBase::printPosAnisotropic(const pMesh mesh, const string baseName)
  {
    string name0 = baseName + "0.pos";
    MAdGmshOutput(mesh, this, name0.c_str(), OD_ANISO_SF_AXIS0);
    string name1 = baseName + "1.pos";
    MAdGmshOutput(mesh, this, name1.c_str(), OD_ANISO_SF_AXIS1);
    string name2 = baseName + "2.pos";
    MAdGmshOutput(mesh, this, name2.c_str(), OD_ANISO_SF_AXIS2);
  }

  // -------------------------------------------------------------------

}
