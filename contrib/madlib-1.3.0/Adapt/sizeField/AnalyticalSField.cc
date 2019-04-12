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

#include "AnalyticalSField.h"
#include "MAdTimeManager.h"
#include "IsoMeshSize.h"
#include "AnisoMeshSize.h"
#include "MathUtils.h"
#include "MAdMessage.h"
#include "MAdStringFieldEvaluator.h"
#include "MAdDefines.h"

#include <math.h>
#include <stdio.h>
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

namespace MAd {

  // -------------------------------------------------------------------
  // Constructors / destructors
  // -------------------------------------------------------------------
  AnalyticalSField::AnalyticalSField():
    SizeFieldBase(), sFct(NULL), isotropic(true), evaluator(NULL)
  {
    setSize("1.");
  }

  // -------------------------------------------------------------------
  AnalyticalSField::AnalyticalSField(std::string h):
    SizeFieldBase(), sFct(NULL), isotropic(true), evaluator(NULL)
  {
    setSize(h);
  }

  // -------------------------------------------------------------------
  AnalyticalSField::AnalyticalSField(std::vector<std::string> _h,
                                     std::vector<std::string> _e0,
                                     std::vector<std::string> _e1,
                                     std::vector<std::string> _e2):
    SizeFieldBase(), sFct(NULL), isotropic(false), evaluator(NULL)
  {
    setSize(_h,_e0,_e1,_e2);
  }

  // -------------------------------------------------------------------
  AnalyticalSField::AnalyticalSField(sizeFunction f): 
    SizeFieldBase(), sFct(f), isotropic(false),
    evaluator(NULL)
  {}

  // -------------------------------------------------------------------
  AnalyticalSField::~AnalyticalSField()
  {
    if (evaluator) delete evaluator;
  }

  // -------------------------------------------------------------------
  void AnalyticalSField::describe() const {
  
    cout << "\nDescribing analytical size field: \n\n";
  
    cout << "  Orientation:   \t";
    if (isotropic) cout << "Isotropic\n\n";
    else           cout << "Anisotropic\n\n";

    cout << "  Representation:\t";
    if (sFct) cout << "Size Function\n\n";
    else {
      if (isotropic) {
        cout << "String:\n";
        cout << "   Size:  " << h0 << "\n\n";
      }
      else {
        cout << "Strings:\n";

        cout << "   - Sizes:\n";
        cout << "     * " << h0 << "\n";
        cout << "     * " << h1 << "\n";
        cout << "     * " << h2 << "\n";

        cout << "   - Vectors:\n";
        cout << "     * " << e0[0] <<"\t" << e0[1] <<"\t" << e0[2] <<"\n";
        cout << "     * " << e1[0] <<"\t" << e1[1] <<"\t" << e1[2] <<"\n";
        cout << "     * " << e2[0] <<"\t" << e2[1] <<"\t" << e2[2] <<"\n";

        cout << "\n";
      }
    }

  }

  // -------------------------------------------------------------------
  // Size imposition
  // -------------------------------------------------------------------
  void AnalyticalSField::setSize(const std::string h)
  {
    if (sFct) throw;

    isotropic = true;

    h0 = h;
    h1 = h;
    h2 = h;
    e0.resize(3);
    e1.resize(3);
    e2.resize(3);
    e0[0] = "1."; e0[1] = "0."; e0[2] = "0.";
    e1[0] = "0."; e1[1] = "1."; e1[2] = "0.";
    e2[0] = "0."; e2[1] = "0."; e2[2] = "1.";

    if(evaluator) delete evaluator;
    evaluator = new MAdStringFieldEvaluator(1,h.c_str());
  }

  // -------------------------------------------------------------------
  void AnalyticalSField::setSize(std::vector<std::string> _h,
                                 std::vector<std::string> _e0,
                                 std::vector<std::string> _e1,
                                 std::vector<std::string> _e2)
  {
    if (sFct) throw;

    isotropic = false;

    h0 = _h[0];
    h1 = _h[1];
    h2 = _h[2];
    e0 = _e0;
    e1 = _e1;
    e2 = _e2;

    if(evaluator) delete evaluator;
    evaluator = new MAdStringFieldEvaluator(12,
                                            h0.c_str(), h1.c_str(), h2.c_str(), 
                                            e0[0].c_str(), e0[1].c_str(), e0[2].c_str(),
                                            e1[0].c_str(), e1[1].c_str(), e1[2].c_str(),
                                            e2[0].c_str(), e2[1].c_str(), e2[2].c_str());
  }

  // -------------------------------------------------------------------
  void AnalyticalSField::setSize(sizeFunction f)
  {
    if (evaluator) throw;
    sFct = f;
  }

  // -------------------------------------------------------------------
  // Evaluate a local size
  // -------------------------------------------------------------------
  pMSize AnalyticalSField::eval(const double xyz[3]) const
  {
    double time = MAdTimeManagerSgl::instance().getTime();

    if(sFct)  return (*sFct)(xyz,time);

    else {
      if (isotropic) {
        double h;
        evaluator->eval(xyz,time,&h);
        return new IsoMeshSize(h);
      }
      else {
        double vals[12];
        evaluator->eval(xyz,time,vals);
        double h[3]        = {vals[0], vals[1], vals[2]};
        double e[3][3]     = { {vals[3], vals[4], vals[5]},
                               {vals[6], vals[7], vals[8]},
                               {vals[9], vals[10], vals[11]} };
        return new AnisoMeshSize(e,h);
      }
    }
    return NULL;
  }

  // -------------------------------------------------------------------
  pMSize AnalyticalSField::getSize(const pVertex pv) const
  {
    double xyz[3];
    V_coord(pv,xyz);
    return eval(xyz);
  }

  // -------------------------------------------------------------------
  pMSize AnalyticalSField::getSizeOnEntity(const pEntity, 
                                           const double xyz[3]) const
  {
    return eval(xyz);
  }

  // -------------------------------------------------------------------
  // Length squared computation
  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  double AnalyticalSField::SF_VV_lengthSq(const pVertex pV0, 
                                          const pVertex pV1) const
  {
    double xyz[2][3];
    V_coord(pV0,xyz[0]);
    V_coord(pV1,xyz[1]);

    pMSize pS[2];
    pS[0] = getSize(pV0);
    pS[1] = getSize(pV1);

    double lSq = SF_XYZ_lengthSq(xyz[0],xyz[1],pS[0],pS[1]);
  
    if ( pS[0] ) delete pS[0];
    if ( pS[1] ) delete pS[1];

    return lSq;
  }

  // -------------------------------------------------------------------
  double AnalyticalSField::SF_XYZ_lengthSq(const double xyz0[3], 
                                           const double xyz1[3], 
                                           const pMSize pS0, 
                                           const pMSize pS1) const
  {
    if( pS0 )
      {
        double e[3];
        diffVec(xyz0,xyz1,e);
        double lenSq0 = pS0->normSq(e);
        if ( pS1 )
          {
            double lenSq1 = pS1->normSq(e);
            return sqrt(lenSq0*lenSq1);
          }
        else return lenSq0;
      }
    else {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,"No size defined");
    }
    return 0.;
  }

  // -------------------------------------------------------------------
  // Area squared computation
  // -------------------------------------------------------------------

  double AnalyticalSField::SF_F_areaSq(const pFace face) const
  {
    double area = 0.;

    double xyz[3][3];
    F_coordP1(face,xyz);
  
    void * temp = 0;
    pPList fVerts = F_vertices(face,1);
    while( pVertex pV = (pVertex)PList_next(fVerts,&temp) )
      {
        pMSize pS = getSize(pV);
        area += SF_XYZ_areaSq(xyz,pS,0);
        if (pS) delete pS;
      }
    PList_delete(fVerts);
  
    area /= F_numVertices(face);

    return area;
  }

  // -------------------------------------------------------------------
  double AnalyticalSField::SF_XYZ_areaSq(const double fxyz[3][3], 
                                         const pMSize pS, 
                                         const double norDir[3]) const
  {
    if( !pS ) {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,"No size defined");
    }

    // get the two first edges
    double e01[3],e02[3];
    diffVec(fxyz[1],fxyz[0],e01);
    diffVec(fxyz[2],fxyz[0],e02);
    
    double nor[3];
    crossProd(e01,e02,nor);

    double l1SqInv = 1. / pS->lengthSqInDir(e01);
    double l2SqInv = 1. / pS->lengthSqInDir(e02);

    if( norDir && dotProd(norDir,nor) < MAdTOL ) return 0.;

    double areaSq = 0.25 * dotProd(nor,nor) * l1SqInv * l2SqInv;
    if( areaSq < MAdTOL ) return 0.;

    return areaSq;
  }

  // -------------------------------------------------------------------
  // Volume computation
  // -------------------------------------------------------------------

  double AnalyticalSField::SF_R_volume(const pRegion region) const
  {
    double vol = 0.;

    double xyz[4][3];
    R_coordP1(region,xyz);

    pPList rVerts = R_vertices(region);
    void * temp = 0;
    while( pVertex pV = (pVertex)PList_next(rVerts,&temp) )
      {
        pMSize pS = getSize(pV);
        vol += SF_XYZ_volume(xyz,pS);
        if (pS) delete pS;
      }
    PList_delete(rVerts);

    vol /= R_numVertices(region);

    return vol;
  }

  // -------------------------------------------------------------------
  double AnalyticalSField::SF_XYZ_volume(const double xyz[4][3], 
                                         const pMSize pS) const
  {
    if( !pS ) {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,"No size defined");
    }

    double physVol = R_XYZ_volume(xyz);
  
    return ( physVol / (pS->size(0)*pS->size(1)*pS->size(2)) );
  }

  // -------------------------------------------------------------------
  // Center of edge computation
  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  double AnalyticalSField::SF_E_center(const pEdge edge, double center[3], 
                                       double * reducSq, pMSize * cSize) const
  {
    return SF_VV_center(E_vertex(edge,0),E_vertex(edge,1),center,reducSq,cSize);
  }

  // -------------------------------------------------------------------
  double AnalyticalSField::SF_VV_center(const pVertex v0, const pVertex v1,
                                        double center[3], double * reducSq, 
                                        pMSize * cSize) const
  {
    double xyz[2][3];
    V_coord(v0,xyz[0]);
    V_coord(v1,xyz[1]);

    pMSize pS[2];
    pS[0] = getSize(v0);
    pS[1] = getSize(v1);

    double cParam = SF_XYZ_center(xyz,pS,center,reducSq,cSize);

    if ( pS[0] ) delete pS[0];
    if ( pS[1] ) delete pS[1];

    return cParam;
  }

  // -------------------------------------------------------------------
  void AnalyticalSField::scale(double fact)
  {
    printf("Not implemented (AnalyticalSField::scale)\n");
    throw;  
  }

  // -------------------------------------------------------------------

}
