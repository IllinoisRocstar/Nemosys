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

#ifndef _H_BACKGROUNDSF
#define _H_BACKGROUNDSF

#include "SizeFieldBase.h"

namespace MAd {

  // -------------------------------------------------------------------
  class BackgroundSF: public SizeFieldBase
  {
  public:

    BackgroundSF(std::string _name="");
    ~BackgroundSF();

    void loadData(std::string fileName);

  public:

    sFieldType getType() const { return BACKGROUNDSFIELD; }

    void scale(double) {throw;}

    // get the size at a location (allocate space!)
    pMSize getSize(const pVertex) const;
    pMSize getSizeOnEntity(const pEntity, const double[3]) const { throw; return 0; }

    // edge length (squared)
    double SF_VV_lengthSq(const pVertex, const pVertex) const { throw; return -1.; }
    double SF_XYZ_lengthSq(const double[3], const double[3],
                           const pMSize, const pMSize=NULL) const { throw; return -1.; }

    // face area (squared)
    double SF_F_areaSq(const pFace) const { throw; return -1.; }
    double SF_XYZ_areaSq(const double[3][3], const pMSize,
                         const double[3]) const { throw; return -1.; }

    // region volume
    double SF_R_volume(const pRegion) const { throw; return -1.; }
    double SF_XYZ_volume(const double[4][3], const pMSize) const { throw; return -1.; }

    // center and its associated size
    double SF_E_center(const pEdge, double[3], double * reducSq, pMSize *) const { throw; return -1.; }
    double SF_VV_center(const pVertex, const pVertex,
                        double[3], double * reducSq, pMSize *) const { throw; return -1.; }

  private:

    void setSize(pEntity, pMSize);
    void setSize(pEntity, double);

    void redistributeToVertices();

    pGModel bgModel;
    pMesh bgMesh;
    pMeshDataId pMSizeFieldId;
  };

}

// -------------------------------------------------------------------
#endif
