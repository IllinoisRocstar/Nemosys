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

#ifndef _H_DISCRETESF
#define _H_DISCRETESF

#include "madlib_export.h"

#include "SizeFieldBase.h"

#include <string>

namespace MAd {

  // -------------------------------------------------------------------
  enum DiscreteSFType {
    UNKNOWN_DSFTYPE,
    VERTEX_P1_DSFTYPE
  };

  // -------------------------------------------------------------------
  class MADLIB_EXPORT DiscreteSF : public SizeFieldBase
  {
  public:

    DiscreteSF(pMesh, std::string name="");
    ~DiscreteSF();

    sFieldType getType() const { return DISCRETESFIELD; }
    virtual DiscreteSFType discretization() const = 0;
    pMesh getMesh() { return mesh; }

    // delete sizes
    virtual void cleanUp() = 0;
    void deleteSize(pEntity);

    // Intersect with another size field
    virtual void intersect(const pSField) = 0;

    // smooth the size field
    virtual void smooth(double) = 0;

    // set a size at all vertices 
    virtual void setCurrentSize() = 0;
    virtual void setCurvatureSize(bool aniso, 
                                  double alpha=2., // prescribe 2*PI*alpha edges around a circle
                                  double hMin=1.e-4) = 0; // minimal size
    virtual void setAllVSizes(pMSize) = 0;
    virtual void setAllVSizes(double[3][3], double[3]) = 0;
    virtual void setAllVSizes(double) = 0;
    virtual void scale(double) = 0;

    // set the size at a vertex 
    void setSize(pEntity, double[3][3], double[3]);
    void setSize(pEntity, pMSize);
    void setSize(pEntity, double);
    void setSize(int, double);

    // get the size at a location (allocate space!)
    virtual pMSize getSize(const pVertex) const = 0;
    virtual pMSize getSizeOnEntity(const pEntity, const double[3]) const = 0;

    // get the size at a vertex (do not allocate space)
    virtual const pMSize findSize(const pVertex) const = 0;
    virtual pMSize findSize(const pVertex) = 0;

    // edge length (squared)
    virtual double SF_VV_lengthSq(const pVertex, const pVertex) const = 0;
    virtual double SF_XYZ_lengthSq(const double[3], const double[3],
                                   const pMSize, const pMSize=NULL) const = 0;

    // face area (squared)
    virtual double SF_F_areaSq(const pFace) const = 0;
    virtual double SF_XYZ_areaSq(const double[3][3], const pMSize,
                                 const double[3]) const = 0;
    
    // region volume
    virtual double SF_R_volume(const pRegion) const = 0;
    virtual double SF_XYZ_volume(const double[4][3], const pMSize) const = 0;

    // center and its associated size
    virtual double SF_E_center(const pEdge, double[3], double *reducSq, pMSize *) const = 0;
    virtual double SF_VV_center(const pVertex, const pVertex,
                                double[3], double *reducSq, pMSize *) const = 0;

  protected:

    pMesh mesh;
    int dim;
    pMeshDataId pMSizeFieldId;

  };

}

// -------------------------------------------------------------------

#endif
