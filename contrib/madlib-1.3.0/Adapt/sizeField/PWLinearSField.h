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

#ifndef _H_PWLINEARSFIELD
#define _H_PWLINEARSFIELD

#include "madlib_export.h"

#include "DiscreteSF.h"
#ifdef PARALLEL
#include "MeshDataBaseComm.h"
#endif

#include <set>

namespace MAd {

  // -------------------------------------------------------------------
  class MADLIB_EXPORT PWLSField : public DiscreteSF
  {
  public:

    PWLSField(pMesh, std::string name="");
    ~PWLSField();

    DiscreteSFType discretization() const { return VERTEX_P1_DSFTYPE; }

    // delete sizes
    void cleanUp();

    // Intersect with another size field
    void intersect(const pSField);

    // smooth the size field
    void smooth(double);

    // set a size at all vertices 
    void setCurrentSize();
    void setCurvatureSize(bool aniso, 
                          double alpha=2., // prescribe 2*PI*alpha edges around a circle
                          double hMin=1.e-4); // minimal size
    void setAllVSizes(pMSize);
    void setAllVSizes(double[3][3], double[3]);
    void setAllVSizes(double);
    void scale(double);

    // get the size at a location (allocate space!)
    pMSize getSize(const pVertex) const;
    pMSize getSizeOnEntity(const pEntity, const double[3]) const;

    // get the size at a vertex (do not allocate space)
    const pMSize findSize(const pVertex) const;
    pMSize findSize(const pVertex);

    // edge length (squared)
    double SF_VV_lengthSq(const pVertex, const pVertex) const;
    double SF_XYZ_lengthSq(const double[3], const double[3],
                           const pMSize, const pMSize=NULL) const;

    // face area (squared)
    double SF_F_areaSq(const pFace) const;
    double SF_XYZ_areaSq(const double[3][3], const pMSize,
                         const double[3]) const;

    // region volume
    double SF_R_volume(const pRegion) const;
    double SF_XYZ_volume(const double[4][3], const pMSize) const;

    // center and its associated size
    double SF_E_center(const pEdge, double[3], double * reducSq, pMSize *) const;
    double SF_VV_center(const pVertex, const pVertex,
                        double[3], double * reducSq, pMSize *) const;

  protected:

    pMSize getSizeOnEdge(const pEdge, const double[3]) const;
    pMSize getSizeOnEdgeParam(const pEdge, const double) const;
    pMSize getSizeOnFace(const pFace, const double[3]) const;
    pMSize getSizeOnFaceParam(const pFace, const double[2]) const;
    pMSize getSizeOnRegion(const pRegion, const double[3]) const;
    pMSize getSizeOnRegionParam(const pRegion, const double[3]) const;

    void smoothOnEdge(const pEdge, double, std::set<pEdge>*);
  };

  // -------------------------------------------------------------------
#ifdef PARALLEL
  class MADLIB_EXPORT PWLSFieldDE : public MDB_DataExchanger
  {
  public :
    
    PWLSField *field;
    
    PWLSFieldDE(PWLSField *f);
    virtual ~PWLSFieldDE();
    
    virtual void * sendData (pEntity pe,    // in
                             int iProcDest, // in
                             int &_size );
    
    virtual void receiveData (pEntity pe,      //in
                              int iProcSender, //in
                              void *buf );
    
    virtual void deleteExternalData(pEntity pe) const;
  };
#endif

}

// -------------------------------------------------------------------

#endif
