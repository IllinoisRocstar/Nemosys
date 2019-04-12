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

#ifndef _H_SIZEFIELDBASE
#define _H_SIZEFIELDBASE

#include "MeshDataBaseInterface.h"

#include <string>

namespace MAd {

  // -------------------------------------------------------------------
  enum sFieldType {
    UNKNOWNSFIELDTYPE,
    NULLSFIELD,
    DISCRETESFIELD,
    ANALYTICALSFIELD,
    LOCALSFIELD,
    BACKGROUNDSFIELD
  };

  // -------------------------------------------------------------------
  typedef class SizeFieldBase * pSField;
  typedef class MeshSizeBase  * pMSize;

  // -------------------------------------------------------------------
  class SizeFieldBase 
  {
  public:
  
    SizeFieldBase(std::string _name=""): name(_name) {};
    virtual ~SizeFieldBase() {};
  
  public:
  
    virtual sFieldType getType() const = 0 ;
    std::string getName() const { return name; }
  
    virtual void scale(double) = 0;
  
  public:

    // edge length (squared)
    virtual double SF_E_lengthSq(const pEdge) const;
    virtual double SF_VV_lengthSq(const pVertex,const pVertex) const = 0;
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
    virtual double SF_E_center(const pEdge, double[3], double * reducSq, pMSize *) const = 0;
    virtual double SF_VV_center(const pVertex, const pVertex,
                                double[3], double * reducSq, pMSize *) const = 0;
    virtual double SF_XYZ_center(const double[2][3], const pMSize[2],
                                 double[3], double * reducSq, pMSize *) const;

  public:

    // get the size at a location
    // --- Memory is allocated ! ---
    virtual pMSize getSize(const pVertex) const = 0;
    virtual pMSize getSizeOnEntity(const pEntity, 
                                   const double[3]) const = 0;

  public:

    // visualisation of the isotropic size in a .pos file
    void printPosIsotropic  (const pMesh mesh, const std::string name);
    // visualisation of the anisotropic size in a .pos file
    void printPosAnisotropic(const pMesh mesh, const std::string baseName);

  private:

    std::string name;

  };

}

// -------------------------------------------------------------------

#endif
