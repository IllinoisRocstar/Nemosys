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

#ifndef _H_MOBILEOBJECT
#define _H_MOBILEOBJECT

#include "MeshDataBaseInterface.h"
#include "AdaptInterface.h"
#include "LocalSizeField.h"

#include <set>
#include <list>
#include <utility>
#include <vector>
#include <string>

namespace MAd {

  class MAdElasticityOp;
  class MAdStringFieldEvaluator;

  // -------------------------------------------------------------------
  // GCTODO: rewrite it simpler
  struct vDisplacement
  {
    vDisplacement(const vDisplacement&);
    vDisplacement(pVertex, double[3]);
    void scale(double);
    pVertex pv;
    double dxyz[3];    // displacement
  };

  struct vDisplacementLess
  {
    bool operator() (const vDisplacement&, const vDisplacement&) const;
  };

  // ----------------------------------------------------------------------
  enum velocityFormulation {
    NOFORMULATION,
    PARSED,
    TRANSLATION,
    FULLRANDOM,
    RANDOM_TRANSLATION,
    ISOMETRY,
    ROTATION,
    SHEAR_YZ,
    SHEAR_ZX,
    SHEAR_XY,
    BENCH2,
    CROSS_EXPANSION_WITH_ROTATION
  };

  // ----------------------------------------------------------------------
  class mobileObject {

  public:

    mobileObject(const pMesh m,std::string _name="");
    mobileObject(const mobileObject & mob);
    ~mobileObject();

    void setName(std::string nm) {name = nm;}
    void setDxKinematics(std::vector<std::string> _Vstr);
    void setVKinematics(velocityFormulation type, double V[3], 
                        double C[3], std::vector<std::string> _Vstr);
    void addLocalSField(LocalSizeField* lsf);
    void addGEntity(int type, int tag);
    void reAddVertices();

    void computePrescribedDisplacement (double t, double dt);

    void describe (std::ostream& out=std::cout) const;
    const std::set<pVertex> getVertices() const {return vertices;}
    std::set<vDisplacement,vDisplacementLess> getPrescribedDisplacement () const {return prescribedDisplacement;}
    std::set<LocalSizeField* > getSizes() const { return sizes; }

  private:

    void addVerticesOnGEntity(int type, int tag);
    void clearDisplacement();

    void randomVelocity(double* V);
    void velocityFunction(double xyz[3], double t, double* V);

  private:

    const pMesh mesh;
    std::string name;
    std::list<std::pair<int,int> > geomEntities; // list of (type,tag) 
    std::set<pVertex> vertices;

    // Parameters describing the kinematics
    std::string prescribeType;
    // for position ...
    MAdStringFieldEvaluator* DxParsed; // displacement relative to and evaluated on initial position
    // ... or velocity
    velocityFormulation velType;
    double Cxyz[3], Vxyz[3];
    MAdStringFieldEvaluator* VParsed; // velocity evaluated on current position

    //   // Parameters describing the size around
    std::set<LocalSizeField* > sizes;

    std::set<vDisplacement,vDisplacementLess> prescribedDisplacement;
  };

  // ----------------------------------------------------------------------
  class mobileObjectSet {

  public:

    mobileObjectSet();
    ~mobileObjectSet();

    // larger possible motion by trials without volume nodes repositioning
    int partlyMove(vertexMoveOp& vMoveOp, double t, double dt, double * part);

    // move objects and reposition volume nodes. Needs the elastic operator. 
    void setupElasticRepositioning(pMesh mesh, double t, double dt,
                                   double qualThr=0.,double chi=-1,
                                   bool meshIsCavity=true, int cavityThickness=3);

    // Advances the relocation as far as possible
    // Returns:
    //    - 0: no relocation possible
    //    - 1: advanced but not to the final position
    //    - 2: reached the full prescribed relocation
    int reposition(double * ratio);

    bool moveAndReposition(pMesh mesh, double t, double dt, bool subAdaptation,
                           double qualityThreshold, double chi, bool meshIsCavity,
                           int cavityThickness, MeshAdapter * ma);

    // empty or fill the set
    void clear();
    void insert(mobileObject*);

    std::set<vDisplacement,vDisplacementLess> getPrescribedDisplacement () const {return prescribedDisplacement;}
    std::set<mobileObject*> getObjects() const {return mobSet;}

    void describe(std::ostream& out=std::cout) const;
 
  private:

    // compute the objects displacement in a time interval
    void computePrescribedDisplacement (double t, double dt) ;

  private:

    std::set<mobileObject*> mobSet;
    std::set<vDisplacement,vDisplacementLess> prescribedDisplacement;

    MAdElasticityOp * elasticOp;
  };

  // ----------------------------------------------------------------------

}

#endif
