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

#ifndef _H_MADELASTICITYOP
#define _H_MADELASTICITYOP

#include "MAdMatrix.h"
#include "VertexMoveOp.h"
#include "MAdLinearSystem.h"

#include <map>
#include <set>

/*
  This operator uses an analogy with an elastic medium in order to compute 
  the displacements of the vertices resulting from a prescribed displacement
  of the boundaries.

  The steps to make a computation are the following:
  
    1. Choose a set of nodal displacements by calling 'addDirichlet'
        --> build part of 'dirichlet'
    2. Build a cavity around the nodes for which you prescribed a displacement
        --> build 'cavity'
    3. Impose homogenous Dirichlet boudnary conditions to the boundaries 
       of the cavity with 'setDirichletBC()' (except the nodes with a 
       prescribed displacement)
        --> finish building 'dirichlet'
    4. Compute the displacements with 'compute'
        --> build 'relocations'
        --> clear 'cavity' and 'dirichlet'
        --> 'computed' set to 'true'
    5. Apply a part of the displacement to the mesh with 'advance' 
       (as much as you can)
    6. Once nodes are fully relocated, call 'clear'
        --> clear all
        --> computed set to 'false'
*/

namespace MAd {

  // -------------------------------------------------------------------
  class MAdElasticityOp {

  public:

    MAdElasticityOp(pMesh m);
    MAdElasticityOp(const MAdElasticityOp&);
    virtual ~MAdElasticityOp();

  public:

    void clear();

    int meshDim() const { return dim; }

    void setMaterials(double _E, double _nu);
    void setStiffnessAlterationCoef(double _chi);
    void setCavityEqualMesh(bool _m, int thickness)
    {
      cavityEqualMesh = _m;
      cavityThickness = thickness;
    }
    void setQualityThreshold (double qual) { qualityThreshold = qual; }

    void buildCavity();

    void setDirichletBC();
    void setHomogDirichletBC();
    void addDirichlet(pVertex pv, smallVector& dxyz);
    void delDirichlet(pVertex pv);

    int  compute();
    bool relocationsComputed() const { return computed; }
    void delRelocation(pVertex pv);

    // Advances the relocation as far as possible
    // Returns:
    //    - 0: no relocation possible
    //    - 1: advanced but not to the final position
    //    - 2: reached the full prescribed relocation
    int  advance(double * ratio, double tolerance=MAdTOL);
    void forceRelocation();

    void removeVertex(pVertex);
    void addVertexOnEdge(pVertex, pEdge);

    // debugging functions
    void printRelocations(std::ostream& out) const;
    void printDirichlet(std::ostream& out) const;

  private:

    int generateVertexIds();
    void collectBCVertices(std::set<pVertex> * bcVerts);
  
    void elementMatrix (pEntity pe, smallMatrix& m) const;
//     void element2DMatrix (pFace pr, smallMatrix& m) const;
//     void element3DMatrix (pRegion pr, smallMatrix& m) const;
    void addToMatrix(const pEntity, int, MAdLinearSystemDef *);
    void allocateMatrix(int, MAdLinearSystemDef *);

    void setVertexResult(int, pVertex, const MAdLinearSystemDef *);
    void setResults(int, const MAdLinearSystemDef *);

    bool checkArea(const pFace, double ratio);
    bool checkAreas(double ratio);
    bool checkVolume(const pRegion, double ratio);
    bool checkVolumes(double ratio);
    bool checkFQuality(const pFace pf, double ratio);
    bool checkRQuality(const pRegion pr, double ratio);
    bool checkQualities(double ratio);
    void relocate(double ratio);

  protected:

    pMesh mesh;
    int dim;

    double E, nu;
    
    // the nodes are relocated until this (low) quality is reached, could be 0.
    double qualityThreshold;

    // exponent of the alteration of stiffness based on elements size:
    // chi=0 => no alteration
    // chi increases => smaller elements stiffer
    double chi;

    bool computed;

    std::map<pVertex,int> localIds;

    bool cavityEqualMesh;
    int cavityThickness;
    std::set<pEntity> cavity;

    std::map<pVertex,smallVector > dirichlet;
    std::map<pVertex,smallVector > relocations;

  };

  // -------------------------------------------------------------------

}

#endif
