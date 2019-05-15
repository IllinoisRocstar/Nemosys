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

#ifndef _H_MADOPERATORBASE
#define _H_MADOPERATORBASE

#include "MAdOperations.h"
#include "DiscreteSF.h"
#include "ElementStatistics.h"
#include "MeshQualityManager.h"
#include "Constraint.h"
#include "History.h"
#include "MAdDefines.h"

#include <string>

namespace MAd {

  // -------------------------------------------------------------------
  typedef class MAdOperatorBase * pMAdOperator;

  // -------------------------------------------------------------------
  class MAdOperatorBase
  {
  public:

    MAdOperatorBase();
    MAdOperatorBase(const MAdOperatorBase &);
    MAdOperatorBase(pMesh, DiscreteSF *);
    virtual ~MAdOperatorBase();

  public:

    virtual operationType type() const = 0;

    void setSizeField(DiscreteSF *);

    double getWorstShape() const;
    double getMinLenSq()   const;
    double getMaxLenSq()   const;

    // get a list of all elements that will be modified
    virtual void getCavity(pPList *) const = 0;
    void exportCavity(std::string) const;

    // check if the operation can be performed and evaluate it
    bool evaluate(double *);

    // apply the operation
    virtual void apply() = 0;
  
  private:

    // --- Checks and evaluations ---
    // ! Supposed to be called in that order !

    // checks the operation regarding the constraints on mesh entities 
    // and geometric entities
    virtual bool checkConstraints() const = 0;

    // checks compatibility of the operation with the geometric model
    virtual bool checkGeometry()   = 0;

    // checks validity of the resulting elements and evaluates their 
    // shapes (saves the worst)
    virtual bool evaluateShapes() = 0;

    // evaluates the resulting minimal and maximal edge lengths
    virtual void evaluateLengths() const = 0;

    // ------------------------------

  protected:

    pMesh mesh;
    DiscreteSF * sizeField;

    // quality evaluator
    MeshQualityManager& mqm;

    // storage for the results of the evaluation
    mutable ElementStatistics * results;

    // mesh dimension
    int dim;

  };

  // -------------------------------------------------------------------
  inline MAdOperatorBase::MAdOperatorBase():
    mesh(NULL), sizeField(NULL), mqm(MeshQualityManagerSgl::instance()), 
    results(NULL), dim(0)
  {}

  // -------------------------------------------------------------------
  inline MAdOperatorBase::MAdOperatorBase(const MAdOperatorBase & _op):
    mesh(_op.mesh), sizeField(_op.sizeField),
    mqm(MeshQualityManagerSgl::instance())
  {
    results = new ElementStatistics(*(_op.results));
    dim = _op.dim;
  } 

  // -------------------------------------------------------------------
  inline MAdOperatorBase::MAdOperatorBase(pMesh _mesh, DiscreteSF * _sf):
    mesh(_mesh), sizeField(_sf), mqm(MeshQualityManagerSgl::instance())
  { 
    results = new ElementStatistics(); 
    dim = M_dim(mesh);
  }

  // -------------------------------------------------------------------
  inline MAdOperatorBase::~MAdOperatorBase() 
  {
    if (results) delete results;
  }

  // -------------------------------------------------------------------
  inline void MAdOperatorBase::setSizeField(DiscreteSF * _sf)
  {
    sizeField = _sf;
  }

  // -------------------------------------------------------------------
  inline double MAdOperatorBase::getWorstShape() const
  { return results->getWorstShape(); }

  // -------------------------------------------------------------------
  inline double MAdOperatorBase::getMinLenSq() const
  { return results->getMinLenSq(); }

  // -------------------------------------------------------------------
  inline double MAdOperatorBase::getMaxLenSq() const
  { return results->getMaxLenSq(); }

  // -------------------------------------------------------------------

}

#endif
