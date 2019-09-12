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

#ifndef _H_MESHQUALITYMANAGER
#define _H_MESHQUALITYMANAGER

#include "SizeFieldBase.h"
#include "ElementEvaluatorBase.h"
#include "MeshDataBaseInterface.h"
#include "MAdSingleton.h"

#include <iostream>

namespace MAd {

  // -------------------------------------------------------------------
  class MeshQualityManager {

  public:
  
    MeshQualityManager():evaluator(NULL), histogram(NULL), histogramAvg(NULL){};
    ~MeshQualityManager() {};
  
    void initialize(pMesh m, DiscreteSF * sf, evaluationType type);
    void setMesh(pMesh m);
    void finalize();

    // -------------------------------------------------------------------

  public:

    const pMeshDataId getShapeId() const { return shapeId; }
    const elementEvaluatorBase * getElementEvaluator() const { return evaluator; }

    // evaluate missing shapes
    void evaluateAllShapes() const;

    // evaluate all shapes, replace existing ones
    void evaluateAndReplaceAllShapes() const;

    // compute and return the shape of the face (no storage for faces)
    // return 0 if the face is not acceptable, 1 otherwise
    int getShape(pFace pf, double normal[3], double * result) const;
    int getShapeWithDisp(pFace pf, double normal[3], 
                         double disp[3][3], double * result) const;

    // get the shape of the region, compute and store the shape if missing
    // return 0 if the region is not acceptable, 1 otherwise
    int getShape(pRegion pr, double * result) const;
    int getShapeWithDisp(pRegion pr, double disp[4][3], double * result) const;

    // delete all shapes
    void clearAllShapes() const;

    // delete the shape of this entity
    void clearShape(pEntity pe) const;

    // delete all shapes of the elements neighbouring the entity
    void clearNeighbourShapes(pVertex pv) const;
    void clearNeighbourShapes(pEdge pe) const;

    // evaluate the worst shape of all elements surrounding one or several entities
    int V_worstShape(pVertex v, double* result) const;
    int E_worstShape(pEdge e,   double* result) const;
    int F_worstShape(pFace f,   double* result) const;

    // evaluate the worst shape in a list of elements
    int FList_worstShape(pPList faces,   double * result) const;
    int RList_worstShape(pPList regions, double * result) const;

    // check that the attached shapes are still correct (debugging function)
    bool checkAttachedShapes() const;

  private:

    // attach a computed shape to the entity
    void attachShape(pEntity pe, double shape) const;

    // get the shape of an entity if it exists, return false otherwise
    bool getAttachedShape(const pEntity pe, double* result) const;

  public:

    // Evaluate statistics about elements sizes
    void evaluateSizes();

    // Evaluate statistics about elements shapes
    void evaluateShapes();

    // Evaluate all statistics
    void evaluateStatistics();

    // get statistics ( evaluate it first )
    double getMeanShape()  const { return meanShape; }
    double getWorstShape() const { return worstShape; }
    double getMinSize()    const { return minAbsoluteSize; }
    double getMaxSize()    const { return maxAbsoluteSize; }
    void printStatistics(std::ostream& out) const;
  
    // -------------------------------------------------------------------

  private:

    pMesh mesh;
    int dim;
    DiscreteSF * sizeField;

    elementEvaluatorBase* evaluator;

    pMeshDataId shapeId;

    // --- statistics ---
    // ------------------
    int nbElem;
    double worstShape;
    double meanShape;
    double minAbsoluteSize; // Area or volume of the smallest element in usual space (no metric)
    double maxAbsoluteSize; // Area or volume of the biggest  element in usual space (no metric)
    double sizesSum;
    int notAcpt;  // number of non acceptable elements
    int* histogram; // counters of elements with different qualities
    double* histogramAvg; // proportion of elements with different qualities

  };

  // -------------------------------------------------------------------

  typedef MAdSingleton<MeshQualityManager> MeshQualityManagerSgl;

}

#endif
