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

#include "MAdOperatorBase.h"
#include "MAdOutput.h"

using std::string;

namespace MAd {

  // -------------------------------------------------------------------
  bool MAdOperatorBase::evaluate(double * worstShape)
  {
    results->reset();
    *worstShape = mqm.getElementEvaluator()->worstShapeEver();

    History& history = HistorySgl::instance();

    // --- check 1: constraints ---
    if ( !checkConstraints() ) {
      history.add((int)type(),OP_CHECKCONSTRAINTS,0);
      return false;
    } else {
      history.add((int)type(),OP_CHECKCONSTRAINTS,1);
    }

    // --- check 2: geometric model ---
    if ( !checkGeometry() ) {
      history.add((int)type(),OP_CHECKGEOMETRY,0);
      return false;
    } else {
      history.add((int)type(),OP_CHECKGEOMETRY,1);
    }

    // --- check 3: elements validity and shapes ---
    if ( !evaluateShapes() ) {
      history.add((int)type(),OP_CHECKSHAPES,0);
      return false;
    } else {
      history.add((int)type(),OP_CHECKSHAPES,1);
    }

    // --- edge lengths ---
    evaluateLengths();
  
    *worstShape = results->getWorstShape();

    return true;
  }

  // -------------------------------------------------------------------
  void MAdOperatorBase::exportCavity(string filename) const
  {
    pPList cavity = PList_new();
    getCavity(&cavity);
    if ( dim == 3 ) {
      void * temp = NULL;
      pEntity pe;
      while ( ( pe = PList_next(cavity,&temp) ) ) {
        if ( EN_type(pe) != 3 ) continue;
        pRegion pr = (pRegion) pe;
        pPList rFaces = R_faces(pr);
        void * temp2 = NULL;
        pFace face;
        while ( ( face = (pFace)PList_next(rFaces,&temp2) ) ) {
          if ( F_whatInType(face) == 2 ) PList_appUnique(cavity,(pEntity)face);
        }
        PList_delete(rFaces);
      }
    }
    printPosEntities(cavity,filename,OD_DIMENSION,sizeField);
    PList_delete(cavity);
  }

  // -------------------------------------------------------------------

}
