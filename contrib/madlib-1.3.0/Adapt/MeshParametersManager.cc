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

#include "MeshParametersManager.h"
#include "MAdMessage.h"

// standart C/C++
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

namespace MAd {

  // -------------------------------------------------------------------
  MeshParametersManager::MeshParametersManager():
    lowerLengthSqBound(-1.), upperLengthSqBound(-1.), 
    sliverTetQualityBound(-1.), sliverTriQualityBound(-1.), 
    sliverPermInESplit(false), sliverPermInECollapse(false),
    sliverLowerLengthSqBound(-1.),sliverUpperLengthSqBound(-1.),
    noSwapQuality(-1.),swapMinImproveRatio(-1.), bigLength(1.e14)
  {}

  // -------------------------------------------------------------------
  void MeshParametersManager::initialize()
  {
  }

  // -------------------------------------------------------------------
  void MeshParametersManager::finalize()
  {
  }

  // -------------------------------------------------------------------
  // -------------------------------------------------------------------

  void MeshParametersManager::setLowerLengthSqBound(double lSq)
  {
    lowerLengthSqBound = lSq;
  }

  // -------------------------------------------------------------------

  void MeshParametersManager::setUpperLengthSqBound(double lSq)
  {
    upperLengthSqBound = lSq;
  }

  // -------------------------------------------------------------------
  double MeshParametersManager::getSliverBound(int dim) const
  {
    if (dim==3) return getSliverTetBound();
    return getSliverTriBound();
  }

  // -------------------------------------------------------------------
  double MeshParametersManager::getSliverBound(const pMesh mesh) const
  {
    return getSliverBound(M_dim(mesh));
  }

  // -------------------------------------------------------------------

  void MeshParametersManager::setSliverTriBound(double bound)
  {
    sliverTriQualityBound = bound;
  }

  // -------------------------------------------------------------------

  void MeshParametersManager::setSliverTetBound(double bound)
  {
    sliverTetQualityBound = bound;
  }

  // -------------------------------------------------------------------
  // -------------------------------------------------------------------

  void MeshParametersManager::setSliverPermissionInESplit(bool perm, 
                                                          double bound)
  {
    sliverPermInESplit = perm;
    sliverUpperLengthSqBound = bound;
    if ( sliverUpperLengthSqBound > 0. && 
         sliverUpperLengthSqBound < upperLengthSqBound ) {
      MAdMsgSgl::instance().warning(__LINE__,__FILE__,
                                    "Incoherent parameters: upper length square bound for edge splits creating slivers (%f) is lower than the global upper legnth square bound (%f)",
                                    sliverUpperLengthSqBound,upperLengthSqBound);
    }
  }

  // -------------------------------------------------------------------

  void MeshParametersManager::setSliverPermissionInECollapse(bool perm, 
                                                             double bound)
  {
    sliverPermInECollapse = perm;
    sliverLowerLengthSqBound = bound;
    if ( sliverLowerLengthSqBound > 0. && 
         sliverLowerLengthSqBound > lowerLengthSqBound ) {
      MAdMsgSgl::instance().warning(__LINE__,__FILE__,
                                    "Incoherent parameters: lower length square bound for edge collapses creating slivers (%f) is bigger than the global lower legnth square bound (%f)",
                                    sliverLowerLengthSqBound,lowerLengthSqBound);
    }
  }

  // -------------------------------------------------------------------
  // -------------------------------------------------------------------

  void MeshParametersManager::diagnostics() const
  {
    cout << "   *** Parameters for the mesh adaptation ***   \n";

    cout << "\n";

    cout << "- Bounds of the edge length squared: \n";
    cout << "\t Lower: " << lowerLengthSqBound << "\n";
    cout << "\t Upper: " << upperLengthSqBound << "\n";

    cout << "\n";

    cout << "- Quality bound for the sliver triangles: " 
         << sliverTriQualityBound << "\n";

    cout << "- Quality bound for the sliver tetrahedra: " 
         << sliverTetQualityBound << "\n";

    cout << "\n";

    cout << "Edge splits can create slivers: ";
    if (sliverPermInESplit) {
      cout << "yes: ";
      if ( sliverUpperLengthSqBound < 0) cout << "for any edge.";
      else cout << "if the edge is longer than "<< sliverUpperLengthSqBound<<" (square of adimensional length)";
    }
    else cout << "no";
    cout << "\n";
    cout << "Edge collapses can create slivers: ";
    if (sliverPermInECollapse) {
      cout << "yes: ";
      if ( sliverLowerLengthSqBound < 0) cout << "for any edge.";
      else cout << "the edge is shorter than "<< sliverLowerLengthSqBound<<" (square of adimensional length)";
    }
    else cout << "no";
    cout << "\n";

    cout << "\n";
  }

  // -------------------------------------------------------------------

}
