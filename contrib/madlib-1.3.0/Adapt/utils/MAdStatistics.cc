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

#include "MAdStatistics.h"

using std::ostream;

namespace MAd {

  // -------------------------------------------------------------------
  void MAdStatistics::initialize()
  {
    t_eSplits    = 0.;
    t_eCollapses = 0.;
    t_eSwaps     = 0.;
    t_rSlivers   = 0.;
    t_fSlivers   = 0.;

    num_eSplits    = 0;
    num_eCollapses = 0;
    num_eSwaps     = 0;

    numInfLoops = 0;
  }

  // -------------------------------------------------------------------
  void MAdStatistics::finalize()
  {
  }

  // -------------------------------------------------------------------
  void MAdStatistics::print(ostream& out) const
  {
    out << "\n*** Statistics about local mesh modifications *** \n\n";

    out << "Time spent in the different operators:\n";
    double t_Total = t_eSplits + t_eCollapses + t_eSwaps + t_rSlivers + t_fSlivers;
    out << "Edge splits\tEdge collapses\tEdge swaps\tSlivers\tTotal\n"
        << t_eSplits <<"\t"<< t_eCollapses <<"\t"<< t_eSwaps <<"\t"
        << t_rSlivers+t_fSlivers <<"\t"<< t_Total <<"\n\n";

    out << "Number of elementary operations (excluding sliver handling):\n";
    out << "Edge splits\tEdge collapses\tEdge swaps\n"
        << num_eSplits <<"\t"<< num_eCollapses <<"\t"<< num_eSwaps << "\n\n";

    out << "Number of infinite loops:\t" << numInfLoops << "\n\n";
  }

  // -------------------------------------------------------------------

}
