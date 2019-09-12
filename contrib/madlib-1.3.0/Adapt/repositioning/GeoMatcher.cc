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

#include "GeoMatcher.h"
#include "CallbackManager.h"
#include "MathUtils.h"
#include "ModelInterface.h"
#include "MAdMessage.h"
#include "AdaptInterface.h"

namespace MAd {

  // -------------------------------------------------------------------
  void GeoMatcherCBFunction (pPList before, pPList after, void * data,
                             operationType type , pEntity ppp) 
  {
#ifdef _HAVE_GMSH_
    geoMatcher * gm = (geoMatcher *)(data);

    switch (type) {
    case MAd_ESPLIT: { 

      if ( gm->relocationsComputed() ) {
        MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                    "Applying an edge split during a geometry snapping procedure");
      }

      // if new point is on a line or surface, you need to find its 
      // projection on the geometrical entity and impose its new 
      // location as a Dirichlet BC in the elastic model

      pGEntity pGE = V_whatIn( (pVertex)ppp );
      int dimBnd = GEN_type(pGE);
      if ( dimBnd < 3 ) {

        double xyz[3]; // current position of the new point
        V_coord( (pVertex)ppp, xyz );
//         printVec(xyz,"Original location");

        double u,v;
        if ( !V_params((pVertex)ppp,&u,&v) ) {
          V_info((pVertex)ppp);
          MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                      "Not a parametric point");
        }

        double tgt[3]; // the target location
        switch (dimBnd) {
        case 2: { // on a surface
          GF_xyz((pGFace)pGE, u, v, tgt);
          break;
        }
        case 1: { // on a line
          GE_xyz((pGEdge)pGE, u, tgt);
          break;
        }
        default: {
          MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                      "Unexpected geometric entity dimension for a new vertex: %d",
                                      dimBnd);
        }
        }
//         printVec(tgt,"Target location");

        // compute the difference with the current position
        double dx[3];
        diffVec(tgt,xyz,dx);
//         printVec(dx,"DISPLACEMENT");
      
        // impose the motion as a Dirichlet BC in the elastic model
        smallVector dxVec(3);
        for (int i=0; i<3; i++) dxVec(i) = dx[i];
        gm->addDirichlet( (pVertex)ppp, dxVec );
      }
      break;
    } 
    case MAd_ECOLLAPSE: {
      gm->delDirichlet( (pVertex)ppp );
      gm->delRelocation( (pVertex)ppp );
      break;
    }
    case MAd_FSWAP:
    case MAd_ESWAP: {
      break;
    }
    case MAd_RREMOVE: {
      throw; // some nodes could become parametric: to be checked
      void * temp = NULL;
      while ( pEntity pE = PList_next(before,&temp) ) {
        if ( EN_type(pE) == 0 ) {
          gm->delDirichlet( (pVertex)pE );
        }
      }
      break;
    }
    default: {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                  "Callback function not implemented for mesh modification %d",
                                  type);
    }
    }
#endif
  }

  // -------------------------------------------------------------------
  geoMatcher::geoMatcher(pMesh _mesh, MeshAdapter * _ma):
    MAdElasticityOp(_mesh), force(false), strictChecking(false),
    adapter(_ma), sliverFOp(NULL), sliverROp(NULL)
  {
    CallBackManagerSgl::instance().registerCallBack(GeoMatcherCBFunction,this);
  }

  // -------------------------------------------------------------------
  geoMatcher::geoMatcher(const geoMatcher& gm):
    MAdElasticityOp(gm), force(gm.force), strictChecking(gm.strictChecking),
    adapter(gm.adapter), sliverFOp(gm.sliverFOp), sliverROp(gm.sliverROp)
  {
    CallBackManagerSgl::instance().registerCallBack(GeoMatcherCBFunction,this);
  }

  // -------------------------------------------------------------------
  geoMatcher::~geoMatcher()
  {
    CallBackManagerSgl::instance().unregisterCallBack(GeoMatcherCBFunction,this);
  }
  
  // -------------------------------------------------------------------
  bool geoMatcher::snap()
  {
    bool res = true;

    // --- prepare the node motion ---
    buildCavity();
    setDirichletBC();
    compute();
    
    // --- perform the node motion ---
    if ( force ) forceRelocation();
    else
      {
        double ratio = 0.;
        int achieved = -1;
        while ( achieved != 2 )
          {
            achieved = advance(&ratio,1.e-12);
          
            if ( achieved == 0 ) {
              MAdMsgSgl::instance().warning(-1,__FILE__,
                                            "Could not advance snapping, achieved: %d, ratio: %f",
                                            achieved,ratio);
            }
            else {
              MAdMsgSgl::instance().info(-1,__FILE__,
                                         "Advance snapping, achieved: %d, ratio: %f",
                                         achieved,ratio);
            }
        
            if ( achieved <= 1 ) {

              // Apply the operations that do not add vertices or manage their 
              // target locations (complicated!)
              bool slivR, slivF;
              if ( sliverROp ) slivR = sliverROp->getNewVertexPermission();
              if ( sliverFOp ) slivF = sliverFOp->getNewVertexPermission();
              if ( sliverROp ) sliverROp->newVertexPermission(false);
              if ( sliverFOp ) sliverFOp->newVertexPermission(false);

              int nbBefore, nbAfter;
              if (M_dim(mesh) == 3) sliverROp->removeSliverRegions(&nbBefore,&nbAfter);
              else                  sliverFOp->removeSliverFaces  (&nbBefore,&nbAfter);
              std::cout << "Removed slivers " << nbBefore << " -> " << nbAfter <<"\n";

              if ( !( nbBefore - nbAfter ) )
                {
                  if ( !adapter->optimiseElementShape() ) {
                    //                 string filename = outPrefix + "dirichlet.txt";
                    //                 std::ofstream dirichletOut(filename.c_str());
                    //                 geoTracker->printDirichlet(dirichletOut);
                    //                 dirichletOut.close();
                    //                 filename = outPrefix + "relocations.txt";
                    //                 std::ofstream relocationsOut(filename.c_str());
                    //                 printRelocations(relocationsOut);
                    //                 relocationsOut.close();
                    MAdMsgSgl::instance().warning(__LINE__,__FILE__,
                                                  "Could not snap vertices, ratio reached: %f",
                                                  ratio);
                    reportFailure(ratio);
                    break;
                  }
                }
              if ( sliverROp ) sliverROp->newVertexPermission( slivR );
              if ( sliverFOp ) sliverFOp->newVertexPermission( slivF );
            }
          }

        if ( strictChecking && achieved != 2 ) { res = false; }
      }
    
    clear();

    return res;
  }
  
  // -------------------------------------------------------------------
  void geoMatcher::printFailures(std::ostream& out) const
  {
    if ( failures.size() ) {
      out << "No failure reported in the GeoMatcher\n";
    }
    else {
      out << failures.size() << " failures reported in the GeoMatcher:\n\n";
      int i = 0;
      std::list<double>::const_iterator it = failures.begin();
      for (; it != failures.end(); it++)
        {
          out << i << "  ratio: " << *it << "\n";
          i++;
        }
    }
  }

  // -------------------------------------------------------------------

}
