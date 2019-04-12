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

#include "MobileObject.h"
#include "NodalDataManager.h"
#include "VertexMoveOp.h"
#include "MAdElasticityOp.h"
#include "DistanceFunction.h"

#include "MAdStringFieldEvaluator.h"

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
#include <string.h>
using std::set;
//using std::multiset;
using std::list;
using std::pair;
using std::vector;
using std::string;
#include <stdlib.h>

namespace MAd {

  // ----------------------------------------------------------------------
  mobileObject::mobileObject(const pMesh m, string _name):
    mesh(m),name(_name),prescribeType(""),DxParsed(NULL),
    velType(NOFORMULATION),VParsed(NULL)//, Cxyz(NULL), Vxyz(NULL)
  {
#ifdef PARALLEL
    MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                "Mobile objects not supported in parallel");
#endif
  }

  // ----------------------------------------------------------------------
  mobileObject::~mobileObject()
  {
    if (DxParsed) { delete DxParsed; DxParsed = NULL; }
    if (VParsed)  { delete VParsed;  VParsed  = NULL; }
  }

  // ----------------------------------------------------------------------
  void mobileObject::setDxKinematics(vector<string> _str)
  {
    if ( strcmp(prescribeType.c_str(),"") ) {
      cerr << "Error: imposing a position on an object with another imposed movement: "<<prescribeType<<"\n";
      throw;
    }
    prescribeType = "Displacement";

    if (DxParsed) delete DxParsed;
    DxParsed = new MAdStringFieldEvaluator(_str);
  }

  // ----------------------------------------------------------------------
  void mobileObject::setVKinematics(velocityFormulation type, double V[3], 
                                    double C[3], vector<string> _str)
  {
    if ( strcmp(prescribeType.c_str(),"") ) {
      cerr << "Error: imposing a velocity on an object with another imposed movement: "<<prescribeType<<"\n";
      throw;
    }
    prescribeType = "Velocity";
      
    velType = type; 
    if ( V ) {
      for (int i=0; i<3; i++)  Vxyz[i] = V[i];
    }
    if ( C ) {
      for (int i=0; i<3; i++)  Cxyz[i] = C[i];
    }
    if (type == PARSED)  VParsed = new MAdStringFieldEvaluator(_str);
  }

  // ----------------------------------------------------------------------
  void mobileObject::addLocalSField(LocalSizeField* lsf)
  {
    sizes.insert(lsf);
  }

  // ----------------------------------------------------------------------
  void mobileObject::addGEntity(int type, int tag)
  {
    geomEntities.push_back(std::make_pair(type,tag));
    addVerticesOnGEntity(type,tag);
  }

  // ----------------------------------------------------------------------
  void mobileObject::reAddVertices()
  {
    vertices.clear();
    list<pair<int,int> >::const_iterator it = geomEntities.begin();
    list<pair<int,int> >::const_iterator itEnd = geomEntities.end();
    for (; it != itEnd; it++)  addVerticesOnGEntity(it->first, it->second);
  }

  // ----------------------------------------------------------------------
  void mobileObject::addVerticesOnGEntity(int type, int tag)
  {
    switch (type) {
    case 0:
      {
        VIter vit = M_vertexIter(mesh);
        while (pVertex pv = VIter_next(vit))
          {
            pGEntity pg = EN_whatIn((pEntity)pv);
            int pgType = EN_whatInType( (pEntity) pv );
            if(pg)
              if (GEN_tag(pg) == tag && pgType == 0)
                {
                  vertices.insert(pv);
                }
          }
        VIter_delete(vit);
        break;
      }
    case 1:
      {
        EIter eit = M_edgeIter(mesh);
        while (pEdge pe = EIter_next(eit))
          {
            pGEntity pg = EN_whatIn((pEntity)pe);
            int pgType = EN_whatInType( (pEntity) pe );
            if(pg)
              if (GEN_tag(pg) == tag && pgType == 1)
                {
                  vertices.insert(E_vertex (pe,0));
                  vertices.insert(E_vertex (pe,1));
                }
          }
        EIter_delete(eit);
      }
    case 2:
      {
        FIter fit = M_faceIter(mesh);
        while (pFace pf = FIter_next(fit))
          {
            pGEntity pg = EN_whatIn((pEntity)pf);
            int pgType = EN_whatInType( (pEntity) pf );
            if(pg)
              if (GEN_tag(pg) == tag && pgType == 2)
                {
                  vertices.insert(F_vertex (pf,0));
                  vertices.insert(F_vertex (pf,1));
                  vertices.insert(F_vertex (pf,2));
                }
          }
        FIter_delete(fit);
      }
    case 3:
      {
        RIter rit = M_regionIter(mesh);
        while (pRegion pr = RIter_next(rit))
          {
            pGEntity pg = EN_whatIn((pEntity)pr);
            int pgType = EN_whatInType( (pEntity) pr );
            if(pg)
              if (GEN_tag(pg) == tag && pgType == 3)
                {
                  vertices.insert(R_vertex (pr,0));
                  vertices.insert(R_vertex (pr,1));
                  vertices.insert(R_vertex (pr,2));
                  vertices.insert(R_vertex (pr,3));
                }
          }
        RIter_delete(rit);
      }
    }
  }

  // ----------------------------------------------------------------------
  void mobileObject::computePrescribedDisplacement (double t, double dt)
  {
    if ( !(NodalDataManagerSgl::instance().isCoordinates()) ) {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                  "Allocate initial coordinates before making any move (use AdaptInterface::storeInitialCoordinates())");
    }
  
    prescribedDisplacement.clear();
    reAddVertices();
  
    if ( !strcmp(prescribeType.c_str(),"Velocity") ) 
      {
        double V[3]; V[0]=0.; V[1]=0.; V[2]=0.; 
        if (velType == RANDOM_TRANSLATION) 
          randomVelocity(V);
        set<pVertex> :: iterator itV    = vertices.begin();
        set<pVertex> :: iterator itVEnd = vertices.end();
        for ( ; itV != itVEnd ; ++itV)
          {
            double dx[3];
            double xyz[3]; V_coord (*itV,xyz);
            if ( velType != RANDOM_TRANSLATION)
              velocityFunction(xyz,t,V);
            dx[0] = V[0] * dt;
            dx[1] = V[1] * dt;
            dx[2] = V[2] * dt;
          
            vDisplacement disp(*itV,dx);
            prescribedDisplacement.insert(disp);
          }
      }

    else if ( !strcmp(prescribeType.c_str(),"Displacement") ) 
      {
        set<pVertex> :: iterator itV    = vertices.begin();
        set<pVertex> :: iterator itVEnd = vertices.end();
        for ( ; itV != itVEnd ; ++itV)
          {
            // get initial position
            vector<double> xyz0;
            NodalDataManagerSgl::instance().getStoredCoordinates(*itV,xyz0);
          
            // get prescribed displacement (relative to and evaluated on initial position)
            double dx0[3];
            DxParsed->eval(xyz0,t,dx0);
          
            // get current position
            double xyz[3]; V_coord (*itV,xyz);

            // get new position
            double xyznew[3];
            for (int i=0; i<3; i++)  xyznew[i] = xyz0[i] + dx0[i];

            // compute the displacement
            double dx[3];
            dx[0] = xyznew[0] - xyz[0];
            dx[1] = xyznew[1] - xyz[1];
            dx[2] = xyznew[2] - xyz[2];
          
            vDisplacement disp(*itV,dx);
            prescribedDisplacement.insert(disp);
          }
      }

    else {
      cerr<< "Error: no kinematics specified for this object\n";
      throw;
    }
  }

  // ----------------------------------------------------------------------
  void mobileObject::describe (std::ostream& out) const 
  {
    out << "Object \'" << name.c_str() << "\' carries "
        << vertices.size() << " nodes with xyz velocity:\n"
        << "  *  Vx: "<<Vxyz[0]<<"\n"
        << "  *  Vy: "<<Vxyz[1]<<"\n"
        << "  *  Vz: "<<Vxyz[2]<<"\n"
        << " and analytical velocity formulation "<<velType<<"\n\n";
  }

  // ----------------------------------------------------------------------
  // generates velocities between -1 and 1
  void mobileObject::randomVelocity(double* V)
  {
    V[0] = (  (  (double)rand() / ((double)(RAND_MAX)+(double)(1))  ) - 0.5  ) * 2.;
    V[1] = (  (  (double)rand() / ((double)(RAND_MAX)+(double)(1))  ) - 0.5  ) * 2.;
    V[2] = (  (  (double)rand() / ((double)(RAND_MAX)+(double)(1))  ) - 0.5  ) * 2.;
  }

  // ----------------------------------------------------------------------
  void mobileObject::velocityFunction(double* xyz, double t, double* V)
  {
    double x = xyz[0];
    double y = xyz[1];
    double z = xyz[2];

    switch (velType) 
      {
      case NOFORMULATION:
        {
          MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                      "Velocity formulation: NOFORMULATION");
        }
      case PARSED:
        {
          VParsed->eval(xyz,t,V);
          break;
        }
      case TRANSLATION:
        {
          V[0] = Vxyz[0]; V[1] = Vxyz[1]; V[2] = Vxyz[2];
          break;
        }
      case FULLRANDOM:
        {
          randomVelocity(V);
          break;
        }
      case RANDOM_TRANSLATION:
        {
          MAdMsgSgl::instance().warning(__LINE__,__FILE__,
                                        "VelocityFunction with tag RANDOM_TRANSLATION should not be used");
          randomVelocity(V);
          break;
        }
      case ISOMETRY: // expansion 
        {
          V[0] = Vxyz[0] * ( x - Cxyz[0] );
          V[1] = Vxyz[1] * ( y - Cxyz[1] );
          V[2] = Vxyz[2] * ( z - Cxyz[2] );
          break;
        }
      case ROTATION: // rotation
        { 
          double velYZ = Vxyz[0];
          double velZX = Vxyz[1];
          double velXY = Vxyz[2];
          double rYZ = sqrt ( (y - Cxyz[1])*(y - Cxyz[1]) + (z - Cxyz[2])*(z - Cxyz[2]) );
          double rZX = sqrt ( (x - Cxyz[0])*(x - Cxyz[0]) + (z - Cxyz[2])*(z - Cxyz[2]) );
          double rXY = sqrt ( (x - Cxyz[0])*(x - Cxyz[0]) + (y - Cxyz[1])*(y - Cxyz[1]) );
          double angleYZ = atan2((y - Cxyz[1]), (z - Cxyz[2]));
          double angleZX = atan2((z - Cxyz[2]), (x - Cxyz[0]));
          double angleXY = atan2((x - Cxyz[0]), (y - Cxyz[1]));
          V[0] = 0                               - velZX * rZX * sin(angleZX) + velXY * rXY * cos(angleXY) ;
          V[1] = 0 +  velYZ * rYZ * cos(angleYZ)                              - velXY * rXY * sin(angleXY) ;
          V[2] = 0 -  velYZ * rYZ * sin(angleYZ) + velZX * rZX * cos(angleZX)                            ;
          break;	
        }
      case SHEAR_YZ: // shear plane YZ  
        {
          V[0] = Vxyz[0] * ( x - Cxyz[0] );
          V[1] = Vxyz[1] * ( x - Cxyz[0] );
          V[2] = Vxyz[2] * ( x - Cxyz[0] );
          break;	
        }
      case SHEAR_ZX: // shear plane ZX  
        {
          V[0] = Vxyz[0] * ( y - Cxyz[1] );
          V[1] = Vxyz[1] * ( y - Cxyz[1] );
          V[2] = Vxyz[2] * ( y - Cxyz[1] );
          break;	
        }
      case SHEAR_XY: // shear plane XY  
        {
          V[0] = Vxyz[0] * ( z - Cxyz[2] );
          V[1] = Vxyz[1] * ( z - Cxyz[2] );
          V[2] = Vxyz[2] * ( z - Cxyz[2] );
          break;	
        }
      case BENCH2: // bench 2 
        {
          double Xmax0 = 10.;
          double X = Xmax0 + Vxyz[0] * t;
          V[0] = Vxyz[0] * x / X;
          V[1] = 0.0;
          V[2] = 0.0;
          break;
        }
      case CROSS_EXPANSION_WITH_ROTATION: // cross expantion with rotation
        {
          double cx = 0.45;
          double cy = 0.35 ;
          double r = sqrt ((x - cx)*(x - cx) + (y - cy)*(y - cy));
          double tt = atan2((x - cx), (y - cy));
          V[0] = 0 +  4*r * cos(tt) + sin(tt) * .3;
          V[1] = 0 -  4*r * sin(tt) + cos(tt) * .3;
          V[2] = 0;
          break;	
        }
      default:
        {
          MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                      "Unknowm velocity formulation %d",velType);
        }
      }
  }

  // ----------------------------------------------------------------------
  // ----------------------- MOBILE OBJECTS SET ---------------------------
  // ----------------------------------------------------------------------
  mobileObjectSet::mobileObjectSet(): elasticOp(NULL)
  {
  }

  // ----------------------------------------------------------------------
  mobileObjectSet::~mobileObjectSet() 
  {
    clear();
    if (elasticOp) delete elasticOp;
  }

  // ----------------------------------------------------------------------
  void mobileObjectSet::clear()
  {
    mobSet.clear();
  }

  // ----------------------------------------------------------------------
  void mobileObjectSet::insert(mobileObject* mob)
  {
    mobSet.insert(mob);
  }

  // ----------------------------------------------------------------------
  void mobileObjectSet::describe(std::ostream& out) const
  {
    out << "\n--- Mobile objects set description ---\n\n";
    set<mobileObject*>::const_iterator groupIt = mobSet.begin(); 
    set<mobileObject*>::const_iterator groupItEnd = mobSet.end(); 
    for ( ; groupIt != groupItEnd ; ++groupIt) {
      (*groupIt)->describe(out);
    }
  }

  // ----------------------------------------------------------------------
  void mobileObjectSet::computePrescribedDisplacement (double t, double dt)
  {
    prescribedDisplacement.clear();

    set<mobileObject*>::const_iterator groupIt = mobSet.begin(); 
    set<mobileObject*>::const_iterator groupItEnd = mobSet.end(); 
    for ( ; groupIt != groupItEnd ; ++groupIt) {

      // compute the forced displacement for each object
      (*groupIt)->computePrescribedDisplacement (t,dt);

      // gather all displacements
      set<vDisplacement,vDisplacementLess> objDisp = (*groupIt)->getPrescribedDisplacement();
      set<vDisplacement,vDisplacementLess>::const_iterator itObjDx = objDisp.begin();
      for(; itObjDx != objDisp.end(); itObjDx++) {
        vDisplacement vdisp(*itObjDx);
        prescribedDisplacement.insert(vdisp);
      }
    }
  }

  // ----------------------------------------------------------------------
  int mobileObjectSet::partlyMove(vertexMoveOp& vMoveOp, double t, double dt, double * part)
  {
    computePrescribedDisplacement(t, dt);
  
    int ok = 0;
    const int maxNumReduct = 8;
    const double reductFactor = 0.5;
    *part = 1.0;
    cout<<"Attempts to move vertices: [ ";
    for (int i=0; i < maxNumReduct; i++) {
      cout <<"* ";
      if (!vMoveOp.move(prescribedDisplacement,*part))  *part *= reductFactor;
      else { ok = 1; break;}
    }
    cout << "]\n";
    return ok;
  }

  // ----------------------------------------------------------------------
  void mobileObjectSet::setupElasticRepositioning(pMesh mesh, double t, double dt, 
                                                  double qualThr,
                                                  double chi, bool meshIsCavity,
                                                  int cavityThickness)
  {
    if ( mobSet.empty() ) return;

    // setup the elastic operator
    if (elasticOp) { delete elasticOp; elasticOp = NULL; }
    elasticOp = new MAdElasticityOp(mesh);
    if ( chi >= 0. ) elasticOp->setStiffnessAlterationCoef(chi);
    elasticOp->setCavityEqualMesh(meshIsCavity,cavityThickness);
    elasticOp->setQualityThreshold(qualThr);

    // prescribe bcs
    computePrescribedDisplacement(t, dt);
    set<vDisplacement,vDisplacementLess>::const_iterator dIter = prescribedDisplacement.begin();
    set<vDisplacement,vDisplacementLess>::const_iterator dLast = prescribedDisplacement.end();
    for (; dIter != dLast; dIter++) {
      pVertex vert = (*dIter).pv;
      smallVector disp(3);
      for (int i=0; i<3; i++)  disp(i) = (*dIter).dxyz[i];
      elasticOp->addDirichlet(vert,disp);
    }

    elasticOp->buildCavity();
    elasticOp->setDirichletBC();

    // compute elastic repositioning
    elasticOp->compute();
  }

  // ----------------------------------------------------------------------
  int mobileObjectSet::reposition(double * ratio)
  {
    if ( mobSet.empty() ) return true;
    int flag = elasticOp->advance(ratio);
    if ( flag == 2 ) elasticOp->clear();
    return flag;
  }

  // ----------------------------------------------------------------------
  bool mobileObjectSet::moveAndReposition(pMesh mesh, double t, double dt,
                                          bool subAdaptation,
                                          double qualityThreshold,
                                          double chi, bool meshIsCavity,
                                          int cavityThickness,
                                          MeshAdapter * ma)
  {
    setupElasticRepositioning(mesh, t, dt, qualityThreshold, chi, 
                              meshIsCavity, cavityThickness);
    bool joker = false;
    int subIter = 0;
    double ratio = 0.;
    int achieved = -1;
    while ( achieved != 2 ) {
      achieved = reposition(&ratio);
      MAdMsgSgl::instance().info(-1,__FILE__,
                                 "Advanced repositioning, achieved: %d, ratio: %f",
                                 achieved,ratio);

      if ( achieved <= 1 ) {
        if ( !subAdaptation ) {
          MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                      "Could not advance objects, ratio reached: %f",ratio);
        }
        else if ( !ma->removeSlivers() && !ma->optimiseElementShape() ) {
          if ( joker ) {
            MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                        "Could not advance objects, ratio reached: %f",
                                        ratio);
          }
          else {
            MAdMsgSgl::instance().warning(__LINE__,__FILE__,
                                          "Could not maintain quality at threshold value %f when repositioning nodes, ratio reached: %f",
                                          qualityThreshold,ratio);
            joker = true;
            elasticOp->setQualityThreshold(0.);
          }
        }
      }

      subIter++;
    }

    if ( joker )  elasticOp->setQualityThreshold(qualityThreshold);

    return true;
  }

  // ----------------------------------------------------------------------

}
