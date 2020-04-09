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

#include <array>

#include "CallbackManager.h"
#include "DistanceFunction.h"
#include "MAdDefines.h"
#include "MAdOutput.h"
#include "MeshSizeBase.h"

/*
// available in Gmsh repository, not in 2.4.2 version
#if defined(_HAVE_GMSH_)
#include "gmsh/GModel.h"
#include "gmsh/dofManager.h"
#include "gmsh/linearSystemGMM.h"
#include "gmsh/distanceTerm.h"
#endif
*/

namespace MAd {

  // -------------------------------------------------------------------
  void DistFct_CBMoveFct(pVertex pv, double * target, void * data)
  {
    // When a vertex is moved, the distance is cleared because it
    // can be recalculated on-the-fly when needed.
    // The gradients and curvatures are not cleared because they 
    // can be needed elsewhere and cannot be recalculated on-the-fly.

    distanceFunction * dFct = (distanceFunction *) data;
    dFct->clearDistance(pv);
  }

  // -------------------------------------------------------------------
  void DistFct_CBFct(pPList before, pPList after, void * data,
                     operationType type , pEntity ppp)
  {
    distanceFunction * dFct = (distanceFunction *) data;
    
    switch (type) {
    case MAd_ESPLIT: {
      // In the edge split case, we have to interpolate the data at the new node

      double grads[2][3], curv[2];
          
      // find the old edge
      void * tmp = 0;
      pEdge pe = (pEdge)PList_next(before,&tmp);
      double t = E_linearParams(pe,(pVertex)ppp);

      // interpolate distance gradient
      if ( dFct->getGradient( E_vertex((pEdge)pe, 0), grads[0] ) &&
           dFct->getGradient( E_vertex((pEdge)pe, 1), grads[1] )   )
        {
          double * newGrad = new double[3];
          for (int i=0; i<3; i++) newGrad[i] = (1.-t) * grads[0][i] + t * grads[1][i];

          // attach it
          dFct->attachGradient((pVertex)ppp, newGrad);
        }

      // interpolate curvature
      if ( dFct->getCurvature( E_vertex((pEdge)pe, 0), &(curv[0]) ) &&
           dFct->getCurvature( E_vertex((pEdge)pe, 1), &(curv[1]) )    )
        {
          double newCurv = (1.-t) * curv[0] + t * curv[1];
          dFct->attachCurvature((pVertex)ppp, newCurv);
        }

      break;
    } 
    case MAd_ECOLLAPSE: {
      // In the edge collapse case, we have to delete the datas attached to the deleted node
      dFct->clearVertexData((pVertex)ppp);
      break;
    }
    case MAd_ESWAP:
    case MAd_FSWAP: {
      // nothing to be done (no modification at nodes)
      break;
    }
    case MAd_RREMOVE: {
      void * temp = NULL;
      while ( pEntity pE = PList_next(before,&temp) ) {
        if ( EN_type(pE) == 0 ) {
          dFct->clearVertexData((pVertex)pE);
        }
      }
      break;
    }
    case MAd_UNKNOWNOPERATION:
    case MAd_VERTEXMOVE:
    case MAd_FCOLLAPSE:
    case MAd_DESPLTCLPS:
    default: {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                  "Not implemented for mesh modification %d", type);
    }
    }
  }

  // -------------------------------------------------------------------
  distanceFunction::distanceFunction(pMesh m, bool _distToEntities): 
    mesh(m), distToEntities(_distToEntities), nVert(0), nEnt(0), 
    nVE(-1), xyzV(NULL), entToV(NULL)
  {
    distId = MD_newMeshDataId();
    vGradId = MD_newMeshDataId();
    rGradId = MD_newMeshDataId();
    vCurvId = MD_newMeshDataId();

    kdSearch = new SearchTool();

    CallBackManagerSgl::instance().registerCallBackMove(DistFct_CBMoveFct,this);
    CallBackManagerSgl::instance().registerCallBack(DistFct_CBFct,this);
  }

  // -------------------------------------------------------------------
  distanceFunction::~distanceFunction()
  {
    clearDistances();
    MD_deleteMeshDataId(distId);

    clearGradientAtVertices();
    MD_deleteMeshDataId(vGradId);

    clearGradientInElements();
    MD_deleteMeshDataId(rGradId);

    clearCurvature();
    MD_deleteMeshDataId(vCurvId);

    if (xyzV) delete [] xyzV;
    if (entToV) delete [] entToV;
    if (kdSearch) delete kdSearch;

    CallBackManagerSgl::instance().unregisterCallBackMove(DistFct_CBMoveFct,this);
    CallBackManagerSgl::instance().unregisterCallBack(DistFct_CBFct,this);
  }
  
  // -------------------------------------------------------------------
  void distanceFunction::computeTree(const std::set<pVertex>& verts, 
                                     const std::set<pEntity>& ents)
  {
    // Build search tree and vertices coordinates

    kdSearch->reset();

    nVert = verts.size();
    kdSearch->allocatePoints(nVert);

    if ( xyzV ) delete [] xyzV;
    xyzV = new double[3*nVert];

    pvToSearchId.clear();
    vToEnt.clear();

    std::set<pVertex> :: const_iterator itV    = verts.begin();
    std::set<pVertex> :: const_iterator itVEnd = verts.end();
    int k = 0;
    for ( ; itV != itVEnd ; ++itV)
      {
        double coord[3]; V_coord(*itV,coord);
        kdSearch->addPoint(k,coord);

        pvToSearchId[*itV] = k;
        for (int i=0; i<3; i++) xyzV[k*3+i] = coord[i];

        k++;
      }

    kdSearch->allocateTree(nVert);

    // Build wall entities topology

    if ( distToEntities ) {
      nEnt = ents.size();
      if ( nEnt > 0 )
        {
          // allocate entities
          nVE = EN_numVertices(*(ents.begin()));
          if ( nVE > 3 ) {
            MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                        "Distance to entities with more than 3 vertices not implemented");
          }
          if (entToV) delete [] entToV;
          entToV = new int[nVE*nEnt];

          // fill it
          int m = 0;
          std::set<pEntity>::const_iterator itE = ents.begin();
          for (; itE != ents.end(); itE++) {
            pPList eVerts = EN_vertices(*itE);
            for (int i=0; i<nVE; i++) {
              int vId = pvToSearchId[(pVertex)PList_item(eVerts,i)];
              entToV[nVE*m+i] = vId;
              vToEnt.insert(std::make_pair(vId,m));
            }
            PList_delete(eVerts);

            // ensure faces normals are outgoing
            if ( EN_type(*itE) == 2 && M_dim(mesh) == 3 ) {
              if ( F_numRegions((pFace)*itE) != 1 ) {
                MAdMsgSgl::instance().warning(__LINE__,__FILE__,
                                              "Faces with %d region(s), normals will not be outgoing",
                                              F_numRegions((pFace)*itE));
              }
              else {
                pRegion pr = F_region((pFace)*itE,0);
                if ( R_dirUsingFace(pr,(pFace)*itE) != 1 ) {
                  int tmp = entToV[nVE*m+1];
                  entToV[nVE*m+1] = entToV[nVE*m+2];
                  entToV[nVE*m+2] = tmp;
                }
              }
            }

            m++;
          }
        }
      else {
        if (entToV) delete [] entToV;
        entToV = NULL;
      }
    }

    // Clear everything else

    clearDistances();
    clearGradientAtVertices();
    clearGradientInElements();
    clearCurvature();
  }
  
  // -------------------------------------------------------------------
  void distanceFunction::computeAllDistances() const
  {
    pVertex pv;
    VIter vit = M_vertexIter(mesh);
    while ( ( pv = VIter_next(vit) ) )
      {
        computeDistance(pv);
      }
    VIter_delete(vit);
  }
  
  // -------------------------------------------------------------------
  /*
  // not available in MAdLib 1.3.0 because not available in Gmsh 2.4.2
  void distanceFunction::computeAllDistancesEDP(const std::set<pGEntity> fixed) const
  {
// #warning "hack exact distance"
//     {
//     pVertex pv;
//     VIter vit = M_vertexIter(mesh);
//     while ( ( pv = VIter_next(vit) ) )
//       {
//         double xyz[3];
//         V_coord(pv,xyz);
//         double dist = sqrt(dotProd(xyz,xyz))-1.;
//         EN_attachDataDbl((pEntity)pv,distId,dist);
//       }
//     VIter_delete(vit);
//     return;
//     }

#if defined(_HAVE_GMSH_)
#if defined(_HAVE_GMM_)

    // GCFIXME: avoid that!!
    M_writeMsh(mesh,"tmpMesh.msh",2);
    GModel * pModel = new GModel();
    pModel->readMSH("tmpMesh.msh");
    
    
    int dimGmsh = pModel->getNumRegions() ? 3 : 2;  
    std::map<int, std::vector<GEntity*> > groups[4];
    pModel->getPhysicalGroups(groups);


    const int _ELTAG = 1;

    linearSystemGmm<double> *lsys = new linearSystemGmm<double>;
    lsys->setPrec(1.e-15);
    lsys->setGmres(1);
    lsys->setNoisy(1);
    dofManager<double,double> * pAssembler = new dofManager<double,double>(lsys);
    
    SBoundingBox3d bbox = pModel->bounds();
    double L = norm(SVector3(bbox.max(), bbox.min())); 
    double mu = L / 28.;
    simpleFunction<double> DIFF(mu * mu), ONE(1.0);
    dofManager<double, double> myAssembler(lsys);
    distanceTerm distanceT(NULL, 1, &DIFF, &ONE);
    printf("Distance computation with mu = %f ( L = %f )\n",mu,L);
    
    // fixation
    for (std::set<pGEntity>::const_iterator it = fixed.begin(); it != fixed.end(); ++it){
      int tpDim = GEN_type(*it);
      int physical = GEN_tag(*it);
      std::vector<MVertex *> v;
      pModel->getMeshVerticesForPhysicalGroup (tpDim,physical,v);
      printf("Fixing Physical %d, dim: %d, nb vert: %d\n",physical,tpDim,v.size());
      for (int i=0;i<v.size();i++){  
        pAssembler-> fixVertex(v[i] , 0, _ELTAG, 0.);
      }
    }
    
    // numbering
    std::map<int, std::vector<GEntity*> > gents = groups[dimGmsh];
    for (std::map<int, std::vector<GEntity*> >::iterator it = gents.begin(); it != gents.end(); it++) {
      int physical = it->first;
      std::vector<MVertex *> allNodes;
      pModel->getMeshVerticesForPhysicalGroup (dimGmsh,physical,allNodes);
      printf("Numbering Physical %d, dim: %d, nb vert: %d\n",physical,dimGmsh,allNodes.size());
      for (int i=0;i<allNodes.size();i++){  
        pAssembler->numberVertex (allNodes[i], 0, _ELTAG);
      }
    }

    // assembling
    for (std::map<int, std::vector<GEntity*> >::iterator it = gents.begin(); it != gents.end(); it++) {
      int physical = it->first;
      std::vector<GEntity*> ent = it->second;
      printf("Assembling Physical %d, dim: %d, nb ents: %d\n",physical,dimGmsh,ent.size());
      for (int i=0;i<ent.size();i++){
        distanceT.addToMatrix(*pAssembler,ent[i]);
        distanceT.addToRightHandSide(*pAssembler,ent[i]);
      }
    }
    
    // solving
    lsys->systemSolve();


    
    std::vector<GEntity*> ent = gents[99];
    int nbNodes;
    FILE *fsig = fopen ("sigma.pos","w");
    fprintf(fsig,"View \"Sigma\"{\n");
  
    for (int i=0;i<ent.size();i++){
      for (int j=0 ; j< ent[i]->getNumMeshElements() ; j++){
        MElement *e = ent[i]->getMeshElement(j);
        if (e->getDim() == dimGmsh){
          
          if ( dimGmsh == 2 ) { fprintf(fsig,"ST("); }
          if ( dimGmsh == 3 ) { fprintf(fsig,"SS("); }
          
          nbNodes = e->getNumVertices();
        
          for (int k=0;k<nbNodes-1;++k) {
            SPoint3 p = e->getVertex(k)->point();
            fprintf(fsig,"%g,%g,%g,",p.x(),p.y(),p.z());
          }
          for (int k=nbNodes-1;k<nbNodes;++k) {
            SPoint3 p = e->getVertex(k)->point();
            fprintf(fsig,"%g,%g,%g){",p.x(),p.y(),p.z());
          }
          for (int k=0;k<nbNodes-1;++k) {
            fprintf(fsig,"%g,",pAssembler->getDofValue(e->getVertex(k), 0, 1));
          }
          for (int k=nbNodes-1;k<nbNodes;++k) {
            fprintf(fsig,"%g};\n",pAssembler->getDofValue(e->getVertex(k), 0, 1));
          }
        }
      }
    }
    fprintf(fsig,"};\n");
    fclose(fsig);

    // transform and attach the distance
    pVertex pv;
    VIter vit = M_vertexIter(mesh);
    while ( ( pv = VIter_next(vit) ) )
      {
        MVertex * mv = pModel->getMeshVertexByTag( V_id(pv) );
        double value = std::min(1.0-1.e-8, pAssembler->getDofValue (mv, 0, _ELTAG));
        double dist = -mu * log(1. - value);
        EN_attachDataDbl((pEntity)pv,distId,dist);
      }
    VIter_delete(vit);

    lsys->clear();

#else
  MAdMsgSgl::instance().error(__LINE__,__FILE__,
                              "Cannot compute distance by EDP: GMM and Gmsh required");
#endif
#else
  MAdMsgSgl::instance().error(__LINE__,__FILE__,
                              "Cannot compute distance by EDP: GMM and Gmsh required");
#endif
  }
  */

  // -------------------------------------------------------------------
  double distanceFunction::getDistance(const pVertex pv) const
  {
    double dist;
    if ( EN_getDataDbl((pEntity)pv,distId,&dist) ) return dist;
    return computeDistance(pv);
  }

  // -------------------------------------------------------------------
  double distanceFunction::computeDistance(const pVertex pv) const
  {
    double xyz[3]; V_coord(pv,xyz);
    double dist = computeDistance(xyz);
    EN_attachDataDbl((pEntity)pv,distId,dist);
    return dist;
  }

  // -------------------------------------------------------------------
  double distanceFunction::computeDistance(const double xyz[3]) const
  {
    return sqrt(computeDistSq(xyz));
  }

  // -------------------------------------------------------------------
  double distanceFunction::computeDistSq(const double xyz[3]) const
  {
    if ( !distToEntities || nEnt == 0 ) return kdSearch->computeDistanceSq(xyz);

    double minD = MAdBIG;

    // --- compute distance to every entity and take minimum ---

    double xyzE[3][3]; // 1 line = coordinates of 1 vertex
    double vec[3];

    //GC note: could be improved

    // --- compute distance to entities neighboring closest point and take minimum ---
    int closestId;
    kdSearch->computeDistanceSq(xyz, &closestId);
    std::multimap<int,int>::const_iterator bfIt = vToEnt.find(closestId);
    assert ( bfIt != vToEnt.end() );
    while ( bfIt != vToEnt.end() && bfIt->first == closestId )
      {
        int iE = (*bfIt).second;
//     for (int iE=0; iE<nEnt; iE++)
//       {
        double dist;

        for(int iV=0; iV<nVE; iV++) {
          for(int iC=0; iC<3; iC++) xyzE[iV][iC] = xyzV[entToV[iE*nVE+iV]*3+iC];
        }

        if ( nVE == 3 ) {
          dist = distToTriangleSq(xyzE,xyz,vec);
        }
        else {
          MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                      "Distance to a set of edge not implemented yet");
        }

        minD = std::min(minD,dist);

        bfIt++;
      }
        
    return minD;
  }

  // -------------------------------------------------------------------
  void distanceFunction::computeAllDistAndGrad() const
  {
    double testVec[3];
    double xyzE[3][3]; // 1 line = coordinates of 1 vertex
    pVertex pv;
    VIter vit = M_vertexIter(mesh);
    while ( ( pv = VIter_next(vit) ) )
      {
        double dist;
        double * grad = new double[3];

        double xyz[3];
        V_coord(pv,xyz);

        if ( !distToEntities || nEnt == 0 )
          {
            int closestId;
            kdSearch->computeDistanceSq(xyz, &closestId);
            double xyzW[3];
            for (int iC=0; iC<3; iC++) xyzW[iC] = xyzV[3*closestId+iC];
            diffVec(xyz,xyzW,grad);
            dist = sqrt( dotProd(grad,grad) );
            if ( dist <= MAdTOL ) {
              MAdMsgSgl::instance().warning(__LINE__,__FILE__,
                                            "Zero vector found for a distance gradient");
            }
            else {
              normalizeVec(grad,grad);
            }
          }
        else
          {
            double minD = MAdBIG;

            //GC note: efficiency could be improved 

            for (int iE=0; iE<nEnt; iE++) {

              double distSq;

              for(int iV=0; iV<nVE; iV++) {
                for(int iC=0; iC<3; iC++) xyzE[iV][iC] = xyzV[entToV[iE*nVE+iV]*3+iC];
              }
            
              if ( nVE == 3 ) {
                distSq = distToTriangleSq(xyzE,xyz,testVec);
              }
              else {
                MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                            "Distance to a set of edge not implemented yet");
              }
              
              if ( distSq < minD ) {
                minD = distSq;
                for (int i=0; i<3; i++) grad[i] = testVec[i];
              }

            }

            // --- compute distance to entities neighboring closest point and take minimum ---
            // GC note: a part is not implemented yet: look at the neighbor entities!
//             int closestNodeId;
//             int closestEntId;
//             kdSearch->computeDistanceSq(xyz, &closestNodeId);
//             std::multimap<int,int>::const_iterator bfIt = vToEnt.find(closestNodeId);
//             assert ( bfIt != vToEnt.end() );
//             while ( bfIt != vToEnt.end() && bfIt->first == closestNodeId )
//               {
//                 int iE = (*bfIt).second;
//                 //             // --- compute distance to every entity and take minimum ---
//                 //         for (int iE=0; iE<nEnt; iE++)
//                 //           {
//                 double distSq;

//                 for(int iV=0; iV<nVE; iV++) {
//                   for(int iC=0; iC<3; iC++) xyzE[iV][iC] = xyzV[entToV[iE*nVE+iV]*3+iC];
//                 }
            
//                 if ( nVE == 3 ) {
//                   distSq = distToTriangleSq(xyzE,xyz,testVec);
//                 }
//                 else {
//                   MAdMsgSgl::instance().error(__LINE__,__FILE__,
//                                               "Distance to a set of edge not implemented yet");
//                 }
            
//                 if ( distSq < minD ) {
//                   closestEntId = iE;
//                   minD = distSq;
//                   for (int i=0; i<3; i++) grad[i] = testVec[i];
//                 }

//                 bfIt++;
//               }

            dist = sqrt(minD);

            // if the point is lying on the iso-zero, take the mean of the 
            // normals of the neighboring elements. The normals are supposed 
            // to be outgoing.
            if ( minD < MAdTOL ) 
              { 
                if (nVE != 3) throw;
                grad[0] = 0.; grad[1] = 0.; grad[2] = 0.;
                int locId = pvToSearchId.find(pv)->second;
                std::multimap<int,int>::const_iterator bfIt = vToEnt.find(locId);
                assert ( bfIt != vToEnt.end() );
                while ( bfIt != vToEnt.end() && bfIt->first == locId )
                  {
                    int entId = bfIt->second;
                    for(int iV=0; iV<nVE; iV++) {
                      for(int iC=0; iC<3; iC++) xyzE[iV][iC] = xyzV[entToV[entId*nVE+iV]*3+iC];
                    }
                
                    double nor[3];
                    XYZ_F_normal(xyzE,nor);
                    normalizeVec(nor,nor);
                
                    // give more weight to the closest information
                    double cen[3], xyzToC[3];
                    meanRow33(xyzE,cen);
                    diffVec(cen,xyz,xyzToC);
                    double invDistToCenSq = 1. / dotProd(xyzToC,xyzToC);
                
                    for (int i=0; i<3; i++) grad[i] -= nor[i] * invDistToCenSq;
                
                    bfIt++;
                  }
            
                normalizeVec(grad,grad);
              }
            else
              {
                double invD = 1./dist;
                for (int i=0; i<3; i++) grad[i] *= invD;
              }
          }
// #warning "hacked gradients"
// //         for (int i=0; i<3; i++) grad[i] = xyz[i];
//         for (int i=0; i<3; i++) grad[i] = 0.;
//         grad[0] = 1.;
//         dist = xyz[0];
//         normalizeVec(grad,grad);

        EN_attachDataDbl((pEntity)pv,distId,dist);
        void * del;
        if ( EN_getDataPtr((pEntity)pv,vGradId,&del) && del ) delete [] (double*)del;
        EN_attachDataPtr((pEntity)pv,vGradId,grad);
      }
    VIter_delete(vit);
  }

  // -------------------------------------------------------------------
  void distanceFunction::clearDistance(const pVertex pv) const
  {
    EN_deleteData((pEntity)pv,distId);
  }

  // -------------------------------------------------------------------
  void distanceFunction::clearDistances() const
  {
    pVertex pv;
    VIter vit = M_vertexIter(mesh);
    while ( ( pv = VIter_next(vit) ) )
      {
        EN_deleteData((pEntity)pv,distId);
      }
    VIter_delete(vit);
  }

  // -------------------------------------------------------------------
  void distanceFunction::computeGradientInElements() const
  {
    double gsf[4][3] = {
      {-1., -1., -1.}, 
      { 1.,  0.,  0.}, 
      { 0.,  1.,  0.}, 
      { 0.,  0.,  1.}, 
    };

    pRegion pr;
    RIter rit = M_regionIter(mesh);
    while ( ( pr = RIter_next(rit) ) )
      {
        pPList rVerts = R_vertices(pr);
        double dist[4];
        void * temp = NULL;
        pVertex pv;
        int i = 0;
        while ( ( pv = (pVertex)PList_next(rVerts,&temp) ) )
          {
            dist[i] = getDistance(pv);
            i++;
          }
        PList_delete(rVerts);

        double ijac[3][3]/*, detj*/;
        //detj = R_invJacobian(pr,ijac);

        void * del;
        if ( EN_getDataPtr((pEntity)pr,rGradId,&del) && del ) delete [] (double*)del;
        double * grad = new double[3];
        for ( int iC=0; iC<3; iC++)
          {
            grad[iC] = 0.;
            for (int iSF=0; iSF<4; iSF++) {
              grad[iC] += dist[iSF] * ( gsf[iSF][0] * ijac[0][iC] +
                                        gsf[iSF][1] * ijac[1][iC] +
                                        gsf[iSF][2] * ijac[2][iC]   );
            }
          }

        EN_attachDataPtr((pEntity)pr,rGradId,grad);
      }
    RIter_delete(rit);
  }

  // -------------------------------------------------------------------
  void distanceFunction::clearGradientInElements() const
  {
    pRegion pr;
    RIter rit = M_regionIter(mesh);
    while ( ( pr = RIter_next(rit) ) )
      {
        EN_deleteData((pEntity)pr,rGradId);
      }
    RIter_delete(rit);
  }

  // -------------------------------------------------------------------
  void distanceFunction::computeGradientAtVertices()
  {
    computeGradientInElements();

    pVertex pv;
    VIter vit = M_vertexIter(mesh);
    while ( ( pv = VIter_next(vit) ) )
      {
        void * del;
        if ( EN_getDataPtr((pEntity)pv,vGradId,&del) && del ) delete [] (double*)del;
        double * grad = new double[3];
        for (int i=0; i<3; i++) grad[i] = 0.;

        pPList vRegs = V_regions(pv);
        void * temp = NULL;
        void * rGrad = NULL;
        pRegion pr;
        while ( ( pr = (pRegion)PList_next(vRegs,&temp) ) )
          {
            if ( !EN_getDataPtr((pEntity)pr,rGradId,&rGrad) ) throw;
            grad[0] += ((double*)rGrad)[0];
            grad[1] += ((double*)rGrad)[1];
            grad[2] += ((double*)rGrad)[2];
          }

        int nRegs = PList_size(vRegs);
        PList_delete(vRegs);
        double invNRegs = 1. / (double)nRegs;
        for (int i=0; i<3; i++) grad[i] *= invNRegs;
        
        EN_attachDataPtr((pEntity)pv,vGradId,grad);
      }
    VIter_delete(vit);

    clearGradientInElements();
  }

  // -------------------------------------------------------------------
  void distanceFunction::clearGradientAtVertices()
  {
    pVertex pv;
    void *tmp;
    VIter vit = M_vertexIter(mesh);
    while ( ( pv = VIter_next(vit) ) )
      {
        if ( EN_getDataPtr((pEntity)pv,vGradId, &tmp) ) {
          EN_deleteData((pEntity)pv,vGradId);
        }
      }
    VIter_delete(vit);
  }

  // -------------------------------------------------------------------
  bool distanceFunction::getGradient(const pVertex pv, double grad[3]) const
  {
    void * tmp;
    if ( !EN_getDataPtr((pEntity)pv,vGradId,&tmp) ) return false;
    for (int i=0; i<3; i++) grad[i] = ((double*) tmp)[i];
    return true;
  }

  // -------------------------------------------------------------------
  bool distanceFunction::getCurvatureOnEntity(const pEntity entity, 
                                              const double xyz[3],
                                              double *curv) const
  {
    int type = EN_type(entity);
    switch(type) {
    case 0: return getCurvature((pVertex)entity,curv);
    case 1: {
      double u = E_linearParams((pEdge)entity,xyz);
      double c[2];
      if ( !getCurvature(E_vertex((pEdge)entity,0), &(c[0])) ||
           !getCurvature(E_vertex((pEdge)entity,1), &(c[1])) ) return false;
      *curv = (1.-u) * c[0] + u * c[1];
      return true;
    }
    case 2: {
      double u[2];
      F_linearParams((pFace)entity,xyz,u);
      double c[3];
      if ( !getCurvature(F_vertex((pFace)entity,0), &(c[0])) ||
           !getCurvature(F_vertex((pFace)entity,1), &(c[1])) ||
           !getCurvature(F_vertex((pFace)entity,2), &(c[2])) ) return false;
  
      *curv = (1.-u[0]-u[1]) * c[0] + u[0] * c[1] + u[1] * c[2];
      return true;
    }
    case 3: {
      double u[3];
      R_linearParams((pRegion)entity, xyz, u);
      double c[4];
      if ( !getCurvature(R_vertex((pRegion)entity,0), &(c[0])) ||
           !getCurvature(R_vertex((pRegion)entity,1), &(c[1])) ||
           !getCurvature(R_vertex((pRegion)entity,2), &(c[2])) ||
           !getCurvature(R_vertex((pRegion)entity,3), &(c[3])) ) return false;
  
      *curv = (1.-u[0]-u[1]-u[2]) * c[0] + u[0] * c[1] + u[1] * c[2] + u[2] * c[3];
      return true;
    }
    default: throw;
    }
  }

  // -------------------------------------------------------------------
  bool distanceFunction::getGradientOnEntity(const pEntity entity, 
                                             const double xyz[3],
                                             double grad[3]) const
  {
    int type = EN_type(entity);
    switch(type) {
    case 0: return getGradient((pVertex)entity,grad);
    case 1: return getGradientOnEdge((pEdge)entity,xyz,grad);
    case 2: return getGradientOnFace((pFace)entity,xyz,grad);
    case 3: return getGradientOnRegion((pRegion)entity,xyz,grad);
    default: throw;
    }
  }

  // -------------------------------------------------------------------
  bool distanceFunction::getGradientOnEdge(const pEdge edge, 
                                           const double xyz[3],
                                           double grad[3]) const
  {
    double u = E_linearParams(edge,xyz);
    return getGradientOnEdgeParam(edge,u,grad);
  }

  // -------------------------------------------------------------------
  bool distanceFunction::getGradientOnEdgeParam(const pEdge pE, 
                                                const double u,
                                                double grad[3]) const
  {
    double g0[3], g1[3];
    if ( !getGradient(E_vertex(pE,0), g0) ||
         !getGradient(E_vertex(pE,1), g1) ) return false;
  
    for (int i=0; i<3; i++) grad[i] = (1.-u) * g0[i] + u * g1[i];
    return true;
  }

  // -------------------------------------------------------------------
  bool distanceFunction::getGradientOnFace(const pFace face, 
                                           const double xyz[3],
                                           double grad[3]) const
  {
    double u[2];
    F_linearParams(face,xyz,u);
    return getGradientOnFaceParam(face,u,grad);
  }

  // -------------------------------------------------------------------
  bool distanceFunction::getGradientOnFaceParam(const pFace face, 
                                                const double u[2],
                                                double grad[3]) const
  {
    double g0[3], g1[3], g2[3];
    if ( !getGradient(F_vertex(face,0), g0) ||
         !getGradient(F_vertex(face,1), g1) ||
         !getGradient(F_vertex(face,2), g2) ) return false;
  
    for (int i=0; i<3; i++) grad[i] = (1.-u[0]-u[1]) * g0[i] + u[0] * g1[i] + u[1] * g2[i];
    return true;
  }

  // -------------------------------------------------------------------
  bool distanceFunction::getGradientOnRegion(const pRegion region, 
                                             const double xyz[3],
                                             double grad[3]) const
  {
    double u[3];
    R_linearParams(region, xyz, u);
    return getGradientOnRegionParam(region, u, grad);
  }

  // -------------------------------------------------------------------
  bool distanceFunction::getGradientOnRegionParam(const pRegion region, 
                                                  const double u[3],
                                                  double grad[3]) const
  {
    double g0[3], g1[3], g2[3], g3[3];
    if ( !getGradient(R_vertex(region,0), g0) ||
         !getGradient(R_vertex(region,1), g1) ||
         !getGradient(R_vertex(region,2), g2) ||
         !getGradient(R_vertex(region,3), g3) ) return false;
  
    for (int i=0; i<3; i++) {
      grad[i] = (1.-u[0]-u[1]-u[2]) * g0[i] + u[0] * g1[i] + u[1] * g2[i] + u[2] * g3[i];
    }
    return true;
  }

  // -------------------------------------------------------------------
  void distanceFunction::attachGradient(pVertex pv, double grad[3]) const
  {
    EN_attachDataPtr((pEntity)pv,vGradId,grad);
  }

  // -------------------------------------------------------------------
  void distanceFunction::outputDistance(const char * fn) const
  {
    MAdAttachedNodalDataOutput(mesh, fn, distId);
  }

  // -------------------------------------------------------------------
  void distanceFunction::outputGradAtVertices(const char * fn) const
  {
    MAdAttachedNodalDataVecOutput(mesh, fn, vGradId);
  }

  // -------------------------------------------------------------------
  void distanceFunction::clearVertexData(pVertex pv) const
  {
    EN_deleteData((pEntity)pv, vGradId);
    EN_deleteData((pEntity)pv, distId);
    EN_deleteData((pEntity)pv, vCurvId);
  }

  // -------------------------------------------------------------------
  void distanceFunction::computeCurvature(const std::set<pRegion>& regs)
  {
    int nR = regs.size();        // number of regions involved
    std::vector<int> r2v(nR*4);  // regions vertices
    std::map<pFace,pRegion> bFaces;  // boundary faces (to 'regs') and 
    //                                  the corresponding region
    std::map<pVertex,int> v2i;   // vertex pointers to local ids

    // --- fill v2i and r2v ---

    pVertex pv;
    pRegion pr;
    pFace pf;
    std::map<pVertex,int>::iterator v2iIt;
    int vId=0, rId=0;
    std::set<pRegion>::const_iterator rIt = regs.begin();
    for (; rIt != regs.end(); rIt++) {
      // vertices
      for (int i=0; i<4; i++) {
        pv = R_vertex(*rIt,i);
        v2iIt = v2i.find(pv);
        if ( v2iIt == v2i.end() ) {
          v2i[pv] = vId;
          r2v[rId*4+i] = vId;
          vId++;
        }
        else {
          r2v[rId*4+i] = (*v2iIt).second;
        }
      }
      // boundary faces
      for (int i=0; i<4; i++) {
        pf = R_face(*rIt,i);
        pr = F_otherRegion(pf,*rIt);
        if ( !pr || regs.find( pr ) == regs.end() ) {
          bFaces[pf] = *rIt;
        }
      }
      rId++;
    }

    // --- make the reverse mapping of v2i: i2v ---

    int nV = v2i.size();           // number of vertices involved
    std::vector<pVertex> i2v(nV);  // local ids to vertex pointers

    std::map<pVertex,int>::const_iterator v2iIter = v2i.begin();
    for(; v2iIter != v2i.end(); v2iIter++) i2v[v2iIter->second] = v2iIter->first;
    
    // --- compute right hand side and mass matrix ---
    
    double dPhidXi[4][3] = {
      {-1., -1., -1.}, {1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}
    };

    double over12 = 1./12.;
    double over24 = 1./24.;

    std::vector<double> rhs(nV), MLump(nV);
    for (int i=0; i<nV; i++) {
      rhs[i] = 0.;
      MLump[i] = 0.;
    }

    void * tmpP;
    double vGrad[4][3];
    double ijac[3][3], detJ, dPhiDx[4][3];
    int iR = 0;
    for ( rIt = regs.begin(); rIt != regs.end(); rIt++ )
      {
        // get inverse Jacobian and determinant
        detJ = R_invJacobian(*rIt,ijac);

        // --- shape functions gradients and distance gradients ---
        for (int iSF=0; iSF<4; iSF++)
          {
            int iId = r2v[iR*4+iSF];
            assert ( EN_getDataPtr( (pEntity) ( i2v[iId] ), vGradId, &tmpP) );
            for ( int iC=0; iC<3; iC++) {
              dPhiDx[iSF][iC] = ( dPhidXi[iSF][0] * ijac[0][iC] +
                                  dPhidXi[iSF][1] * ijac[1][iC] +
                                  dPhidXi[iSF][2] * ijac[2][iC]   );
              vGrad[iSF][iC] = ((double*)tmpP)[iC];
            }
          }

        // --- rhs and mass matrix ---
        int iId, jSF;
        for (int iSF=0; iSF<4; iSF++) {
          iId = r2v[iR*4+iSF];
          for (jSF=0; jSF<4; jSF++) {
            rhs[iId] -= over24 * detJ * dotProd( dPhiDx[iSF], vGrad[jSF] );
          }
          MLump[iId] += over24 * detJ;
        }

        // --- boundary term ---
        
        double nor[3], detJ;
        pPList fVerts;
        for (int iF=0; iF<4; iF++ )
          {
            pf = R_face(*rIt, iF);
            if ( bFaces.find(pf) != bFaces.end() )
              {
                // we will need the outgoing normal...
                //  ... and the ordered vertices ...
                F_normal(pf,nor);
                normalizeVec(nor,nor);
                if ( R_dirUsingFace(*rIt,pf) ) fVerts = F_vertices(pf, 1);
                else {
                  fVerts = F_vertices(pf, 0);
                  for(int i=0; i<3; i++) nor[i] = -1.*nor[i];
                }
                
                // ... the jacobian of the face ...
                detJ = 2. * F_area(pf);
                
                // ... the gradient of distance at vertices.
                for (int iVF=0; iVF<3; iVF++) {
                  pv = (pVertex)PList_item(fVerts,iVF);
                  EN_getDataPtr( (pEntity)pv, vGradId, &tmpP);
                  for(int i=0; i<3; i++) vGrad[iVF][i] = ((double*)tmpP)[i];
                }

                // compute the term
                for (int iVF=0; iVF<3; iVF++) {
                  pv = (pVertex)PList_item(fVerts,iVF);
                  int iId = v2i.find(pv)->second;
                  for (int jVF=0; jVF<3; jVF++) {
                    if ( iVF == jVF ) {
                      rhs[iId] += over12 * detJ * dotProd( nor, vGrad[jVF] );
                    }
                    else {
                      rhs[iId] += over24 * detJ * dotProd( nor, vGrad[jVF] );
                    }
                  }
                }
                PList_delete(fVerts);
              }
          }

        iR++;
      }
    

    // --- solve the system to find curvatures ---
    for (int iV=0; iV<nV; iV++)
      {
        EN_attachDataDbl((pEntity)(i2v[iV]),vCurvId, fabs( rhs[iV] / MLump[iV] ));
      }
  }

  // -------------------------------------------------------------------
  void distanceFunction::computeGradientAndCurvature(const std::set<pRegion>& regs)
  {
    int nR = regs.size();        // number of regions involved
    std::vector<int> r2v(nR*4);  // regions vertices
    std::map<pFace,pRegion> bFaces;  // boundary faces (to 'regs') and 
    //                                  the corresponding region
    std::map<pVertex,int> v2i;   // vertex pointers to local ids

    // --- fill v2i and r2v ---

    pVertex pv;
    pRegion pr;
    pFace pf;
    std::map<pVertex,int>::iterator v2iIt;
    int vId=0, rId=0;
    std::set<pRegion>::const_iterator rIt = regs.begin();
    for (; rIt != regs.end(); rIt++) {
      // vertices
      for (int i=0; i<4; i++) {
        pv = R_vertex(*rIt,i);
        v2iIt = v2i.find(pv);
        if ( v2iIt == v2i.end() ) {
          v2i[pv] = vId;
          r2v[rId*4+i] = vId;
          vId++;
        }
        else {
          r2v[rId*4+i] = (*v2iIt).second;
        }
      }
      // boundary faces
      for (int i=0; i<4; i++) {
        pf = R_face(*rIt,i);
        pr = F_otherRegion(pf,*rIt);
        if ( !pr || regs.find( pr ) == regs.end() ) {
          bFaces[pf] = *rIt;
        }
      }
      rId++;
    }

    // --- make the reverse mapping of v2i: i2v ---

    int nV = v2i.size();           // number of vertices involved
    std::vector<pVertex> i2v(nV);  // local ids to vertex pointers

    std::map<pVertex,int>::const_iterator v2iIter = v2i.begin();
    for(; v2iIter != v2i.end(); v2iIter++) i2v[v2iIter->second] = v2iIter->first;
    
    // --- compute right hand side and mass matrix ---
    
    double dPhidXi[4][3] = {
      {-1., -1., -1.}, {1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}
    };

    double sixth      = 1./6.;
    double oneOver60  = 1./60.;
    double oneOver120 = 1./120.;

    std::vector<double> rhs(nV), MLump(nV);
    std::vector<double> vGrad(nV*3);  // distance gradients computed on the vertices
//     int v_nR[nV];        // number of involved regions around vertices
    for (int i=0; i<nV; i++) {
      rhs[i] = 0.;
      MLump[i] = 0.;
//       v_nR[i] = 0;
      for (int j=0; j<3; j++) vGrad[3*i+j] = 0.;
    }


    double dist[4], ijac[3][3], detJ, dPhiDx[4][3];
    double rGrad[3];
    int iR = 0;
    for ( rIt = regs.begin(); rIt != regs.end(); rIt++ )
      {
        // get distances
        for (int iV = 0; iV<4; iV++) dist[iV] = getDistance(R_vertex(*rIt,iV));
        
        // get inverse Jacobian and determinant
        detJ = R_invJacobian(*rIt,ijac);

        // --- gradients ---
        rGrad[0] = 0.; rGrad[1] = 0.; rGrad[2] = 0.;
        for (int iSF=0; iSF<4; iSF++)
          {
            for ( int iC=0; iC<3; iC++) {
              dPhiDx[iSF][iC] = ( dPhidXi[iSF][0] * ijac[0][iC] +
                                  dPhidXi[iSF][1] * ijac[1][iC] +
                                  dPhidXi[iSF][2] * ijac[2][iC]   );
              rGrad[iC] += dist[iSF] * dPhiDx[iSF][iC];
            }
          }
        if ( fabs(rGrad[0]) > 1.e-12 ||
             fabs(rGrad[1]) > 1.e-12 || 
             fabs(rGrad[2]) > 1.e-12 )   normalizeVec(rGrad,rGrad);
        for (int iSF=0; iSF<4; iSF++)
          {
            for ( int iC=0; iC<3; iC++) {
              vGrad[( r2v[iR*4+iSF] )*3+iC] += rGrad[iC];
            }
//             v_nR[ r2v[iR*4+iSF] ] ++;
          }

        // --- rhs and mass matrix ---
        int iId;
        double detjO60  = detJ * oneOver60;
        double detjO120 = detJ * oneOver120;
        for (int iSF=0; iSF<4; iSF++) {
          iId = r2v[iR*4+iSF];
          rhs[iId] -= sixth * detJ * (dPhiDx[iSF][0] * rGrad[0] +
                                      dPhiDx[iSF][1] * rGrad[1] +
                                      dPhiDx[iSF][2] * rGrad[2]   );
          
          for (int j=0; j<4; j++) {
            if ( iSF == j ) {
              MLump[ iId ]  += detjO60;
            }
            else {
              MLump[ iId ]  += detjO120;
            }
          }
        }

        // --- boundary term ---
        
        double nor[3], term;
        for (int iF=0; iF<4; iF++ )
          {
            pf = R_face(*rIt, iF);
            if ( bFaces.find(pf) != bFaces.end() )
              {
                // we will need the outgoing normal...
                //  ... and the ordered vertices ...
                F_normal(pf,nor);
                normalizeVec(nor,nor);
                pPList fVerts;
                if ( R_dirUsingFace(*rIt,pf) ) fVerts = F_vertices(pf, 1);
                else {
                  fVerts = F_vertices(pf, 0);
                  for(int i=0; i<3; i++) nor[i] = -1.*nor[i];
                }
                
                // ... the jacobian of the face ...
                double detJ = 2. * F_area(pf);

                // compute the term
// #warning "debug: analytical grad"
//                 if ( nor[0] > 0.5 ) {
//                   printVec(rGrad,"rGrad");
//                   rGrad[0]=2;
//                   rGrad[1]=0.; rGrad[2]=0.;
//                 }
//                 if ( nor[0] < -0.5 ) {
//                   rGrad[0]=0;
//                   rGrad[1]=0.; rGrad[2]=0.;
//                 }
//                 else {
//                   rGrad[1]=0.; rGrad[2]=0.;
//                 }

                term = sixth * dotProd( nor, rGrad ) * detJ;
                for (int iVF=0; iVF<3; iVF++) {
                  rhs  [ v2i[(pVertex)PList_item(fVerts,iVF)] ] += term;
                }
              }
          }

        iR++;
      }
    
//     fprintf (frg,"};\n");
//     fclose (frg);
// #warning "debug output rhs"
//     FILE *f = fopen ("testrhs.pos", "w");
//     fprintf (f,"View\" mesh \" {\n");
//     iR = 0;
//     for ( rIt = regs.begin(); rIt != regs.end(); rIt++ )
//       {
//         // get the coordinates
//         double xyz[4][3];
//         R_coordP1(*rIt, xyz);
        
//         // get the data at nodes
//         double data[4];
//         pPList verts = R_vertices(pr);
//         void* tmp = 0;
//         int i = 0;
//         while ( pEntity pent = PList_next(verts,&tmp) ) {
//           //          data[i] = rhs2(v2i[(pVertex)pent]);
//           data[i] = rhs[r2v[iR*4+i]];
//           i++;
//         }
//         PList_delete(verts);

//         // write an element
//         fprintf (f,"SS(%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g,%g};\n",
//                  xyz[0][0],xyz[0][1],xyz[0][2],
//                  xyz[1][0],xyz[1][1],xyz[1][2],
//                  xyz[2][0],xyz[2][1],xyz[2][2],
//                  xyz[3][0],xyz[3][1],xyz[3][2],
//                  data[0],data[1],data[2],data[3]);

//         iR++;
//       }
//     fprintf (f,"};\n");
//     fclose (f);
    

    // --- solve the system to find curvatures ---
    for (int iV=0; iV<nV; iV++)
      {
        EN_attachDataDbl((pEntity)(i2v[iV]),vCurvId, fabs( rhs[iV] / MLump[iV] ));
      }

    // --- compute gradients at vertices ---
    void * del;
    for ( int iV=0; iV<nV; iV++ )
      {
        pv = i2v[iV];
        if ( EN_getDataPtr((pEntity)pv,vGradId,&del) && del ) delete [] (double*)del;

        double * grad = new double[3];
        for (int iC=0; iC<3; iC++) grad[iC] = vGrad[iV*3+iC]; // / (double)(v_nR[iV]);
        if ( fabs(grad[0]) > 1.e-12 ||
             fabs(grad[1]) > 1.e-12 || 
             fabs(grad[2]) > 1.e-12 ) normalizeVec(grad,grad);

        EN_attachDataPtr((pEntity)pv,vGradId,grad);
      }
  }
  // -------------------------------------------------------------------
  //GC: to debug
  void distanceFunction::computeGradientAndCurvature2D(const std::set<pFace>& faces)
  {

    int nF = faces.size();       // number of faces involved
    std::vector<int> f2v(nF*3);  // face vertices
    std::map<pEdge,pFace> bEdges;  // boundary edges (to 'faces') and 
    //                                  the corresponding face
    std::map<pVertex,int> v2i;   // vertex pointers to local ids

    // --- fill v2i and r2v ---

    pVertex pv;
    pFace pf;
    pEdge pe;
    std::map<pVertex,int>::iterator v2iIt;
    int vId=0, fId=0;
    std::set<pFace>::const_iterator fIt = faces.begin();
    for (; fIt != faces.end(); fIt++) {
      // vertices
      for (int i=0; i<3; i++) {
        pv = F_vertex(*fIt,i);
        v2iIt = v2i.find(pv);
        if ( v2iIt == v2i.end() ) {
          v2i[pv] = vId;
          f2v[fId*3+i] = vId;
          vId++;
        }
        else {
          f2v[fId*3+i] = (*v2iIt).second;
        }
      }
      // boundary edges
      for (int i=0; i<3; i++) {
        pe = F_edge(*fIt,i);
        pf = E_otherFace(pe,*fIt);
        if ( !pf || faces.find( pf ) == faces.end() ) {
          bEdges[pe] = *fIt;
        }
      }
      fId++;
    }

    // --- make the reverse mapping of v2i: i2v ---

    int nV = v2i.size();           // number of vertices involved
    std::vector<pVertex> i2v(nV);  // local ids to vertex pointers

    std::map<pVertex,int>::const_iterator v2iIter = v2i.begin();
    for(; v2iIter != v2i.end(); v2iIter++) i2v[v2iIter->second] = v2iIter->first;
    
    // --- compute right hand side and mass matrix ---
    
    double gradShpFct[4][3] =
      { {-1., -1., -1.},
        { 1.,  0.,  0.},
        { 0.,  1.,  0.},
        { 0.,  0.,  1.} };
    if ( M_dim(mesh) < 3 ) gradShpFct[0][2] = 0.;

    double sixth      = 1./6.;
    double oneOver60  = 1./60.;
    double oneOver120 = 1./120.;

    std::vector<double> rhs(nV), MLump(nV);
    std::vector<double> vGrad(nV*3);  // distance gradients computed on the vertices
    std::vector<int> v_nF(nV);        // number of involved regions around vertices
    for (int i=0; i<nV; i++) {
      rhs[i] = 0.;
      MLump[i] = 0.;
      v_nF[i] = 0;
      for (int j=0; j<3; j++) vGrad[3*i+j] = 0.;
    }

    constexpr int nbNodes = 3;
    double xyz[nbNodes][3];
    double dist[nbNodes], invjac[3][3], jac[3][3], detJ;
    double fGrad[3];
    //double Grads[4][3];
    int iF = 0;
    for (fIt = faces.begin(); fIt != faces.end(); fIt++) {
         F_coordP1(*fIt, xyz);

        // get distances
        for (int iV = 0; iV<3; iV++) dist[iV] = getDistance(F_vertex(*fIt,iV));
        
        // get inverse Jacobian and determinant
        jac[0][0] = jac[0][1] = jac[0][2] = 0.;
        jac[1][0] = jac[1][1] = jac[1][2] = 0.;
        jac[2][0] = jac[2][1] = jac[2][2] = 0.;
        for(int i = 0; i < nbNodes; i++) {
          jac[0][0] += xyz[i][0] * gradShpFct[i][0]; jac[0][1] += xyz[i][1] * gradShpFct[i][0]; jac[0][2] += xyz[i][2] * gradShpFct[i][0];
          jac[1][0] += xyz[i][0] * gradShpFct[i][1]; jac[1][1] += xyz[i][1] * gradShpFct[i][1]; jac[1][2] += xyz[i][2] * gradShpFct[i][1];
          jac[2][0] += xyz[i][0] * gradShpFct[i][2]; jac[2][1] += xyz[i][1] * gradShpFct[i][2]; jac[2][2] += xyz[i][2] * gradShpFct[i][2];
        }
        
        if ( M_dim(mesh) < 3 ) {
          for (int i=0;i<2;i++) { jac[i][2] = jac[2][i] = 0.; }
          jac[2][2] = 1.;
        }

        if ( M_dim(mesh) == 3 ) detJ = fabs ( inverseMat (jac, invjac) );
        else {
          double itmp[2][2]; 
          double tmp[2][2]; for (int i=0;i<2;i++) for (int j=0;j<2;j++) tmp[i][j]=jac[i][j];
          detJ = fabs ( inverseMat22(tmp,itmp) );
          for (int i=0;i<2;i++) for (int j=0;j<2;j++) invjac[i][j]=itmp[i][j];
          for (int i=0;i<2;i++) { invjac[i][2] = invjac[2][i] = 0.; }
          invjac[2][2] = 1.;
        }

        // --- gradients ---
        fGrad[0] = 0.; fGrad[1] = 0.; fGrad[2] = 0.;
        for (int j=0;j<nbNodes;j++){
          for ( int iC=0; iC<3; iC++) {
 /*
            Grads[j][iC] =
              invjac[iC][0] * gradShpFct[j][0] + 
              invjac[iC][1] * gradShpFct[j][1] + 
              invjac[iC][2] * gradShpFct[j][2];
*/
            fGrad[iC] += dist[j] * gradShpFct[j][iC];
          }
        }
        normalizeVec(fGrad,fGrad);
        for (int iSF=0; iSF<3; iSF++)
          {
            for ( int iC=0; iC<3; iC++) {
              vGrad[( f2v[iF*3+iSF] )*3+iC] += fGrad[iC];
            }
            v_nF[ f2v[iF*3+iSF] ] ++;
          }

        // --- rhs and mass matrix ---
        int iId;
        double detjO60  = detJ * oneOver60;
        double detjO120 = detJ * oneOver120;
        for (int iSF=0; iSF<3; iSF++) {
          iId = f2v[iF*3+iSF];
          rhs[iId] -= sixth * detJ * (gradShpFct[iSF][0] * fGrad[0] +
                                      gradShpFct[iSF][1] * fGrad[1] +
                                      gradShpFct[iSF][2] * fGrad[2]   );
          
          for (int j=0; j<3; j++) {
            if ( iSF == j ) {
              MLump[ iId ]  += detjO60;
            }
            else {
              MLump[ iId ]  += detjO120;
            }
          }
        }

        // --- boundary term ---
        
//         double nor[3], term;
//         for (int iE=0; iE<3; iE++ )
//           {
//             pe = F_edge(*fIt, iE);
//             if ( bEdges.find(pe) != bEdges.end() )
//               {
//                 // we will need the outgoing normal...
//                 //  ... and the ordered vertices ...
//                 // GCTODO
//                 F_normal(pf,nor);
//                 normalizeVec(nor,nor);
//                 pPList fVerts;
//                 if ( R_dirUsingFace(*rIt,pf) ) fVerts = F_vertices(pf, 1);
//                 else {
//                   fVerts = F_vertices(pf, 0);
//                   for(int i=0; i<3; i++) nor[i] = -1.*nor[i];
//                 }
                
//                 // ... the jacobian of the face ...
//                 double detJ = 2. * F_area(pf);

//                 // compute the term
// // #warning "debug: analytical grad"
// //                 if ( nor[0] > 0.5 ) {
// //                   printVec(rGrad,"rGrad");
// //                   rGrad[0]=2;
// //                   rGrad[1]=0.; rGrad[2]=0.;
// //                 }
// //                 if ( nor[0] < -0.5 ) {
// //                   rGrad[0]=0;
// //                   rGrad[1]=0.; rGrad[2]=0.;
// //                 }
// //                 else {
// //                   rGrad[1]=0.; rGrad[2]=0.;
// //                 }

//                 term = sixth * dotProd( nor, rGrad ) * detJ;
//                 for (int iVF=0; iVF<3; iVF++) {
//                   rhs  [ v2i[(pVertex)PList_item(fVerts,iVF)] ] += term;
//                 }
//               }
//           }

        iF++;
      }
    
//     fprintf (frg,"};\n");
//     fclose (frg);
// #warning "debug output rhs"
//     FILE *f = fopen ("testrhs.pos", "w");
//     fprintf (f,"View\" mesh \" {\n");
//     iR = 0;
//     for ( rIt = regs.begin(); rIt != regs.end(); rIt++ )
//       {
//         // get the coordinates
//         double xyz[4][3];
//         R_coordP1(*rIt, xyz);
        
//         // get the data at nodes
//         double data[4];
//         pPList verts = R_vertices(pr);
//         void* tmp = 0;
//         int i = 0;
//         while ( pEntity pent = PList_next(verts,&tmp) ) {
//           //          data[i] = rhs2(v2i[(pVertex)pent]);
//           data[i] = rhs[r2v[iR*4+i]];
//           i++;
//         }
//         PList_delete(verts);

//         // write an element
//         fprintf (f,"SS(%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g,%g};\n",
//                  xyz[0][0],xyz[0][1],xyz[0][2],
//                  xyz[1][0],xyz[1][1],xyz[1][2],
//                  xyz[2][0],xyz[2][1],xyz[2][2],
//                  xyz[3][0],xyz[3][1],xyz[3][2],
//                  data[0],data[1],data[2],data[3]);

//         iR++;
//       }
//     fprintf (f,"};\n");
//     fclose (f);
    

    // --- solve the system to find curvatures ---
    for (int iV=0; iV<nV; iV++)
      {
        EN_attachDataDbl((pEntity)(i2v[iV]),vCurvId, fabs( rhs[iV] / MLump[iV] ));
      }

    // --- compute gradients at vertices ---
    void * del;
    for ( int iV=0; iV<nV; iV++ )
      {
        pv = i2v[iV];
        if ( EN_getDataPtr((pEntity)pv,vGradId,&del) && del ) delete [] (double*)del;

        double * grad = new double[3];
        for (int iC=0; iC<3; iC++) grad[iC] = vGrad[iV*3+iC] / (double)(v_nF[iV]);
        
        EN_attachDataPtr((pEntity)pv,vGradId,grad);
      }
  }

  // -------------------------------------------------------------------
  void distanceFunction::smoothCurvature(double maxGrad) const
  {
    pEdge pe;
    pVertex v0, v1;
    double curv0, curv1, maxDiff, diff;
    //double xyz[2][3];
    
    bool mod = true;
    while ( mod )
      {
        mod = false;

        EIter eit = M_edgeIter(mesh);
        while ( ( pe = EIter_next(eit) ) )
          {
            v0 = E_vertex(pe,0);
            v1 = E_vertex(pe,1);

            if ( EN_getDataDbl((pEntity)v0,vCurvId,&curv0 ) &&
                 EN_getDataDbl((pEntity)v1,vCurvId,&curv1 ) ) 
              {
                maxDiff = maxGrad * E_length(pe);
                diff = fabs(curv1 - curv0);
                
                if ( diff > 1.000001 * maxDiff ) {
                  if ( curv1 > curv0 ) {
                    curv0 = curv1 - maxDiff;
                    EN_attachDataDbl((pEntity)v0,vCurvId, curv0);
                  }
                  else {
                    curv1 = curv0 - maxDiff;
                    EN_attachDataDbl((pEntity)v1,vCurvId, curv1);
                  }
                  mod = true;
                }
              }
          }
        EIter_delete(eit);
      }
  }

  // -------------------------------------------------------------------
  void distanceFunction::limitCurvature(double maxCurv) const
  {
    double locCurv;
    pVertex pv;
    VIter vit = M_vertexIter(mesh);
    while ( ( pv = VIter_next(vit) ) )
      {
        if ( EN_getDataDbl((pEntity)pv,vCurvId,&locCurv ) ) {
          if ( locCurv > maxCurv ) {
            EN_attachDataDbl((pEntity)pv,vCurvId, maxCurv);
          }
        }
      }
    VIter_delete(vit);
  }

  // -------------------------------------------------------------------
  void distanceFunction::smoothCurvatureDummy(int nbSmoothings) const
  {
    pVertex pv, pv2;
    double curv, weights, locCurv;
    pPList vEdges;
    pEdge edge;

    for (int iS=0; iS < nbSmoothings; iS++)
      {
        VIter vit = M_vertexIter(mesh);
        while ( ( pv = VIter_next(vit) ) )
          {
            if ( EN_getDataDbl((pEntity)pv,vCurvId,&locCurv ) ) 
              {
                weights = 0.;
                curv = 0.;
                vEdges = V_edges(pv);
                void * tmp = NULL;
                while ( ( edge = (pEdge)PList_next(vEdges,&tmp) ) ) {
                  double tmpCurv;
                  pv2 = E_otherVertex(edge,pv);
                  if ( EN_getDataDbl((pEntity)pv2,vCurvId,&tmpCurv ) )
                    {
                      double invDistSq = 1. / E_length(edge);
                      curv += tmpCurv * invDistSq;
                      weights += invDistSq;
                    }
                }
                
                if ( weights < MAdTOL ) curv = locCurv;
                else {
                  curv /= weights;
                  curv = 0.5 * ( curv + locCurv );
                }
                
                EN_attachDataDbl((pEntity)pv,vCurvId, curv);
              }
          }
        VIter_delete(vit);
      }
  }

  // -------------------------------------------------------------------
  bool distanceFunction::getCurvature(const pVertex pv, double *c) const
  {
    if ( EN_getDataDbl((pEntity)pv, vCurvId, c) ) return true;
    return false;
  }

  // -------------------------------------------------------------------
  void distanceFunction::attachCurvature(pVertex pv, double curv) const
  {
    EN_attachDataDbl((pEntity)pv,vCurvId,curv);
  }

  // -------------------------------------------------------------------
  void distanceFunction::clearCurvature() const
  {
    pVertex pv;
    VIter vit = M_vertexIter(mesh);
    while ( ( pv = VIter_next(vit) ) )
      {
        EN_deleteData((pEntity)pv,vCurvId);
      }
    VIter_delete(vit);
  }

  // -------------------------------------------------------------------
  void distanceFunction::outputCurvature(const char * fn) const
  {
    MAdAttachedNodalDataOutput(mesh, fn, vCurvId);
  }

  // -------------------------------------------------------------------
}
