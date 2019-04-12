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

#include "MAdElasticityOp.h"
#include "MAdResourceManager.h"
#include "OperatorTools.h"
#include "MathUtils.h"
#include "MAdOutput.h"
#include "CallbackManager.h"
#include "MAdMessage.h"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
using std::make_pair;
using std::map;
using std::set;

namespace MAd {

  // -------------------------------------------------------------------
  void ElasticityOpCBFunction (pPList before, pPList after, void * data,
                               operationType type , pEntity ppp) 
  {
    MAdElasticityOp * elast = (MAdElasticityOp *)(data);
    if ( !(elast->relocationsComputed()) ) return;

    switch (type) {
    case MAd_ESPLIT: {
      
      // find the old edge
      void * temp = NULL;
      pEntity pE = PList_next(before,&temp);

      // take the interpolated displacements between the extremities 
      // of the edge for the new vertex
      elast->addVertexOnEdge( (pVertex)ppp, (pEdge)pE );

      break;
    } 
    case MAd_ECOLLAPSE: {
      elast->removeVertex( (pVertex)ppp );
      break;
    }
    case MAd_FSWAP:
    case MAd_ESWAP: {
      break;
    }
    case MAd_RREMOVE: {
      void * temp = NULL;
      while ( pEntity pE = PList_next(before,&temp) ) {
        if ( EN_type(pE) == 0 ) {
          elast->removeVertex( (pVertex)pE );
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
  }

  // -------------------------------------------------------------------
  MAdElasticityOp::MAdElasticityOp(pMesh m):
    mesh(m), E(1.), nu(0.45), qualityThreshold(0.), 
    chi(1.0), computed(false), cavityEqualMesh(false), cavityThickness(3)
  {
    dim = M_dim(mesh);
    CallBackManagerSgl::instance().registerCallBack(ElasticityOpCBFunction,this);
  }

  // -------------------------------------------------------------------
  MAdElasticityOp::MAdElasticityOp(const MAdElasticityOp& eop):
    mesh(eop.mesh), dim(eop.dim), 
    E(eop.E), nu(eop.nu), chi(eop.chi),
    computed(eop.computed),
    cavityEqualMesh(eop.cavityEqualMesh),
    cavityThickness(eop.cavityThickness)
  {
    localIds = eop.localIds;
    dirichlet = eop.dirichlet;
    CallBackManagerSgl::instance().registerCallBack(ElasticityOpCBFunction,this);
  }

  // -------------------------------------------------------------------
  MAdElasticityOp::~MAdElasticityOp()
  {
    CallBackManagerSgl::instance().unregisterCallBack(ElasticityOpCBFunction,this);
  }

  // -------------------------------------------------------------------
  void MAdElasticityOp::setMaterials(double _E, double _nu)
  {
    E = _E;  nu = _nu;
  }

  // -------------------------------------------------------------------
  void MAdElasticityOp::setStiffnessAlterationCoef(double _chi)
  {
    chi = _chi;
  }

  // -------------------------------------------------------------------
  void MAdElasticityOp::clear()
  {
    localIds.clear();
    cavity.clear();
    dirichlet.clear();
    relocations.clear();
    computed = false;
  }

  // -------------------------------------------------------------------
  void MAdElasticityOp::buildCavity()
  {
    computed = false;

    if ( cavityEqualMesh ) return;

    cavity.clear();

    if ( cavityThickness < 1 ) return;

    if (dim == 3)
      {
        set<pVertex> innerNodes; // nodes inside the cavity and dirichlet
        set<pRegion> prevLayer;

        // --- first layer ---
        pRegion pr;
        pVertex pv;
        map<pVertex,smallVector>::const_iterator dIter = dirichlet.begin();
        for (; dIter != dirichlet.end(); dIter++) {
          pv = (*dIter).first;
          pPList vRegs = V_regions(pv);
          void * temp = NULL;
          while ( ( pr = (pRegion)PList_next(vRegs,&temp) ) ) {
            cavity.insert((pEntity)pr);
            prevLayer.insert(pr);
          }
          PList_delete(vRegs);
          innerNodes.insert(pv);
        }

        // --- next layers ---
        for (int iL=1; iL < cavityThickness; iL++) {
          set<pRegion> curLayer;
          set<pRegion>::const_iterator rIter = prevLayer.begin();
          for (; rIter != prevLayer.end(); rIter++ ) {
            pPList rVerts = R_vertices(*rIter);
            void * temp = NULL;
            while ( ( pv = (pVertex)PList_next(rVerts,&temp) ) ) {
              if ( innerNodes.find(pv) == innerNodes.end() ) {
                pPList vRegs = V_regions(pv);
                void * temp = NULL;
                while ( ( pr = (pRegion)PList_next(vRegs,&temp) ) ) {
                  if ( cavity.find((pEntity)pr) == cavity.end() ) {
                    cavity.insert((pEntity)pr);
                    curLayer.insert(pr);
                  }
                }
                PList_delete(vRegs);
                innerNodes.insert(pv);
              }
            }
            PList_delete(rVerts);
          }
          prevLayer.clear();
          prevLayer.insert(curLayer.begin(),curLayer.end());
          curLayer.clear();
        }
      }
    else 
      {
        set<pVertex> innerNodes; // nodes inside the cavity and dirichlet
        set<pFace>   prevLayer;

        // --- first layer ---
        pFace pf;
        pVertex pv;
        map<pVertex,smallVector>::const_iterator dIter = dirichlet.begin();
        for (; dIter != dirichlet.end(); dIter++) {
          pv = (*dIter).first;
          pPList vFaces = V_faces(pv);
          void * temp = NULL;
          while ( ( pf = (pFace)PList_next(vFaces,&temp) ) ) {
            cavity.insert((pEntity)pf);
            prevLayer.insert(pf);
          }
          PList_delete(vFaces);
          innerNodes.insert(pv);
        }

        // --- next layers ---
        for (int iL=1; iL < cavityThickness; iL++) {
          set<pFace> curLayer;
          set<pFace>::const_iterator fIter = prevLayer.begin();
          for (; fIter != prevLayer.end(); fIter++ ) {
            pPList fVerts = F_vertices(*fIter,1);
            void * temp = NULL;
            while ( ( pv = (pVertex)PList_next(fVerts,&temp) ) ) {
              if ( innerNodes.find(pv) == innerNodes.end() ) {
                pPList vFaces = V_faces(pv);
                void * temp = NULL;
                while ( ( pf = (pFace)PList_next(vFaces,&temp) ) ) {
                  if ( cavity.find((pEntity)pf) == cavity.end() ) {
                    cavity.insert((pEntity)pf);
                    curLayer.insert(pf);
                  }
                }
                PList_delete(vFaces);
                innerNodes.insert(pv);
              }
            }
            PList_delete(fVerts);
          }
          prevLayer.clear();
          prevLayer.insert(curLayer.begin(),curLayer.end());
          curLayer.clear();
        }
      }

    printf("Built a cavity with %d elements, dirichlet size: %d\n",cavity.size(),dirichlet.size());

  }

  // -------------------------------------------------------------------
  void MAdElasticityOp::setHomogDirichletBC()
  {
    dirichlet.clear();
    setDirichletBC();
  }

  // -------------------------------------------------------------------
  void MAdElasticityOp::setDirichletBC()
  {
    computed = false;
    smallVector zero(3);  zero.set_all(0.);

    set<pVertex> verts;
    collectBCVertices(&verts);

    set<pVertex>::const_iterator vIter = verts.begin();
    set<pVertex>::const_iterator vLast = verts.end();
    for (; vIter != vLast; vIter++) {
      addDirichlet(*vIter,zero);
    }
  }

  // -------------------------------------------------------------------
  void MAdElasticityOp::collectBCVertices(set<pVertex> * bcVerts)
  {
    if ( cavityEqualMesh ) 
      {
        if ( dim == 3 )
          {
            FIter fit = M_faceIter(mesh);
            while ( pFace pf = FIter_next(fit) ) {
              if ( EN_whatInType((pEntity)pf) == 2 ) {
                pPList fV = F_vertices(pf,1);
                void* tmp=0;
                while ( pVertex pv = (pVertex)PList_next(fV,&tmp) ) {
                  (*bcVerts).insert(pv);
                }
                PList_delete(fV);
              }
            }
            FIter_delete(fit);
          }
        else
          {
            EIter eit = M_edgeIter(mesh);
            while ( pEdge pe = EIter_next(eit) ) {
              if ( EN_whatInType((pEntity)pe) == 1 ) {
                (*bcVerts).insert( E_vertex(pe,0) );
                (*bcVerts).insert( E_vertex(pe,1) );
              }
            }
            EIter_delete(eit);
          }
      }
    else
      {
        if ( dim == 3 )
          {
            set<pEntity>::const_iterator cavIter = cavity.begin();
            for (; cavIter != cavity.end(); cavIter++) {
              pPList rFaces = R_faces((pRegion) *cavIter);
              void * temp = NULL;
              pFace pf;
              while ( ( pf = (pFace)PList_next(rFaces,&temp) ) ) {
                pRegion otherR = F_otherRegion(pf,(pRegion) *cavIter);
                if ( !otherR || cavity.find((pEntity)otherR) == cavity.end() ) {
                  pPList fVerts = F_vertices(pf,1);
                  void * temp2 = NULL;
                  pVertex pv;
                  while ( ( pv = (pVertex)PList_next(fVerts,&temp2) ) ) {
                    (*bcVerts).insert(pv);
                  }
                  PList_delete(fVerts);
                }
              }
              PList_delete(rFaces);
            }
          }
        else
          {
            set<pEntity>::const_iterator cavIter = cavity.begin();
            for (; cavIter != cavity.end(); cavIter++) {
              pPList fEdges = F_edges((pFace) *cavIter);
              void * temp = NULL;
              pEdge pe;
              while ( ( pe = (pEdge)PList_next(fEdges,&temp) ) ) {
                pFace otherF = E_otherFace(pe,(pFace) *cavIter);
                if ( !otherF || cavity.find((pEntity)otherF) == cavity.end() ) {
                  (*bcVerts).insert( E_vertex(pe,0) );
                  (*bcVerts).insert( E_vertex(pe,1) );
                }
              }
              PList_delete(fEdges);
            }
          }
      }
  }

  // -------------------------------------------------------------------
  void MAdElasticityOp::addDirichlet(pVertex pv, smallVector& dxyz)
  {
    map<pVertex,smallVector>::iterator dIter = dirichlet.find(pv);
    if ( dIter != dirichlet.end() ) {
      smallVector dxyz1 = (*dIter).second;
      dxyz1.add(dxyz);
      dirichlet.erase(dIter);
      dirichlet.insert(make_pair(pv,dxyz1));
    }
    else {
      dirichlet.insert(make_pair(pv,dxyz));
    }
  }

  // -------------------------------------------------------------------
  void MAdElasticityOp::delDirichlet(pVertex pv)
  {
    map<pVertex,smallVector>::iterator dIter = dirichlet.find(pv);
    if ( dIter != dirichlet.end() ) {
      dirichlet.erase(dIter);
    }
  }

  // -------------------------------------------------------------------
  void MAdElasticityOp::delRelocation(pVertex pv)
  {
    map<pVertex,smallVector>::iterator dIter = relocations.find(pv);
    if ( dIter != relocations.end() ) {
      relocations.erase(dIter);
    }
  }

  // -------------------------------------------------------------------
  // Generates table of vertex pointers to local ids and returns number
  // of vertices in the cavity.
  int MAdElasticityOp::generateVertexIds()
  {
    computed = false;
    localIds.clear();

    int id = 0;

    if (cavityEqualMesh) {
      VIter vit = M_vertexIter(mesh);
      while ( pVertex pv = VIter_next(vit) ) {
        localIds[pv] = id;
        id++;
      }
      VIter_delete(vit);
    }
    else {
      set<pVertex> verts;
      set<pEntity>::const_iterator cavIter = cavity.begin();
      for (; cavIter != cavity.end(); cavIter++) {
        pPList enVerts;
        if ( dim == 3 ) enVerts = R_vertices( (pRegion) *cavIter);
        else            enVerts = F_vertices( (pFace) *cavIter, 1);
        void * temp = NULL;
        while ( pVertex pv = (pVertex)PList_next(enVerts,&temp) ) {
          if ( verts.find(pv) == verts.end() ) {
            localIds[pv] = id;
            id++;
            verts.insert(pv);
          }
        }
        PList_delete(enVerts);
      }
    }
    return id;
  }

  // -------------------------------------------------------------------
  void MAdElasticityOp::addToMatrix(const pEntity pe,
                                    int nbVert,
                                    MAdLinearSystemDef * lsys)
  {
    int nodes;
    if ( dim == 3 ) nodes = R_numVertices((pRegion) pe);
    else            nodes = F_numVertices((pFace) pe);
    int locSyze = nodes*3;

    smallMatrix locK(locSyze,locSyze);
    elementMatrix(pe, locK);
//     if ( dim == 3 ) element3DMatrix((pRegion)pe, locK);
//     else            element2DMatrix((pFace)  pe, locK);

    pPList vertsI,vertsJ;
    if ( dim == 3 ) { vertsI = R_vertices((pRegion)pe); vertsJ = R_vertices((pRegion)pe); }
    else            { vertsI = F_vertices((pFace)pe,1); vertsJ = F_vertices((pFace)pe,1); }

    void* tmpI = 0;
    int i = 0;
    while ( pVertex pvI = (pVertex)PList_next(vertsI,&tmpI) ) {

      int locIdI = localIds[pvI];
      
      void* tmpJ = 0;
      int j = 0;
      while ( pVertex pvJ = (pVertex)PList_next(vertsJ,&tmpJ) ) {
        
        int locIdJ = localIds[pvJ];
        
        for (int k=0; k<3; k++) {
          for (int l=0; l<3; l++) {
            lsys->addToMatrix(nbVert*k + locIdI, nbVert*l+locIdJ, locK(k*nodes+i,l*nodes+j));
          }
        }
        j++;
      }
      i++;
    }
    PList_delete(vertsI);
    PList_delete(vertsJ);
  }

  // -------------------------------------------------------------------
  void MAdElasticityOp::allocateMatrix(int nbVert, MAdLinearSystemDef * lsys)
  {
    if (cavityEqualMesh) {
      VIter vit2 = M_vertexIter(mesh);
      while ( pVertex pv2 = VIter_next(vit2) ) {
        int locId = localIds[pv2];
        int nbNeighbors = V_numEdges(pv2);
        for (int iDir=0; iDir<3; iDir++) {
          int row = nbVert*iDir + locId;
          lsys->set_nnz(row,(nbNeighbors+1)*3);// A node is coupled with all its neighbors, and don't forget himself
        }
      }
      VIter_delete(vit2);
    }
    else {
      map<pVertex,set<pVertex> > conn;

      if ( dim == 3 )
        {
          set<pEntity>::const_iterator cIter = cavity.begin();
          for (; cIter != cavity.end(); cIter++) {
            pPList rEdges = R_edges((pRegion)*cIter);
            void * temp = NULL;
            pEdge edge;
            pVertex v0, v1;
            while ( ( edge = (pEdge)PList_next(rEdges,&temp)) ) {
              v0 = E_vertex(edge,0);
              v1 = E_vertex(edge,1);
              conn[v0].insert(v1);
              conn[v1].insert(v0);
            }
            PList_delete(rEdges);
          }
        }
      else
        {
          set<pEntity>::const_iterator cIter = cavity.begin();
          for (; cIter != cavity.end(); cIter++) {
            pPList fEdges = F_edges((pFace)*cIter);
            void * temp = NULL;
            pEdge edge;
            pVertex v0, v1;
            while ( ( edge = (pEdge)PList_next(fEdges,&temp)) ) {
              v0 = E_vertex(edge,0);
              v1 = E_vertex(edge,1);
              conn[v0].insert(v1);
              conn[v1].insert(v0);
            }
            PList_delete(fEdges);
          }
        }
        
      pVertex pv;
      map<pVertex,set<pVertex> >::const_iterator conIter = conn.begin();
      for (; conIter != conn.end(); conIter++) {
        pv = conIter->first;
        int locId = localIds[pv];
        int nbNeighbors = conIter->second.size();
        for (int iDir=0; iDir<3; iDir++) {
          int row = nbVert*iDir + locId;
          lsys->set_nnz(row,(nbNeighbors+1)*3);// A node is coupled with all its neighbors, and don't forget himself
        }
      }
    }

    lsys->allocate_matrix();
  }

  // -------------------------------------------------------------------
  void MAdElasticityOp::addVertexOnEdge(pVertex pv, pEdge pe)
  {
    if ( computed == true )
      {
        double t = E_linearParams(pe,pv);

        smallVector dx0, dx1;
        for (int i=0; i<2; i++)
          {
            std::map<pVertex,smallVector >::iterator iter = relocations.find(pv);
            if ( iter != relocations.end() ) {
              if (i==0) dx0 = (*iter).second;
              else      dx1 = (*iter).second;
            }
            else return;
          }
        smallVector newDx(3);
        for (int i=0; i<3; i++) newDx(i) = (1.-t) * dx0(i) + t * dx1(i);
        relocations[pv] = newDx;
      }
  }

  // -------------------------------------------------------------------
  void MAdElasticityOp::removeVertex(pVertex pv)
  {
    if ( computed == true ) 
      {
        std::map<pVertex,smallVector >::iterator iter = relocations.find(pv);
        if ( iter != relocations.end() ) {
          relocations.erase(iter);
        }
      }
  }

  // -------------------------------------------------------------------
  void MAdElasticityOp::setVertexResult(int nbVert, pVertex pv, 
                                        const MAdLinearSystemDef * lsys)
  {
    // get the vertex id
    map<pVertex,int>::const_iterator itId = localIds.find(pv);
    int id = (*itId).second;

    // get the dx
    smallVector dX(3);
    for (int iDir=0; iDir<3; iDir++) {
      int row = nbVert*iDir + id;
      dX(iDir) = lsys->getFromSolution(row);
    }
    //    if ( dim < 3 ) dX(2) = 0.;

    relocations[pv] = dX;
  }
  
  // -------------------------------------------------------------------
  void MAdElasticityOp::setResults(int nbVert,
                                   const MAdLinearSystemDef * lsys)
  {
    relocations.clear();

    if (cavityEqualMesh) {
      pVertex pv;
      VIter vit = M_vertexIter(mesh);
      while ( ( pv = VIter_next(vit) ) ) {
        setVertexResult(nbVert,pv,lsys);
      }
      VIter_delete(vit);
    }
    else {
      set<pVertex> vertices;
      if ( dim == 3 )
        {
          set<pEntity>::const_iterator cIter = cavity.begin();
          for (; cIter != cavity.end(); cIter++) {
            pPList rVerts = R_vertices((pRegion)*cIter);
            void * temp = NULL;
            pVertex pv;
            while ( ( pv = (pVertex)PList_next(rVerts,&temp) ) ) {
              vertices.insert(pv);
            }
            PList_delete(rVerts);
          }
        }
      else
        {
          set<pEntity>::const_iterator cIter = cavity.begin();
          for (; cIter != cavity.end(); cIter++) {
            pPList fVerts = F_vertices((pFace)*cIter,1);
            void * temp = NULL;
            pVertex pv;
            while ( ( pv = (pVertex)PList_next(fVerts,&temp) ) ) {
              vertices.insert(pv);
            }
            PList_delete(fVerts);
          }
        }

      set<pVertex>::const_iterator vIter = vertices.begin();
      for (; vIter != vertices.end(); vIter++) {
        setVertexResult(nbVert,*vIter,lsys);
      }
    }
  }

  // -------------------------------------------------------------------
  void MAdElasticityOp::printRelocations(std::ostream& out) const
  {
    out << "Printing relocation: "<<relocations.size()<<" vertices relocated\n";
    std::map<pVertex,smallVector>::const_iterator iter = relocations.begin();
    for(; iter != relocations.end(); iter++)
      {
        pVertex pv = (*iter).first;
        int id = V_id(pv);
        pGEntity ge = V_whatIn(pv);
        int gdim = GEN_type(ge);
        int gtag = GEN_tag(ge);
        smallVector dX = (*iter).second;
        out <<"Vertex "<<pv<<" with id "<<id<<" and class ("<<gdim<<","<<gtag<<"):";
        for (int i=0; i<3; i++) out << " " << dX(i);
        out <<"\n";
      }
    
    out << "Printed relocation: "<<relocations.size()<<" vertices relocated\n";
  }

  // -------------------------------------------------------------------
  void MAdElasticityOp::printDirichlet(std::ostream& out) const
  {
    out << "Printing dirichlet BCs: "<<dirichlet.size()<<" vertices\n";
    std::map<pVertex,smallVector>::const_iterator iter = dirichlet.begin();
    for(; iter != dirichlet.end(); iter++)
      {
        pVertex pv = (*iter).first;
        int id = V_id(pv);
        pGEntity ge = V_whatIn(pv);
        int gdim = GEN_type(ge);
        int gtag = GEN_tag(ge);
        smallVector dX = (*iter).second;
        out <<"Vertex "<<pv<<" with id "<<id<<" and class ("<<gdim<<","<<gtag<<"):";
        for (int i=0; i<3; i++) out << " " << dX(i);
        out <<"\n";
      }
    
    out << "Printed dirichlet BCs: "<<dirichlet.size()<<" vertices\n";
  }

  // -------------------------------------------------------------------
  int MAdElasticityOp::compute()
  {
    MAdResourceManager& tm = MAdResourceManagerSgl::instance();

    relocations.clear();

    // --------------------------------------------
    // If there is no element in the cavity, just 
    // prescribed node relocations (dirichlet) ...
    // --------------------------------------------
    if ( !cavityEqualMesh && cavity.empty() ) {
      MAdMsgSgl::instance().info(__LINE__,__FILE__,
                                 "Empty cavity for elastic computation");
      relocations = dirichlet;
    }

    // --------------------------------------------------
    // ... if an elastic computation has to be performed
    // --------------------------------------------------
    else {
      
      MeshQualityManagerSgl::instance().evaluateSizes(); // required for detJ0

      MAdResourceManagerSgl::instance().printMemoryUsage("Before elastic system allocation");

      // --- produce the mapping pVertex -> localId ---
      int nbVert = generateVertexIds();

      // --- create the linear system ---
      double t_sys0 = tm.getTime();
      int sysSize = 3 * nbVert;
      MAdLinearSystemDef * lsys = new MAdLinearSystemDef();
      lsys->allocate(sysSize);

      lsys->setSolver(CG);
      lsys->setFillIn(4);
      lsys->setEps(1.e-12);

      lsys->setPrec(1.e-12); // for gmm only
      lsys->setNoisy(1);    // for gmm only

      // --- allocate the matrix ---
      allocateMatrix(nbVert,lsys); // for PETSc only

      lsys->zeroMatrix();
      lsys->zeroRightHandSide();

      double dt_sys = tm.getTime() - t_sys0;
      MAdMsgSgl::instance().info(-1,__FILE__,
                                 "Created the system in %f seconds",dt_sys);
      MAdResourceManagerSgl::instance().printMemoryUsage("Elastic system allocated");

      // --- assemble the matrix ---
      double t_ass0 = tm.getTime();
      if (cavityEqualMesh) {
        if ( dim == 3 ) {
          RIter rit = M_regionIter(mesh);
          while ( pRegion pr = RIter_next(rit) ) addToMatrix((pEntity)pr,nbVert,lsys);
          RIter_delete(rit);
        }
        else {
          FIter fit = M_faceIter(mesh);
          while ( pFace pf = FIter_next(fit) ) addToMatrix((pEntity)pf,nbVert,lsys);
          FIter_delete(fit);
        }
      }
      else {
        set<pEntity>::const_iterator cavIter = cavity.begin();
        for (; cavIter != cavity.end(); cavIter++) addToMatrix(*cavIter,nbVert,lsys);
      }
      double dt_ass = tm.getTime() - t_ass0;
      MAdMsgSgl::instance().info(-1,__FILE__,
                                 "Assembled the matrix in %f seconds",dt_ass);
      MAdResourceManagerSgl::instance().printMemoryUsage("Matrix assembled");

      // --- apply bc (ugly but it works hum) ---
      double t_bc0 = tm.getTime();
      map<pVertex,smallVector>::iterator dIter = dirichlet.begin();
      map<pVertex,smallVector>::iterator dLast = dirichlet.end();
      for (; dIter != dLast; dIter++) {

        pVertex pv       = (*dIter).first;
        smallVector dxyz = (*dIter).second;

        // get the vertex id
        map<pVertex,int>::const_iterator itId = localIds.find(pv);
        int id = (*itId).second;
    
        for (int iDir=0; iDir<3; iDir++) {
          int row = nbVert*iDir + id;
          const double BIG = /*lsys->getFromMatrix(row,row) *  */1.e12;
          lsys->addToMatrix(row,row,BIG);
          lsys->addToRightHandSide(row,BIG * dxyz(iDir));
        }
    
      }
      if ( dim == 2 ) {
        MAdMsgSgl::instance().warning(__LINE__,__FILE__,"Fixing Z coordinates");
        for ( int iV=0; iV<nbVert; iV++) {
          for (int iD=2; iD<3; iD++){
            int row = nbVert*iD + iV;
            lsys->addToMatrix(row,row,1.e12);
          }
        }
      }

      double dt_bc = tm.getTime() - t_bc0;
      MAdMsgSgl::instance().info(-1,__FILE__,
                                 "Applied boundary conditions in %f seconds",dt_bc);
      MAdResourceManagerSgl::instance().printMemoryUsage("BC applied");

      // --- reorder matrix ---
      lsys->reorder();
      MAdResourceManagerSgl::instance().printMemoryUsage("Matrix reordered");

      // --- solve the system ---
      double t_slv0 = tm.getTime();
      int converged = lsys->systemSolve();
      double dt_slv = tm.getTime() - t_slv0;
      MAdMsgSgl::instance().info(-1,__FILE__,
                                 "Solved the elastic system (convergence flag: %d) in %f seconds",
                                 converged,dt_slv);
      MAdResourceManagerSgl::instance().printMemoryUsage("System solved");


      // --- fill relocations with the results ---
      setResults(nbVert,lsys);

      localIds.clear();
      delete lsys;
    }

    dirichlet.clear();
    cavity.clear();

    computed = true;

    return 1;
  }


  // -------------------------------------------------------------------
  bool MAdElasticityOp::checkArea(const pFace pf, double ratio)
  {
    bool inCavity = false;

    double fxyz[3][3];
    F_coordP1(pf,fxyz);
    double oriNorm[3];
    XYZ_F_normal(fxyz,oriNorm);

    pPList fVerts = F_vertices(pf,1);
    void * temp = NULL;
    int i = 0;
    pVertex pv;
    while ( ( pv = (pVertex)PList_next(fVerts,&temp) ) ) {
      std::map<pVertex,smallVector>::const_iterator iter = relocations.find(pv);
      if ( iter != relocations.end() ) {
        inCavity = true;
        for (int j=0; j<3; j++) fxyz[i][j] += ratio * (*iter).second(j);
      }
      i++;
    }
    PList_delete(fVerts);

    // Do the check only if one of the nodes moved
    if ( inCavity ) {
      double tgtNorm[3];
      XYZ_F_normal(fxyz,tgtNorm);
      double prod = oriNorm[0] * tgtNorm[0] + oriNorm[1] * tgtNorm[1] + oriNorm[2] * tgtNorm[2];
      if ( prod <= 0. ) return false;
    }

    return true;
  }

  // -------------------------------------------------------------------
  // Cannot use 'cavity': it is hard to maintain if mesh modifications
  // occur during advancement so it is cleared in 'compute()'.
  bool MAdElasticityOp::checkAreas(double ratio)
  {
    pFace pf;
    FIter fit = M_faceIter(mesh);
    while ( ( pf = FIter_next(fit) ) ) {
      if ( !checkArea(pf,ratio) ) {
        FIter_delete(fit);
        return false;
      }
    }
    FIter_delete(fit);

    return true;
  }

  // -------------------------------------------------------------------
  bool MAdElasticityOp::checkVolume(const pRegion pr, double ratio)
  {
    bool inCavity = false;

    double rxyz[4][3];
    R_coordP1(pr,rxyz);

    pPList rVerts = R_vertices(pr);
    void * temp = NULL;
    int i = 0;
    pVertex pv;
    while ( ( pv = (pVertex)PList_next(rVerts,&temp) ) ) {
      std::map<pVertex,smallVector>::const_iterator iter = relocations.find(pv);
      if ( iter != relocations.end() ) {
        inCavity = true;
        for (int j=0; j<3; j++) rxyz[i][j] += ratio * (*iter).second(j);
      }
      i++;
    }
    PList_delete(rVerts);

    // Do the check only if one of the nodes moved
    if ( inCavity && R_XYZ_volume(rxyz) <= MAdTOL ) return false;

    return true;
  }

  // -------------------------------------------------------------------
  // Cannot use 'cavity': it is hard to maintain if mesh modifications
  // occur during advancement so it is cleared in 'compute()'.
  bool MAdElasticityOp::checkVolumes(double ratio)
  {
    pRegion pr;
    RIter rit = M_regionIter(mesh);
    while ( ( pr = RIter_next(rit) ) ) {
      if ( !checkVolume(pr,ratio) ) {
        RIter_delete(rit);
        return false;
      }
    }
    RIter_delete(rit);

    return true;
  }

  // -------------------------------------------------------------------
  bool MAdElasticityOp::checkFQuality(const pFace pf, double ratio)
  {
    bool inCavity = false;

    double disp[3][3] = {
      {0., 0., 0.},
      {0., 0., 0.},
      {0., 0., 0.} };

    pPList fVerts = F_vertices(pf,1);
    void * temp = NULL;
    int i = 0;
    pVertex pv;
    while ( ( pv = (pVertex)PList_next(fVerts,&temp) ) ) {
      std::map<pVertex,smallVector>::const_iterator iter = relocations.find(pv);
      if ( iter != relocations.end() ) {
        inCavity = true;
        for (int j=0; j<3; j++) disp[i][j] += ratio * (*iter).second(j);
      }
      i++;
    }
    PList_delete(fVerts);

    // Do the check only if one of the nodes moved
    if ( inCavity ) {
      double oriNorm[3];
      F_normal(pf,oriNorm);
      double shape;
      MeshQualityManagerSgl::instance().getShapeWithDisp(pf,oriNorm,disp,&shape);
      if ( shape < qualityThreshold + MAdTOL ) return false;
    }

    return true;
  }

  // -------------------------------------------------------------------
  bool MAdElasticityOp::checkRQuality(const pRegion pr, double ratio)
  {
    bool inCavity = false;

    double disp[4][3] = {
      {0., 0., 0.},
      {0., 0., 0.},
      {0., 0., 0.},
      {0., 0., 0.} };
    
    pPList rVerts = R_vertices(pr);
    void * temp = NULL;
    int i = 0;
    pVertex pv;
    while ( ( pv = (pVertex)PList_next(rVerts,&temp) ) ) {
      std::map<pVertex,smallVector>::const_iterator iter = relocations.find(pv);
      if ( iter != relocations.end() ) {
        inCavity = true;
        for (int j=0; j<3; j++) disp[i][j] += ratio * (*iter).second(j);
      }
      i++;
    }
    PList_delete(rVerts);

    // Do the check only if one of the nodes moved
    if ( inCavity ) {
      double shape;
      MeshQualityManagerSgl::instance().getShapeWithDisp(pr,disp,&shape);
      if ( shape < qualityThreshold + MAdTOL ) return false;
    }

    return true;
  }

  // -------------------------------------------------------------------
  // Cannot use 'cavity': it is hard to maintain if mesh modifications
  // occur during advancement so it is cleared in 'compute()'.
  bool MAdElasticityOp::checkQualities(double ratio)
  {
    if ( dim == 3 )
      {
        pRegion pr;
        RIter rit = M_regionIter(mesh);
        while ( ( pr = RIter_next(rit) ) ) {
          if ( !checkRQuality(pr,ratio) ) {
            RIter_delete(rit);
            return false;
          }
        }
        RIter_delete(rit);
      }
    else
      {
        pFace pf;
        FIter fit = M_faceIter(mesh);
        while ( ( pf = FIter_next(fit) ) ) {
          if ( !checkFQuality(pf,ratio) ) {
            FIter_delete(fit);
            return false;
          }
        }
        FIter_delete(fit);
      }
    return true;
  }

  // -------------------------------------------------------------------
  void MAdElasticityOp::relocate(double ratio)
  {
    pVertex pv;
    double xyz[3];
    
    std::map<pVertex,smallVector>::const_iterator iter = relocations.begin();
    for(; iter != relocations.end(); iter++)
      {
        pv = (*iter).first;
        V_coord(pv,xyz);

        smallVector dX = (*iter).second;

        xyz[0] += ratio * dX(0);
        xyz[1] += ratio * dX(1);
        xyz[2] += ratio * dX(2);
        if ( !V_setPosition(pv,xyz) ) {
          MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                      "Could not relocate a node");
          
        }
      }
  }

  // -------------------------------------------------------------------
  // Advances the relocation as far as possible
  // Returns:
  //    - 0: no relocation possible
  //    - 1: advanced but not to the final position
  //    - 2: reached the full prescribed relocation
  int MAdElasticityOp::advance(double * ratio, double tolerance)
  {
    if ( *ratio >= 1.-tolerance ) {
      *ratio = 1.;
      return 2;
    }

    int achieved = 0; // the result returned

    int subIter = 0;
    const int nbMaxIter = 6;

    double margin = 1. - *ratio; // in [0;1] : the rest of the motion to be performed
    double curRatio = 1.;        // in [0;1] : the part of the margin that will be tested next
    double delta = 0.5;          // in [0;1] : the increment/decrement of ratio at next iteration
    double best = -1.;           // in [0;1] : the best current ratio with positive volumes

    while ( subIter < nbMaxIter )
      {
        double totalAdvancement = curRatio * margin + *ratio;
        bool check = checkQualities( curRatio * margin );
//         bool check;
//         if ( dim == 3 ) check = checkVolumes( curRatio * margin );
//         else            check = checkAreas  ( curRatio * margin );
        if ( !check ) {
          curRatio -= delta;
        }
        else {
          best = curRatio;
          if ( totalAdvancement + tolerance >= 1. ) { achieved = 2; break; }
          else                                      { achieved = 1; curRatio += delta; }
        }

        subIter++;
        delta *= 0.5;
//         printf("Ratio: %f, margin: %f, achieved: %d, subIter: %d, best: %f, curRatio: %f, delta: %f\n",
//                *ratio, margin, achieved, subIter, best, curRatio, delta);
      }

    double moveRatio = 0.;
    if ( achieved == 2 ) {
      moveRatio = ( 1. - *ratio );
      relocate(moveRatio);
      *ratio = 1.;
    }
    else if ( achieved == 1 ) {
      moveRatio = best * margin;
      relocate(moveRatio);
      *ratio += moveRatio;
    }
    
    return achieved;
  }

  // -------------------------------------------------------------------
  // Force the relocation even if it leads to negative volumes
  void MAdElasticityOp::forceRelocation()
  {
    relocate(1.);
  }

  // -------------------------------------------------------------------
  void MAdElasticityOp::elementMatrix(pEntity pe, smallMatrix& m) const
  {
    int nbNodes;
    double xyz[4][3];
    if ( dim == 3 ) {
        R_coordP1((pRegion)pe,xyz);
        nbNodes = 4;
    }
    else {
      F_coordP1((pFace)pe,xyz);
      nbNodes = 3;
    }

    double FACT = E / (1 + nu);
    double C11 = FACT * (1 - nu) / (1 - 2 * nu);
    double C12 = FACT * nu / (1 - 2 * nu);
    double C44 = (C11 - C12) / 2;
    const double C[6][6] =
      { {C11, C12, C12,    0,   0,   0}, 
        {C12, C11, C12,    0,   0,   0}, 
        {C12, C12, C11,    0,   0,   0}, 
        {  0,   0,   0,  C44,   0,   0},
        {  0,   0,   0,    0, C44,   0}, 
        {  0,   0,   0,    0,   0, C44} };

    smallMatrix H   (6,6);
    for (int i=0;i<6;i++)
      for (int j=0;j<6;j++)
        H(i,j) = C[i][j];

    smallMatrix B   (6,3*nbNodes);
    smallMatrix BTH (3*nbNodes,6);
    smallMatrix BT  (3*nbNodes,6);

    double gradShpFct[4][3] =
      { {-1., -1., -1.},
        { 1.,  0.,  0.},
        { 0.,  1.,  0.},
        { 0.,  0.,  1.} };

    if ( dim < 3 ) gradShpFct[0][2] = 0.;

    double jac [3][3];
    jac[0][0] = jac[0][1] = jac[0][2] = 0.;
    jac[1][0] = jac[1][1] = jac[1][2] = 0.;
    jac[2][0] = jac[2][1] = jac[2][2] = 0.;
    for(int i = 0; i < nbNodes; i++) {
      jac[0][0] += xyz[i][0] * gradShpFct[i][0]; jac[0][1] += xyz[i][1] * gradShpFct[i][0]; jac[0][2] += xyz[i][2] * gradShpFct[i][0];
      jac[1][0] += xyz[i][0] * gradShpFct[i][1]; jac[1][1] += xyz[i][1] * gradShpFct[i][1]; jac[1][2] += xyz[i][2] * gradShpFct[i][1];
      jac[2][0] += xyz[i][0] * gradShpFct[i][2]; jac[2][1] += xyz[i][1] * gradShpFct[i][2]; jac[2][2] += xyz[i][2] * gradShpFct[i][2];
    }

    if ( dim < 3 ) {
      for (int i=0;i<2;i++) { jac[i][2] = jac[2][i] = 0.; }
      jac[2][2] = 1.;
    }

    double detJ;
    double invjac [3][3];
    if ( dim == 3 ) detJ = fabs ( inverseMat (jac, invjac) );
    else {
      double itmp[2][2]; 
      double tmp[2][2]; for (int i=0;i<2;i++) for (int j=0;j<2;j++) tmp[i][j]=jac[i][j];
      detJ = fabs ( inverseMat22(tmp,itmp) );
      for (int i=0;i<2;i++) for (int j=0;j<2;j++) invjac[i][j]=itmp[i][j];
      for (int i=0;i<2;i++) { invjac[i][2] = invjac[2][i] = 0.; }
      invjac[2][2] = 1.;
    }

    double Grads[4][3];

    B.set_all(0.0);
    BT.set_all(0.0);

    for (int j=0;j<nbNodes;j++){

      Grads[j][0] = invjac[0][0] * gradShpFct[j][0] + invjac[0][1] * gradShpFct[j][1] + invjac[0][2] * gradShpFct[j][2];
      Grads[j][1] = invjac[1][0] * gradShpFct[j][0] + invjac[1][1] * gradShpFct[j][1] + invjac[1][2] * gradShpFct[j][2];
      Grads[j][2] = invjac[2][0] * gradShpFct[j][0] + invjac[2][1] * gradShpFct[j][1] + invjac[2][2] * gradShpFct[j][2];
    
      BT(j,0) = B(0,j) = Grads[j][0];
      BT(j,3) = B(3,j) = Grads[j][1];
      BT(j,4) = B(4,j) = Grads[j][2];
    
      BT(j+nbNodes,1) = B(1,j+nbNodes) = Grads[j][1];
      BT(j+nbNodes,3) = B(3,j+nbNodes) = Grads[j][0];
      BT(j+nbNodes,5) = B(5,j+nbNodes) = Grads[j][2];
    
      BT(j+2*nbNodes,2) = B(2,j+2*nbNodes) = Grads[j][2];
      BT(j+2*nbNodes,4) = B(4,j+2*nbNodes) = Grads[j][0];
      BT(j+2*nbNodes,5) = B(5,j+2*nbNodes) = Grads[j][1];
    }
  
    BTH.set_all(0.0);
    BTH.blas_dgemm (BT,H); 
  
    m.set_all(0.0);
    m.blas_dgemm (BTH,B,detJ,1.0);

    // --- Selective element stiffening ---
    // reference:
    // Stein, Tezduyar, Benney "Mesh Moving Techniques for FSI with Large Displacement", ASME, 2003.

    double detJ0 = MeshQualityManagerSgl::instance().getMaxSize();
    if ( detJ0 <= 0. ) detJ0 = 1.;

    double scaleFactor = pow( detJ0 / detJ, chi);
    if ( scaleFactor > 1.e12 || scaleFactor < 1.e-12 ) {
      MAdMsgSgl::instance().warning(__LINE__,__FILE__,
                                    "Scale factor for element stiffness is %e",
                                    scaleFactor);
    }

    m.scale ( scaleFactor );
  }

  // -------------------------------------------------------------------

}
