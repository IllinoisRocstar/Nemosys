// -------------------------------------------------------------------
// MAdLib - Copyright (C) 2008-2009 Universite catholique de Louvain
//
// See the Copyright.txt and License.txt files for license information. 
// You should have received a copy of these files along with MAdLib. 
// If not, see <http://www.madlib.be/license/>
//
// Please report all bugs and problems to <contrib@madlib.be>
//
// Authors: Jean-Francois Remacle, Gaetan Compere, Cecile Dobrzynski
// -------------------------------------------------------------------

#include "MeshDataBase.h"
#include "MeshDataBaseInterface.h"
#include "MeshDataBaseIO.h"
#include "MeshDataBaseGEntity2Physical.h"
#include <assert.h>
#include "MeshDataBaseCommPeriodic.h"
#include "MshTags.h"

#include <cstdlib>
#include <cstring>
#include <set>
#include <map>
#include <utility>

#ifdef PARALLEL
#include "mpi.h"
#include "autopack.h"
#endif

#include "MeshDataBaseParallelInterface.h"

namespace MAd {

  // -------------------------------------------------------------------
  // ------------- Load a Gmsh mesh ------------------------------------
  // -------------------------------------------------------------------

  void LoadGmshMesh (pMesh m,const char *filename)
  {  
    FILE *fp = fopen(filename, "r");
    if(!fp) {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                  "Unknown File %s",filename);
    }

#ifdef PARALLEL 
    int myRank = 0, nbProc = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &nbProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
#endif

    char String[256];

    if (!m->model) {
      pGModel model = NULL;
      GM_create(&model,"");
      m->model = model; 
    }

    // ----------------------------------------------
    // ------ Reading format
    // ----------------------------------------------

    int version = -1;
    do {
      if(!fgets(String, sizeof(String),fp))
        break;
      if(feof(fp))
        break;
    } while(String[0] != '$' ||
            strncmp(&String[1], "MeshFormat",10));
  
    if(feof(fp)) {
      version = 1;
    }
    else {
      int a, b, c;
      fscanf(fp, "%d %d %d", &a, &b, &c);
      version = a;
    }

    rewind(fp);

#ifdef PARALLEL
    if ( version == 1 && nbProc > 1 ) {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                  "Cannot load Gmsh mesh from file %s: parallel loading not implemented for file format msh1",filename);
    }
#endif

    // -----------------------------------------------
    // ------ Reading elements to find local nodes
    // -----------------------------------------------

    std::set<int> localNodeIds;

    do {
      if(!fgets(String, sizeof(String), fp))
        break;
      if(feof(fp))
        break;
    } while(String[0] != '$' ||
            ( strncmp(&String[1], "ELM",      3) &&
              strncmp(&String[1], "Elements", 8)    ) );

    if(!feof(fp)) {

      int nbElems_world, NbTags, verts[256], Tag;
      int Num, Type, Physical, Elementary, Nbr_Nodes, Partition;
      fscanf(fp, "%d", &nbElems_world);

      for(int i_Element = 0; i_Element < nbElems_world; i_Element++) {
	
        if(version == 1){
          fscanf(fp, "%d %d %d %d %d",
                 &Num, &Type, &Physical, &Elementary, &Nbr_Nodes);
          Partition = -1;
        }
        else{
          fscanf(fp, "%d %d %d", &Num, &Type, &NbTags);
          Elementary = Physical = Partition = -1;
          for(int j = 0; j < NbTags; j++){
            fscanf(fp, "%d", &Tag);	    
            if(j == 0)
              Physical = Tag;
            else if(j == 1)
              Elementary = Tag;
            else if(j == 2)
              Partition = Tag - 1;
            // ignore any other tags for now
          }
          Nbr_Nodes = getNumVerticesForElementTypeMSH(Type);
        }

        for(int j = 0; j < Nbr_Nodes; j++) {
          fscanf(fp, "%d", &verts[j]);
        }

#ifdef PARALLEL 
        if ( nbProc == 1 || Partition == myRank ) {
#endif
          for(int j = 0; j < Nbr_Nodes; j++) {
            localNodeIds.insert(verts[j]);
          }
#ifdef PARALLEL
        }
#endif
      }
    }

    rewind(fp);


    // ----------------------------------------------
    // ------ Reading nodes
    // ----------------------------------------------

    do {
      if(!fgets(String, sizeof(String), fp))
        break;
      if(feof(fp))
        break;
    } while(String[0] != '$' ||
            (strncmp(&String[1], "NOD", 3) &&
             strncmp(&String[1], "NOE", 3) &&
             strncmp(&String[1], "Nodes", 5)  ) );
    
    if(!feof(fp)) {
      int nbNodes_world;
      fscanf(fp, "%d", &nbNodes_world);
      for(int i_Node = 0; i_Node < nbNodes_world; i_Node++) {
        int id;
        double x,y,z;
        fscanf(fp, "%d %lf %lf %lf", &id, &x, &y, &z);
        if ( localNodeIds.find(id) != localNodeIds.end() ) {
          m->add_point(id,x,y,z);
        }
      }
    }

    rewind(fp);

    // ----------------------------------------------
    // ------ Reading parametric nodes
    // ----------------------------------------------

    if ( version == 2 ) {
      do {
        if(!fgets(String, sizeof(String), fp)) break;
        if(feof(fp)) break;
      } while(String[0] != '$' ||
              strncmp(&String[1], "ParametricNodes", 15)  );
    
      if(!feof(fp)) {
        int nbNodes_world;
        fscanf(fp, "%d", &nbNodes_world);
        for(int i_Node = 0; i_Node < nbNodes_world; i_Node++) {
          int id, gDim, gTag;
          double x,y,z;
          fscanf(fp, "%d %lf %lf %lf %d %d", &id, &x, &y, &z, &gDim, &gTag);
          if ( gDim == 3 ) {
            if ( localNodeIds.find(id) != localNodeIds.end() ) {
              m->add_point(id,x,y,z);
            }
          }
          else {
            double u,v;
            if ( gDim == 0 ) {
              u = v = 0.;
            }
            if ( gDim == 1 ) {
              fscanf(fp, "%lf", &u);
              v = 0.;
            } 
            if ( gDim == 2 ) {
              fscanf(fp, "%lf %lf", &u, &v);
            }
            m->add_pointParam(id,x,y,z,u,v);
          }
        }
      }

      rewind(fp);
    }

    // ----------------------------------------------
    // ------ Reading elements
    // ----------------------------------------------

    do {
      if(!fgets(String, sizeof(String), fp))
        break;
      if(feof(fp))
        break;
    } while(String[0] != '$' ||
            ( strncmp(&String[1], "ELM",      3) &&
              strncmp(&String[1], "Elements", 8)    ) );

    if(!feof(fp)) {

      int nbElems_world, NbTags, verts[256], Tag;
      int Num, Type, Physical, Elementary, Nbr_Nodes, Partition;
      fscanf(fp, "%d", &nbElems_world);

      for(int i_Element = 0; i_Element < nbElems_world; i_Element++) {
	
        if(version == 1){
          fscanf(fp, "%d %d %d %d %d",
                 &Num, &Type, &Physical, &Elementary, &Nbr_Nodes);
          Partition = -1;
        }
        else{
          fscanf(fp, "%d %d %d", &Num, &Type, &NbTags);
          Elementary = Physical = Partition = -1;
          for(int j = 0; j < NbTags; j++){
            fscanf(fp, "%d", &Tag);
            if(j == 0)
              Physical = Tag;
            else if(j == 1)
              Elementary = Tag;
            else if(j == 2)
              Partition = Tag - 1;
            // ignore any other tags for now
          }
          Nbr_Nodes = getNumVerticesForElementTypeMSH(Type);
        }

        for(int j = 0; j < Nbr_Nodes; j++) {
          fscanf(fp, "%d", &verts[j]);
        }

#ifdef PARALLEL 
        if ( nbProc == 1 || Partition == myRank ) {
#endif

          pGEntity geom = 0;
	
          switch (Type) {
          case MSH_LIN_2:
            {
              geom = (pGEntity) GM_edgeByTag(m->model,Elementary);
              if ( Physical > 0 ) GEN_setPhysical(geom,GEN_type(geom),Physical);
              m->add_edge(verts[0],verts[1],geom);
            }
            break;
          case MSH_LIN_3:
            {
              geom = (pGEntity) GM_edgeByTag(m->model,Elementary);
              if ( Physical > 0 ) GEN_setPhysical(geom,GEN_type(geom),Physical);
              //m->add_edge(verts[0],verts[2],verts[1],geom);
              m->add_edge(3,geom,verts[0],verts[2],verts[1]);
              //m->add_edge(verts[0],verts[1],ge);
            }
            break;
          case MSH_LIN_4:
            {
              geom = (pGEntity) GM_edgeByTag(m->model,Elementary);
              if ( Physical > 0 ) GEN_setPhysical(geom,GEN_type(geom),Physical);
              m->add_edge(4,geom, verts[0],verts[2],verts[3],verts[1]);
            }
            break;
          case MSH_LIN_5:
            {
              geom = (pGEntity) GM_edgeByTag(m->model,Elementary);
              if ( Physical > 0 ) GEN_setPhysical(geom,GEN_type(geom),Physical);
              m->add_edge(5,geom, verts[0],verts[2],verts[3],verts[4],verts[1]);
            }
            break;
          case MSH_LIN_6:
            {
              geom = (pGEntity) GM_edgeByTag(m->model,Elementary);
              if ( Physical > 0 ) GEN_setPhysical(geom,GEN_type(geom),Physical);
              m->add_edge(6,geom, verts[0],verts[2],verts[3],verts[4],verts[5],verts[1]);
            }
            break;
          case MSH_TRI_3:
            {
              geom = (pGEntity) GM_faceByTag(m->model,Elementary);
              if ( Physical > 0 ) GEN_setPhysical(geom,GEN_type(geom),Physical);
              m->add_triangle(verts[0],verts[1],verts[2],geom); 
            }
            break;
            // 6 nodes triangle
          case MSH_TRI_6:
            {
              geom = (pGEntity) GM_faceByTag(m->model,Elementary);
              if ( Physical > 0 ) GEN_setPhysical(geom,GEN_type(geom),Physical);
              //	    printf ("adding a second order triangle %d %d %d %d %d %d\n",verts[0],verts[1],verts[2],verts[3],verts[4],verts[5]);
              m->add_triangle(2,0,geom,verts[0],verts[3],verts[1],verts[4],verts[2],verts[5]); 
              //m->add_triangle(verts[0],verts[3],verts[1],verts[4],verts[2],verts[5],geom); 
              //m->add_triangle(verts[0],verts[1],verts[2],gf); 
            }
            break;
            // 9 nodes triangle (SERENDIP !)
          case MSH_TRI_9:
            {
              geom = (pGEntity) GM_faceByTag(m->model,Elementary);
              if ( Physical > 0 ) GEN_setPhysical(geom,GEN_type(geom),Physical);
              m->add_triangle(3,0,geom,verts[0],verts[3],verts[4],verts[1],verts[5],verts[6],verts[2],verts[7],verts[8]); 
              //m->add_triangle(verts[0],verts[3],verts[1],verts[4],verts[2],verts[5],geom); 
            }
            break;
            // 10 nodes triangle (LAGRANGE !)
          case MSH_TRI_10:
            {
              geom = (pGEntity) GM_faceByTag(m->model,Elementary);
              if ( Physical > 0 ) GEN_setPhysical(geom,GEN_type(geom),Physical);
              m->add_triangle(3,1,geom,verts[0],verts[3],verts[4],verts[1],verts[5],verts[6],verts[2],verts[7],verts[8],verts[9]); 
            }
            break;
            // 12 nodes triangle (SERENDIP !)
          case MSH_TRI_12:
            {
              geom = (pGEntity) GM_faceByTag(m->model,Elementary);
              if ( Physical > 0 ) GEN_setPhysical(geom,GEN_type(geom),Physical);
              m->add_triangle(4,0,geom,verts[0],verts[3],verts[4],verts[5],verts[1],verts[6],verts[7],verts[8],verts[2],verts[9],verts[10],verts[11]); 
            }
            break;
            // 15 nodes triangle (LAGRANGE !)
          case MSH_TRI_15:
            {
              geom = (pGEntity) GM_faceByTag(m->model,Elementary);
              if ( Physical > 0 ) GEN_setPhysical(geom,GEN_type(geom),Physical);
              m->add_triangle(4,1,geom,verts[0],verts[3],verts[4],verts[5],verts[1],verts[6],verts[7],verts[8],verts[2],verts[9],verts[10],verts[11],verts[12],verts[13],verts[14]); 
            }
            break;
            // 15 nodes triangle (SERENDIP !)
          case MSH_TRI_15I:
            {
              geom = (pGEntity) GM_faceByTag(m->model,Elementary);
              if ( Physical > 0 ) GEN_setPhysical(geom,GEN_type(geom),Physical);
              m->add_triangle(5,0,geom,verts[0],verts[3],verts[4],verts[5],verts[6],verts[1],verts[7],verts[8],verts[9],verts[10],verts[2],verts[11],verts[12],verts[13],verts[14]); 
            }
            break;
          case MSH_QUA_4:
            {
              geom = (pGEntity) GM_faceByTag(m->model,Elementary);
              if ( Physical > 0 ) GEN_setPhysical(geom,GEN_type(geom),Physical);
              m->add_quad(verts[0],verts[1],verts[2],verts[3],geom); 
            }
            break;
            // 9 nodes quad (LAGRANGE)
          case MSH_QUA_9:
            {
              geom = (pGEntity) GM_faceByTag(m->model,Elementary);
              if ( Physical > 0 ) GEN_setPhysical(geom,GEN_type(geom),Physical);
              m->add_quad(2,1,geom,verts[0],verts[4],verts[1],verts[5],verts[2],verts[6],verts[3],verts[7],verts[8]); 
            }
            break;
          case MSH_TET_4:
            {
              geom = (pGEntity) GM_regionByTag(m->model,Elementary);
              if( !geom ) 
                {
                  printf("READ MESH : no classification for  Region\n");
                }
              if ( Physical > 0 ) GEN_setPhysical(geom,GEN_type(geom),Physical);
              m->add_tet(verts[0],verts[1],verts[2],verts[3],geom); 
            }
            break;
          case MSH_PRI_6:
            {
              geom = (pGEntity) GM_regionByTag(m->model,Elementary);
              if ( Physical > 0 ) GEN_setPhysical(geom,GEN_type(geom),Physical);
              m->add_prism(verts[0],verts[1],verts[2],
                           verts[3],verts[4],verts[5],geom); 
            }
            break;
          case MSH_PNT:
            {
              geom = (pGEntity) GM_vertexByTag(m->model,Elementary);
              if ( Physical > 0 ) GEN_setPhysical(geom,GEN_type(geom),Physical);
              MDB_Point *p = m->find_point(verts[0]);
              p->g = geom;
            }
            break;
            // 8 nodes hexahedron
          case MSH_HEX_8:
            {
              geom = (pGEntity) GM_regionByTag(m->model,Elementary);
              if ( Physical > 0 ) GEN_setPhysical(geom,GEN_type(geom),Physical); 
              m->add_hex(verts[0],verts[1],verts[2],verts[3],verts[4],verts[5],verts[6],verts[7],geom); 
            }
            break;
            // 10 node tetrahedron, 2nd order
          case MSH_TET_10:
            {
              geom = (pGEntity) GM_regionByTag(m->model,Elementary);
              if ( Physical > 0 ) GEN_setPhysical(geom,GEN_type(geom),Physical);
              m->add_tet(geom,2,false,verts);
              break;
            }
          case MSH_TET_20:
            {
              geom = (pGEntity) GM_regionByTag(m->model,Elementary);
              if ( Physical > 0 ) GEN_setPhysical(geom,GEN_type(geom),Physical);
              m->add_tet(geom,3,false,verts);
              break;
            }
          case MSH_TET_35:
            {
              geom = (pGEntity) GM_regionByTag(m->model,Elementary);
              if ( Physical > 0 ) GEN_setPhysical(geom,GEN_type(geom),Physical);
              m->add_tet(geom,4,false,verts);
              break;
            }
          case MSH_TET_34:
            {
              geom = (pGEntity) GM_regionByTag(m->model,Elementary);
              if ( Physical > 0 ) GEN_setPhysical(geom,GEN_type(geom),Physical);
              m->add_tet(geom,4,true,verts);
              break;
            }
          default:
            throw;
            break;
          }
          if (geom)
            {
              bool find = false;
              for (std::multimap<int, pGEntity>::iterator it = m->geomFeatures_Tags.lower_bound(Elementary);
                   it != m->geomFeatures_Tags.upper_bound(Elementary);++it)
                if (it->second == geom)find = true;
              if (!find)
                m->geomFeatures_Tags.insert(std::pair<int,pGEntity>(Elementary, geom));
            }
#ifdef PARALLEL
        }
#endif
      }
    }

    m->classify_unclassified_entities();

    m->destroyStandAloneEntities();

    rewind(fp);


    // ----------------------------------------------
    // ------ Reading periodicity properties
    // ----------------------------------------------

    int isPeriodic = 0;
    while(1) {
      do {
        if(!fgets(String, sizeof(String), fp))
          break;
        if(feof(fp))
          break;
      } while(String[0] != '$');
    
      if(feof(fp))
        break;

      else if(!strncmp(&String[1], "PERIODICNOD", 11)) {
        isPeriodic = 1;
      } 
      else if(!strncmp(&String[1], "TRA", 3)) {
        isPeriodic = 2;
      }
      else if (!strncmp(&String[1],"PeriodicNodes",13)) {
        isPeriodic = 2;
      }
      else if (!strncmp(&String[1],"PER", 3)) {
        isPeriodic = 2;
      }

      do {
        if(!fgets(String, sizeof(String), fp))
          throw;
        if(feof(fp))
          throw;
      } while(String[0] != '$');
    }

    rewind(fp);

    // ----------------------------------------------
    // ------ Reading periodic data
    // ----------------------------------------------

    if(isPeriodic) {
      pMeshDataId tagPeriodic = MD_newMeshDataId("PeriodicPoint");
      if(isPeriodic==1) {
        fseek(fp,0,SEEK_SET);
        while(1) {
          do {
            if(!fgets(String, sizeof(String), fp))
              break;
            if(feof(fp))
              break;
          } while(String[0] != '$');
    
          if(feof(fp))
            break;
          if(!strncmp(&String[1], "PERIODICNOD", 11)) {
            int nPeriod,nRefPeriod;
            fscanf(fp, "%d %d", &nPeriod,&nRefPeriod);
            printf("%d Ref Periodic\n",nRefPeriod);
            int *t=new int[nRefPeriod];
            for(int ip = 0; ip < nPeriod; ip++) {
              int Num,v1,v2;
              fscanf(fp, "%d ", &Num);
              std::vector<int> transfo;
              std::vector<int> invtransfo;
              for(int k=0 ; k<nRefPeriod ; k++) {
                fscanf(fp, "%d ", &t[k]);
                transfo.push_back(t[k]);
                invtransfo.push_back(-t[k]);
              } 
              fscanf(fp, "%d %d",&v1, &v2);
              pVertex vt1 = m->find_point(v1);
              pVertex vt2 = m->find_point(v2);
              //printf("point %d %d\n",EN_id((pEntity) vt1),EN_id((pEntity) vt2));
              void *temp_ptr; 
              int isPeriodic = EN_getDataPtr((pEntity) vt1 , tagPeriodic, &temp_ptr);
              if(isPeriodic) {
                std::vector<std::pair<std::vector<int> , pVertex> > *recup = (std::vector<std::pair<std::vector<int> , pVertex> > *) temp_ptr;
                unsigned int i=0;
                for(i=0 ; i<(*recup).size() ; i++) {
                  unsigned int j=0;
                  for(j=0 ; j<(*recup)[i].first.size() ; j++) {
                    if((*recup)[i].first[j] != transfo[j]) break;
                  }
                  if(j==(*recup)[i].first.size() ) break;
                }
                if(i==(*recup).size())
                  (*recup).push_back(std::pair <std::vector<int> , pVertex>(transfo,vt2));
              } else {
                std::vector<std::pair<std::vector<int> , MDB_Point*> > remotePoints;
                (remotePoints).push_back(std::pair <std::vector<int> , pVertex>(transfo,vt2));
                EN_attachDataPtr((pEntity) vt1 , tagPeriodic, 
                                 new std::vector<std::pair<std::vector<int> , pVertex> >(remotePoints));
				   
              }
	   
              isPeriodic = EN_getDataPtr((pEntity) vt2 , tagPeriodic, &temp_ptr);
              if(isPeriodic) {
                std::vector<std::pair<std::vector<int> , pVertex> > *recup = (std::vector<std::pair<std::vector<int> , pVertex> > *) temp_ptr;
                unsigned int i;
                for(i=0 ; i<(*recup).size() ; i++) {
                  unsigned int j;
                  for(j=0 ; j<(*recup)[i].first.size() ; j++) {
                    if((*recup)[i].first[j] != invtransfo[j]) break;
                  }
                  if(j==(*recup)[i].first.size() ) break;
                }
                if(i==(*recup).size())
                  (*recup).push_back(std::pair <std::vector<int> , pVertex>(invtransfo,vt1));
              } else {
                std::vector<std::pair<std::vector<int> , MDB_Point*> > remotePoints;
                (remotePoints).push_back(std::pair <std::vector<int> , pVertex>(invtransfo,vt1));
                EN_attachDataPtr((pEntity) vt2 , tagPeriodic, 
                                 new std::vector<std::pair<std::vector<int> , pVertex> >(remotePoints));	   
              }
            }
            delete []t;
          }  
        }
      } else {

        assert(isPeriodic==2);
        puts("format periodic CENAERO");
        fseek(fp,0,SEEK_SET);
        while(1) {
          do {
            if(!fgets(String, sizeof(String), fp))
              break;
            if(feof(fp))
              break;
          } while(String[0] != '$');
    
          if(feof(fp))
            break;

          if((!strncmp(&String[1],"PER", 3)) || 
             (!strncmp(&String[1],"PeriodicNodes",13))) {
            
            int nPeriod;
            fscanf(fp, "%d ", &nPeriod);

            pMeshDataId tagPeriodic = MD_newMeshDataId("PeriodicVertexID");
            
            for(int ip = 0; ip < nPeriod; ip++) {
              int Num,vtxId;
              fscanf(fp, "%d ", &Num);
              
              std::vector<int> vtcs;

              for (int i=0;i<Num;i++) {
                fscanf(fp, "%d ", &vtxId);
                vtcs.push_back(vtxId);
              }
              
              int trafo;
              fscanf(fp, "%d ",&trafo);
              
              std::vector<int>::const_iterator vIter = vtcs.begin();
              for (;vIter!=vtcs.end();vIter++) {
                
                int vtxId = *vIter;
                pVertex vt1 = m->find_point(vtxId);
                
                if (vt1) {
                  
                  void* tmp_periodic = NULL;
                  int havePeriodic = EN_getDataPtr((pEntity) vt1,tagPeriodic,&tmp_periodic);
                  std::map<int,std::vector<int> > *  connections = NULL;
                  
                  if (havePeriodic) connections = (std::map<int,std::vector<int> > * )  (tmp_periodic);
                  else {
                    connections = new std::map<int,std::vector<int> >;
                    EN_attachDataPtr((pEntity) vt1,tagPeriodic,connections);
                  }
                  
                  // transformations will be written later, using topological classification 
                  
                  std::vector<int>::const_iterator vIter2 = vtcs.begin();
                  for (;vIter2!=vtcs.end();++vIter2) {
                    if (vtxId != *vIter2) connections->insert(std::make_pair(*vIter2,std::vector<int>()));
                  }
                }
              }
            }
          }
          
//           if(!strncmp(&String[1], "PER", 3)) {
//             int nPeriod,nRefPeriod;
//             nRefPeriod = 1;
//             fscanf(fp, "%d ", &nPeriod);
//             printf("%d Ref Periodic\n",nRefPeriod);
//             int *t=new int[nRefPeriod];
//             for(int ip = 0; ip < nPeriod; ip++) {
//               int Num,v1,v2;
//               fscanf(fp, "%d ", &Num);
//               std::vector<int> transfo;
//               std::vector<int> invtransfo;
//               fscanf(fp, "%d %d %d",&v1, &v2,&Num);
//               pVertex vt1 = m->find_point(v1);
//               pVertex vt2 = m->find_point(v2);


//               printf("point %d (%d) %d (%d)\n",EN_id((pEntity) vt1),GEN_tag(vt1->g),EN_id((pEntity) vt2),
//                      GEN_tag(vt2->g));
              
//               for(int k=0 ; k<nRefPeriod ; k++) {
//                 t[k]=-1;
//                 transfo.push_back(t[k]);
//                 invtransfo.push_back(-t[k]);
//               }
//               void *temp_ptr; 
//               int isPeriodic = EN_getDataPtr((pEntity) vt1 , tagPeriodic, &temp_ptr);
//               if(isPeriodic) {
//                 std::vector<std::pair<std::vector<int> , pVertex> > *recup = (std::vector<std::pair<std::vector<int> , pVertex> > *) temp_ptr;
//                 unsigned int i=0;
//                 for(i=0 ; i<(*recup).size() ; i++) {
//                   unsigned int j=0;
//                   for(j=0 ; j<(*recup)[i].first.size() ; j++) {
//                     if((*recup)[i].first[j] != transfo[j]) break;
//                   }
//                   if(j==(*recup)[i].first.size() ) break;
//                 }
//                 if(i==(*recup).size())
//                   (*recup).push_back(std::pair <std::vector<int> , pVertex>(transfo,vt2));
//               } else {
//                 std::vector<std::pair<std::vector<int> , MDB_Point*> > remotePoints;
//                 (remotePoints).push_back(std::pair <std::vector<int> , pVertex>(transfo,vt2));
//                 EN_attachDataPtr((pEntity) vt1 , tagPeriodic, 
//                                  new std::vector<std::pair<std::vector<int> , pVertex> >(remotePoints));
				   
//               }
	   
//               isPeriodic = EN_getDataPtr((pEntity) vt2 , tagPeriodic, &temp_ptr);
//               if(isPeriodic) {
//                 std::vector<std::pair<std::vector<int> , pVertex> > *recup = (std::vector<std::pair<std::vector<int> , pVertex> > *) temp_ptr;
//                 unsigned int i;
//                 for(i=0 ; i<(*recup).size() ; i++) {
//                   unsigned int j;
//                   for(j=0 ; j<(*recup)[i].first.size() ; j++) {
//                     if((*recup)[i].first[j] != invtransfo[j]) break;
//                   }
//                   if(j==(*recup)[i].first.size() ) break;
//                 }
//                 if(i==(*recup).size())
//                   (*recup).push_back(std::pair <std::vector<int> , pVertex>(invtransfo,vt1));
//               } else {
//                 std::vector<std::pair<std::vector<int> , MDB_Point*> > remotePoints;
//                 (remotePoints).push_back(std::pair <std::vector<int> , pVertex>(invtransfo,vt1));
//                 EN_attachDataPtr((pEntity) vt2 , tagPeriodic, 
//                                  new std::vector<std::pair<std::vector<int> , pVertex> >(remotePoints));	   
//               }
//             }
//             delete []t;
//           }
        }
      }  
    }

    fclose(fp);

    
// #ifdef PARALLEL

    // ----------------------------------------------
    // ------ Tagging inter-partition nodes
    // ----------------------------------------------

    pMeshDataId tag = MD_lookupMeshDataId("RemotePoint");
    
    V_createInfoInterface(m,tag);
    E_createInfoInterface(m,tag);
    F_createInfoInterface(m,tag);
  
// #endif

    m->initializeIdData();
  }


  // -------------------------------------------------------------------
  // ------------- Save a Gmsh mesh ------------------------------------
  // -------------------------------------------------------------------

  void SaveGmshMesh (const pMesh m,const char *filename, int version, 
                     bool saveAll, const int * partitionTable)
  {
    if ( version == 1 && partitionTable ) {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                  "Error while saving mesh on file %s: SaveGmshMesh not implemented for msh1 format with partitioning",filename);
    }

    int meshDim = M_dim(m);
    int nbModelVertex = 0;
    pMeshDataId tagPeriodic = MD_lookupMeshDataId("PeriodicPoint");
    // Create a reverse map of physical tags
    GEntity2Physical gentity2phys(m->geomFeatures_Tags);
    FILE *f = fopen(filename, "w");
    if ( !f ) {
      printf("Error: could not open file %s\n",filename);
      throw;
    }
    if ( version != 1 ) {
      fprintf(f, "$MeshFormat\n2.1 0 8\n$EndMeshFormat\n");
    }
    // ------ Writting the nodes ------
    if ( version == 1 || !(m->isParametric()) )
      {
        if ( version == 1 ) { fprintf(f, "$NOD\n");   } 
        else                { fprintf(f, "$Nodes\n"); }
        fprintf(f, "%d\n", m->nbPoints);
        VIter vit = M_vertexIter(m); 
        while (VIter_next(vit)){}
        VIter_reset(vit);
        pVertex pv;  
        int NN = 0;
        while ((pv = VIter_next(vit)))
          {
            NN++;
            if(pv->g)
              {
                int dim = GEN_type(pv->g); 
                if(dim == 0)  nbModelVertex++;
                fprintf(f, "%d %.16g %.16g %.16g\n", pv->iD, pv->X, pv->Y, pv->Z);
              }
            else
              {
                fprintf(f, "%d %.16g %.16g %.16g\n", pv->iD, pv->X, pv->Y, pv->Z);
              }
          }
        if (NN != m->nbPoints) {
          MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                      "Vertex iterator saw %d vertices while there are %d points in the mesh",
                                      NN,m->nbPoints);
        }
        VIter_delete(vit);
        if ( version == 1 ) { fprintf(f, "$ENDNOD\n");   } 
        else                { fprintf(f, "$EndNodes\n"); }
      }
    // ------ Writting the parametric nodes ------
    else
      {
        fprintf(f, "$ParametricNodes\n");
        fprintf(f, "%d\n", m->nbPoints);
        pVertex pv;
        VIter vit = M_vertexIter(m);
        while ((pv = VIter_next(vit)))
          {
            if(pv->g)
              {
                int gdim = GEN_type(pv->g);
                int gtag = GEN_tag(pv->g);
                fprintf(f, "%d %.16g %.16g %.16g %d %d", pv->iD, pv->X, pv->Y, pv->Z, gdim, gtag);
                double u[2]; V_params(pv,&u[0],&u[1]);
                switch (gdim) {
                case 0: { nbModelVertex++; break; }
                case 1: { fprintf(f, " %.16g", u[0]); break; }
                case 2: { fprintf(f, " %.16g %.16g", u[0], u[1]); break; }
                case 3: { break; }
                }
                fprintf(f, "\n");
              }
            else {
              MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                          "Parametric point without classification");
            }
          }
        VIter_delete(vit);
        fprintf(f, "$EndParametricNodes\n");
      }
    // ------ Writting the elements ------
    {
      if ( version == 1 ) { fprintf(f, "$ELM\n");   } 
      else                { fprintf(f, "$Elements\n"); }
    
      bool physic = GM_physical(m->model);

      long int nbClasVerts = 0;
      long int nbClasEdges = 0;
      long int nbClasFaces = 0;
      long int nbClasRegs = 0;
      {
        VIter vit = M_vertexIter(m);
        pVertex pv;
        while ((pv = VIter_next(vit)))
          {
            if( GEN_type(pv->g) == 0 ){
              if ( saveAll ) { nbClasVerts++; continue; }
              if ( physic ) {
                if ( GEN_physTag(pv->g) != 0 ) nbClasVerts++;
              }
              else nbClasVerts++;
            }
          }
        VIter_delete(vit);
      } 
      {
        EIter eit = M_edgeIter(m);
        pEdge pe;  
        while ((pe = EIter_next(eit)))
          {
            if( GEN_type(pe->g) == 1 ){
              if ( saveAll ) { nbClasVerts++; continue; }
              if ( physic ) {
                if ( GEN_physTag(pe->g) != 0 ) nbClasEdges++;
              }
              else nbClasEdges++;
            }
          }
        EIter_delete(eit);
      }  
      {
        FIter fit = M_faceIter(m);
        pFace pf;  
        while ((pf = FIter_next(fit)))
          {
#ifdef PARALLEL
            if ( F_isInterface(pf) ) { nbClasFaces++; continue; }
#endif
            if( GEN_type(pf->g) == 2 ){
              if ( saveAll ) { nbClasVerts++; continue; }
              if ( physic ) {
                if ( GEN_physTag(pf->g) != 0 ) nbClasFaces++;
              }
              else nbClasFaces++;
            }
          }
        FIter_delete(fit);
      }
      {
        RIter rit = M_regionIter(m);
        pRegion pr;
        while ((pr = RIter_next(rit)))
          {
            if( GEN_type(pr->g) == 3 ){
              if ( saveAll ) { nbClasVerts++; continue; }
              if ( physic ) {
                if ( GEN_physTag(pr->g) != 0 ) nbClasRegs++;
              }
              else nbClasRegs++;
            }
          }
        RIter_delete(rit);
      } 

      fprintf(f, "%ld\n", (nbClasVerts + nbClasEdges + nbClasFaces + nbClasRegs) );

      int k = 1;
      // --- Nodes ---
      {
        VIter vit = M_vertexIter(m);
        pVertex pv;  
        while ((pv = VIter_next(vit)))
          {
            if(pv->g)
              {
                int dim = GEN_type(pv->g);
                int nbTags = 3;
                int tag = GEN_tag (pv->g); 
                int phys = 0;
                if ( GEN_type(pv->g) == 0 ) phys = GEN_physTag(pv->g); 
                if ( !saveAll && physic && phys == 0 ) continue;
//                 int phys = gentity2phys.get_first_tag(pv->g);
                int partition = -1;
                if ( partitionTable ) {
                  switch(meshDim) {
                  case 0:
                    partition = partitionTable[pv->iD];
                    break;
                  case 1:
                    for (int i=0; i < V_numEdges(pv); i++) {
                      int ePart = partitionTable[V_edge(pv,i)->iD];
                      partition = std::min(ePart,partition);
                      if (partition < 0) partition = ePart;
                    }
                    break;
                  case 2: {
                    pPList vFaces = V_faces(pv);
                    void *tmp_ptr = 0;
                    pFace pf;
                    while( (pf = (pFace) PList_next(vFaces,&tmp_ptr)) ) {
                      int fPart = partitionTable[pf->iD];
                      partition = std::min(fPart,partition);
                      if (partition < 0) partition = fPart;
                    }
                    PList_delete(vFaces);
                    break;
                  }
                  case 3: {
                    pPList vRegs = V_regions(pv);
                    void *tmp_ptr = 0;
                    pRegion pr;
                    while( (pr = (pRegion) PList_next(vRegs,&tmp_ptr)) ) {
                      int rPart = partitionTable[pr->iD];
                      partition = std::min(rPart,partition);
                      if (partition < 0) partition = rPart;
                    }
                    PList_delete(vRegs);
                    break;
                  }
                  }
                }

                if(dim == 0) {
                  if ( version == 1 ) { fprintf(f, "%d %d %d %d %d %d\n",
                                                k++, 15, phys, tag, 1, pv->iD); }
                  else                { fprintf(f, "%d %d %d %d %d %d %d\n",
                                                k++, 15, nbTags, phys, tag, partition + 1, pv->iD); }
                }
              }
          }
        VIter_delete(vit);    
      }
      // --- Edges ---
      {      
        EIter eit = M_edgeIter(m);
        pEdge pe;  
        while ((pe = EIter_next(eit)))
          {
            if(pe->g)
              {
                int dim = GEN_type(pe->g);
                int nbTags = 3;
                int tag = GEN_tag (pe->g);
                int phys = 0;
                if ( GEN_type(pe->g) == 1 ) phys = GEN_physTag(pe->g);
                if ( !saveAll && physic && phys == 0 ) continue; 
//                 int phys = gentity2phys.get_first_tag(pe->g);
                int partition = -1;
                if ( partitionTable ) {
                  switch(meshDim) {
                  case 1:
                    partition = partitionTable[pe->iD];
                    break;
                  case 2:
                    for (int i=0; i < E_numFaces(pe); i++) {
                      int fPart = partitionTable[E_face(pe,i)->iD];
                      partition = std::min(fPart,partition);
                      if (partition < 0) partition = fPart;
                    }
                    break;
                  case 3: {
                    pPList eRegs = E_regions(pe);
                    void *tmp_ptr = 0;
                    pRegion pr;
                    while( (pr = (pRegion) PList_next(eRegs,&tmp_ptr)) ) {
                      int rPart = partitionTable[pr->iD];
                      partition = std::min(rPart,partition);
                      if (partition < 0) partition = rPart;
                    }
                    PList_delete(eRegs);
                    break;
                  }
                  }
                }
	      
                if(dim ==1) {
                  if ( version == 1 ) { fprintf(f, "%d %d %d %d %d %d %d\n", 
                                                k++, 1, phys, tag, 2, 
                                                pe->p1->iD, pe->p2->iD); }
                  else                { fprintf(f, "%d %d %d %d %d %d %d %d\n",
                                                k++, 1, nbTags, phys, tag, partition + 1, 
                                                pe->p1->iD, pe->p2->iD); }
                }
              }
          }
        EIter_delete(eit);
      }
      // --- Faces ---  
      {
        FIter fit = M_faceIter(m);
        pFace pf;  
        while ((pf = FIter_next(fit)))
          {
            if(pf->g)
              {
                MDB_Point *nod[4];
                int dim = GEN_type(pf->g);
                int nbTags = 3;
                int tag = GEN_tag (pf->g); 
                int phys = 0;
                if ( GEN_type(pf->g) == 2 ) phys = GEN_physTag(pf->g);
#ifdef PARALLEL
                if ( !saveAll && physic && phys == 0 && !F_isInterface(pf) ) continue;
#else
                if ( !saveAll && physic && phys == 0 ) continue;
#endif
//                 int phys = gentity2phys.get_first_tag(pf->g);
                int partition = -1;
                if ( partitionTable ) {
                  switch(meshDim) {
                  case 2:
                    partition = partitionTable[pf->iD];
                    break;
                  case 3:
                    for (int i=0; i < F_numRegions(pf); i++) {
                      int rPart = partitionTable[F_region(pf,i)->iD];
                      partition = std::min(rPart,partition);
                      if (partition < 0) partition = rPart;
                    }
                    break;
                  }
                }
	      
                pf->getNodes(nod);
#ifdef PARALLEL
                if(dim == 2 || F_isInterface(pf)) {
#else
                if(dim == 2) {
#endif
                  if (pf->getNbEdges() == 4) { //quad
                    if ( version == 1 ) { fprintf(f, "%d %d %d %d %d %d %d %d %d\n",
                                                  k++, 3, phys, tag, 4,
                                                  nod[0]->iD, nod[1]->iD, nod[2]->iD, nod[3]->iD); }
                    else                { fprintf(f, "%d %d %d %d %d %d %d %d %d %d\n",
                                                  k++, 3, nbTags, phys, tag, partition + 1, 
                                                  nod[0]->iD, nod[1]->iD, nod[2]->iD, nod[3]->iD); }
                  }
                  else { //triangle
                    if ( version == 1 ) { fprintf(f, "%d %d %d %d %d %d %d %d\n",
                                                  k++, 2, phys, tag, 3,
                                                  nod[0]->iD, nod[1]->iD, nod[2]->iD); }
                    else                { fprintf(f, "%d %d %d %d %d %d %d %d %d\n",
                                                  k++, 2, nbTags, phys, tag, partition + 1, 
                                                  nod[0]->iD, nod[1]->iD, nod[2]->iD); }
                  }
                }
              }
          }
        FIter_delete(fit);
      } 
      {
        RIter rit = M_regionIter(m);
        pRegion pr;  
        while ((pr = RIter_next(rit)))
          {
            if(pr->g)
              {
                MDB_Point *nod[8];
                //int dim = GEN_type(pr->g);
                int nbTags = 3;
                int tag = GEN_tag (pr->g);
                int phys = 0;
                if ( GEN_type(pr->g) == 3 ) phys = GEN_physTag(pr->g);
                if ( !saveAll && physic && phys == 0 ) continue;
//                 int phys = gentity2phys.get_first_tag(pr->g);
                int partition = 0;
                if ( partitionTable ) {
                  partition = partitionTable[pr->iD];
                }
	      
                int numVer =  pr->getNbVertex();
                pPList ll = R_vertices(pr);
                for(int i=0;i<numVer;i++) nod[i] = (pVertex)PList_item(ll, i);
                PList_delete(ll);
              
                if(numVer==4) {
                  if ( version == 1 ) { fprintf(f, "%d %d %d %d %d %d %d %d %d\n", 
                                                k++, 4, phys, tag, 4,
                                                nod[0]->iD, nod[1]->iD, nod[2]->iD, nod[3]->iD); }
                  else                { fprintf(f, "%d %d %d %d %d %d %d %d %d %d\n",
                                                k++, 4, nbTags, phys, tag, partition + 1, 
                                                nod[0]->iD, nod[1]->iD, nod[2]->iD, nod[3]->iD); }
                }
                else if(numVer==8) {
                  if ( version == 1 ) { fprintf(f, "%d %d %d %d %d %d %d %d %d %d %d %d %d\n", 
                                                k++, 5, phys, tag, 8,
                                                nod[0]->iD, nod[1]->iD, nod[2]->iD, nod[3]->iD,
                                                nod[4]->iD, nod[5]->iD, nod[6]->iD, nod[7]->iD); }
                  else                { fprintf(f, "%d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
                                                k++, 5, nbTags, phys, tag, partition + 1, 
                                                nod[0]->iD, nod[1]->iD, nod[2]->iD, nod[3]->iD,
                                                nod[4]->iD, nod[5]->iD, nod[6]->iD, nod[7]->iD); }
                }
                else if(numVer==6) { 
                  if ( version == 1 ) { fprintf(f, "%d %d %d %d %d %d %d %d %d %d %d\n", 
                                                k++, 6, phys, tag, 6,
                                                nod[0]->iD, nod[1]->iD, nod[2]->iD, nod[3]->iD,
                                                nod[4]->iD, nod[5]->iD); }
                  else                { fprintf(f, "%d %d %d %d %d %d %d %d %d %d %d %d\n",
                                                k++, 6, nbTags, phys, tag, partition + 1, 
                                                nod[0]->iD, nod[1]->iD, nod[2]->iD, nod[3]->iD,
                                                nod[4]->iD, nod[5]->iD); }
                }
              }
          }
        RIter_delete(rit);
      }
      if ( version == 1 ) { fprintf(f, "$ENDELM\n");   } 
      else                { fprintf(f, "$EndElements\n"); }
      {
        int nperiodic = 0,nref = 0;;
        VIter vit = M_vertexIter(m);
        pVertex pv;  
        while ((pv = VIter_next(vit)))
          {
            void *temp_ptr; 
            int isPeriodic = EN_getDataPtr((pEntity) pv , tagPeriodic, &temp_ptr);
            if(isPeriodic) {
              std::vector<std::pair<std::vector<int> , pVertex> > *recup = (std::vector<std::pair<std::vector<int> , pVertex> > *) temp_ptr;
              if(!nperiodic) {
                nref = (*recup)[0].first.size();
              }
              nperiodic+=(*recup).size();
            }
          }
        VIter_delete(vit);  
        if(nperiodic) {
          fprintf(f, "$PERIODICNOD\n");
          fprintf(f, "%d %d\n", nperiodic,nref);
          int num = 1;
          vit = M_vertexIter(m);
          while ((pv = VIter_next(vit)))
            {
              void *temp_ptr; 
              int isPeriodic = EN_getDataPtr((pEntity) pv , tagPeriodic, &temp_ptr);
              if(isPeriodic) {
	        std::vector<std::pair<std::vector<int> , pVertex> > *recup = (std::vector<std::pair<std::vector<int> , pVertex> > *) temp_ptr;
		for(unsigned int k=0 ; k<(*recup).size() ; k++) {
	          fprintf(f, "%d", num++);
	          for(unsigned int j=0 ; j<(*recup)[k].first.size() ; j++) fprintf(f, " %d", (*recup)[k].first[j]);
		  fprintf(f, " %d %d \n", EN_id((pEntity) pv),EN_id((pEntity) (*recup)[k].second));
	        }
		//delete recup;
		//EN_deleteData((pEntity) pv , tagPeriodic);
              }
            }
          VIter_delete(vit);
          fprintf(f, "$ENDPERIODICNOD\n");  
        }
      }  
    }
    fclose(f);
    //MD_deleteMeshDataId(tagPeriodic);

  }



  void SaveGmshMeshPer (const pMesh m,const char *filename,MDB_DataExchangerPeriodic &deperiodic, int version)
  {
    int nbModelVertex = 0;
    pMeshDataId tagPeriodic = MD_lookupMeshDataId("PeriodicPoint");
    // Create a reverse map of physical tags
    GEntity2Physical gentity2phys(m->geomFeatures_Tags);
    FILE *f = fopen(filename, "w");
    {
      fprintf(f, "$NOD\n");
      fprintf(f, "%d\n", m->nbPoints);
      VIter vit = M_vertexIter(m); 
      while (VIter_next(vit)){}
      VIter_reset(vit);
      pVertex pv;  
      int NN = 0;
      while ((pv = VIter_next(vit)))
        {
          if(pv->g)
            {
              NN++;
              int dim = GEN_type(pv->g); 
              if (pv->deleted)printf("ouuch\n");
              if(dim == 0)
                nbModelVertex++;
              void *temp_ptr; 
              int isPeriodic = EN_getDataPtr((pEntity) pv , tagPeriodic, &temp_ptr);
              if(pv->X==2.25 && pv->Y==0.25) printf("%e %e %e : %d\n",pv->X,pv->Y,pv->Z,pv->iD);
              if(pv->X==2.75 && pv->Y==0.25) printf("%e %e %e : %d\n",pv->X,pv->Y,pv->Z,pv->iD);
              if(isPeriodic) {
     		        std::vector<std::pair<std::vector<int> , pVertex> > *recup = (
                  std::vector<std::pair<std::vector<int> , pVertex> > *) temp_ptr;
                          
                // the use of allocation is done due to ms visual c++ compiler 
                // that due not support c99 standard (it uses the c95 et c++98 standards)
                int* ntrans = new int[deperiodic.nbRefPeriodic()];
                for(int kk = 0 ; kk<deperiodic.nbRefPeriodic() ; kk++) {
                  ntrans[kk] = 0;
                }
                int nump = EN_id((pEntity) pv);
                pVertex pnump = pv;
                for(unsigned int k=0 ; k<(*recup).size() ; k++) {
                  if ( nump > EN_id((pEntity) (*recup)[k].second) ) {
                    nump  = EN_id((pEntity) (*recup)[k].second);
                    pnump =  (*recup)[k].second;
                    for(unsigned int j=0 ; j<(*recup)[k].first.size() ; j++) ntrans[j] = (-1) * (*recup)[k].first[j];
                  }	        		
     		        }
                if (nump == EN_id((pEntity) pv)) {  
                  fprintf(f, "%d %.16g %.16g %.16g\n", pv->iD, pv->X, pv->Y, pv->Z);					
                } else {
                  double X,Y,Z;
                  X = pnump->X;
                  Y = pnump->Y;
                  Z = pnump->Z; 
                  for(int kk = 0 ; kk<deperiodic.nbRefPeriodic() ; kk++) {
                    if(ntrans[kk]){
                      int inv1 = (ntrans[kk] < 0) ? 1 : 0;
                      for(int nb = 0 ; nb<abs(ntrans[kk]) ; nb++) { //nb de fois qu'il faut appliquer la transfo
                        deperiodic.fperiodic(inv1,X,Y,Z,kk+1,&X,&Y,&Z);
                      }  
                    }
                  }
                  fprintf(f, "%d %.16g %.16g %.16g\n", pv->iD, X, Y, Z);    
                  if((fabs(X-pv->X) > 1e-12) || (fabs(Y-pv->Y) > 1e-12) || (fabs(Z-pv->Z) > 1e-12)  ) {
                    printf("point %d : %e %e %e \n",EN_id((pEntity) pv),pv->X, pv->Y, pv->Z);
                    printf("-- corresp %d : %e %e %e\n",nump,pnump->X,pnump->Y,pnump->Z);
                    printf("*end* %e %e %e\n",X,Y,Z);  			
                  }
                }
                delete [] ntrans;
              } else {
                fprintf(f, "%d %.16g %.16g %.16g\n", pv->iD, pv->X, pv->Y, pv->Z);
              }
            }
          else
            {
              NN++;
              fprintf(f, "%d %.16g %.16g %.16g\n", pv->iD, pv->X, pv->Y, pv->Z);
            }
          // 	   throw;
        }
      if (NN != m->nbPoints) 
        {
          printf("%d != %d\n",NN,m->nbPoints);
          throw;
        }
      VIter_delete(vit);    
      fprintf(f, "$ENDNOD\n");
    }
    {
      fprintf(f, "$ELM\n");
    
      int nbClasEdges = 0;
      int nbClasFaces = 0;
      {
        EIter eit = M_edgeIter(m);
        pEdge pe;  
        while ((pe = EIter_next(eit)))
          {
            int dim = GEN_type(pe->g); 
            if(dim == 1)nbClasEdges++;
          }
        EIter_delete(eit);
      }    
      {
        FIter fit = M_faceIter(m);
        pFace pf;  
        while ((pf = FIter_next(fit)))
          {
            int dim = GEN_type(pf->g); 
            if(dim == 2)nbClasFaces++;
            else if(F_numRegions(pf)==1) nbClasFaces++;
          }
        FIter_delete(fit);
      } 

      fprintf(f, "%ld\n", (long int)(nbClasEdges + nbModelVertex + nbClasFaces + m->nbTets) );

      int k = 1;
      {
        VIter vit = M_vertexIter(m);
        pVertex pv;  
        while ((pv = VIter_next(vit)))
          {
            if(pv->g)
              {
                int dim = GEN_type(pv->g); 
                int tag = GEN_tag (pv->g); 
                int phys = gentity2phys.get_first_tag(pv->g);
	      
                if(dim == 0)
                  fprintf(f, "%d %d %d %d %d %d\n",
                          k++, 15, phys,tag, 1,pv->iD);
              }
          }
        VIter_delete(vit);    
      }
      {      
        EIter eit = M_edgeIter(m);
        pEdge pe;  
        while ((pe = EIter_next(eit)))
          {
            if(pe->g)
              {
                int dim = GEN_type(pe->g); 
                int tag = GEN_tag (pe->g); 
                int phys = gentity2phys.get_first_tag(pe->g);
	      
                if(dim ==1)
                  fprintf(f, "%d %d %d %d %d %d %d\n", k++, 1, phys,tag, 2,pe->p1->iD, pe->p2->iD);
              }
          }
        EIter_delete(eit);
      }
      {
        FIter fit = M_faceIter(m);
        pFace pf;  
        while ((pf = FIter_next(fit)))
          {
            if(pf->g)
              {
                MDB_Point *nod[4];
                int dim = GEN_type(pf->g); 
                int tag = GEN_tag (pf->g); 
                int phys = gentity2phys.get_first_tag(pf->g);
	      
                pf->getNodes(nod);
                if(dim == 2 || F_numRegions(pf)==1){
                  if (0 && nod[3]) //quad
                    fprintf(f, "%d %d %d %d %d %d %d %d %d\n",
                            k++, 3, phys,tag, 4,
                            nod[0]->iD, nod[1]->iD, nod[2]->iD, nod[3]->iD);
                  else //triangle
                    fprintf(f, "%d %d %d %d %d %d %d %d\n",
                            k++, 2, phys,tag, 3,
                            nod[0]->iD, nod[1]->iD, nod[2]->iD);

                }
              }
          }
        FIter_delete(fit);
      } 
      {
        RIter rit = M_regionIter(m);
        pRegion pr;  
        while ((pr = RIter_next(rit)))
          {
            if(pr->g)
              {
                MDB_Point *nod[8];
                //int dim = GEN_type(pr->g); 
                int tag = GEN_tag (pr->g);  
                int phys = gentity2phys.get_first_tag(pr->g);
	      
                int numVer =  pr->getNbVertex();
                pPList ll = R_vertices(pr);
                for(int i=0;i<numVer;i++) nod[i] = (pVertex)PList_item(ll, i);
                PList_delete(ll);
              
                if(numVer==4) fprintf(f, "%d %d %d %d %d %d %d %d %d\n", k++, 4, phys,tag, 4,
                                      nod[0]->iD, nod[1]->iD, nod[2]->iD, nod[3]->iD);
                else if(numVer==8) fprintf(f, "%d %d %d %d %d %d %d %d %d %d %d %d %d\n", k++, 5, phys,tag, 5,
                                           nod[0]->iD, nod[1]->iD, nod[2]->iD, nod[3]->iD, nod[4]->iD, nod[5]->iD, nod[6]->iD, nod[7]->iD);
              }
          }
        RIter_delete(rit);
      }
      fprintf(f, "$ENDELM\n");
      {
        int nperiodic = 0,nref = 0;;
        VIter vit = M_vertexIter(m);
        pVertex pv;  
        while ((pv = VIter_next(vit)))
          {
            void *temp_ptr; 
            int isPeriodic = EN_getDataPtr((pEntity) pv , tagPeriodic, &temp_ptr);
            if(isPeriodic) {
              std::vector<std::pair<std::vector<int> , pVertex> > *recup = (std::vector<std::pair<std::vector<int> , pVertex> > *) temp_ptr;
              if(!nperiodic) {
                nref = (*recup)[0].first.size();
              }
              nperiodic+=(*recup).size();
            }
          }
        VIter_delete(vit);  
        if(nperiodic) {
          fprintf(f, "$PERIODICNOD\n");
          fprintf(f, "%d %d\n", nperiodic,nref);
          int num = 1;
          vit = M_vertexIter(m);
          while ((pv = VIter_next(vit)))
            {
              void *temp_ptr; 
              int isPeriodic = EN_getDataPtr((pEntity) pv , tagPeriodic, &temp_ptr);
              if(isPeriodic) {
	        std::vector<std::pair<std::vector<int> , pVertex> > *recup = (std::vector<std::pair<std::vector<int> , pVertex> > *) temp_ptr;
		for(unsigned int k=0 ; k<(*recup).size() ; k++) {
	          fprintf(f, "%d", num++);
	          for(unsigned int j=0 ; j<(*recup)[k].first.size() ; j++) fprintf(f, " %d", (*recup)[k].first[j]);
		  fprintf(f, " %d %d \n", EN_id((pEntity) pv),EN_id((pEntity) (*recup)[k].second));
	        }
		//delete recup;
		//EN_deleteData((pEntity) pv , tagPeriodic);
              }
            }
          VIter_delete(vit);
          fprintf(f, "$ENDPERIODICNOD\n");  
        }
      }  
    }
    fclose(f);
    //MD_deleteMeshDataId(tagPeriodic);

  }

}



  /*
// -------------------------------------------------------------------
// Old loader (no parallel)
// -------------------------------------------------------------------

void LoadGmshMesh (pMesh m,const char *filename, int version)
{  
  FILE *fp = fopen(filename, "r");
  int  isperiodic=0;
  pMeshDataId tagPeriodic = MD_newMeshDataId("PeriodicPoint");
  if(!fp)
    {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                  "Unknown File %s",filename);
    }
  char String[256];
  if (!m->model)
    m->model = new NullModel; 

  while(1) {
    do {
      if(!fgets(String, sizeof(String), fp))
	break;
      if(feof(fp))
        break;
    } while(String[0] != '$');
    
    if(feof(fp))
      break;
    
    if(!strncmp(&String[1], "MeshFormat", 10)) {
      int a, b, c;
      fscanf(fp, "%d %d %d", &a, &b, &c);
      version = a;      
    }
    
    if(!strncmp(&String[1], "NOD", 3) ||
       !strncmp(&String[1], "NOE", 3) ||
       !strncmp(&String[1], "Nodes", 5)) {
      
      int Nbr_Nodes,Num;
      double x,y,z;
      fscanf(fp, "%d", &Nbr_Nodes);
      for(int i_Node = 0; i_Node < Nbr_Nodes; i_Node++) {
        fscanf(fp, "%d %lf %lf %lf", &Num, &x, &y, &z);
        m->add_point (Num,x,y,z);
      }
    }
    
    // ELEMENTS
    
    else if(!strncmp(&String[1], "ELM", 3) ||
	    !strncmp(&String[1], "Elements", 8)) {
      int Nbr_Elements, NbTags, verts[256],Tag;
      int Num, Type, Physical, Elementary, Nbr_Nodes, Partition;
      fscanf(fp, "%d", &Nbr_Elements);

      for(int i_Element = 0; i_Element < Nbr_Elements; i_Element++) {
	
	if(version == 1){
	  fscanf(fp, "%d %d %d %d %d",
		 &Num, &Type, &Physical, &Elementary, &Nbr_Nodes);
	  Partition = 1;
	}
	else{
	  fscanf(fp, "%d %d %d", &Num, &Type, &NbTags);
	  Elementary = Physical = Partition = 1;
	  for(int j = 0; j < NbTags; j++){
	    fscanf(fp, "%d", &Tag);	    
	    if(j == 0)
	      Physical = Tag;
	    else if(j == 1)
	      Elementary = Tag;
	    else if(j == 2)
	      Partition = Tag;
	    // ignore any other tags for now
	  }
	  Nbr_Nodes = getNumVerticesForElementTypeMSH(Type);
	}

        for(int j = 0; j < Nbr_Nodes; j++)
          fscanf(fp, "%d", &verts[j]);

	GEntity *geom = 0;
	
        switch (Type) {
        case 1:
	  {
	    geom = m->model->edgeByTag(Elementary);
	    m->add_edge(verts[0],verts[1],geom);
	  }
	  break;
        case 8:
	  {
	    geom = m->model->edgeByTag(Elementary);
	    //m->add_edge(verts[0],verts[2],verts[1],geom);
	    m->add_edge(3,geom,verts[0],verts[2],verts[1]);
	    //m->add_edge(verts[0],verts[1],ge);
	  }
	  break;
        case 26:
	  {
	    geom = m->model->edgeByTag(Elementary);
	    m->add_edge(4,geom, verts[0],verts[2],verts[3],verts[1]);
	  }
	  break;
        case 27:
	  {
	    geom = m->model->edgeByTag(Elementary);
	    m->add_edge(5,geom, verts[0],verts[2],verts[3],verts[4],verts[1]);
	  }
	  break;
        case 28:
	  {
	    geom = m->model->edgeByTag(Elementary);
	    m->add_edge(6,geom, verts[0],verts[2],verts[3],verts[4],verts[5],verts[1]);
	  }
	  break;
        case 2:
	  {
	    geom = m->model->faceByTag(Elementary);
	    m->add_triangle(verts[0],verts[1],verts[2],geom); 
	  }
	  break;
// 6 nodes triangle
        case MSH_TRI_6:
	  {
	    geom = m->model->faceByTag(Elementary);
	    //	    printf ("adding a second order triangle %d %d %d %d %d %d\n",verts[0],verts[1],verts[2],verts[3],verts[4],verts[5]);
	    m->add_triangle(2,0,geom,verts[0],verts[3],verts[1],verts[4],verts[2],verts[5]); 
	    //m->add_triangle(verts[0],verts[3],verts[1],verts[4],verts[2],verts[5],geom); 
	    //m->add_triangle(verts[0],verts[1],verts[2],gf); 
	  }
	  break;
// 9 nodes triangle (SERENDIP !)
        case MSH_TRI_9:
	  {
	    geom = m->model->faceByTag(Elementary);
	    m->add_triangle(3,0,geom,verts[0],verts[3],verts[4],verts[1],verts[5],verts[6],verts[2],verts[7],verts[8]); 
	    //m->add_triangle(verts[0],verts[3],verts[1],verts[4],verts[2],verts[5],geom); 
	  }
	  break;
// 10 nodes triangle (LAGRANGE !)
        case MSH_TRI_10:
	  {
	    geom = m->model->faceByTag(Elementary);
	    m->add_triangle(3,1,geom,verts[0],verts[3],verts[4],verts[1],verts[5],verts[6],verts[2],verts[7],verts[8],verts[9]); 
	  }
	  break;
// 12 nodes triangle (SERENDIP !)
        case MSH_TRI_12:
	  {
	    geom = m->model->faceByTag(Elementary);
	    m->add_triangle(4,0,geom,verts[0],verts[3],verts[4],verts[5],verts[1],verts[6],verts[7],verts[8],verts[2],verts[9],verts[10],verts[11]); 
	  }
	  break;
// 15 nodes triangle (LAGRANGE !)
        case MSH_TRI_15:
	  {
	    geom = m->model->faceByTag(Elementary);
	    m->add_triangle(4,1,geom,verts[0],verts[3],verts[4],verts[5],verts[1],verts[6],verts[7],verts[8],verts[2],verts[9],verts[10],verts[11],verts[12],verts[13],verts[14]); 
	  }
	  break;
// 15 nodes triangle (SERENDIP !)
        case MSH_TRI_15I:
	  {
	    geom = m->model->faceByTag(Elementary);
	    m->add_triangle(5,0,geom,verts[0],verts[3],verts[4],verts[5],verts[6],verts[1],verts[7],verts[8],verts[9],verts[10],verts[2],verts[11],verts[12],verts[13],verts[14]); 
	  }
	  break;
        case 3:
	  {
	    geom = m->model->faceByTag(Elementary);
	    m->add_quad(verts[0],verts[1],verts[2],verts[3],geom); 
	  }
	  break;
        case 4:
	  {
	    geom = m->model->regionByTag(Elementary);
	    m->add_tet(verts[0],verts[1],verts[2],verts[3],geom); 
	  }
	  break;
        case 6:
          {
            geom = m->model->regionByTag(Elementary);
            m->add_prism(verts[0],verts[1],verts[2],
                         verts[3],verts[4],verts[5],geom); 
          }
          break;
        case 15:
	  {
	    geom = m->model->vertexByTag(Elementary);
	    MDB_Point *p = m->find_point(verts[0]);
	    p->g = geom;
	  }
	  break;
// 8 nodes hexahedron
        case 5:
	  {
	    geom = m->model->regionByTag(Elementary); 
	    m->add_hex(verts[0],verts[1],verts[2],verts[3],verts[4],verts[5],verts[6],verts[7],geom); 
	  }
	  break;
          // 10 node tetrahedron, 2nd order
        case MSH_TET_10:
          {
            geom = m->model->regionByTag(Elementary);
            m->add_tet(geom,2,false,verts);
            break;
          }
        default:
	  throw;
          break;
        }
	if (geom)
	  {
	    bool find = false;
	    for (std::multimap<int, pGEntity>::iterator it = m->geomFeatures_Tags.lower_bound(Physical);
		 it != m->geomFeatures_Tags.upper_bound(Physical);++it)
	      if (it->second == geom)find = true;
	    if (!find)
	      m->geomFeatures_Tags.insert(std::pair<int,pGEntity>(Physical, geom));
	  }
      }
    }
    else if(!strncmp(&String[1], "PERIODICNOD", 11)) {
      isperiodic = 1;
    } 
    else if(!strncmp(&String[1], "TRA", 3)) {
      isperiodic = 2;
    }

    do {
      if(!fgets(String, sizeof(String), fp))
	throw;
      if(feof(fp))
	throw;
    } while(String[0] != '$');    
  }
  
  m->classify_unclassified_entities();

  m->destroyStandAloneEntities();

  if(isperiodic) {
    if(isperiodic==1) {
      fseek(fp,0,SEEK_SET);
      while(1) {
        do {
          if(!fgets(String, sizeof(String), fp))
	    break;
          if(feof(fp))
            break;
        } while(String[0] != '$');
    
        if(feof(fp))
          break;
        if(!strncmp(&String[1], "PERIODICNOD", 11)) {
          int nPeriod,nRefPeriod;
          fscanf(fp, "%d %d", &nPeriod,&nRefPeriod);
	  printf("%d Ref Periodic\n",nRefPeriod);
	  int *t=new int[nRefPeriod];
         for(int ip = 0; ip < nPeriod; ip++) {
	   int Num,v1,v2;
           fscanf(fp, "%d ", &Num);
	   std::vector<int> transfo;
	   std::vector<int> invtransfo;
	   for(int k=0 ; k<nRefPeriod ; k++) {
	     fscanf(fp, "%d ", &t[k]);
	     transfo.push_back(t[k]);
	     invtransfo.push_back(-t[k]);
	   } 
	   fscanf(fp, "%d %d",&v1, &v2);
	   pVertex vt1 = m->find_point(v1);
	   pVertex vt2 = m->find_point(v2);
//printf("point %d %d\n",EN_id((pEntity) vt1),EN_id((pEntity) vt2));
	   void *temp_ptr; 
           int isPeriodic = EN_getDataPtr((pEntity) vt1 , tagPeriodic, &temp_ptr);
	   if(isPeriodic) {
	      std::vector<std::pair<std::vector<int> , pVertex> > *recup = (std::vector<std::pair<std::vector<int> , pVertex> > *) temp_ptr;
              unsigned int i=0;
	      for(i=0 ; i<(*recup).size() ; i++) {
	        unsigned int j=0;
		for(j=0 ; j<(*recup)[i].first.size() ; j++) {
		  if((*recup)[i].first[j] != transfo[j]) break;
		}
		if(j==(*recup)[i].first.size() ) break;
	      }
	      if(i==(*recup).size())
	        (*recup).push_back(std::pair <std::vector<int> , pVertex>(transfo,vt2));
	   } else {
             std::vector<std::pair<std::vector<int> , MDB_Point*> > remotePoints;
	     (remotePoints).push_back(std::pair <std::vector<int> , pVertex>(transfo,vt2));
	      EN_attachDataPtr((pEntity) vt1 , tagPeriodic, 
	    			new std::vector<std::pair<std::vector<int> , pVertex> >(remotePoints));
				   
	   }
	   
           isPeriodic = EN_getDataPtr((pEntity) vt2 , tagPeriodic, &temp_ptr);
	   if(isPeriodic) {
	      std::vector<std::pair<std::vector<int> , pVertex> > *recup = (std::vector<std::pair<std::vector<int> , pVertex> > *) temp_ptr;
              unsigned int i;
	      for(i=0 ; i<(*recup).size() ; i++) {
	        unsigned int j;
		for(j=0 ; j<(*recup)[i].first.size() ; j++) {
		  if((*recup)[i].first[j] != invtransfo[j]) break;
		}
		if(j==(*recup)[i].first.size() ) break;
	      }
	      if(i==(*recup).size())
                (*recup).push_back(std::pair <std::vector<int> , pVertex>(invtransfo,vt1));
	   } else {
             std::vector<std::pair<std::vector<int> , MDB_Point*> > remotePoints;
	     (remotePoints).push_back(std::pair <std::vector<int> , pVertex>(invtransfo,vt1));
	      EN_attachDataPtr((pEntity) vt2 , tagPeriodic, 
	    			new std::vector<std::pair<std::vector<int> , pVertex> >(remotePoints));	   
	    }
          }
          delete []t;
	}  
      }
    } else {
      assert(isperiodic==2);
      puts("format periodic CENAERO");
      fseek(fp,0,SEEK_SET);
      while(1) {
        do {
          if(!fgets(String, sizeof(String), fp))
	    break;
          if(feof(fp))
            break;
        } while(String[0] != '$');
    
        if(feof(fp))
          break;
        if(!strncmp(&String[1], "PER", 3)) {
          int nPeriod,nRefPeriod;
	  nRefPeriod = 1;
          fscanf(fp, "%d ", &nPeriod);
	  printf("%d Ref Periodic\n",nRefPeriod);
	  int *t=new int[nRefPeriod];
          for(int ip = 0; ip < nPeriod; ip++) {
	   int Num,v1,v2;
           fscanf(fp, "%d ", &Num);
	   std::vector<int> transfo;
	   std::vector<int> invtransfo;
	   fscanf(fp, "%d %d %d",&v1, &v2,&Num);
	   pVertex vt1 = m->find_point(v1);
	   pVertex vt2 = m->find_point(v2);


	printf("point %d (%d) %d (%d)\n",EN_id((pEntity) vt1),GEN_tag(vt1->g),EN_id((pEntity) vt2),
				GEN_tag(vt2->g));
	if(GEN_tag(vt1->g)==132 || GEN_tag(vt1->g)==130 || GEN_tag(vt1->g)== 118 || GEN_tag(vt1->g)== 117) {
	  for(int k=0 ; k<nRefPeriod ; k++) {
             t[k]=-1;
	     transfo.push_back(t[k]);
	     invtransfo.push_back(-t[k]);
	   } 
 
	} else {
	  assert(GEN_tag(vt1->g)==348 || GEN_tag(vt1->g)==31 
	     || GEN_tag(vt1->g)== 93 || GEN_tag(vt1->g)== 345
	     || GEN_tag(vt1->g)== 29 || GEN_tag(vt1->g)== 94
	     || GEN_tag(vt1->g)== 16 || GEN_tag(vt1->g)== 66
	     || GEN_tag(vt1->g)== 13 || GEN_tag(vt1->g)== 74
	     || GEN_tag(vt1->g)== 115 || GEN_tag(vt1->g)== 351
	     || GEN_tag(vt1->g)== 116);
	  for(int k=0 ; k<nRefPeriod ; k++) {
             t[k]=1;
	     transfo.push_back(t[k]);
	     invtransfo.push_back(-t[k]);
	   } 

	}
	   void *temp_ptr; 
           int isPeriodic = EN_getDataPtr((pEntity) vt1 , tagPeriodic, &temp_ptr);
	   if(isPeriodic) {
	      std::vector<std::pair<std::vector<int> , pVertex> > *recup = (std::vector<std::pair<std::vector<int> , pVertex> > *) temp_ptr;
              unsigned int i=0;
	      for(i=0 ; i<(*recup).size() ; i++) {
	        unsigned int j=0;
		for(j=0 ; j<(*recup)[i].first.size() ; j++) {
		  if((*recup)[i].first[j] != transfo[j]) break;
		}
		if(j==(*recup)[i].first.size() ) break;
	      }
	      if(i==(*recup).size())
	        (*recup).push_back(std::pair <std::vector<int> , pVertex>(transfo,vt2));
	   } else {
             std::vector<std::pair<std::vector<int> , MDB_Point*> > remotePoints;
	     (remotePoints).push_back(std::pair <std::vector<int> , pVertex>(transfo,vt2));
	      EN_attachDataPtr((pEntity) vt1 , tagPeriodic, 
	    			new std::vector<std::pair<std::vector<int> , pVertex> >(remotePoints));
				   
	   }
	   
           isPeriodic = EN_getDataPtr((pEntity) vt2 , tagPeriodic, &temp_ptr);
	   if(isPeriodic) {
	      std::vector<std::pair<std::vector<int> , pVertex> > *recup = (std::vector<std::pair<std::vector<int> , pVertex> > *) temp_ptr;
              unsigned int i;
	      for(i=0 ; i<(*recup).size() ; i++) {
	        unsigned int j;
		for(j=0 ; j<(*recup)[i].first.size() ; j++) {
		  if((*recup)[i].first[j] != invtransfo[j]) break;
		}
		if(j==(*recup)[i].first.size() ) break;
	      }
	      if(i==(*recup).size())
                (*recup).push_back(std::pair <std::vector<int> , pVertex>(invtransfo,vt1));
	   } else {
             std::vector<std::pair<std::vector<int> , MDB_Point*> > remotePoints;
	     (remotePoints).push_back(std::pair <std::vector<int> , pVertex>(invtransfo,vt1));
	      EN_attachDataPtr((pEntity) vt2 , tagPeriodic, 
	    			new std::vector<std::pair<std::vector<int> , pVertex> >(remotePoints));	   
	    }
          }
          delete []t;
	}
 
       }
    }  
  } 

  // compute parametric coordinates of vertices on
//   VIter vit = M_vertexIter(m); 
//   pVertex pv;  
//   while (pv = VIter_next(vit))
//     {
//       int vType=V_whatInType(pv);
//       if(vType==1)
// 	{
// 	  double xyz[3];
// 	  V_coord(pv, xyz);
// 	  pGEdge pge = (pGEdge) EN_whatIn(pv);
// 	  double par = pge->parFromPoint(SPoint3(xyz[0],xyz[1],xyz[2])); 
// 	  printf("parametric coord is %12.5E\n",par);
// 	  P_setParam1(pv, par);
// 	}
//     }	
  fclose(fp);
}



*/

