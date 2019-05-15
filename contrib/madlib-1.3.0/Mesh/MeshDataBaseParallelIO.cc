// -------------------------------------------------------------------
// MAdLib - Copyright (C) 2008-2009 Universite catholique de Louvain
//
// See the Copyright.txt and License.txt files for license information. 
// You should have received a copy of these files along with MAdLib. 
// If not, see <http://www.madlib.be/license/>
//
// Please report all bugs and problems to <contrib@madlib.be>
//
// Authors: Jean-Francois Remacle, Cecile Dobrzynski, Gaetan Compere
// -------------------------------------------------------------------

#include "MeshDataBase.h"
#include "MeshDataBaseInterface.h"
#include "MeshDataBaseIO.h"
#include "MeshDataBaseGEntity2Physical.h"
#ifndef _WIN_
#include <unistd.h>
#endif
#include "MshTags.h"

#ifdef PARALLEL
#include "mpi.h"
#include "autopack.h"
#include "fcntl.h"
#include "MeshDataBaseParallelInterface.h"
#endif

#include <cstdlib>
#include <cstring>
// #include <set>
// #include <utility>

extern int getNumVerticesForElementTypeMSH(int type);

namespace MAd {

#ifdef PARALLEL
  // -------------------------------------------------------------------
  // -------------------- SaveGmshMeshParallel -------------------------
  // -------------------------------------------------------------------

  // Save a mesh (for parallel only) with format msh1 or msh2 \ingroup parallel
  // GCRMK: now vertex ids are unique: we can use the ids directly and simplify this function
  void SaveGmshMeshParallel (const pMesh mesh, const char *filename, int version)
  {
    // ***************** update ID interface ***************
    pMeshDataId tagData = MD_lookupMeshDataId("RemotePoint");
    V_createInfoInterface(mesh, tagData);   
    MPI_Barrier(MPI_COMM_WORLD);

    int nbModelVertex = 0;
    int myrank,nproc;

    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if(myrank==0) {
      FILE * nop =  fopen(filename,"w");
      fclose(nop);
    }
    MPI_Barrier(MPI_COMM_WORLD); 

    int npt = 0;
  
    VIter vit = M_vertexIter(mesh); 
    while (VIter_next(vit)){}
    VIter_reset(vit);
    pVertex pv;  
    int NN = 0,nssint = 0;;
    while ((pv = VIter_next(vit)))
      { 
        if(pv->g)
          {
            NN++;
            void *temp_ptr; 
            int isInterface = EN_getDataPtr((pEntity) pv , tagData, &temp_ptr);    
            int j =0,size= 0;
            if(isInterface) {
              const std::vector<std::pair<int , MDB_Point*> > *recup = 
                (const std::vector<std::pair<int , MDB_Point*> > *) temp_ptr;
              size = (*recup).size() ;
              for(j=0 ; j<size; j++) {
                int remoteID  = (*recup)[j].first;
                assert(remoteID!=myrank);
                if(remoteID < myrank) break;
              }
            }
            if(j == size) {
              nssint++;
            } 
          } else {
          throw;
        }
      }
    if (NN != mesh->nbPoints) 
      {
        printf("%d != %d\n",NN,mesh->nbPoints);
        throw;
      }
    VIter_delete(vit); 
    int *tab=new int[nproc];
    int sendnum = nssint; 
    if(myrank) {
      MPI_Send(&sendnum,1,MPI_INT,0,myrank,MPI_COMM_WORLD);   
    } else {
      MPI_Status status;
      tab[0] = sendnum;
      for(int i=1 ; i<nproc ; i++) {
        MPI_Recv(&tab[i],1,MPI_INT,i,i,MPI_COMM_WORLD,&status);
      }
    }
  
    MPI_Bcast(&tab[0],nproc,MPI_INT,0,MPI_COMM_WORLD);
  
    int *tabmax=new int[nproc];
    sendnum = mesh->maxId; 
    if(myrank) {
      MPI_Send(&sendnum,1,MPI_INT,0,myrank,MPI_COMM_WORLD);   
    } else {
      MPI_Status status;
      tabmax[0] = sendnum;
      for(int i=1 ; i<nproc ; i++) {
        MPI_Recv(&tabmax[i],1,MPI_INT,i,i,MPI_COMM_WORLD,&status);
      }
    }
  
    MPI_Bcast(&tabmax[0],nproc,MPI_INT,0,MPI_COMM_WORLD);
  
    // --------------------------------------------------
    // Defines a control point for a sequential writting:
    //    Proc n+1 does not continue while proc n is 
    //    not at the next control point
    int send = 0,recv = 0;
    if(!myrank) recv = 1;
    MPI_Status status;
    while(!recv) {  
      if(myrank!=(nproc-1)) MPI_Send(&send,1,MPI_INT,(myrank+1)%nproc,(myrank+1)%nproc,MPI_COMM_WORLD);
      MPI_Recv(&recv,1,MPI_INT,(myrank+nproc-1)%nproc,myrank,MPI_COMM_WORLD,&status);
    }
    int nseek = recv;
    int nplace = nseek;
    if (!myrank) nplace = 0;
    // --------------------------------------------------

    int nop = open(filename,O_WRONLY);
  

    // ----------------------------------------------------
    // Write format (for msh2 only)
    // ----------------------------------------------------
  
    if ( !myrank && version != 1 ) {
      char format[256];
      sprintf(format, "$MeshFormat\n2 0 8\n$EndMeshFormat\n");
      int size = strlen(format);
      write(nop,format,size);
      nplace += size;
    }

    // ----------------------------------------------------
    // Write partitionning of nodes (for msh1 only)
    // ----------------------------------------------------

    if ( version == 1 ) {
 
      if(!myrank) {
        write(nop,"$PARALLELNOD\n",sizeof("$PARALLELNOD"));
        nplace += sizeof("$PARALLELNOD");
        char nb[256];
        int nbtot = 0;
        for(int i=0 ; i<nproc ; i++) nbtot +=tab[i];
        sprintf(nb,"%d\n",nbtot);
        int size = strlen(nb);
        nplace += size;
        write(nop,nb,size);
      }
      //    int IdGlobal = 0;
      if(myrank) {
        lseek(nop,nseek,0);
        //      for(int j=0 ; j<myrank; j++) IdGlobal += tabmax[j];
      }  
      vit = M_vertexIter(mesh); 
      while (VIter_next(vit)){}
      VIter_reset(vit);  
      NN = 0;
      while ((pv = VIter_next(vit)))
        { 
          if(pv->g)
            {
              NN++;
              if (pv->deleted)printf("ouuch\n");
              void *temp_ptr; 
              int isInterface = EN_getDataPtr((pEntity) pv , tagData, &temp_ptr);    
              int j =0,size= 0/*,nId = IdGlobal*/;
              if(isInterface) {
                const std::vector<std::pair<int , MDB_Point*> > *recup = (const std::vector<std::pair<int , MDB_Point*> > *) temp_ptr;
                size = (*recup).size() ;
                for(j=0 ; j<size; j++) {
                  int remoteID  = (*recup)[j].first;
                  if(remoteID < myrank) break;
                }
              }
              if(j == size) {
                nssint++;
  
                char nb[256];
                if(size) sprintf(nb,"%d %d %d ",pv->iD /*+ nId*/, size + 1,myrank);
                else     sprintf(nb,"%d %d %d \n",pv->iD /*+ nId*/, size + 1,myrank);
                for(j=0 ; j<size; j++) {
                  const std::vector<std::pair<int , MDB_Point*> > *recup = (const std::vector<std::pair<int , MDB_Point*> > *) temp_ptr;
                  char nbtmp[256];
                  if(j == (size-1) )  sprintf(nbtmp,"%d \n",(*recup)[j].first); 
                  else                sprintf(nbtmp,"%d ",(*recup)[j].first); 
                  int sizetmp = strlen(nbtmp);  
                  strncat(nb,nbtmp,sizetmp);    
                } 
  
                int size = strlen(nb);
                nplace += size;
                write(nop,nb,size);
              } 
            } else {
            throw;
          }
        }
      if (NN != mesh->nbPoints) 
        {
          printf("%d != %d\n",NN,mesh->nbPoints);
          throw;
        }
      VIter_delete(vit);    
  
      MPI_Send(&nplace,1,MPI_INT,(myrank+1)%nproc,(myrank+1)%nproc,MPI_COMM_WORLD);
  
      if(!myrank) MPI_Recv(&npt,1,MPI_INT,(myrank+nproc-1)%nproc,myrank,MPI_COMM_WORLD,&status);
    
      MPI_Barrier(MPI_COMM_WORLD);

      send = 0,recv = 0;
      if(!myrank) recv = 1;
  
      while(!recv) {  
        if(myrank!=(nproc-1)) MPI_Send(&send,1,MPI_INT,(myrank+1)%nproc,(myrank+1)%nproc,MPI_COMM_WORLD);
        MPI_Recv(&recv,1,MPI_INT,(myrank+nproc-1)%nproc,myrank,MPI_COMM_WORLD,&status);
      }
  
      nseek = recv;
      nplace = nseek;
  
      if(!myrank) {
        nplace = npt;
        lseek(nop,nplace,0);
        write(nop,"$ENDPARALLELNOD\n",sizeof("$ENDPARALLELNOD"));
        nplace += sizeof("$ENDPARALLELNOD");
      }

    }

    // ----------------------------------------------------
    // Write nodes
    // ----------------------------------------------------

    if(!myrank) {
      char nb[256];
      if ( version == 1 ) sprintf(nb,"$NOD\n");
      else                sprintf(nb,"$Nodes\n");
      int size = strlen(nb);
      nplace += size;
      write(nop,nb,size);
    }
    if(!myrank) {    
      char nb[256];
      int nbtot = 0;
      for(int i=0 ; i<nproc ; i++) nbtot +=tab[i];
      sprintf(nb,"%d\n",nbtot);
      int size = strlen(nb);
      nplace += size;
      write(nop,nb,size);
    }
    //  int IdGlobal = 0;
    if(myrank) {
      lseek(nop,nseek,0);
      //    for(int j=0 ; j<myrank; j++) IdGlobal += tabmax[j];
    }  
    vit = M_vertexIter(mesh); 
    while (VIter_next(vit)){}
    VIter_reset(vit);  
    NN = 0;
    while ((pv = VIter_next(vit)))
      { 
        if(pv->g)
          {
            NN++;
            int dim = GEN_type(pv->g); 

            if (pv->deleted)printf("ouuch\n");
            void *temp_ptr; 
            int isInterface = EN_getDataPtr((pEntity) pv , tagData, &temp_ptr);    
            int j =0,size= 0/*,nId = IdGlobal*/;
            if(isInterface) {
              const std::vector<std::pair<int , MDB_Point*> > *recup = 
                (const std::vector<std::pair<int , MDB_Point*> > *) temp_ptr;
              size = (*recup).size() ;
              for(j=0 ; j<size; j++) {
                int remoteID  = (*recup)[j].first;
                if(remoteID < myrank) break;
              }
            }
            if(j == size) {
              nssint++;
              if(dim == 0)
                nbModelVertex++;
  
              char nb[256];
              sprintf(nb,"%d %g %g %g\n",pv->iD /*+ nId*/, pv->X, pv->Y, pv->Z);
              int size = strlen(nb);
              nplace += size;
              write(nop,nb,size);
            } 
          } else {
          throw;
        }
      }
    if (NN != mesh->nbPoints) 
      {
        printf("%d != %d\n",NN,mesh->nbPoints);
        throw;
      }
    VIter_delete(vit);    
  
    // --------------------------------------------------
    // End of the control point
    MPI_Send(&nplace,1,MPI_INT,(myrank+1)%nproc,(myrank+1)%nproc,MPI_COMM_WORLD);
  
    if(!myrank) MPI_Recv(&npt,1,MPI_INT,(myrank+nproc-1)%nproc,myrank,MPI_COMM_WORLD,&status);
    
    MPI_Barrier(MPI_COMM_WORLD);
    // --------------------------------------------------

    //  pMeshDataId tagGlob = MD_lookupMeshDataId("IdGlobal");
    //  UpdateIDGlobal(mesh,IdGlobal);    
    //  MPI_Barrier(MPI_COMM_WORLD);

    pMeshDataId tagEdge = MD_lookupMeshDataId("WriteEdge"); 
    E_createInfoInterface(mesh,tagEdge);    
    MPI_Barrier(MPI_COMM_WORLD);

    pMeshDataId tagFaceInterface = MD_lookupMeshDataId("FaceInterface"); 
    F_createInfoInterface(mesh,tagFaceInterface);
  
    // *****************ecriture des elt******************* 
    int nbClasEdges = 0;
    int nbClasFaces = 0;
 
    pMeshDataId tagFace = MD_newMeshDataId("WriteFace"); 
 
    {
      EIter eit = M_edgeIter(mesh);
      pEdge pe;  
      while ((pe = EIter_next(eit)))
        {
          int dim = GEN_type(pe->g); 
          if(dim == 1) {
            void * tmpptr;
            int is = EN_getDataPtr((pEntity) pe , tagEdge,&tmpptr);
            if(!is){
              nbClasEdges++;
            } else {
              std::map<int,pEdge> *recup = (std::map<int,pEdge> *) tmpptr;
              int minproc = myrank;
              for( std::map<int,pEdge>::const_iterator iter=(*recup).begin() ; 
                   iter!=(*recup).end() ; iter++){
                minproc = std::min(minproc,(*iter).first);
              }
              if(minproc==myrank) {
                nbClasEdges++;
              } 
            }
          }   
        }
      EIter_delete(eit);
    }    
    {
      FIter fit = M_faceIter(mesh);
      pFace pface;  
      while ((pface = FIter_next(fit)))
        {
          int dim = GEN_type(pface->g);   
          if(dim == 2){
            void * list;
            int is = EN_getDataPtr((pEntity) pface , tagFaceInterface,&list);
            if(!is){
              nbClasFaces++;
              EN_attachDataInt((pEntity) pface,tagFace,1);
            } else {
              std::vector<int> *recup = (std::vector<int> *) list;
              int minproc = myrank;
              for(unsigned int i=0 ; i<(*recup).size() ; i++){
                minproc = std::min(minproc,(*recup)[i]);
              }
              if(minproc==myrank) {
                nbClasFaces++;
                EN_attachDataInt((pEntity) pface,tagFace,1);
              } 
            }
          }
        }
      FIter_delete(fit);
    } 

    sendnum = nbClasEdges + nbModelVertex + nbClasFaces + mesh->nbTets;

    // the first proc collects the number of elements from the others...
    if(myrank) {
      MPI_Send(&sendnum,1,MPI_INT,0,myrank,MPI_COMM_WORLD);   
    } else {
      MPI_Status status;
      tab[0] = sendnum;
      for(int i=1 ; i<nproc ; i++) {
        MPI_Recv(&tab[i],1,MPI_INT,i,i,MPI_COMM_WORLD,&status);
      }
    }
  
    // printf("Before:  Proc %d: tab:\t%d\t%d\n",myrank,tab[0],tab[1]);

    // ... and sends its results to everyone
    MPI_Bcast(&tab[0],nproc,MPI_INT,0,MPI_COMM_WORLD);   

    //  printf("After:   Proc %d: tab:\t%d\t%d\n",myrank,tab[0],tab[1]);

  
    // --------------------------------------------------
    // Defines a control point for a sequential writting:
    //    Proc n+1 does not continue while proc n is 
    //    not at the next control point
    send = 0,recv = 0;
    if(!myrank) recv = 1;
  
    while(!recv) {  
      if(myrank!=(nproc-1)) MPI_Send(&send,1,MPI_INT,(myrank+1)%nproc,(myrank+1)%nproc,MPI_COMM_WORLD);
      MPI_Recv(&recv,1,MPI_INT,(myrank+nproc-1)%nproc,myrank,MPI_COMM_WORLD,&status);
    }
  
    nseek = recv;
    nplace = nseek;
    // --------------------------------------------------
  
    if(!myrank) {
      nplace = npt;
      lseek(nop,nplace,0);
      char nb[256];
      if ( version == 1 ) sprintf(nb,"$ENDNOD\n");
      else                sprintf(nb,"$EndNodes\n");
      int size = strlen(nb);
      write(nop,nb,size);
      nplace += size;
    }

    // ----------------------------------------------------
    // Write partitionning of elements (for msh1 only)
    // ----------------------------------------------------

    if ( version == 1 ) {

      if(!myrank) {
        write(nop,"$PARALLELELM\n",sizeof("$PARALLELELM"));
        nplace += sizeof("$PARALLELELM");
        char nb[256];
        int nbtot = 0;
        for(int i=0 ; i<nproc ; i++) nbtot +=tab[i];
        sprintf(nb,"%d\n",nbtot);
        int size = strlen(nb);
        nplace += size;
        write(nop,nb,size);
      }
      if(myrank) {
        lseek(nop,nseek,0);
      }

      int k = 1;
      if(myrank) {
        for(int i=0 ; i<myrank ; i++) k+= tab[i];
      }  
      {
        VIter vit = M_vertexIter(mesh);
        pVertex pv;  
        while ((pv = VIter_next(vit)))
          {
            int dim = GEN_type(pv->g); 
            if(dim == 0) {
              char nb[256];
              void *temp_ptr; 
              int isInterface = EN_getDataPtr((pEntity) pv , tagData, &temp_ptr);
              if(isInterface) {
                const std::vector<std::pair<int, MDB_Point*> > *recup = (const std::vector<std::pair<int, MDB_Point*> > *) temp_ptr;
                int minproc = nproc + 10;
                for(unsigned int j=0 ; j<(*recup).size() ; j++) minproc = std::min(minproc,(*recup)[j].first);
                minproc = std::min(myrank,minproc);
                if(minproc!=myrank) continue;    
                int size = (*recup).size();
                for(int j=0 ; j<size ; j++) {
                  sprintf(nb,"%d %d %d ",k++, size + 1,myrank);
                  char nbtmp[256];
                  if(j == (size-1) )  sprintf(nbtmp,"%d \n",(*recup)[j].first); 
                  else                sprintf(nbtmp,"%d ",(*recup)[j].first); 
                  int sizetmp = strlen(nbtmp); 
                  strncat(nb,nbtmp,sizetmp);
                }  
              } else {
                sprintf(nb,"%d %d %d \n",k++, 1,myrank);
              }
              int size = strlen(nb);
              nplace += size;
              write(nop,nb,size);
            }
          }
        VIter_delete(vit);   
      }
      { 
        EIter eit = M_edgeIter(mesh);
        pEdge pe;  
        while ((pe = EIter_next(eit)))
          {
            int dim = GEN_type(pe->g); 
            if(dim == 1) {
              void* tmpptr;
              int isMarked = EN_getDataPtr((pEntity) pe,tagEdge,&tmpptr);
              char nb[256];
              if(!isMarked) {
                sprintf(nb,"%d %d %d \n", k++, 1,myrank);
              } else {
                std::map<int,pEdge> *recup = (std::map<int,pEdge> *) tmpptr;
                int minproc = myrank;
                for(std::map<int,pEdge>::const_iterator iter=(*recup).begin() ;
                    iter!=(*recup).end() ; iter++){
                  minproc = std::min(minproc,(*iter).first);
                }
                if(minproc!=myrank) continue;
                sprintf(nb,"%d %d %d ", k++, (*recup).size() + 1,myrank);
                for(std::map<int,pEdge>::const_iterator iter=(*recup).begin() ;
                    iter!=(*recup).end() ; iter++) {
                  char nbtmp[256];
                  if((*iter).first == ((*recup).size()-1) ){
                    sprintf(nbtmp,"%d \n",(*iter).first);
                  } else{
                    sprintf(nbtmp,"%d ",(*iter).first);
                  }
                  int sizetmp = strlen(nbtmp); 
                  strncat(nb,nbtmp,sizetmp);  
                }
              }
              int size = strlen(nb);
              nplace += size;
              write(nop,nb,size);
            }  
          }
        EIter_delete(eit);
      }
      {
        FIter fit = M_faceIter(mesh);
        pFace pf;  
        while ((pf = FIter_next(fit)))
          {
            int ntmp;
            int isMarked = EN_getDataInt((pEntity) pf,tagFace,&ntmp);
            if(!isMarked) continue;
            int dim = GEN_type(pf->g); 
            if(dim == 2) {
              void* tmpptr;
              int isInt = EN_getDataPtr((pEntity) pf,tagFaceInterface,&tmpptr);
              char nb[256];
              if(!isInt) {
                sprintf(nb, "%d %d %d \n", k++, 1,myrank);
              } else {
                std::vector<int> *recup = (std::vector<int> *) tmpptr;
                int size = (*recup).size();    
                sprintf(nb,"%d %d %d ", k++, size + 1,myrank);
                for(int j=0 ; j<size ; j++) {
                  char nbtmp[256];
                  if(j == (size-1) )  sprintf(nbtmp,"%d \n",(*recup)[j]); 
                  else                sprintf(nbtmp,"%d ",(*recup)[j]); 
                  int sizetmp = strlen(nbtmp); 
                  strncat(nb,nbtmp,sizetmp);
                }  
              }  
              int size = strlen(nb);
              nplace += size;
              write(nop,nb,size);
            }  
          }
        FIter_delete(fit);
      } 
      {
        RIter rit = M_regionIter(mesh);
        pRegion pr;  
        while ((pr = RIter_next(rit))) {
          char nb[256];
          sprintf(nb,"%d %d %d \n", k++,1,myrank );
          int size = strlen(nb);
          nplace += size;
          write(nop,nb,size);
        }
        RIter_delete(rit);
      }
 
      // --------------------------------------------------
      // End of the control point
      MPI_Send(&nplace,1,MPI_INT,(myrank+1)%nproc,(myrank+1)%nproc,MPI_COMM_WORLD);
    
      if(!myrank) MPI_Recv(&npt,1,MPI_INT,(myrank+nproc-1)%nproc,myrank,MPI_COMM_WORLD,&status);
    
      MPI_Barrier(MPI_COMM_WORLD);
      // --------------------------------------------------

      // --------------------------------------------------
      // Defines a control point for a sequential writting:
      //    Proc n+1 does not continue while proc n is 
      //    not at the next control point
      send = 0,recv = 0;
      if(!myrank) recv = 1;
    
      while(!recv) {  
        if(myrank!=(nproc-1)) MPI_Send(&send,1,MPI_INT,(myrank+1)%nproc,(myrank+1)%nproc,MPI_COMM_WORLD);
        MPI_Recv(&recv,1,MPI_INT,(myrank+nproc-1)%nproc,myrank,MPI_COMM_WORLD,&status);
      }
    
      nseek = recv;
      nplace = nseek;
      // --------------------------------------------------
       
      if(!myrank) {
        nplace = npt;
        lseek(nop,nplace,0);
        write(nop,"$ENDPARALLELELM\n",sizeof("$ENDPARALLELELM"));
        nplace += sizeof("$ENDPARALLELELM");
      }

    }

    // ----------------------------------------------------
    // Write elements
    // ----------------------------------------------------

    if(!myrank) {
      char nb[256];
      if ( version == 1 ) { sprintf(nb,"$ELM\n"); }
      else                { sprintf(nb,"$Elements\n"); }
      int size = strlen(nb);
      write(nop,nb,size);
      nplace += size;
    }
    if(!myrank) {
      char nb[256];
      int nbtot = 0;
      for(int i=0 ; i<nproc ; i++) nbtot += tab[i];
      sprintf(nb,"%d\n",nbtot);
      int size = strlen(nb);
      nplace += size;
      write(nop,nb,size);
    }
    if(myrank) {
      lseek(nop,nseek,0);
    }

    int k = 1;
    // Create a reverse map of physical tags
    GEntity2Physical gentity2phys(mesh->geomFeatures_Tags);
    if(myrank) {
      for(int i=0 ; i<myrank ; i++) k+= tab[i];
    }
    // --- Nodes ---
    {
      VIter vit = M_vertexIter(mesh);
      pVertex pv;  
      while ((pv = VIter_next(vit)))
        {
          int dim = GEN_type(pv->g); 
          int tag = GEN_tag (pv->g); 
          int phys = gentity2phys.get_first_tag(pv->g);
          if(dim == 0) {
            char nb[256];
            //int nId = IdGlobal;
            void *temp_ptr; 
            int isInterface = EN_getDataPtr((pEntity) pv, tagData, &temp_ptr);
            if(isInterface) {
              const std::vector<std::pair<int , MDB_Point*> > *recup = (const std::vector<std::pair<int , MDB_Point*> > *) temp_ptr;
              int minproc = nproc + 10;
              for(unsigned int j=0 ; j<(*recup).size() ; j++) minproc = std::min(minproc,(*recup)[j].first);
              minproc = std::min(myrank,minproc);
              if(minproc!=myrank) continue;    
            }
            if ( version == 1 ) {
              sprintf(nb,"%d %d %d %d %d %d\n",
                      k++, 15, phys, tag, 1, pv->iD /*+ nId*/);
            } else {
              int nbTags = 3;
              sprintf(nb,"%d %d %d %d %d %d %d\n",
                      k++, 15, nbTags, phys, tag, myrank+1, pv->iD /*+ nId*/);
            }
            int size = strlen(nb);
            nplace += size;
            write(nop,nb,size);
          }
        }
      VIter_delete(vit);   
    }
    // --- Edges ---
    { 
      EIter eit = M_edgeIter(mesh);
      pEdge pe;  
      while ((pe = EIter_next(eit)))
        {
          int dim = GEN_type(pe->g); 
          int tag = GEN_tag (pe->g); 
          int phys = gentity2phys.get_first_tag(pe->g);
          if(dim == 1) {
            char nb[256];
            pVertex p1 = pe->p1;
            pVertex p2 = pe->p2;
            //int nId1 = pe->p1->iD /*+ IdGlobal*/,nId2 = pe->p2->iD /*+ IdGlobal*/;
            //void *temp_ptr1,*temp_ptr2; 
            //int isInterface1 = EN_getDataPtr((pEntity) p1 , tagData, &temp_ptr1);
            //int isInterface2 = EN_getDataPtr((pEntity) p2 , tagData, &temp_ptr2);
            //if(isInterface1) {
            //  int isGlob;
            //  int is = EN_getDataInt((pEntity) p1, tagGlob, &isGlob);
            //  assert(is);
            //  nId1 = --isGlob;
            //}
            //if(isInterface2) {
            //  int isGlob;
            //  int is = EN_getDataInt((pEntity) p2, tagGlob, &isGlob);
            //  assert(is);
            //  nId2 = --isGlob;
            //}
            void * tmpptr;
            int isMarked = EN_getDataPtr((pEntity) pe, tagEdge, &tmpptr);
            if(!isMarked) {
              if ( version == 1 ) {
                sprintf(nb,"%d %d %d %d %d %d %d\n", 
                        k++, 1, phys, tag, 2,p1->iD, p2->iD /* nId1, nId2*/);
              } else {
                int nbTags = 3;
                sprintf(nb,"%d %d %d %d %d %d %d %d\n",
                        k++, 1, nbTags, phys, tag, myrank+1,p1->iD, p2->iD /* nId1, nId2*/);
              }
            } else {
              const std::map<int,pEdge> *recup = (const std::map<int,pEdge> *) tmpptr;
              int minproc = myrank;
              for(std::map<int,pEdge>::const_iterator iter=(*recup).begin() ;
                  iter!=(*recup).end() ; iter++){
                minproc = std::min(minproc,(*iter).first);
              }
              if(minproc!=myrank) continue;
              if ( version == 1 ) {
                sprintf(nb,"%d %d %d %d %d %d %d\n", 
                        k++, 1, phys, tag, 2, p1->iD, p2->iD /* nId1, nId2*/);
              } else {
                int nbTags = 3;
                sprintf(nb,"%d %d %d %d %d %d %d %d\n",
                        k++, 1, nbTags, phys, tag, myrank+1,p1->iD, p2->iD /* nId1, nId2*/);
              }
            }
            int size = strlen(nb);
            nplace += size;
            write(nop,nb,size);
          }  
        }
      EIter_delete(eit);
    }
    // --- Faces ---
    {
      FIter fit = M_faceIter(mesh);
      pFace pf;  
      while ((pf = FIter_next(fit)))
        {
          int ntmp;
          int isMarked = EN_getDataInt((pEntity) pf,tagFace,&ntmp);
          if(!isMarked) continue;
          MDB_Point *nod[3];
          int dim = GEN_type(pf->g); 
          int tag = GEN_tag (pf->g); 
          int phys = gentity2phys.get_first_tag(pf->g);
          pf->getNodes(nod);
          if(dim == 2) {
            char nb[256];
            //int nId[3];
            pf->getNodes(nod);
            //void *temp_ptr1;
            //for(int j=0 ; j<3 ; j++){ 
            //  int isInterface = EN_getDataPtr((pEntity) nod[j], tagData, &temp_ptr1);
            //  if(isInterface) {
            //    int isGlob;
            //    int is = EN_getDataInt((pEntity) nod[j], tagGlob, &isGlob);
            //    assert(is);
            //    nId[j] = --isGlob;
            //  }
            //  else nId[j] = nod[j]->iD + IdGlobal;
            //}
            if ( version == 1 ) {
              sprintf(nb, "%d %d %d %d %d %d %d %d\n",
                      k++, 2, phys,tag, 3, nod[0]->iD, nod[1]->iD,nod[2]->iD/*nId[0],  nId[1],  nId[2]*/);
            } else {
              int nbTags = 3;
              sprintf(nb, "%d %d %d %d %d %d %d %d %d\n",
                      k++, 2, nbTags, phys, tag, myrank+1, nod[0]->iD, nod[1]->iD,nod[2]->iD/* nId[0],  nId[1],  nId[2]*/);
            }
            int size = strlen(nb);
            nplace += size;
            write(nop,nb,size);
          }  
        }
      FIter_delete(fit);
    }
    // --- Regions ---
    {
      RIter rit = M_regionIter(mesh);
      pRegion pr;  
      while ((pr = RIter_next(rit)))
        {
          MDB_Point *nod[4];
          int tag = GEN_tag (pr->g);  
          int phys = gentity2phys.get_first_tag(pr->g);
          pPList ll = R_vertices(pr);       
          nod[0] =  (pVertex)PList_item(ll, 0);
          nod[1] =  (pVertex)PList_item(ll, 1);
          nod[2] =  (pVertex)PList_item(ll, 2);
          nod[3] =  (pVertex)PList_item(ll, 3);
          PList_delete(ll);
          char nb[256];
          // int nId[4];
          // void *temp_ptr1;
          // for(int j=0 ; j<4 ; j++){ 
          //  int isInterface = EN_getDataPtr((pEntity) nod[j] , tagData, &temp_ptr1);
          //   if(isInterface){
          //     int isGlob;
          //     int is = EN_getDataInt((pEntity) nod[j] , tagGlob , &isGlob);
          //     assert(is); 
          //     nId[j] = --isGlob;
          //   }  
          //   else nId[j] = nod[j]->iD + IdGlobal;
          // }
          if ( version == 1 ) {
            sprintf(nb, "%d %d %d %d %d %d %d %d %d\n",
                    k++, 4, phys, tag, 4, nod[0]->iD, nod[1]->iD,nod[2]->iD,nod[3]->iD/* nId[0], nId[1], nId[2], nId[3]*/);
          } else {
            int nbTags = 3;
            sprintf(nb, "%d %d %d %d %d %d %d %d %d %d\n",
                    k++, 4, nbTags, phys, tag, myrank+1, nod[0]->iD, nod[1]->iD,nod[2]->iD,nod[3]->iD/*nId[0], nId[1], nId[2], nId[3]*/);
          }
          int size = strlen(nb);
          nplace += size;
          write(nop,nb,size);
        }
      RIter_delete(rit);
    }
    
    // --------------------------------------------------
    // End of the control point
    MPI_Send(&nplace,1,MPI_INT,(myrank+1)%nproc,(myrank+1)%nproc,MPI_COMM_WORLD);
  
    if(!myrank) MPI_Recv(&npt,1,MPI_INT,(myrank+nproc-1)%nproc,myrank,MPI_COMM_WORLD,&status);
    
    MPI_Barrier(MPI_COMM_WORLD);  
    // --------------------------------------------------  

    if(!myrank) {
      nplace = npt;
      lseek(nop,nplace,0);
      char nb[256];
      if ( version == 1 ) { sprintf(nb,"$ENDELM\n"); }
      else                { sprintf(nb,"$EndElements\n"); }
      int size = strlen(nb);
      write(nop,nb,size);
      nplace += size;
    }    
        
    close(nop);
  
    //vit = M_vertexIter(mesh); 
    //while ((pv=VIter_next(vit))){
    //  int isGlob;
    //  int is = EN_getDataInt((pEntity) pv , tagGlob , &isGlob);
    //  if(is) EN_deleteData((pEntity) pv,tagGlob); 
    //}
    //VIter_reset(vit);  
    //  MD_deleteMeshDataId(tagGlob);

    EIter eit = M_edgeIter(mesh);
    pEdge pe; 
    while ((pe=EIter_next(eit))){
      void* is;
      int isMarked = EN_getDataPtr((pEntity) pe , tagEdge , &is);
      if(isMarked) EN_deleteData((pEntity) pe,tagEdge); 
    }
    EIter_reset(eit);  
    MD_deleteMeshDataId(tagEdge);
  
    FIter fit = M_faceIter(mesh); 
    pFace pface;
    while ((pface=FIter_next(fit))){
      int is;
      int isMarked = EN_getDataInt((pEntity) pface , tagFace , &is);
      if(isMarked) EN_deleteData((pEntity) pface,tagFace); 
      void* isp;
      isMarked = EN_getDataPtr((pEntity) pface , tagFaceInterface , &isp);
      if(isMarked) EN_deleteData((pEntity) pface,tagFaceInterface); 
    }
    FIter_reset(fit);  
    MD_deleteMeshDataId(tagFace);  
    MD_deleteMeshDataId(tagFaceInterface); 
    delete []tab;
    delete []tabmax;
  }

#endif  

}

// -------------------------------------------------------------------
// ---------------- LoadGmshParallelOld ------------------------------
// -------------------------------------------------------------------
/*
void LoadGmshParallelOld (pMesh m,const char *filename,const int numproc,int version)
{  
  FILE *fp = fopen(filename, "r");
  if(!fp)
    {
      Msg(MDB_FATAL,"Unknown File %s\n",filename);
    }
  char String[256];
  //  if (!m->model)
  m->model = new NullModel; 
  int  Nbr_local = 0;
  int *tabelt = NULL;
#ifdef PARALLEL 
  pMeshDataId tagData = MD_lookupMeshDataId("RemotePoint");
#endif 
  int nread = 0;    
  // Interfaces NODES
  while(1) {
    do {
      if(!fgets(String, sizeof(String),fp))
        break;
      if(feof(fp))
        break;
    } while(String[0] != '$' || 
            strncmp(&String[1], "PARALLELNOD",3));
   
    if(feof(fp) || (nread==2)){ 
      break;
    }
    else {   
      if(!strncmp(&String[1], "PARALLELNOD", 11)) {
        nread ++;
        int Nbr_NodesInterfaces,Num,Nb;
        fscanf(fp, "%d", &Nbr_NodesInterfaces);
        for(int i_Node = 0; i_Node < Nbr_NodesInterfaces; i_Node++) {
          fscanf(fp, "%d %d", &Num, &Nb);
          if(Nb == 1) {
            int nproc;
            fscanf(fp,"%d",&nproc);
            if(nproc == numproc) {
              Nbr_local++;
              m->add_point(Num,0.,0.,0.);
            }
          } else {
            int *listproc=new int[Nb];
            int ismine = 0;
            for(int j=0 ; j < Nb ; j++) {
              int nproc;
              fscanf(fp,"%d",&nproc);
              listproc[j] = nproc;
              if(nproc == numproc) {
                Nbr_local++;
                ismine++;
                m->add_point(Num,0.,0.,0.);
              }
            }
#ifdef PARALLEL
            if(ismine) {
              std::vector<std::pair<int , MDB_Point*> > remotePoints;
              for(int j=0 ; j < Nb ; j++) {
                if(listproc[j]!=numproc) {
                  (remotePoints).push_back(std::pair <int , MDB_Point*>(listproc[j],NULL));
                }  
              }
              // attach remotePoints
              assert(remotePoints.size()==(unsigned int)(Nb-1));
              MDB_Point * vt = m->find_point(Num);
              EN_attachDataPtr((pEntity) vt , tagData, 
                               new std::vector<std::pair<int , MDB_Point*> >(remotePoints));
            }
#endif    
            delete []listproc;
          }
        }
        printf("%d Nodes in this partition\n",Nbr_local);
        if(!Nbr_local) Nbr_local++;
      } else if(!strncmp(&String[1], "PARALLELELM", 11)) {
        nread++;
        int Nbr_Elt_local = 0,Num,Nb,Nbr_Elt;
        fscanf(fp, "%d", &Nbr_Elt);
        tabelt = (int *) calloc(Nbr_Elt,sizeof(int));
        for(int i = 0; i < Nbr_Elt; i++) {
          fscanf(fp, "%d %d", &Num, &Nb);
          tabelt[i] = -1;
          if(Nb == 1) {
            int nproc;
            fscanf(fp,"%d",&nproc);
            if(nproc == numproc) {
              Nbr_Elt_local++;
              tabelt[i] = Num;
            }
          } else {
            int *listproc=new int[Nb];
            for(int j=0 ; j < Nb ; j++) {
              int nproc;
              fscanf(fp,"%d",&nproc);
              listproc[j] = nproc;
              if(nproc == numproc) {
                Nbr_Elt_local++;
                tabelt[i] = Num;
              }
            }  
            delete []listproc;
          }
        }
        printf("%d Elt in this partition\n",Nbr_Elt_local);
      }
    }   
  }

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
    
    
    // NODES
    if(!strncmp(&String[1], "NOD", 3) ||
       !strncmp(&String[1], "NOE", 3) ||
       !strncmp(&String[1], "Nodes", 5)) {
      
      int Nbr_Nodes,Num,ncurc = 0;
      double x,y,z;
      fscanf(fp, "%d", &Nbr_Nodes);
      for(int i_Node = 0; i_Node < Nbr_Nodes; i_Node++) {
        fscanf(fp, "%d %lf %lf %lf", &Num, &x, &y, &z);
  if(Nbr_local) {
    MDB_Point * p = m->find_point(Num);
          if(p) {
      ncurc++;
      p->X = x;
      p->Y = y;
      p->Z = z;
    }  
  } else {
    m->add_point(Num,x,y,z);
  }
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
        //  int myrank;
  //      MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

        for(int j = 0; j < Nbr_Nodes; j++)
          fscanf(fp, "%d", &verts[j]);

  GEntity *geom = 0;
        switch (Type) {

          // -------------------------------------------------------------------
          // linear elements
          // -------------------------------------------------------------------
          
        case MSH_LIN_2:
    {
      if(Nbr_local) {
        if(tabelt[i_Element]>=0) {
          if(tabelt[i_Element]!=Num) printf(" %d tab %d num %d\n",i_Element,tabelt[i_Element],Num);
    assert(tabelt[i_Element]==Num);
                for (int k=0;k<2;k++) assert(m->find_point(verts[k]));
                geom = m->model->edgeByTag(Elementary);
          m->add_edge(verts[0],verts[1],geom);  
        }       
      } else {
        GEdge* geom = m->model->edgeByTag(Elementary);
              m->add_edge(verts[0],verts[1],geom);  
      }  
    }
    break;
        case MSH_LIN_3:
    {
      if(Nbr_local) {
        if(tabelt[i_Element]>=0) {
          if(tabelt[i_Element]!=Num) printf(" %d tab %d num %d\n",i_Element,tabelt[i_Element],Num);
    assert(tabelt[i_Element]==Num);
                assert(m->find_point(verts[0]) && m->find_point(verts[1]) && m->find_point(verts[2]));
                geom = m->model->edgeByTag(Elementary);
                m->add_edge(3,geom,verts[0],verts[2],verts[1]);
        }       
      } else {
        geom = m->model->edgeByTag(Elementary);
        m->add_edge(3,geom,verts[0],verts[2],verts[1]);
      }  
    }
    break;
          
        case MSH_LIN_4:
    {
      if(Nbr_local) {
        if(tabelt[i_Element]>=0) {
          if(tabelt[i_Element]!=Num) printf(" %d tab %d num %d\n",i_Element,tabelt[i_Element],Num);
    assert(tabelt[i_Element]==Num);
                assert(m->find_point(verts[0]) && m->find_point(verts[1]));
                geom = m->model->edgeByTag(Elementary);
        m->add_edge(4,geom, verts[0],verts[2],verts[3],verts[1]);
        }       
      } else {
        geom = m->model->edgeByTag(Elementary);
          m->add_edge(4,geom, verts[0],verts[2],verts[3],verts[1]);
      }  
    }
    break;
        case MSH_LIN_5:
    {
      if(Nbr_local) {
        if(tabelt[i_Element]>=0) {
          if(tabelt[i_Element]!=Num) printf(" %d tab %d num %d\n",i_Element,tabelt[i_Element],Num);
    assert(tabelt[i_Element]==Num);
                assert(m->find_point(verts[0]) && m->find_point(verts[1]));
                geom = m->model->edgeByTag(Elementary);
        m->add_edge(5,geom, verts[0],verts[2],verts[3],verts[4],verts[1]);
        }       
      } else {
        geom = m->model->edgeByTag(Elementary);
          m->add_edge(5,geom, verts[0],verts[2],verts[3],verts[4],verts[1]);
      }  
    }
    break;
        case MSH_LIN_6:
    {
      if(Nbr_local) {
        if(tabelt[i_Element]>=0) {
          if(tabelt[i_Element]!=Num) printf(" %d tab %d num %d\n",i_Element,tabelt[i_Element],Num);
    assert(tabelt[i_Element]==Num);
                assert(m->find_point(verts[0]) && m->find_point(verts[1]));
                geom = m->model->edgeByTag(Elementary);
        m->add_edge(6,geom, verts[0],verts[2],verts[3],verts[4],verts[5],verts[1]);
        }       
      } else {
        geom = m->model->edgeByTag(Elementary);
          m->add_edge(6,geom, verts[0],verts[2],verts[3],verts[4],verts[5],verts[1]);
      }  
    }
    break;
          
          // -------------------------------------------------------------------
          // TRIANGLES
          // -------------------------------------------------------------------
          
        case MSH_TRI_3:
    {
      if(Nbr_local) {
        if(tabelt[i_Element]>=0) {
          assert(tabelt[i_Element]==Num);
          assert(m->find_point(verts[0]) && m->find_point(verts[1]) && m->find_point(verts[2]));
                 geom = m->model->faceByTag(Elementary);
          m->add_triangle(verts[0],verts[1],verts[2],geom);   
        }
      } else {
        geom = m->model->faceByTag(Elementary);
        m->add_triangle(verts[0],verts[1],verts[2],geom);
      }   
    }
    break;

          // order 2 
          
        case MSH_TRI_6:
    { 
      if(Nbr_local) {
        if(tabelt[i_Element]>=0) {
          assert(tabelt[i_Element]==Num);
          assert(m->find_point(verts[0]) && m->find_point(verts[1]) && m->find_point(verts[2]));
                geom = m->model->faceByTag(Elementary);
        m->add_triangle(2,0,geom,verts[0],verts[3],verts[1],verts[4],verts[2],verts[5]); 
        }
      } else {
        geom = m->model->faceByTag(Elementary);
          m->add_triangle(2,0,geom,verts[0],verts[3],verts[1],verts[4],verts[2],verts[5]); 
      }   
    }
    break;

          // order 3
          
        case MSH_TRI_9:
          {       
            if (Nbr_local) {
              
        if(tabelt[i_Element]>=0) {
          assert(tabelt[i_Element]==Num);
                for (int i=0;i<9;i++) assert(m->find_point(verts[i]));
                geom = m->model->faceByTag(Elementary);
                m->add_triangle(3,false,geom,verts[0],verts[3],verts[4],verts[1],verts[5],verts[6],verts[2],verts[7],verts[8]); 
        }
      }
            else {
              geom = m->model->faceByTag(Elementary);
              m->add_triangle(3,false,geom,verts[0],verts[3],verts[4],verts[1],verts[5],verts[6],verts[2],verts[7],verts[8]); 
      }  
            break;
          }

          // Serendipity order 4
          
        case MSH_TRI_12:
          {
            if (Nbr_local) {
              
        if(tabelt[i_Element]>=0) {
          assert(tabelt[i_Element]==Num);
                for (int i=0;i<12;i++) assert(m->find_point(verts[i]));
                geom = m->model->faceByTag(Elementary);
                m->add_triangle(4,false,geom,verts[0],verts[3],verts[4],verts[5],verts[1],verts[6],verts[7],verts[8],verts[2],verts[9],verts[10],verts[11]); 
        }
      }
            else {
              geom = m->model->faceByTag(Elementary);
              m->add_triangle(4,false,geom,verts[0],verts[3],verts[4],verts[5],verts[1],verts[6],verts[7],verts[8],verts[2],verts[9],verts[10],verts[11]); 
      }  
            break;
          }

          // Complete 
          
        case MSH_TRI_15:
          {
            if (Nbr_local) {
  
        if(tabelt[i_Element]>=0) {
          assert(tabelt[i_Element]==Num);
                for (int i=0;i<24;i++) assert(m->find_point(verts[i]));
                geom = m->model->faceByTag(Elementary);
                m->add_triangle(4,true,geom,verts[0],verts[3],verts[4],verts[5],verts[6],verts[1],verts[7],verts[8],verts[9],verts[10],verts[2],verts[11],verts[12],verts[13],verts[14]); 
        }
      }
            else {
              geom = m->model->faceByTag(Elementary);
              m->add_triangle(4,true,geom,verts[0],verts[3],verts[4],verts[5],verts[6],verts[1],verts[7],verts[8],verts[9],verts[10],verts[2],verts[11],verts[12],verts[13],verts[14]); 
      }
            break;
          }
          
          // -------------------------------------------------------------------
          // Tetrahedra
          // -------------------------------------------------------------------
          
          // linear
          
        case MSH_TET_4:
          {
            
            if(Nbr_local) {
              if(tabelt[i_Element]>=0) {
                assert(tabelt[i_Element]==Num);
                assert(m->find_point(verts[0]) && m->find_point(verts[1]) && m->find_point(verts[2]) &&m->find_point(verts[3]));
                geom = m->model->regionByTag(Elementary);
                m->add_tet(verts[0],verts[1],verts[2],verts[3],geom); 
        }
            } else {
              geom = m->model->regionByTag(Elementary);
              m->add_tet(verts[0],verts[1],verts[2],verts[3],geom); 
            }  
          }
          break;

          // quadratic
          
        case MSH_TET_10:
          {
            if (Nbr_local) {
  
        if(tabelt[i_Element]>=0) {
          assert(tabelt[i_Element]==Num);
                for (int i=0;i<10;i++) assert(m->find_point(verts[i]));
                geom = m->model->regionByTag(Elementary);
                m->add_tet(geom,2,false,verts);
        }
      }
            else {
              geom = m->model->regionByTag(Elementary);
              m->add_tet(geom,2,false,verts);
      }
            break;
          }
          // -------------------------------------------------------------------
          // Hexahedra
          // -------------------------------------------------------------------
        
        case MSH_HEX_8:
          {
            
            if(Nbr_local) { cout << "not yet implemented!!!" << endl; throw; }
            geom = m->model->regionByTag(Elementary);
            m->add_hex(verts[0],verts[1],verts[2],verts[3],
                       verts[4],verts[5],verts[6],verts[7],geom); 
          }
          break;

          // -------------------------------------------------------------------
          // Prisms
          // -------------------------------------------------------------------
          
        case MSH_PRI_6:
          {
            
            if(Nbr_local) { cout << "not yet implemented!!!" << endl; throw; }
            geom = m->model->regionByTag(Elementary);
            m->add_prism(verts[0],verts[1],verts[2],
                         verts[3],verts[4],verts[5],geom); 
          }
          break;
          
          // -------------------------------------------------------------------
          // node 
          // -------------------------------------------------------------------
          
        case 15:
    {
      if(Nbr_local) {
        MDB_Point *p = m->find_point(verts[0]);
        if(p) {
          assert(tabelt[i_Element]==Num);
          GVertex *gv = m->model->vertexByTag(Elementary);
    p->g = gv;
        }
      } else {
        geom = m->model->vertexByTag(Elementary);
        MDB_Point *p = m->find_point(verts[0]);
        p->g = geom;
      }  
    }
    break;
        default:
    throw;
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
    
    do {
      if(!fgets(String, sizeof(String), fp))
  throw;
      if(feof(fp))
  throw;
    } while(String[0] != '$');    
  }

  m->classify_unclassified_entities();
    
#ifdef PARALLEL   
//   linkRemotePoints(m);
  pMeshDataId tagVertex = MD_lookupMeshDataId("RemotePoint");
  V_createInfoInterface(m, tagVertex);
#endif
}
*/

