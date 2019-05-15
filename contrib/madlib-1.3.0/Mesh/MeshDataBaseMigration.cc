// -------------------------------------------------------------------
// MAdLib - Copyright (C) 2008-2009 Universite catholique de Louvain
//
// See the Copyright.txt and License.txt files for license information. 
// You should have received a copy of these files along with MAdLib. 
// If not, see <http://www.madlib.be/license/>
//
// Please report all bugs and problems to <contrib@madlib.be>
//
// Authors: Josue Barboza
// -------------------------------------------------------------------

#ifdef PARALLEL
#include "mpi.h"
#include "autopack.h"
#include "MeshDataBaseParallelInterface.h"
#include "MeshDataBaseLoadBalancing.h"
#include "MeshDataBase.h"
#include "assert.h"
#include <iostream>

namespace MAd {

  typedef std::set< int >                     SetOfDestination_type;
  // ------------------------------------------------------------------- 
  // -------------------------------------------------------------------
  struct coor_comm2
  {
    double   X,Y,Z;
    int tag,dim;
    int id;
  }; 

  // -------------------------------------------------------------------
  void MarkVertex( pVertex pv,  pMeshDataId tagD ,int d_proc )
  {
    void *temp_ptr;
    SetOfDestination_type *recup = NULL;
    int is = EN_getDataPtr((pEntity) pv, tagD, &temp_ptr);    
    if( !is ){
      recup = new SetOfDestination_type;
      EN_attachDataPtr((pEntity) pv , tagD, recup);
    }else
      { 
        recup = reinterpret_cast< SetOfDestination_type *> (temp_ptr);     
      }
    recup->insert( d_proc );
  }
  // -------------------------------------------------------------------
  void MarkEdge( pEdge pe,  pMeshDataId tagD ,int d_proc )
  {
    void *temp_ptr;
    SetOfDestination_type *recup = NULL;
    int is = EN_getDataPtr((pEntity) pe, tagD, &temp_ptr);    
    if( !is ){
      recup = new SetOfDestination_type;
      EN_attachDataPtr((pEntity) pe, tagD, recup);
    }else { 
      recup = reinterpret_cast<SetOfDestination_type *> (temp_ptr);     
    }
    recup->insert(d_proc);
    pVertex  nod[2];
    nod[0] = E_vertex(pe,0);
    nod[1] = E_vertex(pe,1);
    MarkVertex( nod[0], tagD, d_proc);
    MarkVertex( nod[1], tagD, d_proc);
  }
  // -------------------------------------------------------------------
  void MarkFace( pFace pf,  pMeshDataId tagD ,int d_proc )
  {
    void *temp_ptr;
    SetOfDestination_type *recup = NULL;
    int is = EN_getDataPtr((pEntity) pf, tagD, &temp_ptr);    
    if( !is ){
      recup = new SetOfDestination_type;
      EN_attachDataPtr((pEntity) pf, tagD, recup);
    }else { 
      recup = reinterpret_cast<SetOfDestination_type *> (temp_ptr);     
    }
    recup->insert(d_proc);
    for( int i = 0; i<(pf->getNbEdges()); i++) {
      MarkEdge( F_edge( pf,i ), tagD, d_proc );
    }
  }
  // -------------------------------------------------------------------
  void MarkRegion( pRegion pr,  pMeshDataId tagD ,int d_proc )
  {
    void *temp_ptr;
    SetOfDestination_type *recup = NULL;
    int is = EN_getDataPtr((pEntity) pr, tagD, &temp_ptr);    
    if( !is ){
      recup = new SetOfDestination_type;
      EN_attachDataPtr((pEntity) pr, tagD, recup);
    }else { 
      recup = reinterpret_cast<SetOfDestination_type *> (temp_ptr);     
    }
    recup->insert(d_proc);
    for( int i = 0; i<(pr->getNbFace()); i++) {
      MarkFace( R_face( pr,i ), tagD, d_proc );
    }
  }

  // -------------------------------------------------------------------
  void MarkEltSubEntities( pMesh mesh, pMeshDataId tagElt )
  {
    pMeshDataId tagDest = MD_lookupMeshDataId("tagDestinations");
    MD_deleteMeshDataId(tagDest);
  
    int dim = (mesh->tets.empty()) ? 2 : 3;
    if( dim == 2 ){
      FIter fit = M_faceIter(mesh);
      pFace pface;  
      while ((pface = FIter_next(fit))) {
        int dest; 
        int migre = EN_getDataInt((pEntity) pface, tagElt, &dest);
        if(!migre) continue;    
        MarkFace( pface, tagDest ,dest-1 );      
      }
      FIter_delete(fit);    
    }
    else {
      if(dim!=3) throw;
      RIter rit = M_regionIter(mesh);
      pRegion pr;  
      while ((pr = RIter_next(rit))) {
        int dest; 
        int migre = EN_getDataInt((pEntity) pr, tagElt, &dest);
        if(!migre) continue;  
        MarkRegion( pr, tagDest ,dest-1 );         
      }    
      RIter_delete(rit);
    }
  }
  // -------------------------------------------------------------------
  void * VertexAndDataPackaging( pVertex pv, int d, MDB_DataExchanger &de)
  {
    int sizebuf;
    void *msg = de.sendData ((pEntity) pv, d, sizebuf );
    void *buf = AP_alloc(d,de.tag(),sizeof(coor_comm2)+sizebuf);

    char *cast = reinterpret_cast< char *> (buf);  
    coor_comm2 castbuf; 
    castbuf.X     = P_x(pv);
    castbuf.Y     = P_y(pv);
    castbuf.Z     = P_z(pv);
    pGEntity  pg  = EN_whatIn(pv);
    castbuf.tag   = GEN_tag(pg);
    castbuf.dim   = GEN_type(pg);
    castbuf.id = EN_id( pv );
    memcpy(&cast[0],&castbuf,sizeof(coor_comm2));           
    memcpy(&cast[sizeof(coor_comm2)],msg,sizebuf);          
    free(msg);
    return buf;
  }
  // -------------------------------------------------------------------
  void VertexAndDataUnpackaging(  pMesh mesh, void *msg, int from, 
                                  MDB_DataExchanger &de)
  {
    char * castbuf = reinterpret_cast< char *> ( msg );
    coor_comm2 * coorcom = reinterpret_cast< coor_comm2 *> (castbuf);
    pGEntity pg = NULL;
    int dim  = coorcom->dim;
    int tag  = coorcom->tag;
    if(dim == 0) {
      pg = (pGEntity) GM_vertexByTag(mesh->model,tag);
    } else if(dim==1) {
      pg = (pGEntity) GM_edgeByTag(mesh->model,tag);     
    } else if(dim==2) {
      pg = (pGEntity) GM_faceByTag(mesh->model,tag);
    } else if(dim==3) {
      pg = (pGEntity) GM_regionByTag(mesh->model,tag);
    } 
    pVertex pnew =  mesh->find_point(coorcom->id);
    if( !pnew ){
      
      //       double  X = coorcom->X;
      //       double  Y = coorcom->Y;
      //       double  Z = coorcom->Z;
      //       pnew = M_createVP(mesh,X,Y,Z,coorcom->id,pg);
      
      double  XYZ[3] = {coorcom->X,coorcom->Y,coorcom->Z};
      pnew = M_createV2(mesh,XYZ,coorcom->id,pg);
    }
    de.receiveData (pnew,from, &castbuf[sizeof(coor_comm2)]);
  
  }
  // -------------------------------------------------------------------
  void MigrateVerticesAndData( pMesh mesh, pMeshDataId tagDest, 
                               MDB_DataExchanger &de)
  {
    int nproc,myrank;
  
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    int *sendcounts = new int[nproc];
    for( int i=0; i<nproc; i++){
      sendcounts[i]= 0;
    }
    VIter vit = M_vertexIter(mesh);
    pVertex pv;  
    while ((pv = VIter_next(vit))) {
      void *temp_ptr;
      int migre = EN_getDataPtr((pEntity) pv, tagDest, &temp_ptr);
      if(!migre) continue;
      const SetOfDestination_type *recup =
        reinterpret_cast< const SetOfDestination_type* > (temp_ptr);
      assert( (recup->size())>0 );
      SetOfDestination_type::const_iterator itd = recup->begin();
      for(; itd!=recup->end(); itd ++) {
        int dest = *itd;
        assert(dest!=myrank);
        void *buf = VertexAndDataPackaging( pv, dest, de );    
        AP_send(buf);
        sendcounts[dest]++; 
      } 
    }    
    VIter_delete(vit);
  

    AP_check_sends(AP_NOFLAGS);
    AP_reduce_nsends(sendcounts);
  
    /*receive Vertices and data*/
    int message=0;
    int count;
    while (!AP_recv_count(&count) || message<count) {
      void *msg;
      int  from;
      int  tag;
      int  size;
      int  recv;
      recv = AP_recv(MPI_ANY_SOURCE, de.tag(), AP_BLOCKING|AP_DROPOUT,
                     &msg, &size, &from, &tag);
      if(recv) {
        message++;
        VertexAndDataUnpackaging( mesh, msg, from, de);
        AP_free(msg);    
      }
    }
    AP_check_sends(AP_WAITALL);
    MPI_Barrier(MPI_COMM_WORLD);
    delete [] sendcounts;  
  }
  // -------------------------------------------------------------------
  void * EdgeAndDataPackaging( pEdge pe, int d, MDB_DataExchanger &de)
  {
    int sizebuf;
    void *msg = de.sendData ((pEntity) pe, d, sizebuf );
    int castbuf[4];
    int sizeOfcastbuf = sizeof(int) *4;
    void *buf = AP_alloc(d,de.tag(), sizeOfcastbuf +sizebuf );
    char *cast = reinterpret_cast< char *> (buf);  
    pGEntity pg = EN_whatIn( pe );
    castbuf[0] = GEN_tag( pg );
    castbuf[1] = GEN_type( pg ); 
    castbuf[2] = EN_id( E_vertex( pe, 0 ) );
    castbuf[3] = EN_id( E_vertex( pe, 1 ) ); 
    memcpy( &cast[0], castbuf, sizeOfcastbuf );             
    memcpy( &cast[sizeOfcastbuf],msg,sizebuf);
    free(msg);     
    return buf;      
  }
  // -------------------------------------------------------------------
  void EdgeAndDataUnpackaging(  pMesh mesh, void *msg, int from, 
                                MDB_DataExchanger &de)
  {
    char* castbuf = reinterpret_cast<char*> (msg);
    int * Edge_comm = reinterpret_cast<int*> (castbuf);
    int tag = *(Edge_comm++);
    int dim = *(Edge_comm++);
    pGEntity pg = NULL;
    if(dim==2) {
      pg = (pGEntity) GM_faceByTag(mesh->model,tag);
    } else if(dim==3) {
      pg = (pGEntity) GM_regionByTag(mesh->model,tag);
    }
  
    pVertex pv[2];
    pv[0]= mesh->find_point( *(Edge_comm++) );
    pv[1]= mesh->find_point( *(Edge_comm++) );


    pEdge pe =  E_exist(pv[0],pv[1]);
    if( !pe ) {
      pe = M_createE(mesh,pv[0],pv[1],pg);
    }
    assert(pe);
    de.receiveData (pe,from, &castbuf[sizeof(int)*4]);
  }
  // -------------------------------------------------------------------
  // -------------------------------------------------------------------
  void MigrateEdgesAndData(pMesh mesh, pMeshDataId tagDest,
                           MDB_DataExchanger &de)
  {
    int nproc,myrank;
  
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    int *sendcounts = new int[nproc];
    for( int i=0; i<nproc; i++){
      sendcounts[i]= 0;
    }
    EIter eit = M_edgeIter(mesh);
    pEdge pe;  
    while ((pe = EIter_next(eit))) {
      void *temp_ptr;
      int migre = EN_getDataPtr((pEntity) pe, tagDest, &temp_ptr);
      if(!migre) continue;
      const SetOfDestination_type *recup =
        reinterpret_cast< const SetOfDestination_type* > (temp_ptr);
      assert( (recup->size())>0 );
      SetOfDestination_type::const_iterator itd = recup->begin();
      for(; itd!=recup->end(); itd ++) {
        int dest = *itd;
        assert(dest!=myrank);
        void *buf = EdgeAndDataPackaging( pe, dest, de );    
        AP_send(buf);
        sendcounts[dest]++; 
      } 
    }    
    EIter_delete(eit);
  

    AP_check_sends(AP_NOFLAGS);
    AP_reduce_nsends(sendcounts);
  
    /*receive Edges and data*/
    int message=0;
    int count;
    while (!AP_recv_count(&count) || message<count) {
      void *msg;
      int  from;
      int  tag;
      int  size;
      int  recv;
      recv = AP_recv(MPI_ANY_SOURCE, de.tag(), AP_BLOCKING|AP_DROPOUT,
                     &msg, &size, &from, &tag);
      if(recv) {
        message++;
        EdgeAndDataUnpackaging( mesh, msg, from, de);
        AP_free(msg);    
      }
    }
    AP_check_sends(AP_WAITALL);
    MPI_Barrier(MPI_COMM_WORLD);
    delete [] sendcounts;  
  }
  // -------------------------------------------------------------------
  // -------------------------------------------------------------------
  void * FaceAndDataPackaging( pFace pf, int d, MDB_DataExchanger &de )
  {
    int sizebuf;
    void *msg = de.sendData ((pEntity) pf, d, sizebuf );
    int nV = F_numVertices( pf );
    int nComm = nV+3;
    int *castbuf = new int[nComm];
    int sizeOfcastbuf = sizeof(int) * nComm;
    void *buf = AP_alloc(d,de.tag(), sizeOfcastbuf +sizebuf );
    char *cast = reinterpret_cast< char *> (buf);  
    pGEntity pg = EN_whatIn( pf );
    castbuf[0] = GEN_tag( pg );
    castbuf[1] = GEN_type( pg ); 
    castbuf[2] = nV;
    for( int i = 0; i< nV; i++) {
      int vId = EN_id( F_vertex( pf, i ) );
      castbuf[3+i]= vId ;
    }   
    memcpy( &cast[0], castbuf, sizeOfcastbuf );             
    memcpy( &cast[sizeOfcastbuf],msg,sizebuf );
    delete [] castbuf;
    free(msg);     
    return buf;    
  }
  // -------------------------------------------------------------------
  void FaceAndDataUnpackaging(  pMesh mesh, void *msg, int from, 
                                MDB_DataExchanger &de )
  {
    /*! \TODO: only triangle is taken into account. Should be extended to other*/
    char* castbuf = reinterpret_cast< char* >( msg );
    int* pf_com = reinterpret_cast< int* > (castbuf);
    int tag = pf_com[0];
    int dim = pf_com[1];
    int nbV = pf_com[2];
    if( nbV != 3 ) {
      std::cout<<"Unpackaging Face that is diff. form triangle is not implemented: "
               <<nbV<<std::endl;
      assert(0);
    }
    pGEntity pg = NULL;
    if(dim==2) {
      pg = (pGEntity) GM_faceByTag(mesh->model,tag);
    } else if(dim==3) {
      pg = (pGEntity) GM_regionByTag(mesh->model,tag);
    }
  
    pVertex *ListOfpVertex = new pVertex[nbV];
    for( int i = 0; i< nbV; i++ ) {
      ListOfpVertex[i] = mesh->find_point( pf_com[3+i] );
    }
    pEdge    pe[3];
    pe[0] =  E_exist(ListOfpVertex[0],ListOfpVertex[1]);
    assert(pe[0]);
    pe[1] =  E_exist(ListOfpVertex[1],ListOfpVertex[2]);
    assert(pe[1]);
    pe[2] =  E_exist(ListOfpVertex[0],ListOfpVertex[2]);
    assert(pe[2]);

    // pFace pface =  F_exist(2,pe[0],pe[1],pe[2],0);
    pFace pface =  F_exist(pe[0],pe[1],pe[2],0);
    if( !pface ) {
      pface = M_createF(mesh,3,pe,pg);
    }
    assert(pface);
    de.receiveData (pface,from, &castbuf[sizeof(int)*(nbV+3)]);

    delete []  ListOfpVertex ;
  }
  // -------------------------------------------------------------------
  // -------------------------------------------------------------------
  void MigrateFacesAndData( pMesh mesh, pMeshDataId tagDest, 
                            MDB_DataExchanger &de)
  {
    int nproc,myrank;
  
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    int *sendcounts = new int[nproc];
    for( int i=0; i<nproc; i++){
      sendcounts[i]= 0;
    }
    FIter fit = M_faceIter(mesh);
    pFace pf;  
    while ((pf = FIter_next(fit))) {
      void *temp_ptr;
      int migre = EN_getDataPtr((pEntity) pf, tagDest, &temp_ptr);
      if(!migre) continue;
      const SetOfDestination_type *recup =
        reinterpret_cast< const SetOfDestination_type* > (temp_ptr);
      assert( (recup->size())>0 );
      SetOfDestination_type::const_iterator itd = recup->begin();
      for(; itd!=recup->end(); itd ++) {
        int dest = *itd;
        assert(dest!=myrank);
        void *buf = FaceAndDataPackaging( pf, dest, de );    
        AP_send(buf);
        sendcounts[dest]++; 
      }
    }    
    FIter_delete(fit);
  

    AP_check_sends(AP_NOFLAGS);
    AP_reduce_nsends(sendcounts);
  
    /*receive Faces and data*/
    int message=0;
    int count;
    while (!AP_recv_count(&count) || message<count) {
      void *msg;
      int  from;
      int  tag;
      int  size;
      int  recv;
      recv = AP_recv(MPI_ANY_SOURCE, de.tag(), AP_BLOCKING|AP_DROPOUT,
                     &msg, &size, &from, &tag);
      if(recv) {
        message++;
        FaceAndDataUnpackaging( mesh, msg, from, de);
        AP_free(msg);    
      }
    }
    AP_check_sends(AP_WAITALL);
    MPI_Barrier(MPI_COMM_WORLD);
    delete [] sendcounts;  
  }
  // -------------------------------------------------------------------
  // -------------------------------------------------------------------
  void * RegionAndDataPackaging( pRegion pr, int d, MDB_DataExchanger &de)
  {
    int sizebuf;
    void *msg = de.sendData ((pEntity) pr, d, sizebuf );
    int nV = R_numVertices( pr );
    int nComm = nV+3;
    int *castbuf = new int[nComm];
    int sizeOfcastbuf = sizeof(int) * nComm;
    void *buf = AP_alloc(d,de.tag(), sizeOfcastbuf +sizebuf );
    char *cast = reinterpret_cast<char *> (buf);  
    pGEntity pg = EN_whatIn( pr );
    castbuf[0] = GEN_tag( pg );
    castbuf[1] = GEN_type( pg ); 
    castbuf[2] = nV;
    for( int i = 0; i< nV; i++) {
      int vId = EN_id( R_vertex( pr, i ) );
      castbuf[3+i]= vId ;
    } 
    memcpy( &cast[0], castbuf, sizeOfcastbuf );             
    memcpy( &cast[sizeOfcastbuf],msg,sizebuf);
    delete [] castbuf;
    free(msg);     
    return buf;  
  }
  // -------------------------------------------------------------------
  void RegionAndDataUnpackaging( pMesh mesh, void *msg, int from, 
                                 MDB_DataExchanger &de)
  {
    /*! \TODO: only tethra is taken into account. Should be extended to other*/
    char* castbuf = reinterpret_cast<char*> (msg);
    int * Region_comm = reinterpret_cast<int*> (castbuf);
    int tag = *(Region_comm++);
    int dim = *(Region_comm++);
    int nbVertices = *(Region_comm++);
    if( nbVertices != 4 ) {
      std::cout<<"Unpackaging Region that is diff. form tethra is not implemented"
               <<nbVertices<<std::endl;
      assert(0);
    }
    pGEntity pg = NULL;
    if(dim==2) {
      pg = (pGEntity) GM_faceByTag(mesh->model,tag);
    } else if(dim==3) {
      pg = (pGEntity) GM_regionByTag(mesh->model,tag);
    }
  
    pVertex *ListOfpVertex = new pVertex[nbVertices];
    for( int i = 0; i< nbVertices; i++ ) {
      ListOfpVertex[i] = mesh->find_point( *(Region_comm++) );
    }


    pFace    pface[4];
    
    //     pface[0] =  F_exist(2,ListOfpVertex[0],ListOfpVertex[1],ListOfpVertex[2],0);
    //     assert(pface[0]);
    //     pface[1] =  F_exist(2,ListOfpVertex[0],ListOfpVertex[1],ListOfpVertex[3],0);
    //     assert(pface[1]);
    //     pface[2] =  F_exist(2,ListOfpVertex[1],ListOfpVertex[2],ListOfpVertex[3],0);
    //     assert(pface[2]); 
    //     pface[3] =  F_exist(2,ListOfpVertex[0],ListOfpVertex[2],ListOfpVertex[3],0);
    //     assert(pface[3]);
    
    pface[0] =  F_exist(ListOfpVertex[0],ListOfpVertex[1],ListOfpVertex[2],0);
    assert(pface[0]);
    pface[1] =  F_exist(ListOfpVertex[0],ListOfpVertex[1],ListOfpVertex[3],0);
    assert(pface[1]);
    pface[2] =  F_exist(ListOfpVertex[1],ListOfpVertex[2],ListOfpVertex[3],0);
    assert(pface[2]);
    pface[3] =  F_exist(ListOfpVertex[0],ListOfpVertex[2],ListOfpVertex[3],0);
    assert(pface[3]);
    
   
    pRegion pr = M_createR(mesh,4,pface,pg);
    assert(pr);
    de.receiveData (pr,from, &castbuf[sizeof(int)*(nbVertices+3)]);

    delete []  ListOfpVertex ;
  }
  // -------------------------------------------------------------------
  // -------------------------------------------------------------------
  void MigrateRegionsAndData( pMesh mesh, pMeshDataId tagDest, 
                              MDB_DataExchanger &de)
  {
    int nproc,myrank;
  
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    int *sendcounts = new int[nproc];
    for( int i=0; i<nproc; i++){
      sendcounts[i]= 0;
    }
    RIter rit = M_regionIter(mesh);
    pRegion pr;  
    while ((pr = RIter_next(rit))) {
      void *temp_ptr;
      int migre = EN_getDataPtr((pEntity) pr, tagDest, &temp_ptr);
      if(!migre) continue;
      const SetOfDestination_type *recup =
        reinterpret_cast< const SetOfDestination_type*> (temp_ptr);
      assert( (recup->size())>0 );
      int dest = *(recup->begin());
      assert(dest!=myrank);
      void *buf = RegionAndDataPackaging( pr, dest, de );    
      AP_send(buf);
      sendcounts[dest]++; 
    }    
    RIter_delete(rit);
  

    AP_check_sends(AP_NOFLAGS);
    AP_reduce_nsends(sendcounts);
  
    /*receive Regions and data*/
    int message=0;
    int count;
    while (!AP_recv_count(&count) || message<count) {
      void *msg;
      int  from;
      int  tag;
      int  size;
      int  recv;
      recv = AP_recv(MPI_ANY_SOURCE, de.tag(), AP_BLOCKING|AP_DROPOUT,
                     &msg, &size, &from, &tag);
      if(recv) {
        message++;
        RegionAndDataUnpackaging( mesh, msg, from, de);
        AP_free(msg);    
      }
    }
    AP_check_sends(AP_WAITALL);
    MPI_Barrier(MPI_COMM_WORLD);
    delete [] sendcounts;
  }
  // -------------------------------------------------------------------
  // ------------------------------------------------------------------- 
  // ------------------------------------------------------------------- 
  void MigrateEntitiesAndData(pMesh mesh, pMeshDataId tagDest, 
                              MDB_DataExchanger &de)
  {
    MigrateVerticesAndData( mesh, tagDest, de);
    MigrateEdgesAndData( mesh, tagDest, de);
    MigrateFacesAndData( mesh, tagDest, de);
    MigrateRegionsAndData( mesh, tagDest, de);
  }

}

#endif

