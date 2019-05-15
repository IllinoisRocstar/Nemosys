// -------------------------------------------------------------------
// MAdLib - Copyright (C) 2008-2009 Universite catholique de Louvain
//
// See the Copyright.txt and License.txt files for license information. 
// You should have received a copy of these files along with MAdLib. 
// If not, see <http://www.madlib.be/license/>
//
// Please report all bugs and problems to <contrib@madlib.be>
//
// Authors: Jean-Francois Remacle, Gaetan Compere, Koen Hillewaert
// -------------------------------------------------------------------

#include "MeshDataBase.h"
#include "MeshDataBaseInterface.h"
#include "MeshDataBaseGEntity2Physical.h"
#include "MAdMessage.h"
#include "MAdDefines.h"
#include "MathUtils.h"

#include <set>
using std::set;
using std::cout;

#ifdef PARALLEL
#include "mpi.h"
#endif

namespace MAd {

  MDB_Mesh::MDB_Mesh (int _MAXX)
    :  parametric(false), shrinked(false), maxId(_MAXX), idIncrement(1),
       nbPoints(0), nbEdges(0), nbTriangles(0), nbQuads(0), nbTets(0), 
       nbHexes(0), nbPrisms(0)
  {
    initializeIdData();
  }

  void MDB_Mesh::initializeIdData( )
  {
#ifdef PARALLEL
    int myRank = 0, nbProc = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &nbProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    int maxIdGlob = -1;
    MPI_Allreduce(&maxId,&maxIdGlob,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    maxId = maxIdGlob + myRank;
    idIncrement = nbProc;
    MPI_Barrier(MPI_COMM_WORLD);
#else
    idIncrement = 1;
#endif
  }

  void MDB_Mesh::deleteTheHighOrderPoint(std::set < MDB_Point *, PointLessThan >::iterator it){
    MDB_Point *pt;
    if (it != points.end())
      {
        pt = *it;
        points.erase(it);
        nbPoints--;
      }
    else throw;
  }

  void MDB_Point::getFaces(MDB_ListFace &t) const
  {
    t.clear();
    MDB_VectorE::const_iterator it = edges.begin();
    MDB_VectorE::const_iterator ite = edges.end();
    while(it != ite) {
      int NF = (*it)->numfaces();
      for(int i = 0; i < NF; ++i) {
        MDB_Face *tt = (*it)->faces(i);
        if(tt) {
          MDB_ListFace::iterator tit = t.begin();
          MDB_ListFace::iterator tite = t.end();
          int found = 0;
          while(tit != tite) {
            if(tt == *tit)
              found = 1;
            ++tit;
          }
          if(!found)
            t.push_back(tt);
        }
      }
      ++it;
    }
  }

  MDB_Point * MDB_Mesh::add_point(int num, double x, double y, double z, pGEntity pg)
  {
    int newId = num;
    if ( newId < 0 ) {
      maxId += idIncrement;
      newId = maxId;
    }

    MDB_Point *pp = new MDB_Point(newId, x, y, z);
    pp->g = pg;
    points.insert(pp);
    maxId = (maxId < newId) ? newId : maxId;
    nbPoints++;
    return pp;
  }

  MDB_PointParam * MDB_Mesh::add_pointParam(int num , double x, double y, double z, 
                                            double u, double v, pGEntity pg)
  {
    int newId = num;
    if ( newId < 0 ) {
      maxId += idIncrement;
      newId = maxId;
    }

    MDB_PointParam *pp = new MDB_PointParam(newId, x, y, z, u, v);
    pp->g = pg;
    points.insert(pp);
    maxId = (maxId < newId) ? newId : maxId;
    nbPoints++;
    parametric = true;
    return pp;
  }

  MDB_Point *MDB_Mesh::find_point(int p)
  {
    MDB_Point P(p);
    std::set < MDB_Point *, PointLessThan >::iterator it = points.find(&P);
    if(it != points.end()) return (MDB_Point *) (*it);
    return 0;

    //   MDB_PointParam VP(p);
    //   it = points.find(&VP);
    //   if(it != points.end()) return (MDB_Point *) (*it);
    //   else return 0;
  }

  MDB_Edge *MDB_Mesh::find_edge(int num1, int num2)
  {
    MDB_Point *p = find_point(num1);
    // MS
    if (!p)  // TODO: Use MAdLib's message class.
       std::cout << "Point does not exist." << std::endl;
    // MS End
    MDB_VectorE::iterator eit = p->edges.begin();
    while(eit != p->edges.end()) {
      if((*eit)->p1 == p && (*eit)->p2->iD == num2)
        return (*eit);
      if((*eit)->p2 == p && (*eit)->p1->iD == num2)
        return (*eit);
      ++eit;
    }
    return 0;
  }

  MDB_Edge *MDB_Mesh::find_edge(MDB_Point * p1, MDB_Point * p2) const
  {
    const int nbe1 = p1->edges.size();
    const int nbe2 = p2->edges.size();
    if (nbe1 < nbe2)
      {
        for (int i=0;i<nbe1;i++)
          {
            MDB_Edge *e = p1->edges [i];
            if (e->p1 == p1 && e->p2 == p2)return e;
            if (e->p2 == p1 && e->p1 == p2)return e;
          }
      }
    else
      {
        for (int i=0;i<nbe2;i++)
          {
            MDB_Edge *e = p2->edges [i];
            if (e->p1 == p1 && e->p2 == p2)return e;
            if (e->p2 == p1 && e->p1 == p2)return e;
          }
      }
    return 0;
  }


  MDB_Edge *MDB_Quad::find_edge(const MDB_Point * p1, const MDB_Point * p2) const
  {
    if((e1->p1->iD == p1->iD && e1->p2->iD == p2->iD) ||
       (e1->p1->iD == p2->iD && e1->p2->iD == p1->iD))
      return e1;
    if((e2->p1->iD == p1->iD && e2->p2->iD == p2->iD) ||
       (e2->p1->iD == p2->iD && e2->p2->iD == p1->iD))
      return e2;
    if((e3->p1->iD == p1->iD && e3->p2->iD == p2->iD)||
       (e3->p1->iD == p2->iD && e3->p2->iD == p1->iD))
      return e3;
    if((e4->p1->iD == p1->iD && e4->p2->iD == p2->iD)||
       (e4->p1->iD == p2->iD && e4->p2->iD == p1->iD))
      return e4;
    return 0; 
  }

  MDB_Edge *MDB_Triangle::find_edge(const MDB_Point * p1, const MDB_Point * p2) const
  {
    if((e1->p1->iD == p1->iD && e1->p2->iD == p2->iD) ||
       (e1->p1->iD == p2->iD && e1->p2->iD == p1->iD))
      return e1;
    if((e2->p1->iD == p1->iD && e2->p2->iD == p2->iD) ||
       (e2->p1->iD == p2->iD && e2->p2->iD == p1->iD))
      return e2;
    if((e3->p1->iD == p1->iD && e3->p2->iD == p2->iD)||
       (e3->p1->iD == p2->iD && e3->p2->iD == p1->iD))
      return e3;
    return 0;
  }

  MDB_Edge *MDB_Mesh::find_clone_edge(MDB_Edge * edge) const
  {
    double xyz[3];
    double eXYZ[2][3];
    MDB_Point * P = edge->p1;
    eXYZ[0][0] = P->X; eXYZ[0][1] = P->Y; eXYZ[0][2] = P->Z;
    P = edge->p2;
    eXYZ[1][0] = P->X; eXYZ[1][1] = P->Y; eXYZ[1][2] = P->Z;
    
    MDB_ListE::const_iterator eit = edges.begin();
    for(; eit != edges.end(); ++eit) {

      if ( (*eit) == edge ) continue; 

      P = (*eit)->p1;
      xyz[0] = P->X; xyz[1] = P->Y; xyz[2] = P->Z;
      if ( distanceSq(xyz,eXYZ[0]) < MAdTOL ) {
        P = (*eit)->p2;
        xyz[0] = P->X; xyz[1] = P->Y; xyz[2] = P->Z;
        if ( distanceSq(xyz,eXYZ[1]) < MAdTOL ) {
          return (*eit);
        }
      }

      P = (*eit)->p2;
      xyz[0] = P->X; xyz[1] = P->Y; xyz[2] = P->Z;
      if ( distanceSq(xyz,eXYZ[0]) < MAdTOL ) {
        P = (*eit)->p1;
        xyz[0] = P->X; xyz[1] = P->Y; xyz[2] = P->Z;
        if ( distanceSq(xyz,eXYZ[1]) < MAdTOL ) {
          return (*eit);
        }
      }
    }
    return NULL;
  }

  MDB_Triangle *MDB_Mesh::find_triangle(MDB_Edge * e1, MDB_Edge * e2,
                                        MDB_Edge * e3)
  {
    int i;
    for(i = 0; i < e1->numfaces(); i++) {
      MDB_Face *t = e1->faces(i);
      if(t->getNbEdges()!=3)
        continue;
      MDB_Edge *o1 = t->getEdge(0);
      MDB_Edge *o2 = t->getEdge(1);
      MDB_Edge *o3 = t->getEdge(2);
      if((o1 == e1 && o2 == e2 && o3 == e3) ||
         (o1 == e1 && o2 == e3 && o3 == e2) ||
         (o1 == e2 && o2 == e1 && o3 == e3) ||
         (o1 == e2 && o2 == e3 && o3 == e1) ||
         (o1 == e3 && o2 == e1 && o3 == e2) ||
         (o1 == e3 && o2 == e2 && o3 == e1))
        return (MDB_Triangle *)t;
    }

    return 0;
  }


  MDB_Quad *MDB_Mesh::find_quad(MDB_Edge * e1, MDB_Edge * e2,
                                MDB_Edge * e3, MDB_Edge * e4)
  {
    MDB_Edge *eds[4] = {e1,e2,e3,e4};
    std::sort(eds,eds+4);
    for(int i = 0; i < e1->numfaces(); i++) {
      MDB_Face *q = e1->faces(i);
      if(q->getNbEdges()!=4)
        continue;
      MDB_Edge *edsb[4] = {q->getEdge(0),q->getEdge(1),q->getEdge(2),q->getEdge(3)};
      std::sort(edsb,edsb+4);
      if (edsb[0] == eds[0] && edsb[1] == eds[1] && edsb[2] == eds[2] && edsb[3] == eds[3])
        return (MDB_Quad*)q;
    }
    return 0;
  }

  MDB_Edge *MDB_Mesh::add_edge(int p1, int p2, pGEntity pg)
  {
    MDB_Edge *efound = find_edge(p1, p2);
    if(efound)
      {
        if (efound->g == 0 || ( pg && GEN_type(pg) <  GEN_type(efound->g )))
          {
            efound->g = pg;
          }
        return efound;
      }

    MDB_Point *pp1 = find_point(p1);
    MDB_Point *pp2 = find_point(p2);
    if(!pp1 || !pp2)
      throw;
    return add_edge (pp1,pp2,pg);
  }

  /*MDB_Edge *MDB_Mesh::add_edge(int p1, int p2, int p3, pGEntity pg)
    {
    MDB_Edge *efound = find_edge(p1, p3);

    if(efound)
    return efound;
    MDB_Point *pp1 = find_point(p1);
    MDB_Point *pp2 = find_point(p2);
    MDB_Point *pp3 = find_point(p3);

    if(!pp1 || !pp2 || !pp3)
    throw;
    return add_edge (pp1,pp2,pp3,pg);
    }*/

  MDB_Edge *MDB_Mesh::add_edge(int numberOfPoints, pGEntity pg, int p1, ...)
  {
    va_list ap;
    // va_start(ap,numberOfPoints);
    va_start(ap,p1);

    int p2=0;int p3=0;int p4=0;int p5=0;
    if (numberOfPoints>=3){
      p2 = va_arg(ap,int);
      if (numberOfPoints>=4){
        p3 = va_arg(ap,int);
        if (numberOfPoints>=5){
          p4 = va_arg(ap,int);
          if (numberOfPoints>=6){
            p5 = va_arg(ap,int);
          }
        }
      }
    }
    int pend = va_arg(ap,int);
    va_end(ap);

    MDB_Edge *efound = find_edge(p1, pend);

    if(efound)
      return efound;
    MDB_Point *pp1 = find_point(p1);
    MDB_Point *pp2=NULL;
    MDB_Point *pp3=NULL;
    MDB_Point *pp4=NULL;
    MDB_Point *pp5=NULL;// = find_point(p2);
    MDB_Point *ppend = find_point(pend);
    switch(numberOfPoints){
    case 3:
      pp2 = find_point(p2);
      assert(pp1 && pp2 && ppend);
      return add_edge (3,pg,pp1,pp2,ppend);
    case 4:
      pp2 = find_point(p2);
      pp3 = find_point(p3);
      assert(pp1 && pp2 && pp3 && ppend);
      return add_edge (4,pg,pp1,pp2,pp3,ppend);
    case 5:
      pp2 = find_point(p2);
      pp3 = find_point(p3);
      pp4 = find_point(p4);
      assert(pp1 && pp2 && pp3 && pp4 && ppend);
      return add_edge (5,pg,pp1,pp2,pp3,pp4,ppend);
    case 6:
      pp2 = find_point(p2);
      pp3 = find_point(p3);
      pp4 = find_point(p4);
      pp5 = find_point(p5);
      assert(pp1 && pp2 && pp3 && pp4 && pp5 && ppend);
      return add_edge (6,pg,pp1,pp2,pp3,pp4,pp5,ppend);
    default:
      MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                  "add_edge only implemented up to order 5 (this edge has %d nodes)",
                                  numberOfPoints);
    }
    throw;
  }

  MDB_Edge *MDB_Mesh::add_edge(MDB_Point *pp1 , MDB_Point *pp2, pGEntity pg)
  {
    MDB_Edge *e = new MDB_Edge(pp1, pp2);
    e->g = pg;
    edges.push_back(e);
    nbEdges++;
    return e;
  }

  MDB_Edge *MDB_Mesh::add_edge(MDB_Point* pp1,MDB_Point *pp2,pGEntity pg,int order,MDB_Point** pp)
  {
    //   int numberOfPoints = order + 1;

    MDB_Edge *e = NULL;
    switch(order){
    case 1:
      e = new MDB_Edge(pp1,pp2);
      break;
    case 2:
      e = new MDB_EdgeP2(pp1, pp[0], pp2);
      break;
    case 3:
      e = new MDB_EdgeP3(pp1, pp[0], pp[1], pp2);
      break;
    case 4:
      e = new MDB_EdgeP4(pp1, pp[0], pp[1], pp[2], pp2);
      break;
    case 5:
      e = new MDB_EdgeP5(pp1, pp[0], pp[1], pp[2], pp[3] , pp2);
      break;
    default:
      MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                  "add_edge only implemented up to order 5");
    }
    e->g = pg;

    for (int i=0;i<order-1;i++) {
      if (pp[i]->g == NULL) pp[i]->g = pg;
      std::set < MDB_Point *, PointLessThan >::iterator it = points.find(pp[i]);
      deleteTheHighOrderPoint(it);
    }

    edges .push_back(e);
    nbEdges++;
    return e;
  }



  MDB_Edge *MDB_Mesh::add_edge(int numberOfPoints, pGEntity pg, MDB_Point *firstPoint, ...)
  {
    va_list ap;
    // va_start(ap,numberOfPoints);
    va_start(ap,firstPoint);

    MDB_Point *pp2 = NULL;
    MDB_Point *pp3 = NULL;
    MDB_Point *pp4 = NULL;
    MDB_Point *pp5 = NULL;

    if (numberOfPoints>=3){
      pp2 = va_arg(ap,MDB_Point*);
      if (numberOfPoints>=4){
        pp3 = va_arg(ap,MDB_Point*);
        if (numberOfPoints>=5){
          pp4 = va_arg(ap,MDB_Point*);
          if (numberOfPoints>=6){
            pp5 = va_arg(ap,MDB_Point*);
          }
        }
      }
    }
    MDB_Point *ppend = va_arg(ap,MDB_Point*);
    va_end(ap);

    MDB_Edge *e = NULL;
    switch(numberOfPoints){
    case 2:
      break;
    case 3:
      e = new MDB_EdgeP2(firstPoint, pp2, ppend);
      break;
    case 4:
      e = new MDB_EdgeP3(firstPoint, pp2, pp3, ppend);
      break;
    case 5:
      e = new MDB_EdgeP4(firstPoint, pp2, pp3, pp4, ppend);
      break;
    case 6:
      e = new MDB_EdgeP5(firstPoint, pp2, pp3, pp4, pp5, ppend);
      break;
    default:
      MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                  "add_edge only implemented up to order 5");
    }
    e->g = pg;
    edges.push_back(e);
    nbEdges++;

    if (numberOfPoints==2) return e;

    std::set < MDB_Point *, PointLessThan >::iterator it1 = points.find(pp2);
    deleteTheHighOrderPoint(it1);

    std::set < MDB_Point *, PointLessThan >::iterator it2;
    std::set < MDB_Point *, PointLessThan >::iterator it3;
    std::set < MDB_Point *, PointLessThan >::iterator it4;
    switch(numberOfPoints){
    case 4:
      it2 = points.find(pp3);
      deleteTheHighOrderPoint(it2);
      break;
    case 5:
      it2 = points.find(pp3);
      it3 = points.find(pp4);
      deleteTheHighOrderPoint(it2);
      deleteTheHighOrderPoint(it3);
      break;
    case 6:
      it2 = points.find(pp3);
      it3 = points.find(pp4);
      it4 = points.find(pp5);
      deleteTheHighOrderPoint(it2);
      deleteTheHighOrderPoint(it3);
      deleteTheHighOrderPoint(it4);
      break;
    }
    return e;
  }

  /*MDB_Edge *MDB_Mesh::add_edge(MDB_Point *pp1 , MDB_Point *pp2, MDB_Point *pp3, pGEntity pg)
    {
    MDB_Edge *e = new MDB_EdgeP2(pp1, pp2,pp3);
    e->g = pg;
    edges.push_back(e);
    nbEdges++;

    std::set < MDB_Point *, PointLessThan >::iterator it = points.find(pp2);
    MDB_Point *pt;
    if (it != points.end())
    {
    pt = *it;
    points.erase(it);
    nbPoints--;
    }
    else throw;
    return e;
    }*/

  MDB_Triangle *MDB_Mesh::add_triangle(int p1, int p2, int p3, pGEntity pg)
  {
    MDB_Edge *e1 = add_edge(p1, p2);
    MDB_Edge *e2 = add_edge(p2, p3);
    MDB_Edge *e3 = add_edge(p3, p1);


    MDB_Triangle *efound = find_triangle(e1, e2, e3);
    if(efound)
      {
        if (efound->g == 0 || ( pg && GEN_type(pg) <  GEN_type(efound->g )))
          {
            efound->g = pg;
          }
        return efound;
      }

    return add_triangle(e1, e2, e3,pg);
  }

  MDB_Triangle *MDB_Mesh::add_triangle(MDB_Point* p1, MDB_Point* p2, 
                                       MDB_Point* p3, pGEntity pg)
  {
    MDB_Edge *e1 = find_edge(p1,p2);
    if (!e1) e1 = add_edge(p1, p2, pg);
    MDB_Edge *e2 = find_edge(p2,p3);
    if (!e2) e2 = add_edge(p2, p3, pg);
    MDB_Edge *e3 = find_edge(p3,p1);
    if (!e3) e3 = add_edge(p3, p1, pg);

    return add_triangle(e1, e2, e3, pg);
  }

  /*MDB_Triangle *MDB_Mesh::add_triangle(int numberOfPoints, pGEntity pg, int p1, int p2, ...)
    {
    va_list ap;
    // va_start(ap,numberOfPoints);
    va_start(ap,p2);
    int p3=0;int p4=0;int p5=0;
    int p6=0;int p7=0;int p8=0;int p9=0;
    int p10=0;int p11=0;int p12=0;int p13=0;
    int p14=0;int p15=0;
    if (numberOfPoints>=6){
    p3 = va_arg(ap,int);
    p4 = va_arg(ap,int);
    p5 = va_arg(ap,int);
    p6 = va_arg(ap,int);
    if (numberOfPoints>=9){
    p7 = va_arg(ap,int);
    p8 = va_arg(ap,int);
    p9 = va_arg(ap,int);
    if (numberOfPoints>=12){
    p10 = va_arg(ap,int);
    p11 = va_arg(ap,int);
    p12 = va_arg(ap,int);
    if (numberOfPoints>=15){
    p13 = va_arg(ap,int);
    p14 = va_arg(ap,int);
    p15 = va_arg(ap,int);
    }
    }
    }
    }
    va_end(ap);

    MDB_Edge *e1,*e2,*e3;
    switch(numberOfPoints){
    case 3:
    e1 = add_edge(2,0,p1, p2);
    e2 = add_edge(2,0,p2, p3);
    e3 = add_edge(2,0,p3, p1);
    break;
    case 6:
    e1 = add_edge(3,0,p1, p2, p3);
    e2 = add_edge(3,0,p3, p4, p5);
    e3 = add_edge(3,0,p5, p6, p1);
    break;
    case 9:
    e1 = add_edge(4,0,p1, p2, p3, p4);
    e2 = add_edge(4,0,p4, p5, p6, p7);
    e3 = add_edge(4,0,p7, p8, p9, p1);
    break;
    case 12:
    e1 = add_edge(5,0,p1, p2, p3, p4, p5);
    e2 = add_edge(5,0,p5, p6, p7, p8, p9);
    e3 = add_edge(5,0,p9, p10,p11,p12,p1);
    break;
    case 15:
    e1 = add_edge(6,0,p1, p2, p3, p4, p5, p6);
    e2 = add_edge(6,0,p6, p7, p8, p9, p10,p11);
    e3 = add_edge(6,0,p11,p12,p13,p14,p15,p1);
    break;
    }

    MDB_Triangle *efound = find_triangle(e1, e2, e3);
    if(efound)
    {
    if (efound->g == 0 || ( pg && GEN_type(pg) <  GEN_type(efound->g )))
    {
    efound->g = pg;
    }
    return efound;
    }
    return add_triangle(e1, e2, e3,pg);
    }*/

  // high order triangle ...
  // template is the following v0 ho[v0->v1] v1 ho[v1->v2] v2 ho[v2->v0]
  MDB_Triangle *MDB_Mesh::add_triangle(int order, bool complete, pGEntity pg, int p1, int p2, ...)
  {
    // int numberOfPoints = complete ? (order+2)*(order+1)/2 : 3*order;
    va_list ap;
    // va_start(ap,numberOfPoints);
    va_start(ap,p2);
    int p3=0;int p4=0;int p5=0;
    int p6=0;int p7=0;int p8=0;int p9=0;
    int p10=0;int p11=0;int p12=0;int p13=0;
    int p14=0;int p15=0;
    if (order>=2){
      p3 = va_arg(ap,int);
      p4 = va_arg(ap,int);
      p5 = va_arg(ap,int);
      p6 = va_arg(ap,int);
      if (order>=3){
        p7 = va_arg(ap,int);
        p8 = va_arg(ap,int);
        p9 = va_arg(ap,int);
        if (order>=4){
          p10 = va_arg(ap,int);
          p11 = va_arg(ap,int);
          p12 = va_arg(ap,int);
          if (order>=5){
            p13 = va_arg(ap,int);
            p14 = va_arg(ap,int);
            p15 = va_arg(ap,int);
          }
        }
      }
    }
    if (complete && order>2){
      switch(order){
      case 3:
        p10 = va_arg(ap,int);
        break;
      case 4:
        p13 = va_arg(ap,int);
        p14 = va_arg(ap,int);
        p15 = va_arg(ap,int);
        break;
      default:
        MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                    "Complete triangles not defined for order %d",
                                    order);
      }
    }
    va_end(ap);

    MDB_Edge *e1 = NULL;
    MDB_Edge *e2 = NULL;
    MDB_Edge *e3 = NULL;

    switch(order){
    case 1:
      e1 = add_edge(2,0,p1, p2);
      e2 = add_edge(2,0,p2, p3);
      e3 = add_edge(2,0,p3, p1);
      break;
    case 2:
      e1 = add_edge(3,0,p1, p2, p3);
      e2 = add_edge(3,0,p3, p4, p5);
      e3 = add_edge(3,0,p5, p6, p1);
      break;
    case 3:
      e1 = add_edge(4,0,p1, p2, p3, p4);
      e2 = add_edge(4,0,p4, p5, p6, p7);
      e3 = add_edge(4,0,p7, p8, p9, p1);
      break;
    case 4:
      e1 = add_edge(5,0,p1, p2, p3, p4, p5);
      e2 = add_edge(5,0,p5, p6, p7, p8, p9);
      e3 = add_edge(5,0,p9, p10,p11,p12,p1);
      break;
    case 5:
      e1 = add_edge(6,0,p1, p2, p3, p4, p5, p6);
      e2 = add_edge(6,0,p6, p7, p8, p9, p10,p11);
      e3 = add_edge(6,0,p11,p12,p13,p14,p15,p1);
      break;
    }
    MDB_Point** innerPoints = NULL;
    if (complete && order>2){
      int nbFacePoints = (order-1) * (order-2)/2;
      innerPoints = new MDB_Point*[nbFacePoints];
      switch (order){
      case 3:
        innerPoints[0] = find_point(p10);
        break;
      case 4:
        innerPoints[0] = find_point(p13);
        innerPoints[1] = find_point(p14);
        innerPoints[2] = find_point(p15);
        break;
      default:
        MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                    "Complete triangles not defined for order %d",
                                    order);
      }
    }

    MDB_Triangle *efound = find_triangle(e1, e2, e3);
    if(efound)
      {
        if (efound->g == 0 || ( pg && GEN_type(pg) <  GEN_type(efound->g )))
          {
            efound->g = pg;
          }
        return efound;
      }
    return add_triangle(e1,e2,e3,pg,order,!complete,innerPoints);
    //return add_triangle(e1, e2, e3,pg);
  }

  // high order quad ...
  // template is the following v0 ho[v0->v1] v1 ho[v1->v2] v2 ho[v2->v3] v3 ho[v3->v0]
  MDB_Quad *MDB_Mesh::add_quad(int order, bool complete, pGEntity pg, int p1, int p2, int p3, ...)
  {
    va_list ap;
    va_start(ap,p3);
    int p4(0),p5(0),p6(0),p7(0),p8(0),p9(0);
    if(order>=1){
      p4 = va_arg(ap,int);
      if (order>=2){
        p5 = va_arg(ap,int);
        p6 = va_arg(ap,int);
        p7 = va_arg(ap,int);
        p8 = va_arg(ap,int);
      }
    }
    if (complete && order>=2){
      switch(order){
      case 2:
        p9 = va_arg(ap,int);
        break;
      default:
        MAdMsgSgl::instance().error(__LINE__,__FILE__, "Complete quad not defined for order %d", order);
      }
    }
    va_end(ap);

    MDB_Edge *e1 = NULL;
    MDB_Edge *e2 = NULL;
    MDB_Edge *e3 = NULL;
    MDB_Edge *e4 = NULL;

    switch(order){
    case 1:
      e1 = add_edge(2,0,p1, p2);
      e2 = add_edge(2,0,p2, p3);
      e3 = add_edge(2,0,p3, p4);
      e4 = add_edge(2,0,p4, p1);
      break;
    case 2:
      e1 = add_edge(3,0,p1, p2, p3);
      e2 = add_edge(3,0,p3, p4, p5);
      e3 = add_edge(3,0,p5, p6, p7);
      e4 = add_edge(3,0,p7, p8, p1);
      break;
      default :
        MAdMsgSgl::instance().error(__LINE__,__FILE__, "Quad not defined for order %d", order);
    }
    MDB_Point** innerPoints = NULL;
    if (complete && order>=2){
      int nbFacePoints = (order-1)*(order-1);
      innerPoints = new MDB_Point*[nbFacePoints];
      switch (order){
      case 2:
        innerPoints[0] = find_point(p9);
        break;
      default:
        MAdMsgSgl::instance().error(__LINE__,__FILE__, "Complete triangles not defined for order %d", order);
      }
    }

    MDB_Quad *efound = find_quad(e1, e2, e3,e4);
    if(efound)
      {
        if (efound->g == 0 || ( pg && GEN_type(pg) <  GEN_type(efound->g )))
          {
            efound->g = pg;
          }
        return efound;
      }
    return add_quad(e1,e2,e3,e4,pg,order,!complete,innerPoints);
  }

  // template for higher order triangle, following gmsh template
  // v0 v1 v2
  // higher order points
  // edges  loop following edge direction E0(V0-V1) E1(V1->V2) E2(V2->V0)
  // face   loop (eta(ksi)), ksi=e0, eta=-e2

  /*MDB_Triangle *MDB_Mesh::add_triangle(pGEntity pg,int order,bool serendip,int* p) {


  int nbEdgePoints = (order-1);
  int nbFacePoints = (order+1) * (order+2)/2 - (serendip ? std::max(0,(order-2)*(order-1)/2) : 0);

  MDB_Point** point = new MDB_Point*[nbFacePoints];

  for (int i=0;i<nbFacePoints;i++) point[i] = find_point(p[i]);

  MDB_Point** pp = point + 3;

  // create edges

  MDB_Edge* e1 = add_edge(point[0],point[1],pg,order,pp);
  pp += nbEdgePoints;

  MDB_Edge* e2 = add_edge(point[1],point[2],pg,order,pp);
  pp += nbEdgePoints;

  MDB_Edge* e3 = add_edge(point[2],point[0],pg,order,pp);
  pp += nbEdgePoints;

  //

  return add_triangle(e1,e2,e3,pg,order,serendip,pp);
  }*/

  MDB_Triangle *MDB_Mesh::add_triangle(MDB_Edge * e1, MDB_Edge * e2,
                                       MDB_Edge * e3, pGEntity pg)
  {
    MDB_Triangle *t = new MDB_Triangle(e1, e2, e3);
    t->g = pg;
    triangles.push_back(t);
    nbTriangles++;
    return t;
  }

  MDB_Triangle *MDB_Mesh::add_triangle(MDB_Edge *e1,MDB_Edge* e2,
                                       MDB_Edge *e3,pGEntity pg,
                                       int order,bool serendip,
                                       MDB_Point** point) {




    MDB_Triangle *t = NULL;
    if (serendip || order <= 2) t = new MDB_Triangle(e1,e2,e3);
    else {
      switch (order) {

      case 3:
        t = new MDB_CompleteTriangle<3>(e1,e2,e3,point);
        break;
      case 4:
        t = new MDB_CompleteTriangle<4>(e1,e2,e3,point);
        break;
        /*    case 5:
              t = new MDB_CompleteTriangle<5>(e1,e2,e3,point);
              break;*/
      default:
        MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                    "Complete triangle of order %d",
                                    order);
      }

      int nbPoints = (order-2)*(order-1)/2;

      for (int i=0;i<nbPoints;i++) {
        std::set<MDB_Point*,PointLessThan>::iterator it = points.find(point[i]);
        deleteTheHighOrderPoint(it);
      }
    }

    if (t) {
      t->g = pg;
      triangles.push_back(t);
      nbTriangles++;
    }
    return t;
  }



  MDB_Quad *MDB_Mesh::add_quad(int p1, int p2, int p3, int p4, pGEntity pg)
  {
    MDB_Edge *e1 = add_edge(p1, p2);
    MDB_Edge *e2 = add_edge(p2, p3);
    MDB_Edge *e3 = add_edge(p3, p4);
    MDB_Edge *e4 = add_edge(p4, p1);

    MDB_Quad *efound = find_quad(e1, e2, e3, e4);
    if(efound)
      {
        if (efound->g == 0 || ( pg && GEN_type(pg) <  GEN_type(efound->g )))
          {
            efound->g = pg;
          }
        return efound;
      }

    return add_quad(e1, e2, e3, e4,pg);
  }

  MDB_Quad *MDB_Mesh::add_quad(MDB_Point* p1, MDB_Point* p2, 
                               MDB_Point* p3, MDB_Point* p4, pGEntity pg)
  {
    MDB_Edge *e1 = find_edge(p1,p2);
    if (!e1) e1 = add_edge(p1, p2, pg);
    MDB_Edge *e2 = find_edge(p2,p3);
    if (!e2) e2 = add_edge(p2, p3, pg);
    MDB_Edge *e3 = find_edge(p3,p4);
    if (!e3) e3 = add_edge(p3, p4, pg);
    MDB_Edge *e4 = find_edge(p4,p1);
    if (!e4) e4 = add_edge(p4, p1, pg);

    return add_quad(e1, e2, e3, e4, pg);
  }

  MDB_Quad *MDB_Mesh::add_quad(MDB_Edge * e1, MDB_Edge * e2,
                               MDB_Edge * e3, MDB_Edge * e4, pGEntity pg)
  {
    MDB_Quad *q = new MDB_Quad(e1, e2, e3, e4);
    q->g = pg;
    quads.push_back(q);
    nbQuads++;
    return q;
  }


  static double VOLUME ( MDB_Tet *t )
  {
    pVertex v1 = R_vertex (t,0);
    pVertex v2 = R_vertex (t,1);
    pVertex v3 = R_vertex (t,2);
    pVertex v4 = R_vertex (t,3);

    double v12[3] = { v2->X - v1->X,  v2->Y - v1->Y,  v2->Z - v1->Z};
    double v13[3] = { v3->X - v1->X,  v3->Y - v1->Y,  v3->Z - v1->Z};
    double v14[3] = { v4->X - v1->X,  v4->Y - v1->Y,  v4->Z - v1->Z};

    double pv[3] = { v12[1]*v13[2]-v12[2]*v13[1],
                     v12[2]*v13[0]-v12[0]*v13[2],
                     v12[0]*v13[1]-v12[1]*v13[0]};
    double v =
      pv[0] * v14[0]
      + pv[1] * v14[1]
      + pv[2] * v14[2];

    return v/6;
  }

  static double VOLUME ( MDB_Hex *h )
  {
    pVertex v1 = R_vertex (h,0);
    pVertex v2 = R_vertex (h,1);
    pVertex v3 = R_vertex (h,2);
    pVertex v4 = R_vertex (h,3);
    pVertex v5 = R_vertex (h,4);
    pVertex v6 = R_vertex (h,5);
    pVertex v7 = R_vertex (h,6);
    pVertex v8 = R_vertex (h,7);

    double v21[3] = { v1->X - v2->X,  v1->Y - v2->Y,  v1->Z - v2->Z};
    double v23[3] = { v3->X - v2->X,  v3->Y - v2->Y,  v3->Z - v2->Z};
    double v24[3] = { v4->X - v2->X,  v4->Y - v2->Y,  v4->Z - v2->Z};
    double v25[3] = { v5->X - v2->X,  v5->Y - v2->Y,  v5->Z - v2->Z};
    double v26[3] = { v6->X - v2->X,  v6->Y - v2->Y,  v6->Z - v2->Z};
    double v27[3] = { v7->X - v2->X,  v7->Y - v2->Y,  v7->Z - v2->Z};
    double v28[3] = { v8->X - v2->X,  v8->Y - v2->Y,  v8->Z - v2->Z};

    double pv248[3] = { v24[1]*v28[2]-v24[2]*v28[1], v24[2]*v28[0]-v24[0]*v28[2], v24[0]*v28[1]-v24[1]*v28[0]};
    double pv258[3] = { v25[1]*v28[2]-v25[2]*v28[1], v25[2]*v28[0]-v25[0]*v28[2], v25[0]*v28[1]-v25[1]*v28[0]};
    double pv278[3] = { v27[1]*v28[2]-v27[2]*v28[1], v27[2]*v28[0]-v27[0]*v28[2], v27[0]*v28[1]-v27[1]*v28[0]};

    double vo1 = pv248[0] * v21[0] + pv248[1] * v21[1] + pv248[2] * v21[2];
    double vo2 = pv258[0] * v21[0] + pv258[1] * v21[1] + pv258[2] * v21[2];
    double vo3 = pv258[0] * v26[0] + pv258[1] * v26[1] + pv258[2] * v26[2];
    double vo4 = pv248[0] * v23[0] + pv248[1] * v23[1] + pv248[2] * v23[2];
    double vo5 = pv278[0] * v23[0] + pv278[1] * v23[1] + pv278[2] * v23[2];
    double vo6 = pv278[0] * v26[0] + pv278[1] * v26[1] + pv278[2] * v26[2];

    //NP The volume of hex was 0, I think there is an error in the tet volume formula.
    //   I fix the problem by changing the following line.
    //   return (vo1+vo2+vo3+vo4+vo5+vo6)/6;
    return (-vo1+vo2-vo3+vo4-vo5+vo6)/6;
  }
  MDB_Quad *MDB_Mesh::add_quad(MDB_Edge *e1,MDB_Edge* e2,
                                       MDB_Edge *e3,MDB_Edge *e4,pGEntity pg,
                                       int order,bool serendip,
                                       MDB_Point** point) {




    MDB_Quad *q = NULL;
    if (serendip || order < 2) q = new MDB_Quad(e1,e2,e3,e4);
    else {
      switch (order) {
        case 2:
          q = new MDB_CompleteQuad<2>(e1,e2,e3,e4,point);
          break;
        default:
          MAdMsgSgl::instance().error(__LINE__,__FILE__, "Complete quad of order %d", order);
      }
      int nbPoints = (order-1)*(order-1);
      for (int i=0;i<nbPoints;i++) {
        std::set<MDB_Point*,PointLessThan>::iterator it = points.find(point[i]);
        deleteTheHighOrderPoint(it);
      }
    }

    if (q) {
      q->g = pg;
      quads.push_back(q);
      nbQuads++;
    }
    return q;
  }


  // -----------------------------------------------------------------------------
  // template tet
  // -----------------------------------------------------------------------------
  //
  // gmsh template
  // V0 V1 V2 V3
  // higher order points
  //
  // edges:
  // ------
  //
  //   E0 : V0 V1
  //   E1 : V1 V2
  //   E2 : V2 V0
  //   E3 : V0 V3
  //   E5 : V2 V3
  //   E4 : V1 V3
  //
  // faces:
  // ------
  //
  // loop (eta(ksi)), ksi e0 (v0-v1) , eta -e2 (v0-v2)
  // orientation such that rotation axis points inward
  // face starts at vertex with corresponding index
  //
  //      |  e0  e1  e2 | v0 v1 v2 |
  //   -----------------------------
  //   F0 | +E0 +E1 +E2 | V0 V1 V2 |
  //   F1 | -E0 +E3 -E4 | V1 V0 V3 |
  //   F2 | +E5 -E3 -E2 | V2 V3 V0 |
  //   F3 | -E5 -E1 +E4 | V3 V2 V1 |
  //
  //
  // volume:
  // -------
  //
  // loop (zeta(eta(ksi))), ksi (V0-V1), eta (V0-V2), zeta (V0-V3)
  //
  // -----------------------------------------------------------------------------



  MDB_Tet *MDB_Mesh::add_tet(pGEntity pg,int order,bool serendip,int* pointID) {

    int nbEdgePoints = order > 1 ? order - 1 : 0;
    int nbFacePoints = order > 2 ? (order-2)*(order-1)/2 : 0;
    int nbVolPoints = ((serendip || order < 4) ? 0 : (order-3)*(order-2)*(order-1)/6);

    int nb = 4 + 6 * nbEdgePoints +  4 * nbFacePoints + nbVolPoints;

    MDB_Point** point = new MDB_Point*[nb];
    for (int i=0;i<nb;i++) point[i] = find_point(pointID[i]);

    MDB_Point** pp = point + 4;

    // create edges
    // to be tested : ordering of higher order nodes ...

    MDB_Edge *e1 = find_edge(point[0],point[1]);
    if (!e1) e1 = add_edge(point[0], point[1],pg,order,pp);
    pp += nbEdgePoints;

    MDB_Edge *e2 = find_edge(point[1],point[2]);
    if (!e2) e2 = add_edge(point[1], point[2],pg,order,pp);
    pp += nbEdgePoints;

    MDB_Edge *e3 = find_edge(point[2],point[0]);
    if (!e3) e3 = add_edge(point[2], point[0],pg,order,pp);
    pp += nbEdgePoints;

    // KH: edges 4,5 and 6 have opposite orientation in gmsh

    MDB_Edge *e4 = find_edge(point[0],point[3]);
    if (!e4) e4 = add_edge(point[3], point[0],pg,order,pp);
    pp += nbEdgePoints;

    // KH: edges 5 and 6 switched wrt gmsh, opposite direction

    MDB_Edge *e6 = find_edge(point[2],point[3]);
    if (!e6) e6 = add_edge(point[3], point[2],pg,order,pp);
    pp += nbEdgePoints;

    MDB_Edge *e5 = find_edge(point[1],point[3]);
    if (!e5) e5 = add_edge(point[3], point[1],pg,order,pp);
    pp += nbEdgePoints;

    // create faces - complete

    // face 1 is face 1 in Gmsh, the latter is oriented as v0-v2-v1

    MDB_Triangle *t1 = find_triangle(e1, e2, e3);
    // if(!t1) t1 = add_triangle(e1, e2, e3,pg,order,false,pp);
    if(!t1) t1 = add_triangle(e3, e2, e1,pg,order,false,pp);
    pp += nbFacePoints;

    // face 2 is the same as in Gmsh

    MDB_Triangle *t2 = find_triangle(e1, e4, e5);
    if(!t2) t2 = add_triangle(e1, e5, e4,pg,order,false,pp);
    pp += nbFacePoints;

    // face 4 is face 3 in gmsh, the latter is oriented as v0-v3-v2

    MDB_Triangle *t4 = find_triangle(e3, e4, e6);
    // if(!t4)t4 = add_triangle (e3, e4, e6,pg,order,false,pp);
    if(!t4) t4 = add_triangle (e4, e6, e3,pg,order,false,pp);
    pp += nbFacePoints;

    // face 3 in mdb is  face 4 in gmsh, the latter is oriented as v3-v1-v2

    MDB_Triangle *t3 = find_triangle(e2, e6, e5);
    // if(!t3) t3 = add_triangle(e2, e6, e5,pg,order,false,pp);
    if(!t3) t3 = add_triangle(e5, e2, e6,pg,order,false,pp);
    pp += nbFacePoints;

    // create tets

    MDB_Tet* t = serendip ? add_tet(t1,t2,t3,t4,pg) : add_tet(t1,t2,t3,t4,pg,order,serendip,pp);

    // clean up

    delete [] point;
    return t;
  }




  MDB_Tet *MDB_Mesh::add_tet(MDB_Triangle *t1,MDB_Triangle *t2,
                             MDB_Triangle *t3,MDB_Triangle *t4,
                             pGEntity pg,int order,bool serendip,MDB_Point** pp) {


    MDB_Tet* t = NULL;

    if (serendip || order <= 2)  {
      t = new MDB_Tet(t1,t2,t3,t4);
    }
    else {
      switch (order){

      case 3:
        t = new MDB_CompleteTet<3> (t1,t2,t3,t4,pp);
        break;
      case 4:
        t = new MDB_CompleteTet<4> (t1,t2,t3,t4,pp);
        break;
      case 5:
        t = new MDB_CompleteTet<5> (t1,t2,t3,t4,pp);
        break;
      default:
        MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                    "Complete tetrahedron of order %d not yet implemented",
                                    order);
        break;
      }

      int nbPoints = std::max(0,(order-3)*(order-2)*(order-1)/6);

      for (int i=0;i<nbPoints;i++) {
        std::set<MDB_Point*,PointLessThan>::iterator it = points.find(pp[i]);
        deleteTheHighOrderPoint(it);
      }
    }

    if (t) {
      t->g = pg;
      tets.push_back(t);
      nbTets++;
    }
    return t;
  }


  MDB_Tet *MDB_Mesh::add_tet(MDB_Triangle *t1,
                             MDB_Triangle *t2,
                             MDB_Triangle *t3,
                             MDB_Triangle *t4,
                             pGEntity pg)
  {
    MDB_Tet *t = new MDB_Tet(t1, t2, t3, t4);

    double v = VOLUME (t);
    if (v < 0)
      {
        t->f1 = t2;
        t->f2 = t1;
        //      v = VOLUME (t);
        //      if (v < 0)
        //	throw;
      }

    t->g = pg;
    tets.push_back(t);
    nbTets++;
    return t;
  }

  MDB_Tet *MDB_Mesh::add_tet(int p1, int p2, int p3, int p4, pGEntity pg)
  {
    MDB_Edge *e1 = add_edge(p1, p2);
    MDB_Edge *e2 = add_edge(p2, p3);
    MDB_Edge *e3 = add_edge(p3, p1);
    MDB_Edge *e4 = add_edge(p1, p4);
    MDB_Edge *e5 = add_edge(p2, p4);
    MDB_Edge *e6 = add_edge(p3, p4);

    MDB_Triangle *t1 = find_triangle(e1, e2, e3);
    if(!t1)
      t1 = add_triangle(e1, e2, e3);
    MDB_Triangle *t2 = find_triangle(e1, e4, e5);
    if(!t2)
      t2 = add_triangle(e1, e4, e5);
    MDB_Triangle *t3 = find_triangle(e2, e6, e5);
    if(!t3)
      t3 = add_triangle(e2, e6, e5);
    MDB_Triangle *t4 = find_triangle(e3, e4, e6);
    if(!t4)
      t4 = add_triangle(e3, e4, e6);
    return add_tet (t1, t2, t3, t4,pg);
  }

  MDB_Tet *MDB_Mesh::add_tet(MDB_Point *p1, MDB_Point *p2,MDB_Point *p3, MDB_Point *p4, pGEntity pg)
  {
    MDB_Edge *e1 = find_edge(p1,p2);
    if (!e1) e1 = add_edge(p1, p2, pg);
    MDB_Edge *e2 = find_edge(p2,p3);
    if (!e2) e2 = add_edge(p2, p3, pg);
    MDB_Edge *e3 = find_edge(p3,p1);
    if (!e3) e3 = add_edge(p3, p1, pg);
    MDB_Edge *e4 = find_edge(p1,p4);
    if (!e4) e4 = add_edge(p1, p4, pg);
    MDB_Edge *e5 = find_edge(p2,p4);
    if (!e5) e5 = add_edge(p2, p4, pg);
    MDB_Edge *e6 = find_edge(p3,p4);
    if (!e6) e6 = add_edge(p3, p4, pg);

    MDB_Triangle *t1 = find_triangle(e1, e2, e3);
    if(!t1)t1 = add_triangle(e1, e2, e3, pg);
    MDB_Triangle *t2 = find_triangle(e1, e4, e5);
    if(!t2)t2 = add_triangle(e1, e4, e5, pg);
    MDB_Triangle *t3 = find_triangle(e2, e6, e5);
    if(!t3)t3 = add_triangle(e2, e6, e5, pg);
    MDB_Triangle *t4 = find_triangle(e3, e4, e6);
    if(!t4)t4 = add_triangle(e3, e4, e6, pg);
    return add_tet (t1, t2, t3, t4, pg);
  }

  void MDB_Mesh::del_tet(MDB_Tet * t)
  {
    t->f1->del(t);
    t->f2->del(t);
    t->f3->del(t);
    t->f4->del(t);
    nbTets--;
    t->deleted = true;
  }


  void MDB_Mesh::del_face(MDB_Face *f)
  {
    int nbE = f->getNbEdges();
    if ( nbE == 3 ) {
      del_triangle ( dynamic_cast<MDB_Triangle *>(f) ); return;
    }
    else if ( nbE == 4 ) {
      del_quad ( dynamic_cast<MDB_Quad *>(f) ); return;
    }
    throw;
  }

  MDB_Hex *MDB_Mesh::add_hex(MDB_Quad *q1, MDB_Quad *q2, MDB_Quad *q3, MDB_Quad *q4, MDB_Quad *q5, MDB_Quad *q6, pGEntity pg)
  {
    MDB_Hex *h= new MDB_Hex(q1, q2, q3, q4, q5, q6);

    double v = VOLUME (h);
    if (v < 0)
      {
        printf("Hexahedron has a negative volume.");
        throw;
      }

    h->g = pg;
    hexes.push_back(h);
    nbHexes++;

    return h;
  }

  MDB_Hex *MDB_Mesh::add_hex(int p1, int p2, int p3, int p4,
                             int p5, int p6, int p7, int p8, pGEntity pg)
  {
    MDB_Edge *e1 = add_edge(p1, p2);
    MDB_Edge *e2 = add_edge(p2, p3);
    MDB_Edge *e3 = add_edge(p3, p4);
    MDB_Edge *e4 = add_edge(p4, p1);
    MDB_Edge *e5 = add_edge(p1, p5);
    MDB_Edge *e6 = add_edge(p2, p6);
    MDB_Edge *e7= add_edge(p3, p7);
    MDB_Edge *e8 = add_edge(p4, p8);
    MDB_Edge *e9 = add_edge(p5, p6);
    MDB_Edge *e10 = add_edge(p6, p7);
    MDB_Edge *e11 = add_edge(p7, p8);
    MDB_Edge *e12 = add_edge(p8, p5);

    MDB_Quad *q1 = find_quad(e4, e3, e2, e1);
    if(!q1) q1 = add_quad(e4, e3, e2, e1);
    MDB_Quad *q2 = find_quad(e1, e6, e9, e5);
    if(!q2) q2 = add_quad(e1, e6, e9, e5);
    MDB_Quad *q3 = find_quad(e2, e7, e10, e6);
    if(!q3) q3 = add_quad(e2, e7, e10, e6);
    MDB_Quad *q4 = find_quad(e3, e8, e11, e7);
    if(!q4) q4 = add_quad(e3, e8, e11, e7);
    MDB_Quad *q5 = find_quad(e4, e5, e12, e8);
    if(!q5) q5 = add_quad(e4, e5, e12, e8);
    MDB_Quad *q6 = find_quad(e9, e10, e11, e12);
    if(!q6) q6 = add_quad(e9, e10, e11, e12);
    return add_hex (q1, q2, q3, q4, q5, q6, pg);
  }

  MDB_Hex *MDB_Mesh::add_hex(MDB_Point *p1, MDB_Point *p2,MDB_Point *p3, MDB_Point *p4,
                             MDB_Point *p5, MDB_Point *p6,MDB_Point *p7, MDB_Point *p8, pGEntity pg)
  {
    MDB_Edge *e1 = find_edge(p1,p2);
    if (!e1) e1 = add_edge(p1, p2);
    MDB_Edge *e2 = find_edge(p2,p3);
    if (!e2) e2 = add_edge(p2, p3);
    MDB_Edge *e3 = find_edge(p3,p4);
    if (!e3) e3 = add_edge(p3, p4);
    MDB_Edge *e4 = find_edge(p4,p1);
    if (!e4) e4 = add_edge(p4, p1);
    MDB_Edge *e5 = find_edge(p1,p5);
    if (!e5) e5 = add_edge(p1, p5);
    MDB_Edge *e6 = find_edge(p2,p6);
    if (!e6) e6 = add_edge(p2, p6);
    MDB_Edge *e7 = find_edge(p3,p7);
    if (!e7) e7 = add_edge(p3, p7);
    MDB_Edge *e8 = find_edge(p4,p8);
    if (!e8) e8 = add_edge(p4, p8);
    MDB_Edge *e9 = find_edge(p5,p6);
    if (!e9) e9 = add_edge(p5, p6);
    MDB_Edge *e10 = find_edge(p6,p7);
    if (!e10) e10 = add_edge(p6, p7);
    MDB_Edge *e11 = find_edge(p7,p8);
    if (!e11) e11 = add_edge(p7, p8);
    MDB_Edge *e12 = find_edge(p8,p5);
    if (!e12) e12 = add_edge(p8, p5);

    MDB_Quad *q1 = find_quad(e4, e3, e2, e1);
    if(!q1) q1 = add_quad(e4, e3, e2, e1);
    MDB_Quad *q2 = find_quad(e1, e6, e9, e5);
    if(!q2) q2 = add_quad(e1, e6, e9, e5);
    MDB_Quad *q3 = find_quad(e2, e7, e10, e6);
    if(!q3) q3 = add_quad(e2, e7, e10, e6);
    MDB_Quad *q4 = find_quad(e3, e8, e11, e7);
    if(!q4) q4 = add_quad(e3, e8, e11, e7);
    MDB_Quad *q5 = find_quad(e4, e5, e12, e8);
    if(!q5) q5 = add_quad(e4, e5, e12, e8);
    MDB_Quad *q6 = find_quad(e9, e10, e11, e12);
    if(!q6) q6 = add_quad(e9, e10, e11, e12);
    return add_hex (q1, q2, q3, q4, q5, q6, pg);
  }




  MDB_Prism *MDB_Mesh::add_prism(MDB_Triangle *t1,MDB_Triangle *t2,
                                 MDB_Quad *t3,MDB_Quad *t4,MDB_Quad *t5, pGEntity pg)
  {
    MDB_Prism *p = new MDB_Prism(t1, t2, t3, t4, t5);

    //SL has not implemented the check on volume...
    //double v = VOLUME (p);
    //if (v < 0)
    //  {
    //   p->f1 = t2;
    //   p->f2 = t1;
    //  //      v = VOLUME (p);
    //  //      if (v < 0)
    //  //	throw;
    //  }

    p->g = pg;
    prisms.push_back(p);
    nbPrisms++;
    return p;
  }

  MDB_Prism *MDB_Mesh::add_prism(int p1, int p2, int p3,
                                 int p4, int p5, int p6, pGEntity pg)
  {
    MDB_Edge *e1 = add_edge(p1, p2);
    MDB_Edge *e2 = add_edge(p1, p3);
    MDB_Edge *e3 = add_edge(p2, p3);
    MDB_Edge *e4 = add_edge(p4, p5);
    MDB_Edge *e5 = add_edge(p4, p6);
    MDB_Edge *e6 = add_edge(p5, p6);
    MDB_Edge *e7= add_edge(p1, p4);
    MDB_Edge *e8 = add_edge(p2, p5);
    MDB_Edge *e9 = add_edge(p3, p6);

    MDB_Triangle *q1 = find_triangle(e1, e2, e3);
    if(!q1)
      q1 = add_triangle(e1, e2, e3);
    MDB_Triangle *q2 = find_triangle(e4, e6, e5);
    if(!q2)
      q2 = add_triangle(e4, e6, e5);
    MDB_Quad *q3 = find_quad(e7, e4, e8, e1);
    if(!q3)
      q3 = add_quad(e7, e4, e8, e1);
    MDB_Quad *q4 = find_quad(e8, e6, e9, e3);
    if(!q4)
      q4 = add_quad(e8, e6, e9, e3);
    MDB_Quad *q5 = find_quad(e7, e2, e9, e5);
    if(!q5)
      q5 = add_quad(e7, e2, e9, e5);
    return add_prism (q1, q2, q3, q4, q5, pg);
  }

  MDB_Prism *MDB_Mesh::add_prism(MDB_Point *p1, MDB_Point *p2,MDB_Point *p3,
                                 MDB_Point *p4, MDB_Point *p5, MDB_Point *p6, pGEntity pg)
  {
    MDB_Edge *e1 = find_edge(p1,p2);
    if (!e1) e1 = add_edge(p1, p2);
    MDB_Edge *e2 = find_edge(p1,p3);
    if (!e2) e2 = add_edge(p1, p3);
    MDB_Edge *e3 = find_edge(p2,p3);
    if (!e3) e3 = add_edge(p2, p3);
    MDB_Edge *e4 = find_edge(p4,p5);
    if (!e4) e4 = add_edge(p4, p5);
    MDB_Edge *e5 = find_edge(p4,p6);
    if (!e5) e5 = add_edge(p4, p6);
    MDB_Edge *e6 = find_edge(p5,p6);
    if (!e6) e6 = add_edge(p5, p6);
    MDB_Edge *e7 = find_edge(p1,p4);
    if (!e7) e7 = add_edge(p1, p4);
    MDB_Edge *e8 = find_edge(p2,p5);
    if (!e8) e8 = add_edge(p2, p5);
    MDB_Edge *e9 = find_edge(p3,p6);
    if (!e9) e9 = add_edge(p3, p6);

    MDB_Triangle *q1 = find_triangle(e1, e2, e3);
    if(!q1) q1 = add_triangle(e1, e2, e3);
    MDB_Triangle *q2 = find_triangle(e4, e6, e5);
    if(!q2) q2 = add_triangle(e4, e6, e5);
    MDB_Quad *q3 = find_quad(e7, e4, e8, e1);
    if(!q3) q3 = add_quad(e7, e4, e8, e1);
    MDB_Quad *q4 = find_quad(e8, e6, e9, e3);
    if(!q4) q4 = add_quad(e8, e6, e9, e3);
    MDB_Quad *q5 = find_quad(e7, e2, e9, e5);
    if(!q5) q5 = add_quad(e7, e2, e9, e5);
    return add_prism (q1, q2, q3, q4, q5, pg);
  }

  void MDB_Mesh::del_hex(MDB_Hex * h)
  {
    h->f1->del(h);
    h->f2->del(h);
    h->f3->del(h);
    h->f4->del(h);
    h->f5->del(h);
    h->f6->del(h);
    nbHexes--;
    h->deleted = true;
  }

  void MDB_Mesh::del_prism(MDB_Prism * p)
  {
    p->f1->del(p);
    p->f2->del(p);
    p->f3->del(p);
    p->f4->del(p);
    p->f5->del(p);
    nbPrisms--;
    p->deleted = true;
  }

  void MDB_Mesh::del_triangle(MDB_Triangle * t)
  {
    t->e1->del(t);
    t->e2->del(t);
    t->e3->del(t);
    nbTriangles--;
    t->deleted = true;
    if (t->getNbRegions())throw;
  }

  void MDB_Mesh::del_quad(MDB_Quad * q)
  {
    q->e1->del(q);
    q->e2->del(q);
    q->e3->del(q);
    q->e4->del(q);
    nbQuads--;
    q->deleted = true;
    if (q->getNbRegions())throw;
  }

  void MDB_Mesh::del_edge(MDB_Edge * e)
  {
    e->p1->del(e);
    e->p2->del(e);
    nbEdges--;
    e->deleted = true;

    if (e->numfaces())throw;
  }

  void MDB_Mesh::del_point(MDB_Point * p)
  {
    if (p->deleted) { 
      MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                  "point already deleted");
    }
    p->deleted = true;
    nbPoints--;
    if (p->edges.size())throw;
  }


  void MDB_Edge::oppositeof(MDB_Point * oface[2]) const
  {
    oface[0] = oface[1] = 0;
    if(faces(0)) {
      MDB_Point *pts[3];
      faces(0)->getNodes(pts);
      if(pts[0] != p1 && pts[0] != p2)
        oface[0] = pts[0];
      else if(pts[1] != p1 && pts[1] != p2)
        oface[0] = pts[1];
      else
        oface[0] = pts[2];
    }
    if(faces(1)) {
      MDB_Point *pts[3];
      faces(1)->getNodes(pts);
      if(pts[0] != p1 && pts[0] != p2)
        oface[1] = pts[0];
      else if(pts[1] != p1 && pts[1] != p2)
        oface[1] = pts[1];
      else
        oface[1] = pts[2];
    }
  }

  template < class IT > void DESTROOOY(IT beg, IT end)
  {
    while(beg != end) {
      delete *beg;
      beg++;
    }
  }

  MDB_Mesh ::~ MDB_Mesh ()
  {
    DESTROOOY ( points.begin(),points.end());
    DESTROOOY ( edges.begin(),edges.end());
    DESTROOOY ( triangles.begin(),triangles.end());
    DESTROOOY ( quads.begin(),quads.end());
    DESTROOOY ( tets.begin(),tets.end());
    DESTROOOY ( hexes.begin(),hexes.end());
    DESTROOOY ( prisms.begin(),prisms.end());
  }

  MDB_Edge::~MDB_Edge ()
  {
  }


  MDB_Mesh::MDB_Mesh(const MDB_Mesh & other)
  {
    throw;
  }

  void MDB_Mesh :: load ( FILE *f )
  {
    _load (f);
    shrinked = true;
  }

  void MDB_Mesh :: flush ( FILE *f ) const
  {
    _flush (f);
  }

  void MDB_Mesh :: shrink ()
  {
    VIter vit  = M_vertexIter ( this );
    pVertex pv;
    GEntity2Physical gentity2phys(geomFeatures_Tags);
    coords = new double[3*nbPoints];
    ids    = new int   [  nbPoints];
    int i=0;
    nbInfos = 0;

    int SIZE_INFO_PT       = 4;
    int SIZE_INFO_EDGE     = 5;
    int SIZE_INFO_TRIANGLE = 6;
    int SIZE_INFO_TET      = 7;

    while ((pv = VIter_next (vit)))
      {
        int dim = EN_whatInType(pv);
        if (dim == 0) nbInfos += SIZE_INFO_PT;
        ids[i] = pv->iD;
        pv->iD= i+1;
        V_coord ( pv, &coords[3*i++]);
      }
    nbVertex = nbPoints;
    EIter eit  = M_edgeIter ( this );
    pEdge pe;
    while ((pe = EIter_next (eit)))
      {
        int dim = EN_whatInType(pe);
        if (dim == 1) nbInfos += SIZE_INFO_EDGE;
      }
    FIter fit  = M_faceIter ( this );
    pFace pf;
    while ((pf = FIter_next (fit)))
      {
        int dim = EN_whatInType(pf);
        if (dim == 2) nbInfos += SIZE_INFO_TRIANGLE;
      }
    nbInfos += SIZE_INFO_TET * nbTets;

    infos = new int [nbInfos];

    VIter_reset(vit);
    EIter_reset(eit);
    FIter_reset(fit);

    int K = 0;

    while ((pv = VIter_next (vit)))
      {
        int dim = EN_whatInType(pv);
        if (dim == 0)
          {
            pGEntity pg = EN_whatIn(pv);
            int tag = GEN_tag ( pg );
            int phys = gentity2phys.get_first_tag( pg );
            infos[K++] = 15;
            infos[K++] = phys;
            infos[K++] = tag;
            infos[K++] = EN_id(pv);
          }
      }

    while ((pe = EIter_next (eit)))
      {
        int dim = EN_whatInType(pe);
        if (dim == 1)
          {
            pGEntity pg = EN_whatIn(pe);
            int tag = GEN_tag ( pg );
            int phys = gentity2phys.get_first_tag( pg );
            pVertex v1 = E_vertex ( pe, 0);
            pVertex v2 = E_vertex ( pe, 1);
            infos[K++] = 1;
            infos[K++] = phys;
            infos[K++] = tag;
            infos[K++] = EN_id (v1);
            infos[K++] = EN_id (v2);
          }
      }
    while ((pf = FIter_next (fit)))
      {
        int dim = EN_whatInType(pf);
        if (dim == 2)
          {
            pGEntity pg = EN_whatIn(pf);
            int tag = GEN_tag ( pg );
            int phys = gentity2phys.get_first_tag( pg );
            pVertex v1 = F_vertex ( pf, 0);
            pVertex v2 = F_vertex ( pf, 1);
            pVertex v3 = F_vertex ( pf, 2);
            infos[K++] = 2;
            infos[K++] = phys;
            infos[K++] = tag;
            infos[K++] = EN_id (v1);
            infos[K++] = EN_id (v2);
            infos[K++] = EN_id (v3);
          }
      }

    pRegion pr;
    RIter rit  = M_regionIter ( this );
    while ((pr = RIter_next (rit)))
      {
        pGEntity pg = EN_whatIn(pr);
        int tag = GEN_tag ( pg );
        int phys = gentity2phys.get_first_tag( pg );
        pVertex v1 = R_vertex ( pr, 0);
        pVertex v2 = R_vertex ( pr, 1);
        pVertex v3 = R_vertex ( pr, 2);
        pVertex v4 = R_vertex ( pr, 3);
        infos[K++] = 4;
        infos[K++] = phys;
        infos[K++] = tag;
        infos[K++] = EN_id (v1);
        infos[K++] = EN_id (v2);
        infos[K++] = EN_id (v3);
        infos[K++] = EN_id (v4);
      }

    //  printf("%d %d\n",K,nbInfos);

    VIter_delete(vit);
    EIter_delete(eit);
    FIter_delete(fit);
    RIter_delete(rit);

    DESTROOOY ( points.begin(),points.end());
    DESTROOOY ( edges.begin(),edges.end());
    DESTROOOY ( triangles.begin(),triangles.end());
    DESTROOOY ( tets.begin(),tets.end());
    DESTROOOY ( quads.begin(),quads.end());
    DESTROOOY ( hexes.begin(),hexes.end());
    DESTROOOY ( prisms.begin(),prisms.end());

    points.clear();
    edges.clear();
    triangles.clear();
    quads.clear();
    tets.clear();
    hexes.clear();
    prisms.clear();

    nbPoints = nbTriangles = nbEdges = nbTets = nbQuads = nbHexes = nbPrisms = 0;
    maxId=0;

    shrinked = true;
  }

  // MS
  void MDB_Mesh :: classify_grid_boundaries(pGEntity bndGEntity)
  {
    // going through all faces and find boundary ones
    FIter fit = M_faceIter(this);
    pFace pf;
    int nBndFace = 0;
    while ((pf = FIter_next(fit)))
      {
       if (pf->getNbRegions() < 2){
         // applying given entity to the boundary face
         // and all of its edges and all its vertices
         nBndFace++;
         pf->g = bndGEntity;
         pEdge e1 = F_edge (pf,0);
         e1->g = pf->g;
         pEdge e2 = F_edge (pf,1);
         e2->g = pf->g;
         pEdge e3 = F_edge (pf,2);
         e3->g = pf->g;
         pEdge e4 = F_edge (pf,3);
         if (e4 ) e4->g = pf->g;
         //vertices
         pVertex v1 = F_vertex(pf, 0);
         v1->g = pf->g;
         pVertex v2 = F_vertex(pf, 1);
         v2->g = pf->g;
         pVertex v3 = F_vertex(pf, 2);
         v3->g = pf->g;
         for (int iv=3; iv<pf->getNbNodes(); iv++){
           pVertex v4 = F_vertex(pf, iv);
           if (v4) v4->g = pf->g;
         }
       }
      }

    std::cout << "Found " << nBndFace << " boundary faces.\n";
    FIter_delete(fit);

  }

  void MDB_Mesh :: unclassify_grid_boundaries ()
  {
    int nbV=0, nbE=0, nbF=0, nbR=0;

    RIter rit = M_regionIter(this);
    pRegion pr;
    while ((pr = RIter_next(rit)))
      {
        if (pr->g)
          {
            pFace f1 = R_face (pr,0);
            f1->g = pr->g;
            pFace f2 = R_face (pr,1);
            f2->g = pr->g;
            pFace f3 = R_face (pr,2);
            f3->g = pr->g;
            pFace f4 = R_face (pr,3);
            f4->g = pr->g;
            pFace f5 = R_face (pr,4);
            if (f5) f5->g = pr->g;
            pFace f6 = R_face (pr,5);
            if (f6) f6->g = pr->g;
          }
        else
          {
            nbR++;
          }
      }

    RIter_delete(rit);

    FIter fit = M_faceIter(this);
    pFace pf;
    while ((pf = FIter_next(fit)))
      {
        if (pf->g)
          {
            pEdge e1 = F_edge (pf,0);
            e1->g = pf->g;
            pEdge e2 = F_edge (pf,1);
            e2->g = pf->g;
            pEdge e3 = F_edge (pf,2);
            e3->g = pf->g;
            pEdge e4 = F_edge (pf,3);
            if (e4) e4->g = pf->g;
          }
        else
          {
            nbF++;
          }
      }
    FIter_delete(fit);


    EIter eit = M_edgeIter(this);
    pEdge pe;
    while ((pe = EIter_next(eit)))
      {
        if (pe->g)
          {
            pVertex v1 = E_vertex (pe,0);
            v1->g = pe->g;
            pVertex v2 = E_vertex (pe,1);
            v2->g = pe->g;
          }
        else
          {
            nbE++;
          }
      }
    EIter_delete(eit);

    VIter vit = M_vertexIter(this);
    pVertex pv;
    while ((pv = VIter_next(vit)))
      {
        if (!pv->g)
          {
            nbV++;
          }
      }
    VIter_delete(vit);

    if ( nbV || nbE || nbF || nbR ) {
      MAdMsgSgl::instance().info(__LINE__,__FILE__,
                                 "There are unclassified entities: %d, regions, %d faces, %d edges, %d vertices\n",
                                 nbR,nbF,nbE,nbV);
    }
  }

  void MDB_Mesh :: skin_me(MDB_Mesh *skinMesh, std::vector<int>& regIds)
  {
    // going through all faces and find boundary ones
    FIter fit = M_faceIter(this);
    pFace pf;
    int nBndVrtx = 0;
    int nBndFace = 0;
    int nNewRegion = 1;
    std::set<int> bndVrtxIds;
    while ((pf = FIter_next(fit)))
    {
       if (pf->getNbRegions() < 2)
       {
         // applying given entity to the boundary face
         // and all of its edges and all its vertices
         nBndFace++;
//         pf->g = bndGEntity;
//         pEdge e1 = F_edge (pf,0);
//         e1->g = pf->g;
//         pEdge e2 = F_edge (pf,1);
//         e2->g = pf->g;
//         pEdge e3 = F_edge (pf,2);
//         e3->g = pf->g;
//         pEdge e4 = F_edge (pf,3);
//         if (e4 ) e4->g = pf->g;
         MDB_Region *myReg = pf->getRegion(0);
         if (myReg->iD == 0)
         {
           // this means we have  new region after
           // mesh adaptation we number the new ones with
           // negative numbers
           myReg->iD = -nNewRegion++;
         }
         regIds.push_back(myReg->iD);

         //std::cout << "My reg id = " << myReg->iD << std::endl;
         //vertices
         pVertex v1 = F_vertex(pf, 0);
         pVertex v2 = F_vertex(pf, 1);
         pVertex v3 = F_vertex(pf, 2);
         if (pf->getNbNodes()>3)
         {
           std::cerr << "Only triangular elements can be added to skin mesh.\n";
           continue;
         }
         // adding point
         //MAd::pGEntity pGeom = (MAd::pGEntity) MAd::GM_vertexByTag(skinMesh->model, 0);
         if (bndVrtxIds.find(v1->iD) == bndVrtxIds.end())
         {
           bndVrtxIds.insert(v1->iD);
           skinMesh->add_point(v1->iD, v1->X, v1->Y, v1->Z);
           nBndVrtx+=1;
         }
         if (bndVrtxIds.find(v2->iD) == bndVrtxIds.end())
         {
           bndVrtxIds.insert(v2->iD);
           skinMesh->add_point(v2->iD, v2->X, v2->Y, v2->Z);
           nBndVrtx+=1;
         }
         if (bndVrtxIds.find(v3->iD) == bndVrtxIds.end())
         {
           bndVrtxIds.insert(v3->iD);
           skinMesh->add_point(v3->iD, v3->X, v3->Y, v3->Z);
           nBndVrtx+=1;
         }
         // ading triangle
         MAd::pGEntity fGeom = (MAd::pGEntity) MAd::GM_faceByTag(skinMesh->model, 0);
         skinMesh->add_triangle(v1->iD, v2->iD, v3->iD, fGeom);

       }
    }
    skinMesh->classify_unclassified_entities();
    skinMesh->destroyStandAloneEntities();
    std::cout << nBndVrtx << " boundary vertices added to the mesh.\n";
    std::cout << nBndFace << " boundary faces added to the mesh.\n";
    FIter_delete(fit);

    // going through all regions to find corresponding region index for ids
    // since mesh might have updated this needs to be checked everytime
    RIter rit = M_regionIter(this);
    pRegion pr;
    int iReg = 1;
    std::map<int,int> old2new;
    while ( (pr = RIter_next(rit)) )
    {
      //std::cout << pr->iD << " -> " << iReg << std::endl;
      old2new[pr->iD] = iReg++;
    }
    RIter_delete(rit);
    // converting region ids from old to new
    for (int it=0; it<regIds.size(); it++)
      regIds[it] = old2new[regIds[it]];
  }
  // MS END

  void MDB_Mesh :: classify_unclassified_entities ()
  {
    int nbV=0, nbE=0, nbF=0, nbR=0;

    RIter rit = M_regionIter(this);
    pRegion pr;
    while ((pr = RIter_next(rit)))
      {
        if (pr->g)
          {
            pFace f1 = R_face (pr,0);
            if (!f1->g) f1->g = pr->g;
            pFace f2 = R_face (pr,1);
            if (!f2->g) f2->g = pr->g;
            pFace f3 = R_face (pr,2);
            if (!f3->g) f3->g = pr->g;
            pFace f4 = R_face (pr,3);
            if (!f4->g) f4->g = pr->g;
            pFace f5 = R_face (pr,4);
            if (f5 && !f5->g) f5->g = pr->g;
            pFace f6 = R_face (pr,5);
            if (f6 && !f6->g) f6->g = pr->g;
          }
        else
          {
            nbR++;
          }
      }

    RIter_delete(rit);

    FIter fit = M_faceIter(this);
    pFace pf;
    while ((pf = FIter_next(fit)))
      {
        if (pf->g)
          {
            pEdge e1 = F_edge (pf,0);
            if (!e1->g) e1->g = pf->g;
            pEdge e2 = F_edge (pf,1);
            if (!e2->g) e2->g = pf->g;
            pEdge e3 = F_edge (pf,2);
            if (!e3->g) e3->g = pf->g;
            pEdge e4 = F_edge (pf,3);
            if (e4 && !e4->g) e4->g = pf->g;
          }
        else
          {
            nbF++;
          }
      }
    FIter_delete(fit);


    EIter eit = M_edgeIter(this);
    pEdge pe;
    while ((pe = EIter_next(eit)))
      {
        if (pe->g)
          {
            pVertex v1 = E_vertex (pe,0);
            if (!v1->g) v1->g = pe->g;
            pVertex v2 = E_vertex (pe,1);
            if (!v2->g) v2->g = pe->g;
          }
        else
          {
            nbE++;
          }
      }
    EIter_delete(eit);

    VIter vit = M_vertexIter(this);
    pVertex pv;
    while ((pv = VIter_next(vit)))
      {
        if (!pv->g)
          {
            nbV++;
          }
      }
    VIter_delete(vit);

    if ( nbV || nbE || nbF || nbR ) {
      MAdMsgSgl::instance().info(__LINE__,__FILE__,
                                 "There are unclassified entities: %d, regions, %d faces, %d edges, %d vertices\n",
                                 nbR,nbF,nbE,nbV);
    }
  }

  void MDB_Mesh :: destroyStandAloneEntities() {

    // First build a list of entities < mesh dim and not to be deleted
    // = those which are part of an entity of dimension = mesh dim
    set<pVertex> usedV;
    set<pEdge>   usedE;
    set<pFace>   usedF;

    // 3D case
    if ( nbTets + nbHexes + nbPrisms != 0 ) {

      RIter rit = M_regionIter(this);
      pRegion pr;
      while ((pr = RIter_next(rit)))
        {
          for (int i=0; i < pr->getNbFace(); i++) {
            pFace f = R_face (pr,i);
            usedF.insert(f);
          }
          for (int i=0; i < pr->getNbEdge(); i++) {
            pEdge e = R_edge (pr,i);
            usedE.insert(e);
          }
          for (int i=0; i < pr->getNbVertex(); i++) {
            pVertex v = R_vertex (pr,i);
            usedV.insert(v);
          }
        }
      RIter_delete(rit);
    }

    // 2D case
    else if ( nbTriangles + nbQuads != 0 ) {

      FIter fit = M_faceIter(this);
      pFace pf;
      while ((pf = FIter_next(fit)))
        {
          usedF.insert(pf);
          for (int i=0; i < pf->getNbEdges(); i++) {
            pEdge e = F_edge (pf,i);
            usedE.insert(e);
          }
          for (int i=0; i < pf->getNbNodes(); i++) {
            pVertex v = F_vertex (pf,i);
            usedV.insert(v);
          }
        }
      FIter_delete(fit);
    }

    // 1D case
    else if ( nbEdges != 0 ) {

      EIter eit = M_edgeIter(this);
      pEdge pe;
      while ((pe = EIter_next(eit)))
        {
          usedE.insert(pe);
          for (int i=0; i < 2; i++) {
            pVertex v = E_vertex (pe,i);
            usedV.insert(v);
          }
        }
      EIter_delete(eit);
    }

    // 0D case
    else {

      VIter vit = M_vertexIter(this);
      pVertex pv;
      while ((pv = VIter_next(vit)))
        {
          usedV.insert(pv);
        }
      VIter_delete(vit);
    }

    // then delete the others
    int delV = 0;
    int delE = 0;
    int delF = 0;

    FIter fit = M_faceIter(this);
    pFace pf;
    while ((pf = FIter_next(fit)))
      {
        if (usedF.find(pf) == usedF.end() ) {
          del_face(pf);
          delF++;
        }
      }
    FIter_delete(fit);

    EIter eit = M_edgeIter(this);
    pEdge pe;
    while ((pe = EIter_next(eit)))
      {
        if (usedE.find(pe) == usedE.end() ) {
// #warning "debug"
//           printf("Looking for a clone...");
//           pEdge peClone = find_clone_edge(pe);
//           if (peClone) printf("found one\n");
//           else         printf("no clone\n");
//           if (peClone) {
//             printf("This edge will be deleted:\n");
//             E_info(pe);
//             printf("This edge is its clone:\n");
//             E_info(peClone);
//           }
//           if ( (peClone) && ( GEN_type(peClone->g) > GEN_type(pe->g) ) ) {
//             peClone->g = pe->g;
//           }

          del_edge(pe);
          delE++;
        }
      }
    EIter_delete(eit);

    VIter vit = M_vertexIter(this);
    pVertex pv;
    while ((pv = VIter_next(vit)))
      {
        if (usedV.find(pv) == usedV.end() ) {
          del_point(pv);
          delV++;
        }
      }
    VIter_delete(vit);

    if ( delV != 0 || delE != 0 || delF != 0 ) {
      cout <<"WARNING: deleted stand alone entities: "
           << delV <<" vertices, "
           << delE <<" edges, "
           << delF <<" faces\n";
    }
  }

  // void MDB_Mesh :: checkGeomTagsCoherence () const
  // {
  //   std::set<int> volumeTags;
  //   std::multimap<int, pGEntity >::const_iterator tagIt = geomFeatures_Tags.begin();
  //   for (; tagIt != geomFeatures_Tags.end(); tagIt++)
  //     if (GEN_type(tagIt->second) == 3)
  //       volumeTags.insert(tagIt->first);

  //   tagIt = geomFeatures_Tags.begin();
  //   for (; tagIt != geomFeatures_Tags.end(); tagIt++)
  //     if (GEN_type(tagIt->second) == 2)
  //       if ( volumeTags.find(tagIt->first) != volumeTags.end() ) {
  // 	printf("Wrong input geometry: the geometric ID \'%d\' is used for both a volume and a surface\n",tagIt->first);
  // 	exit(0);
  //       }
  // }

  void MDB_Mesh :: expand ()
  {
    initializeIdData();

    pVertex *vs = new pVertex[nbVertex];

    for (int i=0;i<nbVertex;i++)
      {
        //       vs[i] = add_point (i+1,coords[3*i],coords[3*i+1],coords[3*i+2]);
        vs[i] = add_point (ids[i],coords[3*i],coords[3*i+1],coords[3*i+2]);
      }


    nbVertex = 0;
    delete [] coords;

    int i=0;
    while (i<nbInfos)
      {
        int type = infos[i++];
        int physicalTag  = infos[i++];
        int elementaryTag = infos[i++];
        pGEntity g=NULL;

        if (type == 15)
          {
            pVertex p = vs[infos[i++]-1];
            g = GM_entityByTag(model, 0, elementaryTag);
            p->g = g;
          }
        else if (type == 1)
          {
            pVertex v1 = vs[infos[i++]-1];
            pVertex v2 = vs[infos[i++]-1];
            g = GM_entityByTag(model, 1, elementaryTag);
            add_edge(v1,v2,g);
          }
        else if (type == 2)
          {
            int v1 = ids[infos[i++]-1];
            int v2 = ids[infos[i++]-1];
            int v3 = ids[infos[i++]-1];
            g = GM_entityByTag(model, 2, elementaryTag);
            add_triangle(v1,v2,v3,g);
          }
        else if (type == 4)
          {
            pVertex v1 = vs[infos[i++]-1];
            pVertex v2 = vs[infos[i++]-1];
            pVertex v3 = vs[infos[i++]-1];
            pVertex v4 = vs[infos[i++]-1];
            g = GM_entityByTag(model, 3, elementaryTag);
            add_tet(v1,v2,v3,v4,g);
          }
	if (g)
	  {
	    bool find = false;
	    for (std::multimap<int, pGEntity>::iterator it = geomFeatures_Tags.lower_bound(physicalTag);
		 it != geomFeatures_Tags.upper_bound(physicalTag);++it)
	      if (it->second == g)find = true;
	    if (!find)
	      geomFeatures_Tags.insert(std::pair<int,pGEntity>(physicalTag, g));
	  }
      }
    delete [] infos; infos = NULL;
    delete [] ids;   ids   = NULL;
    delete [] vs;    vs    = NULL;
    nbInfos = 0;
    classify_unclassified_entities();
    shrinked = false;
  }

  int MDB_Tet::face_vtx[4][3] = {{0,1,2},{0,1,3},{1,2,3},{0,2,3}};

  // return 1 if the face direction points to outside of the tet
  // return 0                                 inside
  int MDB_Tet::getFaceDir (int n) const
  {
    int tetfacedir [] = {1,0,0,1};
    int corr = tetfacedir[n];
    pFace f  = getFace(n);
    pVertex v1f, v2f;
    v1f = F_vertex (f,0);
    v2f = F_vertex (f,1);

    pVertex v[3];
    switch (n)
      {
      case 0 : //inside
        v[0] = getVertex(0);
        v[1] = getVertex(1);
        v[2] = getVertex(2);
        break;
      case 1 : // outside
        v[0] = getVertex(0);
        v[1] = getVertex(1);
        v[2] = getVertex(3);
        break;
      case 2 : // outside
        v[0] = getVertex(1);
        v[1] = getVertex(2);
        v[2] = getVertex(3);
        break;
      case 3 : //inside
        v[0] = getVertex(0);
        v[1] = getVertex(2);
        v[2] = getVertex(3);
        break;
      default :
        throw;
      }


    for(int j=0;j<3;j++)
      {
        if (v[j] == v1f)
          {
            pVertex vtest = v[(j+1)%3];
            if (vtest == v2f) return 1-corr;
            vtest = v[(j+2)%3];
            if (vtest == v2f) return corr;
            throw;
          }
      }
    throw;
  }

  int MDB_Tet::getFaceOrientation (int n) const {

    pFace   f  = getFace(n);
    pVertex v0 = getVertex(face_vtx[n][0]);
    pVertex v1 = getVertex(face_vtx[n][1]);

    int orientation = 0;

    for (int j=0;j<3;j++) if (F_vertex(f,j) == v0) orientation = j;

    int direction = 1;
    if (v1 == F_vertex(f,(orientation+2)%3)) direction = -1;

    return direction * (orientation + 1);
  }


  int MDB_Hex::getFaceDir (int n) const {  return 0;}

  // -----------------------------------------------------------------------------
  /*! \brief Compute orientation wrt face given by a list of principal nodes \ingroup internal  */
  /*! abs(rot) - 1 = the index of the node in the current element corresponding to node 0 in the remote element */
  /*! sgn(rot)     = 1 if the element is oriented in the same fashion, -1 if inverted */

  int MDB_Face::orientation(MDB_Point* v0,
                            MDB_Point* v1,
                            MDB_Point* v2,
                            MDB_Point* v3) const {

    if ((getNbNodes() == 4) && (v3 == NULL)) return 0;

    int rot = 0;
    for (int j=0;j<getNbNodes();j++) {
      if (getNode(j) == v0) rot = j+1;
    }

    if (getNode(rot%getNbNodes()) != v1)  {
      if (getNode((rot-2+getNbNodes())%getNbNodes()) == v1) rot *= -1;
    }
    return rot;
  }

  // -----------------------------------------------------------------------------
  /*! \brief Realign the edge with the direction given by a set of two nodes \ingroup internal  */
  /*! returns 1 if already aligned, -1 if inversely aligned, 0 if failed  */

  int MDB_Edge::align(MDB_Point* po1,MDB_Point* po2)
  {
    if ((p1 == po1) && (p2 == po2)) return 1;
    if ((p2 == po1) && (p1 == po2)) {
      std::swap(p1,p2);
      swapHONodes();
      return -1;
    }
    return 0;
  }

  // -----------------------------------------------------------------------------
  /*! \brief Realign the triangle with the direction given by a set of three nodes \ingroup internal */
  /*! Returns 0 if failed, otherwise
    abs(rot) = index of node corresponding to v0
    sng(rot) = 1 if same orientation, -1 if inverted orientation */

  int MDB_Triangle::align(MDB_Point* v0,MDB_Point* v1,MDB_Point* v2,MDB_Point* v3) {

    int rot = orientation(v0,v1,v2,v3);

    if (rot != 1 && rot != 0) {

      MDB_Edge* tmpEdge[3];
      for (int i=0;i<3;i++) tmpEdge[i] = getEdge((i+abs(rot)-1)%3);
      
      if (rot < 0) {
        MDB_Edge* tmp = tmpEdge[0];
        tmpEdge[0] = tmpEdge[2];
        tmpEdge[2] = tmp;
      }

      e1 = tmpEdge[0];
      e2 = tmpEdge[1];
      e3 = tmpEdge[2];
    }
    return rot;
  }

  // -----------------------------------------------------------------------------
  /*! \brief Realign the quadrilateral with the direction given by a set of four nodes \ingroup internal */
  /*! Returns 0 if failed, otherwise
    abs(rot) = index of node corresponding to v0
    sng(rot) = 1 if same orientation, -1 if inverted orientation */

  int MDB_Quad::align(MDB_Point* v0,MDB_Point* v1,MDB_Point* v2,MDB_Point* v3) {

    int rot = orientation(v0,v1,v2,v3);
    if (rot != 1 && rot != 0) {

      MDB_Edge* tmpEdge[4];
      for (int i=0;i<getNbEdges();i++) tmpEdge[i] = getEdge((i+abs(rot)-1)%getNbEdges());

      if (rot < 0) {
        MDB_Edge* tmp = tmpEdge[0];
        tmpEdge[0] = tmpEdge[3];
        tmpEdge[3] = tmp;
        tmp = tmpEdge[1];
        tmpEdge[1] = tmpEdge[2];
        tmpEdge[2] = tmp;
      }

      e1 = tmpEdge[0];
      e2 = tmpEdge[1];
      e3 = tmpEdge[2];
      e4 = tmpEdge[3];
    }

    return rot;

  }

} // End of namespace MAd

