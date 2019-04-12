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
// Authors: Jean-Francois Remacle, Gaetan Compere, Koen Hillewaert
// -------------------------------------------------------------------

#ifndef H_MESHDATABASE
#define H_MESHDATABASE
#include <set>
#include <vector>
#include <algorithm>
#include <list>
#include <map>
#include <cmath>
#include <assert.h>
#include <cstdarg>
#include <string.h>
#include <iostream>
#include "ModelInterface.h"
#include "MeshDataBaseIterators.h"
#include "MeshDataBaseAttachable.h"
#include "MeshDataBaseMiniMesh.h"
#include "MAdMessage.h"
#include "MshTags.h"

namespace MAd {

  class MDB_Region;
  class MDB_Tet;
  class MDB_Hex;
  class MDB_Prism;
  class MDB_Quad;
  class MDB_Face;
  class MDB_Edge;
  class MDB_EdgeP2;
  class MDB_Triangle;
  class MDB_Mesh;
  class MDB_Point;

  typedef std::list<MDB_Point *>      MDB_ListP ;
  typedef std::list<MDB_Edge *>       MDB_ListE ;
  typedef std::list<MDB_Triangle *>   MDB_ListF ;
  typedef std::list<MDB_Quad*>        MDB_ListQ ;
  typedef std::list<MDB_Face*>        MDB_ListFace;
  typedef std::list<MDB_Tet *>        MDB_ListT ;
  typedef std::list<MDB_Hex *>        MDB_ListH ;
  typedef std::list<MDB_Prism *>      MDB_ListPr;
  typedef std::list<MDB_Region *>     MDB_ListR ;
  typedef std::vector<MDB_Triangle *> MDB_VectorF ;
  typedef std::vector<MDB_Quad *>     MDB_VectorQ ;
  typedef std::vector<MDB_Face *>     MDB_VectorFace ;
  typedef std::vector<MDB_Edge *>     MDB_VectorE ;
  typedef std::vector<MDB_Tet *>      MDB_VectorT ;
  typedef std::vector<MDB_Hex *>      MDB_VectorH ;
  typedef std::vector<MDB_Prism *>    MDB_VectorPr;

  class MDB_MeshEntity : public mAttachableDataContainer  
  {
  public :
    int iD;
    bool deleted;
    pGEntity g;
    virtual int getDim()   const = 0;
    virtual int getOrder() const = 0;
    virtual int getMshTag() const = 0;
    MDB_MeshEntity ( int _ID = 0, pGEntity _g = 0) : iD(_ID), deleted (false), g(_g) {}    
    virtual ~MDB_MeshEntity() {}
  };

  /*! \brief Vertex */ 

  class MDB_Point : public MDB_MeshEntity
  {
  public:
    // SHOULD BE PRIVATE, 2B DONE
    double X,Y,Z;
    MDB_VectorE edges;

    virtual int getDim()   const {return 0;}
    virtual int getOrder() const {return 0;}
    virtual int getMshTag() const {return MSH_PNT;}
    virtual bool isParametric() const { return false; }
    virtual bool getParams(double * u, double * v) const { return false; }

    inline bool operator < (const MDB_Point & other) const
    {
      return iD < other.iD;
    }
    inline void del(MDB_Edge *e)
    {
      MDB_VectorE::iterator it  = edges.begin();
      MDB_VectorE::iterator ite = edges.end();
      while(it!=ite){
        if(*it == e){
          edges.erase(it);
          break;
        }
        ++it;
      }
    }
    void getFaces(MDB_ListFace &) const; 	
    void getTets(MDB_ListT &) const; 
    void getHexes(MDB_ListH &) const;
  
    //void getHexs     (MDB_ListH &) const;

    MDB_Point(int id, double x=0, double y=0, double z=0) 
      : MDB_MeshEntity(id),X(x),Y(y),Z(z)
    {	    
      //    edges.reserve(12);
    }
    virtual ~MDB_Point() {}
  };

  /*! \brief Parametric vertex */ 

  class MDB_PointParam : public MDB_Point
  {
  public:

    MDB_PointParam(int id, double x=0, double y=0, double z=0, double u=0, double v=0)
      :  MDB_Point(id,x,y,z),U(u),V(v)
    {}
    ~MDB_PointParam() {}

    bool isParametric() const { return true; }
    bool getParams(double * u, double * v) const { *u=U; *v=V; return true; }

  private:

    double U,V;
  };

  /*! \brief Straight edge or base class for higher order edges */ 

  class MDB_Edge: public MDB_MeshEntity
  {
    MDB_VectorFace _faces;
  public:

    MDB_Point *p1,*p2;

    virtual int getDim()   const {return 1;}
    virtual int getOrder() const {return 1;}
    virtual int getMshTag() const {return MSH_LIN_2;}

    inline MDB_Face* faces(int i) const
    {
      if ( _faces.size() == 0 ) return NULL;
      if ( _faces.size() == 1 && i == 1 ) return NULL;
      return _faces [i];
    }

    inline double length() const
    {
      return sqrt((p1->X-p2->X)*(p1->X-p2->X)+(p1->Y-p2->Y)*(p1->Y-p2->Y)+(p1->Z-p2->Z)*(p1->Z-p2->Z));
    }

    inline int numfaces() const
    {
      return _faces.size();
    }

    inline MDB_Point * commonvertex(const MDB_Edge *other) const
    {
      if(p1 == other->p1 || p1 == other->p2) return p1;
      if(p2 == other->p1 || p2 == other->p2) return p2;
      printf("two edges %d %d -- %d %d without common vertex\n",p1->iD,p2->iD, other->p1->iD, other->p2->iD);
      throw;
    }

    inline MDB_Point * othervertex(const MDB_Point *p) const
    {
      if(p1 == p) return p2;
      if(p2 == p) return p1;
      return 0;
    }

    inline void addface(MDB_Face *f)
    {
      _faces.push_back(f);
    }

    inline MDB_Face * otherFace(const MDB_Face *f) const
    {
      if(numfaces()!=2) return NULL;
      if(f == _faces[0]) return _faces[1];
      if(f == _faces[1]) return _faces[0];
      throw;
    }

    inline void del( MDB_Face *t )
    {
      _faces.erase(std::remove_if(_faces.begin(),_faces.end() , std::bind2nd(std::equal_to<MDB_Face*>(), t)) , 
                   _faces.end());
    }
  
    inline void oppositeof(MDB_Point * oface[2]) const; 
  
    virtual int getNbHighOrderPoints () const {return 0;}
    virtual MDB_Point *getHighOrderPoint (int i) {throw;}

    MDB_Edge(MDB_Point *A, MDB_Point *B)
      : p1(A),p2(B)
    {
//       _faces.reserve(2);
//       _quads.reserve(1);
      p1->edges.push_back(this);
      p2->edges.push_back(this);
    }

    virtual void swapHONodes() {}
    virtual int align(MDB_Point* po1,MDB_Point* po2);
    virtual ~MDB_Edge();
  };

  /*! \brief 2nd order equidistant Lagrange edge */

  class MDB_EdgeP2: public MDB_Edge
  {
    MDB_Point *secondOrderPt;
  public :
    virtual int getOrder() const {return 2;}
    virtual int getMshTag() const {return MSH_LIN_3;}
    virtual int getNbHighOrderPoints() const {return 1;}
    virtual MDB_Point *getHighOrderPoint (int i) {assert (i==0);return secondOrderPt;}
    MDB_EdgeP2(MDB_Point *A, MDB_Point *B, MDB_Point *C)
      : MDB_Edge ( A, C ), secondOrderPt(B)
    {
    }
    virtual ~MDB_EdgeP2() {delete secondOrderPt;}
  };

  /*! \brief 3rd order equidistant Lagrange edge */

  class MDB_EdgeP3: public MDB_Edge
  {
    MDB_Point *Pt1,*Pt2;//third order points
  public :
    virtual int getOrder() const {return 3;}
    virtual int getMshTag() const {return MSH_LIN_4;}
    virtual int getNbHighOrderPoints() const {return 2;}
    virtual MDB_Point *getHighOrderPoint (int i) {
      switch (i){
      case 0:
        return Pt1;
      case 1:
        return Pt2;
      default:
        std::cout << "Point " << i << "not defined for third order edges\n"; throw;
      }
    }
    MDB_EdgeP3(MDB_Point *A, MDB_Point *B1, MDB_Point *B2, MDB_Point *C)
      : MDB_Edge ( A, C ), Pt1(B1), Pt2(B2)
    {
    }
    virtual ~MDB_EdgeP3() {delete Pt1;delete Pt2;}

    virtual void swapHONodes() {

      MDB_Point* tmp = Pt2;
      Pt2 = Pt1;
      Pt1 = tmp;
    }
  };

  /*! \brief 4th order equidistant Lagrange edge */

  class MDB_EdgeP4: public MDB_Edge
  {
    MDB_Point *Pt1,*Pt2,*Pt3;//fourth order points
  public :
    virtual int getOrder() const {return 4;}
    virtual int getMshTag() const {return MSH_LIN_5;}
    virtual int getNbHighOrderPoints() const {return 3;}
    virtual MDB_Point *getHighOrderPoint (int i) {
      switch (i){
      case 0:
        return Pt1;
      case 1:
        return Pt2;
      case 2:
        return Pt3;
      default:
        std::cout << "Point " << i << "not defined for fourth order edges\n"; throw;
      }
    }
    MDB_EdgeP4(MDB_Point *A, MDB_Point *B1, MDB_Point *B2, MDB_Point *B3, MDB_Point *C)
      : MDB_Edge ( A, C ), Pt1(B1), Pt2(B2), Pt3(B3)
    {
    }
    virtual ~MDB_EdgeP4() {delete Pt1;delete Pt2;delete Pt3;}

  protected:
  
    virtual void swapHONodes() 
    {
      MDB_Point* tmp = Pt1;
      Pt1 = Pt3;
      Pt3 = tmp;
    }
  };

  /*! \brief 5th order equidistant Lagrange edge */

  class MDB_EdgeP5: public MDB_Edge
  {
    MDB_Point *Pt1,*Pt2,*Pt3, *Pt4;//fifth order points
  public :
    virtual int getOrder() const {return 5;}
    virtual int getMshTag() const {return MSH_LIN_6;}
    virtual int getNbHighOrderPoints() const {return 4;}
    virtual MDB_Point *getHighOrderPoint (int i) {
      switch (i){
      case 0:
        return Pt1;
      case 1:
        return Pt2;
      case 2:
        return Pt3;
      case 3:
        return Pt4;
      default:
        std::cout << "Point " << i << "not defined for fifth order edges\n"; throw;
      }
    }
    MDB_EdgeP5(MDB_Point *A, MDB_Point *B1, MDB_Point *B2, MDB_Point *B3, MDB_Point *B4, MDB_Point *C)
      : MDB_Edge ( A, C ), Pt1(B1), Pt2(B2), Pt3(B3), Pt4(B4)
    {
    }
    virtual ~MDB_EdgeP5() {delete Pt1;delete Pt2;delete Pt3;delete Pt4;}

  protected:
    virtual void swapHONodes() 
    {
      MDB_Point* tmp = Pt1;
      Pt1 = Pt4;
      Pt4 = tmp;

      tmp = Pt2;
      Pt2 = Pt3;
      Pt3 = tmp;
    }  
  };

  class MDB_Face: public MDB_MeshEntity
  {
  protected:
    std::vector<MDB_Region *> _regions;
  public :
    virtual int        getNbEdges ()  const = 0;
    virtual int        getNbNodes ()  const = 0;
    virtual MDB_Edge * getEdge  (int) const = 0;
    virtual int getDim() const {return 2;}
    virtual int getMshTag() const = 0;

    inline MDB_Region* getRegion (int i) const
    {
      if (getNbRegions()==0)return 0;
      if (getNbRegions()==1 && i == 1)return 0;
      return _regions[i];
    }

    inline MDB_Region * opposite_region(MDB_Region *t)
    {
      if(t == getRegion(0))return getRegion(1);
      if(t == getRegion(1))return getRegion(0);
      throw;
    }

    inline int getNbRegions() const 
    {
      return (_regions.size());
    }
  
    inline void addregion(MDB_Region *t)
    {
      _regions.push_back(t);
    }

    inline void del( MDB_Region *t )
    {
      _regions.erase(std::remove_if(_regions.begin(),_regions.end() , std::bind2nd(std::equal_to<MDB_Region*>(), t)) , 
                     _regions.end());
    }

    inline MDB_Edge * commonedge(const MDB_Face *other) const
    {
      for (int i=0;i<getNbEdges ();i++)
        for (int j=0;j<other->getNbEdges ();j++)
          {
            MDB_Edge *e = getEdge(i);
            if (e == other->getEdge(j)) return  e;
          }
      throw;
    }

    inline MDB_Point *getNode(int i) const
    {
      if (!i)
        {
          int nbe = getNbEdges ();
          return getEdge(0)->commonvertex(getEdge(nbe-1));
        }
      return getEdge(i)->commonvertex(getEdge(i-1));        
    }

    inline void getNodes(MDB_Point **n) const
    {
      for (int i=0;i<getNbEdges ();i++)n[i] = getNode(i);
    }
    virtual MDB_Edge * find_edge(const MDB_Point*, const MDB_Point*) const = 0;

    MDB_Face() /*: r1(0), r2(0) */{};
    virtual MDB_Point *getHighOrderPoint (int i) {throw;}
    virtual int getNbHighOrderPoints() const {return 0;}
  
    virtual int align(MDB_Point* v0,MDB_Point* v1,MDB_Point* v2,MDB_Point* v3) = 0;

  protected:
  
    virtual int orientation(MDB_Point* v0,MDB_Point* v1,MDB_Point* v2,MDB_Point* v3) const;
  
  };

  /*! \brief Straight edged triangle and base class for higher order versions */ 

  class MDB_Triangle: public MDB_Face
  {
  public:
    MDB_Edge *e1,*e2,*e3;
    
    virtual int        getOrder()     const {return 1;}
    virtual int        getMshTag()    const {return MSH_TRI_3;}
    virtual int        getNbEdges ()  const {return 3;}
    virtual int        getNbNodes ()  const {return 3;}
    virtual int        align(MDB_Point*,MDB_Point*,MDB_Point*,MDB_Point*);
    virtual MDB_Edge * getEdge  (int i) const 
    {
      switch (i){
      case 0 : return e1;
      case 1 : return e2;
      case 2 : return e3;
      case 3 : return (MDB_Edge *)  NULL;
      }
      throw;
    }
  
    MDB_Triangle(MDB_Edge *A, MDB_Edge *B, MDB_Edge *C)
      : e1(A),e2(B),e3(C)
    {	
      e1->addface(this);
      e2->addface(this);
      e3->addface(this);
    }
    virtual ~MDB_Triangle() {}

    virtual MDB_Edge * find_edge(const MDB_Point*, const MDB_Point*) const;

    virtual int getNbHighOrderPoints() const { return 0;}

  };


  /*! \brief Generic order equidistant Lagrange triangle with complete interpolation */ 
  template <int order>
  class MDB_CompleteTriangle : public MDB_Triangle {

  public:
  
    MDB_CompleteTriangle(MDB_Edge *A, MDB_Edge *B, MDB_Edge *C,MDB_Point** P):
      MDB_Triangle(A,B,C),nbPoints(std::max(0,(order-2)*(order-1)/2)),point(NULL) {
      point = nbPoints > 0 ? new MDB_Point*[nbPoints] : NULL;
      for (int i=0;i<nbPoints;i++) point[i] = P[i];  
    }
    
    ~MDB_CompleteTriangle() {if (point) delete [] point;}
    
    virtual int getOrder()     const {return order;}
    int getMshTag() const {
      switch (nbPoints) {
      case 0: return MSH_TRI_3;
      case 3: return MSH_TRI_6;
      case 6: return MSH_TRI_9;
      case 7: return MSH_TRI_10;
      case 9: return MSH_TRI_12;
      case 12: return MSH_TRI_15;
        //      case 12: return MSH_TRI_15I;
      case 18: return MSH_TRI_21;
      default: return MSH_UNDEFINED_ELEM;
      }
      return MSH_UNDEFINED_ELEM;
    }
    int getNbHighOrderPoints() const {return nbPoints;}
    MDB_Point* getHighOrderPoint(int i) {return (i<nbPoints) ?  point[i] : NULL;}

  protected:
    
    int nbPoints;
    MDB_Point** point;


  protected:


    virtual int align(MDB_Point* v0,MDB_Point* v1,MDB_Point* v2,MDB_Point* v3) {


      if (v3 != NULL) {
        MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                    "Aligning a triangle with a quadrilateral is not allowed");
      }
    
      int rot = MDB_Triangle::align(v0,v1,v2,NULL);


      /* we only store the internal points (the higher order nodes on the edges are stored elsewhere)
         this means that this corresponds to a triangle of order = order - 3 */ 
      rotateHOPoints(rot,point,order-3);
    
      return rot;
    }

  
    /*! rotate a triangular set of points such that the node on position abs(rot) - 1
      is migrated to node 0; swap if rot < 0 */

    void rotateHOPoints(int rot,MDB_Point** op,int recursiveOrder) {
    
      if (recursiveOrder <= 0) return; // only one node needs no further reorientation 
      if (rot   == 1) return; // already aligned 
  
      int nbRecursivePoints = (recursiveOrder+1)*(recursiveOrder+2)/2;
      int nbEdgePoints  = (recursiveOrder - 1);
  
      MDB_Point** np = new MDB_Point*[nbRecursivePoints];
  
      // first rearrange principal vertices
  
      memcpy(np,op,nbRecursivePoints*sizeof(MDB_Point*));
  
      for (int i=0;i<3;i++) np[i] = op[(i + abs(rot) - 1)%3];
      if (rot < 0) std::swap(np[1],np[2]);
  
  
      if (recursiveOrder > 1)  {
    
        // then edge nodes
      
        if (rot > 0) {
          for (int i=0;i<3;i++) {
            MDB_Point** nb = np + 3 + i * nbEdgePoints;
            MDB_Point** ob = op + 3 + (i + rot - 1) % 3 * nbEdgePoints;
            for (int j=0;j<nbEdgePoints;j++) memcpy(nb,ob,nbEdgePoints * sizeof(MDB_Point*));
          }
        }
      
        if (rot < 0) {
          for (int i=0;i<3;i++) {
            MDB_Point** nb = np + 2 + ((abs(rot) - i + 1) % 3 + 1)  * nbEdgePoints;      
            MDB_Point** ob = op + 3 + i * nbEdgePoints;
            for (int j=0;j<nbEdgePoints;j++) *nb-- = *ob++;
          }
        }
      
        rotateHOPoints(rot,np + 3 * (nbEdgePoints+1) , recursiveOrder - 3);
      }
    
      memcpy(op,np,nbRecursivePoints * sizeof(MDB_Point*));
      delete [] np;
    } 

   
  };

  class MDB_Quad : public MDB_Face
  {
  public:
    MDB_Edge *e1,*e2,*e3,*e4;
  
    
    virtual int        getOrder()     const {return 2;}
    virtual int        getMshTag()    const {return MSH_QUA_4;}
    virtual int        getNbEdges ()  const {return 4;}
    virtual int        getNbNodes ()  const {return 4;}
    virtual int        align(MDB_Point*,MDB_Point*,MDB_Point*,MDB_Point*);
    virtual MDB_Edge * getEdge  (int i) const 
    {
      switch (i){
      case 0 : return e1;
      case 1 : return e2;
      case 2 : return e3;
      case 3 : return e4;
      }
      throw;
    }
  
    MDB_Quad(MDB_Edge *A, MDB_Edge *B, MDB_Edge *C, MDB_Edge *D)
      : e1(A),e2(B),e3(C),e4(D)
    {	
      e1->addface(this);
      e2->addface(this);
      e3->addface(this);
      e4->addface(this);
    }
    virtual ~MDB_Quad() {}
    virtual MDB_Edge * find_edge(const MDB_Point*, const MDB_Point*) const;
  };

  /*! \brief Generic order equidistant Lagrange quad with complete interpolation */ 
  template <int order>
  class MDB_CompleteQuad : public MDB_Quad {

  public:
  
    MDB_CompleteQuad(MDB_Edge *A, MDB_Edge *B, MDB_Edge *C,MDB_Edge *D, MDB_Point** P):
      MDB_Quad(A,B,C,D),nbPoints(std::max(0,(order-1)*(order-1))),point(NULL) {
      point = nbPoints > 0 ? new MDB_Point*[nbPoints] : NULL;
      for (int i=0;i<nbPoints;i++) point[i] = P[i];  
    }
    
    ~MDB_CompleteQuad() {if (point) delete [] point;}
    
    virtual int getOrder()     const {return order;}
    int getNbHighOrderPoints() const {return nbPoints;}
    MDB_Point* getHighOrderPoint(int i) {return (i<nbPoints) ?  point[i] : NULL;}

  protected:
    
    int nbPoints;
    MDB_Point** point;

  };


  class MDB_Region: public MDB_MeshEntity
  {
  public:
    virtual MDB_Face *getFace (int i) const = 0;
    virtual MDB_Edge *getEdge (int i) const = 0;
    virtual MDB_Point *getVertex (int i) const = 0;
    virtual MDB_Point *getHighOrderPoint (int i) {throw;}
    virtual int getNbHighOrderPoints() const {return 0;}
    virtual void      getNodes(MDB_Point **n) const = 0;
    virtual int       getNbFace () const = 0;
    virtual int       getNbEdge () const = 0;
    virtual int       getNbVertex () const = 0;
    virtual int       getFaceDir (int) const = 0;
    virtual int       getFaceOrientation(int) const {throw;}
  
  };

  class MDB_Hex: public MDB_Region
  {
  public:
    MDB_Quad  *f1,*f2,*f3,*f4,*f5,*f6;
  
    virtual int getOrder() const {return 3;}
    virtual int getDim()   const {return 3;}
    virtual int getMshTag() const {return MSH_HEX_8;}

    virtual int       getNbFace () const {return 6;}
    virtual int       getNbEdge () const {return 12;}
    virtual int       getNbVertex () const{return 8;}
    virtual MDB_Face *getFace (int i) const
    {
      switch(i){
      case 0:return f1;
      case 1:return f2;
      case 2:return f3;
      case 3:return f4;
      case 4:return f5;
      case 5:return f6;
      }
      throw;
    }
    virtual MDB_Edge *getEdge (int n) const
    {
      switch (n){
      case 0  : return f1->commonedge ( f2 );
      case 1  : return f1->commonedge ( f3 );
      case 2  : return f1->commonedge ( f4 );
      case 3  : return f1->commonedge ( f5 );
      case 4  : return f2->commonedge ( f5 );
      case 5  : return f2->commonedge ( f3 );    
      case 6  : return f3->commonedge ( f4 );    
      case 7  : return f4->commonedge ( f5 );    
      case 8  : return f2->commonedge ( f6 );    
      case 9  : return f3->commonedge ( f6 );    
      case 10 : return f4->commonedge ( f6 );    
      case 11 : return f5->commonedge ( f6 );    
      }
      throw;    
    }

    virtual MDB_Point *getVertex (int n) const /*SL OK*/
    {
      MDB_Edge *e1=0,*e2=0; 
      switch (n){
      case 0 : e1 = getEdge(0);e2 = getEdge(3);break;
      case 1 : e1 = getEdge(0);e2 = getEdge(1);break;
      case 2 : e1 = getEdge(1);e2 = getEdge(2);break;
      case 3 : e1 = getEdge(2);e2 = getEdge(3);break;
      case 4 : e1 = getEdge(4);e2 = getEdge(8);break;
      case 5 : e1 = getEdge(5);e2 = getEdge(9);break;
      case 6 : e1 = getEdge(6);e2 = getEdge(10);break;
      case 7 : e1 = getEdge(7);e2 = getEdge(11);break;
      default: throw;
      }

      return e1->commonvertex ( e2 );
    }

    virtual int getFaceDir (int n) const;


    virtual void getNodes(MDB_Point **n) const
    {
      MDB_Point *o[4];
      f1->getNodes(n);
      f2->getNodes(o);
      for(int i=0;i<4;i++) n[i+4]=o[i];
    }
  
    MDB_Hex(MDB_Quad *A, MDB_Quad *B, MDB_Quad *C, MDB_Quad *D, MDB_Quad *E,  MDB_Quad *F)
      : f1(A),f2(B),f3(C),f4(D),f5(E),f6(F)
    { 
      f1->addregion(this);
      f2->addregion(this);
      f3->addregion(this);
      f4->addregion(this);
      f5->addregion(this);
      f6->addregion(this);
    }
    virtual ~MDB_Hex() {}

  };

  class MDB_Tet: public MDB_Region
  {
  public:
    MDB_Triangle  *f1,*f2,*f3,*f4;
  
    virtual int getOrder() const {return 1;}
    virtual int getDim()   const {return 3;}
    virtual int getMshTag() const {return MSH_TET_4;}

    virtual int       getNbFace () const {return 4;}
    virtual int       getNbEdge () const {return 6;}
    virtual int       getNbVertex () const{return 4;}
    virtual MDB_Face *getFace (int i) const
    {
      switch(i){
      case 0:return f1;
      case 1:return f2;
      case 2:return f3;
      case 3:return f4;
        // NP default NULL added to avoid call to throw 
        //    in classify_unclassified_entities
      default: return (MDB_Face *)  NULL;
      }
      throw;
    }
    virtual MDB_Edge *getEdge (int n) const
    {
      switch (n){
      case 0 : return f1->commonedge ( f2 );
      case 1 : return f1->commonedge ( f3 );
      case 2 : return f1->commonedge ( f4 );
      case 3 : return f2->commonedge ( f4 );
      case 4 : return f2->commonedge ( f3 );
      case 5 : return f3->commonedge ( f4 );    
      }
      throw;    
    }

    virtual MDB_Point *getVertex (int n) const
    {
      MDB_Edge *e1=0,*e2=0; 
      switch (n){
      case 0 : e1 = getEdge(0);e2 = getEdge(2);break;
      case 1 : e1 = getEdge(0);e2 = getEdge(1);break;
      case 2 : e1 = getEdge(1);e2 = getEdge(2);break;
      case 3 : e1 = getEdge(3);e2 = getEdge(5);break;
      default: throw;
      }
      return e1->commonvertex ( e2 );
    }

    // return 1 if the face direction points to outside of the tet
    // return 0                                 inside
    virtual int getFaceDir (int n) const;
    // return number of times we need to rotate the element to make
    // the principal node coincide with that of the tetrahedron (before normal inversion)
    virtual int getFaceOrientation (int n) const;

    virtual void getNodes(MDB_Point **n) const
    {
      MDB_Point *o[3];
      f1->getNodes(n);
      f2->getNodes(o);	  
      n[3] = 0; //for stupid gcc warning
      if(o[0] != n[0] && o[0] != n[1] &&o[0] != n[2])n[3] = o[0];
      if(o[1] != n[0] && o[1] != n[1] &&o[1] != n[2])n[3] = o[1];
      if(o[2] != n[0] && o[2] != n[1] &&o[2] != n[2])n[3] = o[2];
    }
  
    MDB_Tet(MDB_Triangle *A, MDB_Triangle *B, MDB_Triangle *C, MDB_Triangle *D)
      : f1(A),f2(B),f3(C),f4(D)
    {	
      f1->addregion(this);
      f2->addregion(this);
      f3->addregion(this);
      f4->addregion(this);
    }
    virtual ~MDB_Tet() {}

  protected:

    static int face_vtx[4][3];
  
  };

  // 2 be done
  class MDB_Prism: public MDB_Region
  {
  public:
    MDB_Triangle  *f1,*f2;
    MDB_Quad      *f3,*f4,*f5;
  
    virtual int getOrder() const {return 2;}
    virtual int getDim()   const {return 3;} /*SL OK*/
    virtual int getMshTag() const {return MSH_PRI_6;}

    virtual int       getNbFace () const {return 5;} /*SL OK*/ 
    virtual int       getNbEdge () const {return 9;} /*SL OK*/
    virtual int       getNbVertex () const{return 6;} /*SL OK*/
    virtual MDB_Face *getFace (int i) const /*SL OK*/
    {
      switch(i){
      case 0:return f1;
      case 1:return f2;
      case 2:return f3;
      case 3:return f4;
      case 4:return f5;
      case 5:return (MDB_Face *) NULL;
      }
      throw;
    }
    virtual MDB_Edge *getEdge (int n) const /*SL OK*/
    {
      switch (n){
      case 0  : return f1->commonedge ( f3 );
      case 1  : return f1->commonedge ( f5 );
      case 2  : return f1->commonedge ( f4 );
      case 3  : return f2->commonedge ( f3 );
      case 4  : return f2->commonedge ( f5 );
      case 5  : return f2->commonedge ( f4 );    
      case 6  : return f3->commonedge ( f5 );    
      case 7  : return f3->commonedge ( f4 );    
      case 8  : return f4->commonedge ( f5 );       
      }
      throw;    
    }
  
    virtual MDB_Point *getVertex (int n) const /*SL OK*/
    {
      MDB_Edge *e1=0,*e2=0; 
      switch (n){
      case 0 : e1 = getEdge(0);e2 = getEdge(1);break;
      case 1 : e1 = getEdge(0);e2 = getEdge(2);break;
      case 2 : e1 = getEdge(1);e2 = getEdge(2);break;
      case 3 : e1 = getEdge(4);e2 = getEdge(3);break;
      case 4 : e1 = getEdge(3);e2 = getEdge(5);break;
      case 5 : e1 = getEdge(4);e2 = getEdge(5);break;
      default: throw;
      }
      return e1->commonvertex ( e2 );
    }

    virtual int getFaceDir (int n) const
    { std::cerr << "not yet implemented\n"; abort(); return 0;}


    virtual void getNodes(MDB_Point **n) const
    {
      std::cerr << "not yet implemented\n"; abort();
      /*    MDB_Point *o[3];
            f1->getNodes(n);
            f2->getNodes(o);    
            n[3] = 0; //for stupid gcc warning
            if(o[0] != n[0] && o[0] != n[1] &&o[0] != n[2])n[3] = o[0];
            if(o[1] != n[0] && o[1] != n[1] &&o[1] != n[2])n[3] = o[1];
            if(o[2] != n[0] && o[2] != n[1] &&o[2] != n[2])n[3] = o[2];*/
    }
  
    MDB_Prism(MDB_Triangle *A, MDB_Triangle *B, MDB_Quad *C, 
              MDB_Quad *D, MDB_Quad *E) /*SL OK*/
      : f1(A),f2(B),f3(C),f4(D),f5(E)
    { 
      f1->addregion(this);
      f2->addregion(this);
      f3->addregion(this);
      f4->addregion(this);
      f5->addregion(this);
    }
    virtual ~MDB_Prism() {} /*SL OK*/
  };

  /*! \brief Generic order equidistant Lagrange tetrahedron */ 

  template <int order>
  class MDB_CompleteTet : public MDB_Tet {

  public:

    MDB_CompleteTet(MDB_Triangle *A, MDB_Triangle *B,
                    MDB_Triangle *C, MDB_Triangle *D,MDB_Point** pp):
    
      MDB_Tet(A,B,C,D),nbPoints(std::max(0,(order-3)*(order-2)*(order-1)/6))
    {
      point = nbPoints ? new MDB_Point*[nbPoints] : NULL;
      for (size_t i=0;i<(unsigned int)nbPoints;i++) point[i] = pp[i];
    }
    
    int getNbHighOrderPoints() const {return nbPoints;}
    MDB_Point* getHighOrderPoint(int i) {return (i<nbPoints) ? point[i] : NULL;}
    virtual int getOrder() {return order;}
    virtual int getMshTag() const {
      switch (nbPoints) {
      case 0: return MSH_TET_4;
      case 6: return MSH_TET_10;
      case 16: return MSH_TET_20;
      case 31: return MSH_TET_35;
      case 52: return MSH_TET_56;
      case 30: return MSH_TET_34;
      default: return MSH_UNDEFINED_ELEM;
      }
      return MSH_UNDEFINED_ELEM;
    }
    
  protected:
    int nbPoints;
    MDB_Point** point;
  };

  class PointLessThan
  {
  public:
    bool operator()(const MDB_Point* ent1, const MDB_Point* ent2) const
    {
      return *ent1 < *ent2;
    }
  };


  class MDB_Mesh : public MDB_MiniMesh
  {        
    // no copies allowed
    MDB_Mesh(const MDB_Mesh &other);

  private:
    bool parametric;

  public:

    // this is an add on to the classical interface
    // one can build lists of geometrical "features" that
    // encompass more than one geometrical entities
    std::multimap<int, pGEntity > geomFeatures_Tags;
    // sometimes, features have names
    std::map<std::string, int >         geomFeatures_Names;

    bool shrinked;

    // constructor of an empty mesh
    MDB_Mesh(int _MAXX = 1); 
    // the destructor
    ~MDB_Mesh();
    // initialize first vertex id and increments
    void initializeIdData();
    // load the mini mesh into the binary file f, mesh is still empty
    bool isParametric() const { return parametric; }
    void load   ( FILE *f );
    // expand the mini mesh in the bi directional data structure
    // then, delete the mini mesh. Mesh adaptation can be performed
    // only when the mesh is expanded 
    void expand ( );
    // shrink the mini mesh from the bi-directional data structure
    // then, delete the bi-directional data structure. The mesh
    // is "idle", i.e. is shrinked in memory, allowing the solver 
    // to take the whole memory space
    void shrink ( );
    bool isShrinked() const { return shrinked; };
    // flush the mini mesh
    void flush ( FILE *f ) const;
    // data about vertex ids (globally unique ids in parallel)
    int maxId;       // max vertex id
    int idIncrement; // id increment size (=nbProcs)
    // the geometric model on which the mesh entities are classified
    pGModel model; 
    // bounding boxes
    double Min[3],Max[3],LC;
    // the datas
    int nbPoints, nbEdgePoints, nbEdges, nbTriangles, nbQuads, nbTets, nbHexes, nbPrisms;
    std::set<MDB_Point*,PointLessThan>     points; 
    MDB_ListE edges; 
    MDB_ListF triangles; 
    MDB_ListQ quads; 
    MDB_ListT tets; 
    MDB_ListH hexes; 
    MDB_ListPr prisms; 
    // Points operations
    MDB_Point * add_point(int num , double x, double y,double z, pGEntity=0);
    MDB_PointParam * add_pointParam(int num , double x, double y, double z, 
                                    double u, double v, pGEntity=0);
    void del_point(MDB_Point *p);  
    MDB_Point *find_point(int num);
    void deleteTheHighOrderPoint(std::set < MDB_Point *, PointLessThan >::iterator it);
    // Edges operations
    /*! straight edge from point id's */
    MDB_Edge  * add_edge(int p1, int p2, pGEntity=0);
    //  MDB_Edge  * add_edge(int p1, int p2, int p3, pGEntity=0);
    /*! straight edge from points */ 
    MDB_Edge  * add_edge(MDB_Point *p1, MDB_Point *p2, pGEntity=0);
    //  MDB_Edge  * add_edge(MDB_Point *p1, MDB_Point *p2, MDB_Point *p3, pGEntity=0);
    /*! straight/curved edge from point id's */
    MDB_Edge  * add_edge(int numberOfPoints, pGEntity pg, int p1, ...);
    /*! straight/curved edge from points */
    MDB_Edge  * add_edge(int numberOfPoints, pGEntity pg, MDB_Point *firstPoint, ...);
    /*! straight/curved edge from points (gmsh order) */
    MDB_Edge  * add_edge(MDB_Point* p1,MDB_Point*p2,pGEntity,int o,MDB_Point** iP);
    void del_edge(MDB_Edge *e);
    MDB_Edge  *find_edge(int p1, int p2);
    MDB_Edge  *find_edge(MDB_Point *p1, MDB_Point *p2)const;
    MDB_Edge  *find_clone_edge(MDB_Edge *edge)const;
    // Faces operations
    void del_face(MDB_Face *f);
    // Triangles operations
    /*! straight triangle from vertex id's */ 
    MDB_Triangle *add_triangle(int p1, int p2, int p3, pGEntity=0); 
    /*! straight triangle from vertices */ 
    MDB_Triangle *add_triangle(MDB_Point* v0, MDB_Point* v1, MDB_Point* v2, pGEntity=0); 
    /* curvilinear triangle */
    //  MDB_Triangle *add_triangle(int p1, int p2, int p3, int p4, int p5, int p6, pGEntity=0);
    /*! straight/serendipity triangle from edges */
    MDB_Triangle *add_triangle(MDB_Edge *e1,MDB_Edge *e2,MDB_Edge *e3, pGEntity=0);
    /*! straight/complete triangle triangle from edges and internal points */
    MDB_Triangle *add_triangle(MDB_Edge *e1,MDB_Edge *e2,MDB_Edge *e3, pGEntity,int order,
                               bool serendip=true,MDB_Point** iP=NULL);
    /*! straight/complete triangle triangle from point id's*/
    //MDB_Triangle *add_triangle(pGEntity,int order,bool serendip=true,int* pp=NULL);
    /*! straight/serendipity triangle from point id's */ 
    MDB_Triangle *add_triangle(int order, bool complete, pGEntity pg, int p1, int p2, ...);
    void del_triangle(MDB_Triangle *t);
    MDB_Triangle  *find_triangle(MDB_Edge *e1,MDB_Edge *e2,MDB_Edge *e3);
    // Quad operations
    /*! straight quadrangle from vertex id's */ 
    MDB_Quad *add_quad(int p1, int p2, int p3, int p4, pGEntity=0); 
    /*! straight quadrangle from vertices */ 
    MDB_Quad *add_quad(MDB_Point* v0, MDB_Point* v1, MDB_Point* v2, MDB_Point* v3, pGEntity=0); 
    MDB_Quad *add_quad(MDB_Edge *e1,MDB_Edge *e2,MDB_Edge *e3,MDB_Edge *e4, pGEntity=0);
    /*! straight/serendipity quad from point id's */ 
    MDB_Quad *add_quad(int order, bool complete, pGEntity pg, int p1, int p2, int p3, ...);
    /*! straight/complete quad from edges and internal points */
    MDB_Quad *add_quad(MDB_Edge *e1,MDB_Edge *e2,MDB_Edge *e3,MDB_Edge *e4, pGEntity,int order,
                               bool serendip=true,MDB_Point** iP=NULL);
    void del_quad(MDB_Quad *q);
    MDB_Quad  *find_quad(MDB_Edge *e1,MDB_Edge *e2,MDB_Edge *e3,MDB_Edge *e4);
    // Tets operations
    /*! straight tet from point id's */
    MDB_Tet *add_tet(int p1, int p2, int p3,int p4, pGEntity=0);
    /*! straight/serendipity/complete triangle from point id's (in gmsh order) */
    MDB_Tet *add_tet(pGEntity,int order,bool serendip,int* points);
    /*! straight tet from points */ 
    MDB_Tet *add_tet(MDB_Point *p1, MDB_Point *p2,MDB_Point *p3, MDB_Point *p4, pGEntity=0);
    /*! straight/serendipity tetrahedron from faces */
    MDB_Tet *add_tet(MDB_Triangle  *,MDB_Triangle  *,MDB_Triangle  *,MDB_Triangle  *, pGEntity=0);
    /*! straight/complete tetrahedron from faces and points */
    MDB_Tet *add_tet(MDB_Triangle  *,MDB_Triangle  *,MDB_Triangle  *,MDB_Triangle  *, pGEntity,
                     int order,bool serendip=true,MDB_Point** iP=NULL);
  
    /*! straight/complete tet from faces and internal point id's */
    MDB_Tet *find_tet(MDB_Triangle * t1, MDB_Triangle * t2, MDB_Triangle * t3,  MDB_Triangle * t4);
    void del_tet(MDB_Tet *);
    // Hexes operations
    MDB_Hex *add_hex(int p1, int p2, int p3,int p4,int p5, int p6, int p7, int p8, pGEntity=0);
    MDB_Hex *add_hex(MDB_Point *p1, MDB_Point *p2,MDB_Point *p3, MDB_Point *p4,
                     MDB_Point *p5, MDB_Point *p6,MDB_Point *p7, MDB_Point *p8, pGEntity=0);
    MDB_Hex *add_hex(MDB_Quad  *,MDB_Quad  *,MDB_Quad  *,MDB_Quad  *,MDB_Quad  *,MDB_Quad  *, pGEntity=0);
    MDB_Hex *find_hex(MDB_Quad  *,MDB_Quad  *,MDB_Quad  *,MDB_Quad  *,MDB_Quad  *,MDB_Quad  *);
    void del_hex(MDB_Hex *);
    // Prisms operations
    MDB_Prism *add_prism(int p1, int p2, int p3,int p4,int p5, int p6, pGEntity=0);
    MDB_Prism *add_prism(MDB_Point *p1, MDB_Point *p2,MDB_Point *p3, 
                         MDB_Point *p4, MDB_Point *p5, MDB_Point *p6, pGEntity=0);
    MDB_Prism *add_prism(MDB_Triangle  *,MDB_Triangle  *,MDB_Quad  *,MDB_Quad  *,MDB_Quad  *, pGEntity=0);
    MDB_Prism *find_prism(MDB_Triangle  *,MDB_Triangle  *,MDB_Quad  *,MDB_Quad  *,MDB_Quad  *);
    void del_prism(MDB_Prism *);
    // classify entities that are not classified
    void classify_unclassified_entities();
    void destroyStandAloneEntities();

// #ifdef PARALLEL
//     void bdryLinkRemotePoint();
// #endif  
    
    // check if there is no identical geometrical tags between regions and faces
    /*   void checkGeomTagsCoherence () const; */
  };

  typedef std::set<MDB_Point*,PointLessThan> MDB_SetV;

  class MDB_RIter : public MDB_Iter < MDB_ListR , MDB_Region , pGEntity>
  {
  public:
    MDB_RIter (MDB_ListR*_l): MDB_Iter< MDB_ListR , MDB_Region , pGEntity > (_l){}
    MDB_RIter (MDB_ListR*_l,pGEntity _g): MDB_Iter< MDB_ListR , MDB_Region , pGEntity > (_l,_g){}
  };
  class MDB_TIter : public MDB_Iter < MDB_ListT , MDB_Tet , pGEntity>
  {
  public:
    MDB_TIter (MDB_ListT*_l): MDB_Iter< MDB_ListT , MDB_Tet , pGEntity > (_l){}
    MDB_TIter (MDB_ListT*_l,pGEntity _g): MDB_Iter< MDB_ListT , MDB_Tet , pGEntity > (_l,_g){}
  };
  class MDB_HIter : public MDB_Iter < MDB_ListH , MDB_Hex , pGEntity>
  {
  public:
    MDB_HIter (MDB_ListH*_l): MDB_Iter< MDB_ListH , MDB_Hex , pGEntity > (_l){}
    MDB_HIter (MDB_ListH*_l,pGEntity _g): MDB_Iter< MDB_ListH , MDB_Hex , pGEntity > (_l,_g){}
  };
  class MDB_PIter : public MDB_Iter < MDB_ListPr , MDB_Prism , pGEntity>
  {
  public:
    MDB_PIter (MDB_ListPr*_l): MDB_Iter< MDB_ListPr , MDB_Prism , pGEntity > (_l){}
    MDB_PIter (MDB_ListPr*_l,pGEntity _g): MDB_Iter< MDB_ListPr , MDB_Prism , pGEntity > (_l,_g){}
  };
  class MDB_EIter : public MDB_Iter < MDB_ListE , MDB_Edge , pGEntity >
  {
  public:
    MDB_EIter (MDB_ListE*_l): MDB_Iter< MDB_ListE , MDB_Edge, pGEntity > (_l){}
    MDB_EIter (MDB_ListE*_l,pGEntity _g,int _c): MDB_Iter< MDB_ListE , MDB_Edge, pGEntity > (_l,_g,_c){}
  };
  class MDB_FIter : public MDB_Iter < MDB_ListF , MDB_Triangle , pGEntity >
  {
  public:
    MDB_FIter (MDB_ListF*_l): MDB_Iter <MDB_ListF , MDB_Triangle ,pGEntity> (_l){}
    MDB_FIter (MDB_ListF*_l,pGEntity _g,int _c): MDB_Iter <MDB_ListF , MDB_Triangle ,pGEntity> (_l,_g,_c){}
  };
  class MDB_QIter : public MDB_Iter < MDB_ListQ , MDB_Quad , pGEntity>
  {
  public:
    MDB_QIter (MDB_ListQ*_l): MDB_Iter <MDB_ListQ , MDB_Quad , pGEntity> (_l){}
    MDB_QIter (MDB_ListQ*_l,pGEntity _g,int _c): MDB_Iter <MDB_ListQ , MDB_Quad ,pGEntity> (_l,_g,_c){}
  };
  class MDB_VIter : public MDB_Iter < MDB_SetV , MDB_Point ,pGEntity >
  {
  public:
    MDB_VIter (MDB_SetV *_l): MDB_Iter < MDB_SetV , MDB_Point , pGEntity > (_l){}
    MDB_VIter (MDB_SetV *_l,pGEntity _g,int _c): MDB_Iter < MDB_SetV , MDB_Point , pGEntity > (_l,_g,_c){}
  };

  class MDB_FaceIter
  {
  public:
    MDB_FIter itf;
    MDB_QIter itq;

    MDB_Triangle *t;
    MDB_Quad *q;
  public:
    MDB_FaceIter ( MDB_ListF* _t , MDB_ListQ* _q):itf(_t),itq(_q),t(0),q(0){}
    MDB_FaceIter ( MDB_ListF* _t , MDB_ListQ* _q, pGEntity _g,int _c):itf(_t,_g,_c),itq(_q,_g,_c),t(0),q(0){}
    inline MDB_Face * next ()
    {
      t = itf.next();
      if (t) return (MDB_Face*)t;
      q = itq.next();
      return (MDB_Face*)q;
    }
    inline void reset()
    {
      itf.reset();
      itq.reset();
    }
    inline void cleanup()
    {
      itf.cleanup();
      itq.cleanup();
    }

  };

  class MDB_RegionIter
  {
  public:
    MDB_TIter itt;
    MDB_HIter ith;
    MDB_PIter itp;

    MDB_Tet *t;
    MDB_Hex *h;
    MDB_Prism *p;
  public:
    MDB_RegionIter ( MDB_ListT* _t , MDB_ListH* _h, MDB_ListPr* _p ):itt(_t),ith(_h),itp(_p),t(0),h(0),p(0){}
    MDB_RegionIter ( MDB_ListT* _t , MDB_ListH* _h, MDB_ListPr* _p, pGEntity _g):itt(_t,_g),ith(_h,_g),itp(_p,_g),t(0),h(0),p(0){}
    inline MDB_Region * next ()
    {
      t = itt.next();
      if (t) return (MDB_Region*)t;
      h = ith.next();
      if (h) return (MDB_Region*)h;
      p = itp.next();
      return (MDB_Region*)p;
    }
    inline void reset()
    {
      itt.reset();
      ith.reset();
      itp.reset();
    }
    inline void cleanup()
    {
      itt.cleanup();
      ith.cleanup();
      itp.cleanup();
    }

  };

} // End of namespace MAd

#endif






