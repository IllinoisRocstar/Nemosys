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
// Authors: Arnaud Francois, Gaetan Compere, Jean-Francois Remacle
// -------------------------------------------------------------------

#ifndef _H_EDGESWAPCONFIG
#define _H_EDGESWAPCONFIG

namespace MAd
{
  // -------------------------------------------------------------------
  class EdgeSwapTemplate  // 3D case
  {
  public:
    EdgeSwapTemplate() {}
    virtual ~EdgeSwapTemplate() {}
    virtual int getConfig(){ return 0; }       // nb of faces attached to edge
    virtual int nb_triangulations() = 0;       // nb of possible triangulations
    virtual int nb_triangles() = 0;            // nb of different triangles
    virtual int nb_tri_triangulation() = 0;    // nb of triangles in a triangulation
    virtual const int* triangle( int i ) = 0;  // return the triangle i 
    virtual const int* triangulation( int i ) = 0;  // return the triangulation i
  };

  // -------------------------------------------------------------------
  // Null edge swap template
  class EdgeSwap0 : public EdgeSwapTemplate
  {
  public:
    EdgeSwap0() { }
    int getConfig(){ return 0; }
    int nb_triangles(){ return -1; }
    int nb_triangulations(){ return -1; }
    int nb_tri_triangulation(){ return -1; }
    const int* triangle( int i) { return 0; }
    const int* triangulation( int i) { return 0; }
  };

  // -------------------------------------------------------------------
  // Edge swap template with 3 faces connected to the edge
  class EdgeSwap3 : public EdgeSwapTemplate
  {
  public:
    EdgeSwap3() { }
    int getConfig(){ return 3; }
    int nb_triangles(){ return 1; }
    int nb_triangulations(){ return 1; }
    int nb_tri_triangulation(){ return 1; }
    const int* triangle( int i) { return triangles[i]; }
    const int* triangulation( int i) { return triangulations[i]; }
  private:
    static const int triangles[1][3];
    static const int triangulations[1][1];
  };

  // -------------------------------------------------------------------
  // Edge swap template with 4 faces connected to the edge
  class EdgeSwap4 : public EdgeSwapTemplate
  {
  public:
    EdgeSwap4() { }
    int getConfig(){ return 4; }
    int nb_triangles(){ return 4; }
    int nb_triangulations(){ return 2; }
    int nb_tri_triangulation(){ return 2; }
    const int* triangle( int i) { return triangles[i]; }
    const int* triangulation( int i) { return triangulations[i]; }
  private:
    static const int triangles[4][3];
    static const int triangulations[2][2];
  };

  // -------------------------------------------------------------------
  // Edge swap template with 5 faces connected to the edge
  class EdgeSwap5 : public EdgeSwapTemplate
  {
  public:
    EdgeSwap5() { }
    int getConfig(){ return 5; }
    int nb_triangles(){ return 10; }
    int nb_triangulations(){ return 5; }
    int nb_tri_triangulation(){ return 3; }
    const int* triangle( int i) { return triangles[i]; }
    const int* triangulation( int i) { return triangulations[i]; }
  private:
    static const int triangles[10][3];
    static const int triangulations[5][3];
  };

  // -------------------------------------------------------------------
  // Edge swap template with 6 faces connected to the edge
  class EdgeSwap6 : public EdgeSwapTemplate
  {
  public:
    EdgeSwap6() { }
    int getConfig(){ return 6; }
    int nb_triangles(){ return 20; }
    int nb_triangulations(){ return 14; }
    int nb_tri_triangulation(){ return 4; }
    const int* triangle( int i) { return triangles[i]; }
    const int* triangulation( int i) { return triangulations[i]; }
  protected:
    static const int triangles[20][3];
    static const int triangulations[14][4];
  };

  // -------------------------------------------------------------------
  // Edge swap template with 7 faces connected to the edge
  class EdgeSwap7 : public EdgeSwapTemplate
  {
  public:
    EdgeSwap7() { }
    int getConfig(){ return 7; }
    int nb_triangles(){ return 35; }
    int nb_triangulations(){ return 42; }
    int nb_tri_triangulation(){ return 5; }
    const int* triangle( int i) { return triangles[i]; }
    const int* triangulation( int i) { return triangulations[i]; }
  protected:
    static const int triangles[35][3];
    static const int triangulations[42][5];
  };

  // -------------------------------------------------------------------
  // Interface
  class EdgeSwapConfiguration
  {
  public:

    EdgeSwapConfiguration( ) { set(0); }
    EdgeSwapConfiguration( int n ) { set( n ); }
    EdgeSwapConfiguration( const EdgeSwapConfiguration &x ) { set( x.get() ); }
    ~EdgeSwapConfiguration() { }

    // select edge swap template with n faces connected to the edge
    void set( int n );

    // return the nodes of triangle i for the selected configuration
    const int* triangle( int i ) const {  return (*c).triangle(i); }

    // return the node j of triangle i for the selected configuration
    int triangle( int i, int j ) const {  return (*c).triangle(i)[j]; }

    // return the triangle j of triangulation i
    int triangulation( int i, int j ) const {  return (*c).triangulation(i)[j]; }

    int nb_triangles() const { return c->nb_triangles(); }
    int nb_triangulations() const { return c->nb_triangulations(); }
    int nb_tri_triangulation() const { return c->nb_tri_triangulation(); }

    int get() const { return c->getConfig(); }

  private:
    EdgeSwap0         cNull;
    EdgeSwap3         c3;
    EdgeSwap4         c4;
    EdgeSwap5         c5;
    EdgeSwap6         c6;
    EdgeSwap7         c7;
    EdgeSwapTemplate * c;
  };

  // -------------------------------------------------------------------
}
#endif
