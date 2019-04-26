// -------------------------------------------------------------------
// MAdLib - Copyright (C) 2008-2009 Universite catholique de Louvain
//
// See the Copyright.txt and License.txt files for license information. 
// You should have received a copy of these files along with MAdLib. 
// If not, see <http://www.madlib.be/license/>
//
// Please report all bugs and problems to <contrib@madlib.be>
//
// Authors: Jean-Francois Remacle, Gaetan Compere
// -------------------------------------------------------------------

// from Common/
#include "MathUtils.h"
#include "MAdMessage.h"
#include "MAdDefines.h"

#include "ModelInterface.h"
#ifdef _HAVE_GMSH_
#include "GmshModel.h"
#else
#include "NullModel.h"
#endif
#include "GM_Iterators.h"
#include "PGList.h"

#include <string.h>

namespace MAd {

  // -------------------------------------------------------------------
  void GM_create(pGModel* model, std::string name)
  { 
    if (*model) delete (*model);
#ifdef _HAVE_GMSH_
    *model = new GmshModel(name);
#else
    *model = new NullModel(name);
#endif
  }

  // -------------------------------------------------------------------
  void GM_delete(pGModel model)
  { 
    if (model) { delete model; model=NULL; }
  }
  
  // -------------------------------------------------------------------
  enum GeoFileFormat {
    MAD_FORMAT_MSH,
    MAD_FORMAT_GEO,
    MAD_FORMAT_STEP,
    MAD_FORMAT_BREP,
    MAD_FORMAT_IGES,
    MAD_FORMAT_UNKNOWN
  };

  // -------------------------------------------------------------------
  std::vector<std::string> SplitFileName(std::string fileName)
  {
    // returns [path, baseName, extension]
    unsigned int idot = fileName.find_last_of('.');
    unsigned int islash = fileName.find_last_of("/\\");
    if(idot == std::string::npos) idot = -1;
    if(islash == std::string::npos) islash = -1;
    std::vector<std::string> s(3);
    if(idot > 0)
      s[2] = fileName.substr(idot);
    if(islash > 0)
      s[0] = fileName.substr(0, islash + 1);
    s[1] = fileName.substr(s[0].size(), fileName.size() - s[0].size() - s[2].size());
    return s;
  }

  // -------------------------------------------------------------------
  GeoFileFormat guessFormatFromExtension(const std::string fileName)
  {
    std::string ext = SplitFileName(fileName)[2];
    if     ( !strcmp(ext.c_str(),".msh" ) )  return MAD_FORMAT_MSH;
    if     ( !strcmp(ext.c_str(),".geo" ) )  return MAD_FORMAT_GEO;
    if     ( !strcmp(ext.c_str(),".stp" ) )  return MAD_FORMAT_STEP;
    if     ( !strcmp(ext.c_str(),".step") )  return MAD_FORMAT_STEP;
    if     ( !strcmp(ext.c_str(),".brep") )  return MAD_FORMAT_BREP;
    if     ( !strcmp(ext.c_str(),".iges") )  return MAD_FORMAT_IGES;
    return MAD_FORMAT_UNKNOWN;
  }

  // -------------------------------------------------------------------
  int GM_read(pGModel model, const std::string fileName)
  {
    GeoFileFormat format = guessFormatFromExtension(fileName);
    if ( format == MAD_FORMAT_MSH  ) return GM_readFromMSH (model,fileName);
    if ( format == MAD_FORMAT_GEO  ) return GM_readFromGEO (model,fileName);
    if ( format == MAD_FORMAT_STEP ) return GM_readFromSTEP(model,fileName);
    if ( format == MAD_FORMAT_BREP ) return GM_readFromBREP(model,fileName);
    if ( format == MAD_FORMAT_IGES ) return GM_readFromIGES(model,fileName);
    MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                "Unknown geo file format %d",format);
    return 0;
  }

  // -------------------------------------------------------------------
  int GM_readFromMSH(pGModel model, const std::string name)
  {
    return (model->readMSH(name));
  }

  // -------------------------------------------------------------------
  int GM_readFromGEO(pGModel model, const std::string name)
  {
    return (model->readGEO(name));
  }

  // -------------------------------------------------------------------
  int GM_readFromSTEP(pGModel model, const std::string name)
  {
    return (model->readSTEP(name));
  }

  // -------------------------------------------------------------------
  int GM_readFromBREP(pGModel model, const std::string name)
  {
    return (model->readBREP(name));
  }

  // -------------------------------------------------------------------
  int GM_readFromIGES(pGModel model, const std::string name)
  {
    return (model->readIGES(name));
  }

  // -------------------------------------------------------------------
  bool GM_physical(const pGModel model)
  {
    return model->physical();
  }

  // -------------------------------------------------------------------
  pGEntity GM_entityByTag(pGModel model, int dim, int tag)
  {
    return model->getEntityByTag(dim,tag);
  }
  pGRegion GM_regionByTag(const pGModel model, int tag)
  {
    return model->getRegionByTag(tag);
  }
  pGFace GM_faceByTag(const pGModel model, int tag)
  {
    return model->getFaceByTag(tag);
  }
  pGEdge GM_edgeByTag(const pGModel model, int tag)
  {
    return model->getEdgeByTag(tag);
  }
  pGVertex GM_vertexByTag(const pGModel model, int tag)
  {
    return model->getVertexByTag(tag);
  }

  // -------------------------------------------------------------------
  int GM_numVertices(const pGModel model)
  {
    return model->getNumVertices();
  }
  int GM_numEdges(const pGModel model)
  {
    return model->getNumEdges();
  }
  int GM_numFaces(const pGModel model)
  {
    return model->getNumFaces();
  }
  int GM_numRegions(const pGModel model)
  {
    return model->getNumRegions();
  }

  // -------------------------------------------------------------------
  GRIter GM_regionIter(pGModel model)
  {
    return new GM_RegionIterator(model);
  }
  GFIter GM_faceIter(pGModel model)
  {
    return new GM_FaceIterator(model);
  }
  GEIter GM_edgeIter(pGModel model)
  {
    return new GM_EdgeIterator(model);
  }
  GVIter GM_vertexIter(pGModel model)
  {
    return new GM_VertexIterator(model);
  }

  pGRegion GRIter_next(GRIter iter) { return iter->next(); }
  pGFace   GFIter_next(GFIter iter) { return iter->next(); }
  pGEdge   GEIter_next(GEIter iter) { return iter->next(); }
  pGVertex GVIter_next(GVIter iter) { return iter->next(); }

  void GRIter_delete(GRIter iter) { delete (iter); }
  void GFIter_delete(GFIter iter) { delete (iter); }
  void GEIter_delete(GEIter iter) { delete (iter); }
  void GVIter_delete(GVIter iter) { delete (iter); }

  void GRIter_reset(GRIter iter) { iter->reset(); }
  void GFIter_reset(GFIter iter) { iter->reset(); }
  void GEIter_reset(GEIter iter) { iter->reset(); }
  void GVIter_reset(GVIter iter) { iter->reset(); }

  // -------------------------------------------------------------------
  int GEN_tag(const pGEntity ent)
  {
    return ent->tag();
  }

  // -------------------------------------------------------------------
  int GEN_type(const pGEntity ent)
  {
    return ent->dim();
  }

  // -------------------------------------------------------------------
  void GEN_setPhysical(pGEntity ent, int dim, int tag)
  {
    ent->setPhysical(dim,tag);
  }

  // -------------------------------------------------------------------
  int GEN_physTag(const pGEntity ent)
  {
    return ent->pTag();
  }

  // -------------------------------------------------------------------
  int GEN_physDim(const pGEntity ent)
  {
    return ent->pDim();
  }

#ifdef _HAVE_GMSH_
  // -------------------------------------------------------------------
  std::vector<pGEntity> GEN_closure(const pGEntity pGE)
  {
    std::vector<pGEntity> theList;

    int type = GEN_type(pGE);
    switch (type) {
    case 0: break;
    case 1: {
      std::vector<pGVertex> vList = GE_vertices((pGEdge)pGE);
      std::vector<pGVertex>::const_iterator vIter = vList.begin();
      for (; vIter != vList.end(); vIter++) {
        theList.push_back(*vIter);
      }
      break;
    }
    case 2: {
      std::vector<pGEdge> eList = GF_edges((pGFace)pGE);
      std::vector<pGEdge>::const_iterator eIter = eList.begin();
      for (; eIter != eList.end(); eIter++) {
        theList.push_back(*eIter);
      }
      break;
    }
    case 3: {
      std::vector<pGFace> fList = GR_faces((pGRegion)pGE);
      std::vector<pGFace>::const_iterator fIter = fList.begin();
      for (; fIter != fList.end(); fIter++) {
        theList.push_back(*fIter);
      }
      break;
    }
    }

    return theList;
  }

  // -------------------------------------------------------------------
  std::vector<pGFace> GR_faces(const pGRegion pGR)
  {
    return pGR->faces();
  }

  // -------------------------------------------------------------------
  int GF_numRegions(const pGFace f)
  {
    return f->numRegions();
  }

  // -------------------------------------------------------------------
  std::vector<pGEdge> GF_edges(const pGFace pGF)
  {
    return pGF->edges();
  }

  // -------------------------------------------------------------------
  //! Computes the parameter location of xyz in the surface. 
  //! Returns false if xyz is not in the surface.
  bool GF_getParams(const pGFace pGF, const double xyz[3], 
                    double params[2])
  {
    //   SPoint2 sP = pGF->parFromPoint(SPoint3(xyz));
    //   params[0] = sP.x();
    //   params[1] = sP.y();

    pGF->XYZtoUV( xyz[0], xyz[1], xyz[2],
                  params[0], params[1], 1.);

    //    throw;

    // a test should be done to check that xyz was on the surface

    return true;
  }

  // -------------------------------------------------------------------
  //! Computes the coordinates of the point on the surface closest to xyz
  void GF_closestPoint(const pGFace pGF, const double xyz[3],
                       const double initGuess[2], double xyzOnF[3])
  {
    //   GPoint gP = pGF->closestPoint( SPoint3(xyz), initGuess );
    //   xyzOnF[0] = gP.x();
    //   xyzOnF[1] = gP.y();
    //   xyzOnF[2] = gP.z();
    throw;
  }

  // -------------------------------------------------------------------
  void GF_xyz(const pGFace pGF, double u, double v, double xyz[3])
  {
    GPoint gP = pGF->point(u,v);
    xyz[0] = gP.x();
    xyz[1] = gP.y();
    xyz[2] = gP.z();
  }

  // -------------------------------------------------------------------
  //! Gets the curvature of the surface computed as the divergence of its normal
  //! Result is bounded by cMaxBound. NaN curvatures are turned into cMaxBound.
  double GF_curvatureDiv(const pGFace surface, const double u[2], 
                         double cMaxBound)
  {
    SPoint2 param(u[0],u[1]);
    double curv = surface->curvatureDiv(param);
    if ( std::isnan(curv) ) {
      MAdMsgSgl::instance().warning(__LINE__,__FILE__,"NaN curvature");
      curv = cMaxBound;
    }
    return std::min(curv,cMaxBound);
  }

  // -------------------------------------------------------------------
  //! Compute the min and max curvatures and the corresponding directions.
  //! Returns the max curvature.
  //! Min and max curvatures are bounded by cMaxBound. 
  //! NaN curvatures are turned into cMaxBound.
  double GF_curvatures(const pGFace surface, const double u[2],
                       double dirMax[3], double dirMin[3],
                       double *curvMax, double *curvMin, 
                       double cMaxBound)
  {
    SPoint2 param(u[0],u[1]);
    SVector3 dirMaxTmp = SVector3();
    SVector3 dirMinTmp = SVector3();

    surface->curvatures(param, &dirMaxTmp, &dirMinTmp, curvMax, curvMin);

    dirMax[0] = dirMaxTmp.x();
    dirMax[1] = dirMaxTmp.y();
    dirMax[2] = dirMaxTmp.z();
    if ( ( dotProd(dirMax,dirMax) <= MAdTOL ) ||
         ( std::isnan(dirMax[0]) || std::isnan(dirMax[1]) || std::isnan(dirMax[2]) ) )
      {
        MAdMsgSgl::instance().warning(__LINE__,__FILE__,
                                      "NaN direction for maximum curvature");
        dirMax[0] = 0.36436431;
        dirMax[1] = 0.76356436;
        dirMax[2] = 0.96983673;
      }
    
    dirMin[0] = dirMinTmp.x();
    dirMin[1] = dirMinTmp.y();
    dirMin[2] = dirMinTmp.z();
    if ( ( dotProd(dirMin,dirMin) <= MAdTOL ) ||
         ( std::isnan(dirMin[0]) || std::isnan(dirMin[1]) || std::isnan(dirMin[2]) ) )
      {
        MAdMsgSgl::instance().warning(__LINE__,__FILE__,
                                      "Inconsistent direction for minimum curvature (u,v)=(%f,%f), direction: %f, %f, %f, curvature: %f",
                                      u[0],u[1],dirMin[0],dirMin[1],dirMin[2],*curvMin);
        double tmp[3] = { 0.86684859, 0.69576964, 0.39876864 };
        crossProd(tmp,dirMax,dirMin);
      }
    
    if ( std::isnan(*curvMax) ) {
      MAdMsgSgl::instance().warning(__LINE__,__FILE__,"NaN maximum curvature");
      *curvMax = cMaxBound;
    }
    *curvMax = std::min(*curvMax,cMaxBound);
    if ( std::isnan(*curvMin) ) {
      MAdMsgSgl::instance().warning(__LINE__,__FILE__,"NaN minimum curvature");
      *curvMin = cMaxBound;
    }
    *curvMin = std::min(*curvMin,cMaxBound);

    return *curvMax;
  }

  // -------------------------------------------------------------------
  //! Computes the parametric coordinates of the point of an edge on 
  //! a geodesic of the surface. t is the location on the edge ( 0 <= t <= 1 ).
  void GF_centerOnGeodesic(const pGFace face, double t, 
                           const double e[2][2], double c[2])
  {
    SPoint2 pt1(e[0][0],e[0][1]);
    SPoint2 pt2(e[1][0],e[1][1]);
    SPoint2 res = face->geodesic(pt1,pt2,t);
    c[0] = res.x();
    c[1] = res.y();
  }

  // -------------------------------------------------------------------
  std::vector<pGVertex> GE_vertices(const pGEdge pGE)
  {
    return pGE->vertices();
  }

  // -------------------------------------------------------------------
  //! Computes the coordinates of the point on the line closest to xyz
  void GE_closestPoint(const pGEdge pGE, const double xyz[3], 
                       double xyzOnE[3])
  {
    //   GPoint gP = pGE->closestPoint( SPoint3(xyz) );
    //   xyzOnE[0] = gP.x();
    //   xyzOnE[1] = gP.y();
    //   xyzOnE[2] = gP.z();
    throw;
  }

  // -------------------------------------------------------------------
  void GE_xyz(const pGEdge pGE, double u, double xyz[3])
  {
    GPoint gP = pGE->point(u);
    xyz[0] = gP.x();
    xyz[1] = gP.y();
    xyz[2] = gP.z();
  }

  // -------------------------------------------------------------------
  //! Given an edge which is a seam of the face, and a point of 
  //! the edge with parametric coordinate uOnE on the edge, 
  //! find the parametric coordinates of the point on the face.
  //! uClose are the parametric coordinates of a close point of 
  //! face used to determine the direction for the parametrization.
  void GE_reparamOnFace(const pGEdge edge, const pGFace face, 
                        double uOnE, double uOnF[2], double uClose[2])
  {
    //   assert( edge->isSeam(face) );

    SPoint2 pt = edge->reparamOnFace(face, uOnE, 0);
    uOnF[0] = pt.x();
    uOnF[1] = pt.y();

    if ( uClose )
      {
        double dist = 
          ( uClose[0] - pt.x() ) * ( uClose[0] - pt.x() ) +
          ( uClose[1] - pt.y() ) * ( uClose[1] - pt.y() );

        SPoint2 pt2 = edge->reparamOnFace(face, uOnE, 1);
        double dist2 = 
          ( uClose[0] - pt2.x() ) * ( uClose[0] - pt2.x() ) +
          ( uClose[1] - pt2.y() ) * ( uClose[1] - pt2.y() );

        if ( dist2 < dist ) {
          uOnF[0] = pt2.x();
          uOnF[1] = pt2.y();
        }
      }
  }

  // -------------------------------------------------------------------
  //! true if the edge is a seam for the given face.
  bool GE_isSeam(const pGEdge edge, const pGFace face)
  {
    return edge->isSeam(face);
  }

  // -------------------------------------------------------------------
  //! gets the curvature of the line at that point bounded by cMaxBound
  double GE_curvature(const pGEdge line, double u, double cMaxBound)
  {
    double curv = line->curvature(u);
    if ( std::isnan(curv) ) {
      MAdMsgSgl::instance().warning(__LINE__,__FILE__,"NaN curvature");
      curv = cMaxBound;
    }
    return std::min(curv,cMaxBound);
  }

  // -------------------------------------------------------------------
  //! returns a list of the lines including the vertex
  std::vector<pGEdge> GV_edges(const pGVertex pGV)
  {
    return pGV->edges();
  }

  // -------------------------------------------------------------------
  //! Given a geometric vertex which is on the face, find the parametric 
  //! coordinates of the vertex on the face.
  //! uClose are the parametric coordinates of a close point of 
  //! face used to determine the direction for the parametrization.
  void GV_reparamOnFace(const pGVertex vertex, const pGFace face, 
                        double uOnF[2], double uClose[2])
  {
    SPoint2 pt = vertex->reparamOnFace(face, 0);
    uOnF[0] = pt.x();
    uOnF[1] = pt.y();

    if ( uClose )
      {
        double dist = 
          ( uClose[0] - pt.x() ) * ( uClose[0] - pt.x() ) +
          ( uClose[1] - pt.y() ) * ( uClose[1] - pt.y() );

        SPoint2 pt2 = vertex->reparamOnFace(face, 1);
        double dist2 = 
          ( uClose[0] - pt2.x() ) * ( uClose[0] - pt2.x() ) +
          ( uClose[1] - pt2.y() ) * ( uClose[1] - pt2.y() );

        if ( dist2 < dist ) {
          uOnF[0] = pt2.x();
          uOnF[1] = pt2.y();
        }
      }
  }

  // -------------------------------------------------------------------
  //! Given a geometric vertex which is on the edge, find the parametric 
  //! coordinates of the vertex on the edge.
  //! uClose is the parametric coordinate of a close point of 
  //! edge used to determine the direction for the parametrization.
  void GV_reparamOnEdge(const pGVertex vertex, const pGEdge edge, 
                        double * uOnE, double uClose)
  {
    Range<double> range = edge->parBounds(0);

    if ( uClose >= 0. && 
         vertex == edge->getBeginVertex() && 
         vertex == edge->getEndVertex() ) {
      *uOnE = range.low();
      double dist  = ( uClose - *uOnE ) * ( uClose - *uOnE );
      double dist2 = ( uClose - range.high() ) * ( uClose - range.high() );
      if ( dist2 < dist ) *uOnE = range.high();
      return;
    }

    if ( vertex == edge->getBeginVertex() ) { *uOnE = range.low();  return; }
    if ( vertex == edge->getEndVertex() )   { *uOnE = range.high(); return; }

    MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                "GVertex not in the closure of the GEdge");
  }

  // -------------------------------------------------------------------
  //! Return true if the vertex is on a seam of the given face
  bool GV_isOnSeam(const pGVertex vert, const pGFace face)
  {
    return vert->isOnSeam(face);
  }

#else

  // -------------------------------------------------------------------
  void GF_centerOnGeodesic(const pGFace face, double t,
                           const double e[2][2], double c[2])
  {
    MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                "Using geodesics requires Gmsh");
  }

#endif

  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  // PGList functions
  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  PGList * PGList_new()
  {
    return new PGList();
  }

  // -------------------------------------------------------------------
  void PGList_delete(PGList * lst)
  {
    if (lst) delete lst;
    lst = NULL;
  }

  // -------------------------------------------------------------------
  void PGList_clear(PGList * lst)
  {
    lst->clear();
  }

  // -------------------------------------------------------------------
  PGList * PGList_appUnique(PGList * lst, pGEntity ent)
  {
    for ( unsigned int i=0; i<lst->entities.size(); i++ ) {
      if ( lst->entities[i] == ent ) return lst;
    }
    lst->entities.push_back(ent);
    return lst;
  }

  // -------------------------------------------------------------------
  PGList * PGList_appPGListUnique(PGList * lst, PGList * source)
  {
    for ( unsigned int iSrc=0; iSrc<source->entities.size(); iSrc++ ) {
      PGList_appUnique(lst,source->entities[iSrc]);
    }
    return lst;
  }

  // -------------------------------------------------------------------
  PGList * PGList_append(PGList * lst, pGEntity ent)
  {
    lst->entities.push_back(ent);
    return lst;
  }

  // -------------------------------------------------------------------
  int PGList_size(PGList * lst)
  {
    return lst->entities.size();
  }

  // -------------------------------------------------------------------
  pGEntity PGList_item(PGList * lst, int i)
  {
    return lst->entities[i];
  }

  // -------------------------------------------------------------------
  pGEntity PGList_next(PGList * lst, void **restart)
  {
    if( *(int*)(restart) >= (int)lst->entities.size() ) return NULL;
    return lst->entities[(*(int*)(restart))++];
  }

  // -------------------------------------------------------------------
  int PGList_inList(PGList * lst, pGEntity ent)
  {
    for ( unsigned int i=0; i<lst->entities.size(); i++ ) {
      if ( lst->entities[i] == ent ) return 1;
    }
    return 0;
  }

  // -------------------------------------------------------------------
  void PGList_remItem(PGList * lst, pGEntity ent)
  {
    std::vector<pGEntity>::iterator eIter = lst->entities.begin();
    for (; eIter != lst->entities.end() ; eIter++) {
      if ( *eIter == ent ) lst->entities.erase(eIter);
    }
  }

  // -------------------------------------------------------------------


}

