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

#include "MAdOutput.h"
#include "MeanRatioEvaluator.h"
#include "OrientedMeanRatioEvaluator.h"
#include "MAdMessage.h"
#include "MathUtils.h"
#include "LocalSizeField.h"
#include "DistanceFunction.h"
#include "MAdResourceManager.h"
#include "MeshSizeBase.h"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <vector>
using std::vector;
using std::string;
using std::set;

namespace MAd {

  // ----------------------------------------------------------------------
  enum dataType {
    SCALAR,
    VECTORIAL
  };

  dataType getOutputDataType(MAdOutputData t) {
    switch (t) {
    case OD_CONSTANT:
    case OD_MEANRATIO:
    case OD_ORIENTEDMEANRATIO:
    case OD_SIZEFIELD_MEAN:
    case OD_SIZEFIELD_MIN:
    case OD_SIZEFIELD_MAX:
    case OD_DIMENSION:
    case OD_ITERORDER:
    case OD_CURVATURE_DIV:
    case OD_CURVATURE_MAX:
    case OD_CURVATURE_MIN:
    case OD_DISTANCE:
      { return SCALAR; break; }
    case OD_CURVATURE_MAX_VEC:
    case OD_CURVATURE_MIN_VEC:
    case OD_ANISO_SF_AXIS0:
    case OD_ANISO_SF_AXIS1:
    case OD_ANISO_SF_AXIS2:
      { return VECTORIAL; break; }
    }
    throw;
  };

  std::string getOutputName( MAdOutputData type ) {
    switch (type) {
    case OD_CONSTANT: return "constant field";
    case OD_MEANRATIO: return "mean ratio";
    case OD_ORIENTEDMEANRATIO: return "oriented mean ratio";
    case OD_SIZEFIELD_MEAN: return "mean size (among the 3 directions) in the size field";
    case OD_SIZEFIELD_MIN: return "minimum size (among the 3 directions) in the size field";
    case OD_SIZEFIELD_MAX: return "maximum size (among the 3 directions) in the size field";
    case OD_DIMENSION: return "element dimension";
    case OD_ITERORDER: return "element id regarding place in iterator";
    case OD_CURVATURE_DIV: return "divergence of the curvature";
    case OD_CURVATURE_MAX: return "surface maximum curvature (scalar field)";
    case OD_CURVATURE_MIN: return "surface minimum curvature (scalar field)";
    case OD_CURVATURE_MAX_VEC: return "surface maximum curvature (vectorial field)";
    case OD_CURVATURE_MIN_VEC: return "surface minimum curvature (vectorial field)";
    case OD_ANISO_SF_AXIS0: return "anisotropic size field along first direction";
    case OD_ANISO_SF_AXIS1: return "anisotropic size field along second direction";
    case OD_ANISO_SF_AXIS2: return "anisotropic size field along third direction";
    case OD_DISTANCE: return "distance to the walls";
    }
    return "unknown output type";
  };

  // ----------------------------------------------------------------------
  double* getData(MAdOutputData type, const pEntity pe, const pSField sf=NULL, int id=0)
  {
    int dim = EN_type(pe);
    double * result;

    switch (type) {

    case OD_DIMENSION: {
      result = new double[4];
      for (int i=0; i<4; i++) result[i] = (double)dim;
      return result;
    }

    case OD_CONSTANT: {
      result = new double[4];
      for (int i=0; i<4; i++) result[i] = 1.;
      return result;
    }

    case OD_MEANRATIO: {
      if ( !sf ) {
        MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                    "No size field given to compute mean ratio");
      }
      if ( sf->getType() != DISCRETESFIELD ) {
        MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                    "Discrete size field required for quality evaluation");
      }
      meanRatioEvaluator evalu((DiscreteSF*)sf);
      double tmp;
      switch (dim) {
      case 3:
        result = new double[4];
        evalu.R_shape((pRegion)pe,&tmp);
        for (int i=0; i<4; i++) result[i] = pow(tmp,1./3.);
        return result;
      case 2:
        result = new double[3];
        evalu.F_shape((pFace)pe,0,&tmp); 
        for (int i=0; i<3; i++) result[i] = sqrt(tmp);
        return result;
      default:
        cerr << "Error in outputs: getData on an element if dimension inferior to 2\n";
        throw;
      }
    }

    case OD_ORIENTEDMEANRATIO: {
      if ( !sf ) {
        MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                    "No size field given to compute mean ratio");
      }
      if ( sf->getType() != DISCRETESFIELD ) {
        MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                    "Discrete size field required for quality evaluation");
      }
      orientedMeanRatioEvaluator evalu((DiscreteSF*)sf);
      double tmp;
      switch (dim) {
      case 3:
        result = new double[4];
        evalu.R_shape((pRegion)pe,&tmp);
        for (int i=0; i<4; i++) result[i] = tmp;
        return result;
      case 2:
        result = new double[3];
        evalu.F_shape((pFace)pe,0,&tmp); 
        for (int i=0; i<3; i++) result[i] = tmp;
        return result;
      default:
        cerr << "Error in outputs: getData on an element if dimension inferior to 2\n";
        throw;
      }
    }

    case OD_ITERORDER: {
      result = new double[4];
      for (int i=0; i<4; i++) result[i] = (double)id;
      return result;
    }

    case OD_SIZEFIELD_MEAN:
    case OD_SIZEFIELD_MIN:
    case OD_SIZEFIELD_MAX: {
      if ( !sf ) {
        MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                    "No size field given");
      }
      void * temp; 
      pPList verts;
      int k;
      pVertex pv;
      pMSize size;
      double h;
      switch (dim) {
      case 3:
        result = new double[4];
        verts = R_vertices((pRegion)pe);
        temp = 0 ; k=0;
        while ( ( pv = (pVertex)PList_next(verts,&temp) ) )
          {
            size = sf->getSize(pv);
            if ( type == OD_SIZEFIELD_MEAN ) h = size->getMeanLength();
            if ( type == OD_SIZEFIELD_MIN )  h = size->getMinLength();
            if ( type == OD_SIZEFIELD_MAX )  h = size->getMaxLength();
            if (size) delete size;
            result[k] = h;
            k++;
          }
        PList_delete(verts);
        return result;
      case 2:
        result = new double[3];
        verts = F_vertices( (pFace)pe, 1);
        temp = 0 ; k=0;
        while ( ( pv = (pVertex)PList_next(verts,&temp) ) )
          {
            size = sf->getSize(pv);
            if ( type == OD_SIZEFIELD_MEAN ) h = size->getMeanLength();
            if ( type == OD_SIZEFIELD_MIN )  h = size->getMinLength();
            if ( type == OD_SIZEFIELD_MAX )  h = size->getMaxLength();
            if (size) delete size;
            result[k] = h;
            k++;
          }
        PList_delete(verts);
        return result;
      default:
        cerr << "Error in outputs: getData on an element if dimension inferior to 2\n";
        throw;
      }
    }

    case OD_CURVATURE_DIV: {
      if ( dim < 3 )
        {
#ifdef _HAVE_GMSH_
          result = new double[4];
          pGEntity pge = EN_whatIn(pe);
          if ( GEN_type(pge) != 2 ) {
            for (int i=0; i<4; i++) result[i] = -1.;
            return result;
          }
          double u[4][2];
          F_params((pFace)pe,u);
          for (int iV=0; iV<F_numVertices((pFace)pe); iV++)
            {
              result[iV] = GF_curvatureDiv((pGFace)pge, u[iV], 1.e6);
            }
          return result;
#else
          MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                      "Gmsh required for divergence of surface curvature");
#endif
        }
      else
        {
          if ( !sf || ( sf->getType() != LOCALSFIELD ) ) {
            MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                        "A local size field is required to output divergence of the curvature on elements");
          }
          result = new double[4];
          LocalSizeField * sf_cast = (LocalSizeField *) sf;
          pRegion pr = (pRegion) pe;
          for (int iV=0; iV<4; iV++)
            {
              if ( sf_cast->getCurvature(R_vertex(pr,iV), &(result[iV])) ) {}
              else result[iV] = -1.;
            }
          return result;
        }
      return NULL;
    }

    case OD_CURVATURE_MAX: {
#ifdef _HAVE_GMSH_
      if ( dim != 2 ) {
        MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                    "Cannot output curvatures on other entities than faces");
      }
      result = new double[4];
      pGEntity pge = EN_whatIn(pe);
      if ( GEN_type(pge) != 2 ) {
        for (int i=0; i<4; i++) result[i] = -1.;
        return result;
      }
      double u[4][2];
      F_params((pFace)pe,u);
      for (int iV=0; iV<F_numVertices((pFace)pe); iV++)
        {
          double dirMax[3], dirMin[3], curvMax, curvMin;
          GF_curvatures((pGFace)pge, u[iV],
                        dirMax, dirMin, &curvMax, &curvMin, 1.e6);
          result[iV] = curvMax;
        }
      return result;
#else
      MAdMsgSgl::instance().error(__LINE__,__FILE__,"Gmsh required for curvature");
#endif
    }

    case OD_CURVATURE_MIN: {
#ifdef _HAVE_GMSH_
      if ( dim != 2 ) {
        MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                    "Cannot output curvatures on other entities than faces");
      }
      result = new double[4];
      pGEntity pge = EN_whatIn(pe);
      if ( GEN_type(pge) != 2 ) {
        for (int i=0; i<4; i++) result[i] = -1.;
        return result;
      }
      double u[4][2];
      F_params((pFace)pe,u);
      for (int iV=0; iV<F_numVertices((pFace)pe); iV++)
        {
          double dirMax[3], dirMin[3], curvMax, curvMin;
          GF_curvatures((pGFace)pge, u[iV],
                        dirMax, dirMin, &curvMax, &curvMin, 1.e6);
          result[iV] = curvMin;
        }
      return result;
#else
      MAdMsgSgl::instance().error(__LINE__,__FILE__,"Gmsh required for curvature");
#endif
    }

    case OD_CURVATURE_MAX_VEC: {
#ifdef _HAVE_GMSH_
      if ( dim != 2 ) {
        MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                    "Cannot output curvatures on other entities than faces");
      }
      result = new double[9];
      pGEntity pge = EN_whatIn(pe);
      if ( GEN_type(pge) != 2 ) {
        for (int i=0; i<9; i++) result[i] = -1.;
        return result;
      }
      double u[3][2];
      F_params((pFace)pe,u);
      for (int iV=0; iV<3; iV++)
        {
          double dirMax[3], dirMin[3], curvMax, curvMin;
          GF_curvatures((pGFace)pge, u[iV],
                        dirMax, dirMin, &curvMax, &curvMin, 1.e6);
          normalizeVec(dirMax,dirMax);
          for (int iC=0; iC<3; iC++) result[3*iV+iC] = dirMax[iC] * curvMax;
        }
      return result;
#else
      MAdMsgSgl::instance().error(__LINE__,__FILE__,"Gmsh required for curvature");
#endif
    }

    case OD_CURVATURE_MIN_VEC: {
#ifdef _HAVE_GMSH_
      if ( dim != 2 ) {
        MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                    "Cannot output curvatures on other entities than faces");
      }
      result = new double[9];
      pGEntity pge = EN_whatIn(pe);
      if ( GEN_type(pge) != 2 ) {
        for (int i=0; i<9; i++) result[i] = -1.;
        return result;
      }
      double u[3][2];
      F_params((pFace)pe,u);
      for (int iV=0; iV<3; iV++)
        {
          double dirMax[3], dirMin[3], curvMax, curvMin;
          GF_curvatures((pGFace)pge, u[iV],
                        dirMax, dirMin, &curvMax, &curvMin, 1.e6);
          normalizeVec(dirMin,dirMin);
          for (int iC=0; iC<3; iC++) result[3*iV+iC] = dirMin[iC] * curvMin;
        }
      return result;
#else
      MAdMsgSgl::instance().error(__LINE__,__FILE__,"Gmsh required for curvature");
#endif
    }

    case OD_ANISO_SF_AXIS0: 
    case OD_ANISO_SF_AXIS1: 
    case OD_ANISO_SF_AXIS2:
      {
        if ( !sf ) MAdMsgSgl::instance().error(__LINE__,__FILE__,"No size field given");
      
        void * temp; 
        pPList verts;
        int k;
        pVertex pv;
        switch (dim) {
        case 3:
          result = new double[12];
          verts = R_vertices((pRegion)pe);
          temp = 0 ; k=0;
          while ( ( pv = (pVertex)PList_next(verts,&temp) ) )
            {
              pMSize size = sf->getSize(pv);
              double h, dir[3];
              if ( type == OD_ANISO_SF_AXIS0 ) h = size->direction(0,dir);
              if ( type == OD_ANISO_SF_AXIS1 ) h = size->direction(1,dir);
              if ( type == OD_ANISO_SF_AXIS2 ) h = size->direction(2,dir);
              if (size) delete size;
              for (int iC=0; iC<3; iC++) result[3*k+iC] = h * dir[iC];
              k++;
            }
          PList_delete(verts);
          return result;
        case 2:
          result = new double[9];
          verts = F_vertices( (pFace)pe, 1);
          temp = 0 ; k=0;
          while ( ( pv = (pVertex)PList_next(verts,&temp) ) )
            {
              pMSize size = sf->getSize(pv);
              double h, dir[3];
              if ( type == OD_ANISO_SF_AXIS0 ) h = size->direction(0,dir);
              if ( type == OD_ANISO_SF_AXIS1 ) h = size->direction(1,dir);
              if ( type == OD_ANISO_SF_AXIS2 ) h = size->direction(2,dir);
              if (size) delete size;
              for (int iC=0; iC<3; iC++) result[3*k+iC] = h * dir[iC];
              k++;
            }
          PList_delete(verts);
          return result;
        default:
          MAdMsgSgl::instance().error(__LINE__,__FILE__,"Dimension < 2: %d",dim);
        }
      }

    case OD_DISTANCE: {
      if ( !sf || ( sf->getType() != LOCALSFIELD ) ) {
        MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                    "A local size field is required to output distance to walls");
      }
      LocalSizeField * sf_cast = (LocalSizeField *) sf;
      if ( dim == 3 ) {
        result = new double[4];
        pRegion pr = (pRegion) pe;
        for (int iV=0; iV<4; iV++) {
          result[iV] = sf_cast->getDistance(R_vertex(pr,iV));
        }
        return result;
      }
      else {
        assert ( dim == 2 );
        result = new double[3];
        pFace pf = (pFace) pe;
        for (int iV=0; iV<3; iV++) {
          result[iV] = sf_cast->getDistance(F_vertex(pf,iV));
        }
        return result;
      }
      return NULL;
    }

    default:
      MAdMsgSgl::instance().error(__LINE__,__FILE__,"Data type not handled in switch");
    }
    return 0;
  }

  // ----------------------------------------------------------------------
  void writeFaces (const pMesh m, const pSField sf, FILE *f, MAdOutputData type)
  {
    int count = 0;
    FIter fit = M_faceIter(m);
    while (pFace pf = FIter_next(fit))
      {
        double xyz[3][3];
        F_coordP1(pf, xyz);
        double* data = getData(type,(pEntity)pf,sf,count);
        if ( getOutputDataType(type) == SCALAR ) {
          fprintf (f,"ST(%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g};\n",
                   xyz[0][0],xyz[0][1],xyz[0][2],
                   xyz[1][0],xyz[1][1],xyz[1][2],
                   xyz[2][0],xyz[2][1],xyz[2][2],
                   data[0],data[1],data[2]);
        }
        if ( getOutputDataType(type) == VECTORIAL ) {
          fprintf (f,"VT(%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g,%g,%g,%g,%g,%g,%g};\n",
                   xyz[0][0],xyz[0][1],xyz[0][2],
                   xyz[1][0],xyz[1][1],xyz[1][2],
                   xyz[2][0],xyz[2][1],xyz[2][2],
                   data[0],data[1],data[2],
                   data[3],data[4],data[5],
                   data[6],data[7],data[8]);
        }
        delete [] data;
        count++;
      }
    FIter_delete(fit);
  }

  // ----------------------------------------------------------------------
  void writeRegions (const pMesh m, const pSField sf, FILE *f, MAdOutputData type)
  {
    int count = 0;
    RIter rit = M_regionIter(m);
    while (pRegion pr = RIter_next(rit))
      {
        double xyz[4][3];
        R_coordP1(pr, xyz);
        double* data = getData(type,(pEntity)pr,sf,count);
        if ( getOutputDataType(type) == SCALAR ) {
          fprintf (f,"SS(%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g,%g};\n",
                   xyz[0][0],xyz[0][1],xyz[0][2],
                   xyz[1][0],xyz[1][1],xyz[1][2],
                   xyz[2][0],xyz[2][1],xyz[2][2],
                   xyz[3][0],xyz[3][1],xyz[3][2],
                   data[0],data[1],data[2],data[3]);
        }
        if ( getOutputDataType(type) == VECTORIAL ) {
          fprintf (f,"VS(%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g};\n",
                   xyz[0][0],xyz[0][1],xyz[0][2],
                   xyz[1][0],xyz[1][1],xyz[1][2],
                   xyz[2][0],xyz[2][1],xyz[2][2],
                   xyz[3][0],xyz[3][1],xyz[3][2],
                   data[0],data[1],data[2],data[3],
                   data[4],data[5],data[6],data[7],
                   data[8],data[9],data[10],data[11]);
        }
        delete [] data;
        count++;
      }
    RIter_delete(rit);
  }

  // ----------------------------------------------------------------------
  void MAdGmshOutput (const pMesh m, const pSField sf, const char *fn, MAdOutputData type)
  {
    MAdResourceManager& tm   = MAdResourceManagerSgl::instance();
    double t0 = tm.getTime();

    if ( sf ) MAdMsgSgl::instance().info(-1,__FILE__,
                                         "Generating output \'%s\' on file \'%s\' with size field \'%s\'",
                                         getOutputName(type).c_str(), fn, sf->getName().c_str());
    else      MAdMsgSgl::instance().info(-1,__FILE__,
                                         "Generating output \'%s\' on file \'%s\'",
                                         getOutputName(type).c_str(), fn);
  
    FILE *f = fopen (fn, "w");
    if ( !f ) {
      cerr << "Error: could not open file " << fn << endl; throw;
    }

    fprintf (f,"View\" mesh \" {\n");

    if (M_dim(m)==2) writeFaces(m,sf,f,type);
    else             writeRegions(m,sf,f,type);

    fprintf (f,"};\n");
    fclose (f);

    MAdMsgSgl::instance().info(-1,__FILE__,"Output generated in %f seconds",
                               tm.getTime()-t0);
  }

  // ----------------------------------------------------------------------
  void MAdAttachedNodalDataOutput(const pMesh m, const char *fn, pMeshDataId id)
  {
    FILE *f = fopen (fn, "w");
    if ( !f ) {
      cerr << "Error: could not open file " << fn << endl; throw;
    }

    fprintf (f,"View\" mesh \" {\n");

    if (M_dim(m)==2) {
    
      FIter fit = M_faceIter(m);
      while (pFace pf = FIter_next(fit))
        {
          // get the coordinates
          double xyz[3][3];
          F_coordP1(pf, xyz);
        
          // get the data at nodes
          double data[3];
          pPList verts = F_vertices(pf,1);
          void* tmp = 0;
          int i = 0;
          while ( pEntity pv = PList_next(verts,&tmp) ) {
            if ( EN_getDataDbl((pEntity)pv,id,&(data[i])) ) {}
            else data[i] = -10.;
            i++;
          }
          PList_delete(verts);

          // write an element
          fprintf (f,"ST(%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g};\n",
                   xyz[0][0],xyz[0][1],xyz[0][2],
                   xyz[1][0],xyz[1][1],xyz[1][2],
                   xyz[2][0],xyz[2][1],xyz[2][2],
                   data[0],data[1],data[2]);
        }
      FIter_delete(fit);
    }
    else {

      RIter rit = M_regionIter(m);
      while (pRegion pr = RIter_next(rit))
        {
          // get the coordinates
          double xyz[4][3];
          R_coordP1(pr, xyz);
        
          // get the data at nodes
          double data[4];
          pPList verts = R_vertices(pr);
          void* tmp = 0;
          int i = 0;
          while ( pEntity pv = PList_next(verts,&tmp) ) {
            if ( EN_getDataDbl(pv,id,&(data[i])) ) {}
            else data[i] = -10.;
            i++;
          }
          PList_delete(verts);

          // write an element
          fprintf (f,"SS(%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g,%g};\n",
                   xyz[0][0],xyz[0][1],xyz[0][2],
                   xyz[1][0],xyz[1][1],xyz[1][2],
                   xyz[2][0],xyz[2][1],xyz[2][2],
                   xyz[3][0],xyz[3][1],xyz[3][2],
                   data[0],data[1],data[2],data[3]);
        }
      RIter_delete(rit);
  
    }

    fprintf (f,"};\n");
    fclose (f);
  }

  // ----------------------------------------------------------------------
  void MAdAttachedNodalDataVecOutput(const pMesh m, const char *fn, pMeshDataId id)
  {
    FILE *f = fopen (fn, "w");
    if ( !f ) {
      cerr << "Error: could not open file " << fn << endl; throw;
    }

    fprintf (f,"View\" mesh \" {\n");

    if (M_dim(m)==2) {
    
      FIter fit = M_faceIter(m);
      while (pFace pf = FIter_next(fit))
        {
          // get the coordinates
          double xyz[3][3];
          F_coordP1(pf, xyz);
        
          // get the data at nodes
          double data[3][3];
          pPList verts = F_vertices(pf,1);
          void* tmp = 0;
          int i = 0;
          while ( pEntity pv = PList_next(verts,&tmp) ) {
            void * tmpDat = NULL;
            if ( EN_getDataPtr((pEntity)pv,id,&tmpDat) )
              {
                //  vector<double> * vec = (vector<double>*) tmpDat;
                //  for (int j=0; j<3; j++)  data[i][j] = (*vec)[j];
                for (int j=0; j<3; j++)  data[i][j] = ((double*)tmpDat)[j];
              }
            else
              {
                for (int j=0; j<3; j++)  data[i][j] = -10.;
              }
            i++;
          }
          PList_delete(verts);

          // write an element
          fprintf (f,"VT(%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g,%g,%g,%g,%g,%g,%g};\n",
                   xyz[0][0],xyz[0][1],xyz[0][2],
                   xyz[1][0],xyz[1][1],xyz[1][2],
                   xyz[2][0],xyz[2][1],xyz[2][2],
                   data[0][0],data[0][1],data[0][2],
                   data[1][0],data[1][1],data[1][2],
                   data[2][0],data[2][1],data[2][2]);
        }
      FIter_delete(fit);
    }
    else {

      RIter rit = M_regionIter(m);
      while (pRegion pr = RIter_next(rit))
        {
          // get the coordinates
          double xyz[4][3];
          R_coordP1(pr, xyz);
        
          // get the data at nodes
          double data[4][3];
          pPList verts = R_vertices(pr);
          void* tmp = 0;
          int i = 0;
          while ( pEntity pv = PList_next(verts,&tmp) ) {
            void * tmpDat = NULL;
            if ( EN_getDataPtr((pEntity)pv,id,&tmpDat) )
              {
                //  vector<double> * vec = (vector<double>*) tmpDat;
                //  for (int j=0; j<3; j++)  data[i][j] = (*vec)[j];
                for (int j=0; j<3; j++)  data[i][j] = ((double*)tmpDat)[j];
              }
            else
              {
                for (int j=0; j<3; j++)  data[i][j] = -10.;
              }
            i++;
          }
          PList_delete(verts);

          // write an element
          fprintf (f,"VS(%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g};\n",
                   xyz[0][0],xyz[0][1],xyz[0][2],
                   xyz[1][0],xyz[1][1],xyz[1][2],
                   xyz[2][0],xyz[2][1],xyz[2][2],
                   xyz[3][0],xyz[3][1],xyz[3][2],
                   data[0][0],data[0][1],data[0][2],
                   data[1][0],data[1][1],data[1][2],
                   data[2][0],data[2][1],data[2][2],
                   data[3][0],data[3][1],data[3][2]);
        }
      RIter_delete(rit);
  
    }

    fprintf (f,"};\n");
    fclose (f);
  }

  // -------------------------------------------------------------------
  void printPosRegions(const set<pRegion> regs, string fn, MAdOutputData type, 
                       const pSField sf, int id)
  {
    FILE *f = fopen (fn.c_str(), "w");
    fprintf (f,"View\" mesh \" {\n");

    set<pRegion>::const_iterator rIter = regs.begin();
    for (; rIter != regs.end(); rIter++) {
      double* data = getData(type,(pEntity)(*rIter),sf);
      double xyz[4][3];
      R_coordP1(*rIter, xyz);
      fprintf (f,"SS(%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g,%g};\n",
               xyz[0][0],xyz[0][1],xyz[0][2],
               xyz[1][0],xyz[1][1],xyz[1][2],
               xyz[2][0],xyz[2][1],xyz[2][2],
               xyz[3][0],xyz[3][1],xyz[3][2],
               data[0],data[1],data[2],data[3]);
      delete [] data;
    }

    fprintf (f,"};\n");
    fclose (f);
  }

  // -------------------------------------------------------------------
  void printPosEntities(const pPList ents, string fn, MAdOutputData type, 
                        const pSField sf, int id)
  {
    FILE *f = fopen (fn.c_str(), "w");
    fprintf (f,"View\" mesh \" {\n");

    void *iter=0;
    while( pEntity ent = PList_next(ents,&iter) ) {
      double* data = getData(type,ent,sf);
      switch ( EN_type(ent) ) {
      case 0: 
        {
          double xyz[3];
          V_coord((pVertex)ent, xyz);
          fprintf (f,"SP(%g,%g,%g) {%g};\n",
                   xyz[0],xyz[1],xyz[2],
                   data[0]);
          break;
        }
      case 1: 
        {
          double xyz[2][3];
          E_coordP1((pEdge)ent, xyz);
          fprintf (f,"SL(%g,%g,%g,%g,%g,%g) {%g,%g};\n",
                   xyz[0][0],xyz[0][1],xyz[0][2],
                   xyz[1][0],xyz[1][1],xyz[1][2],
                   data[0],data[1]);
          break;
        }
      case 2: 
        {
          double xyz[3][3];
          F_coordP1((pFace)ent, xyz);
          fprintf (f,"ST(%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g};\n",
                   xyz[0][0],xyz[0][1],xyz[0][2],
                   xyz[1][0],xyz[1][1],xyz[1][2],
                   xyz[2][0],xyz[2][1],xyz[2][2],
                   data[0],data[1],data[2]);
          break;
        }
      case 3: 
        {
          double xyz[4][3];
          R_coordP1((pRegion)ent, xyz);
          fprintf (f,"SS(%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g) {%g,%g,%g,%g};\n",
                   xyz[0][0],xyz[0][1],xyz[0][2],
                   xyz[1][0],xyz[1][1],xyz[1][2],
                   xyz[2][0],xyz[2][1],xyz[2][2],
                   xyz[3][0],xyz[3][1],xyz[3][2],
                   data[0],data[1],data[2],data[3]);
          break;
        }
      }
      delete [] data;
    }

    fprintf (f,"};\n");
    fclose (f);
  }

  // -------------------------------------------------------------------

}
