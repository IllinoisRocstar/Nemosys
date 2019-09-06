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

#include "MeshQualityManager.h"
#include "MeanRatioEvaluator.h"
#include "OrientedMeanRatioEvaluator.h"
#include "CallbackManager.h"
#include "MAdDefines.h"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
using std::ostream;
#include <iomanip>
#include <math.h>
#include <stdio.h>

namespace MAd {

  // -------------------------------------------------------------------
  void QualityDeletingCBFunctionMove (pVertex pv, double *, void *) 
  {
    MeshQualityManagerSgl::instance().clearNeighbourShapes(pv);
  }

  // -------------------------------------------------------------------
  void QualityDeletingCBFunction (pPList before, pPList after, void *,
                                  operationType type , pEntity ppp) 
  {
    switch (type) {
    case MAd_ESPLIT: {
      // In the edge split case, we have to delete the shapes associated to the old elements

      // find the old edge
      void *tmp=0;
      pEntity pE = PList_next(before,&tmp);

      // clear all shapes in the neighbour of the edge
      MeshQualityManagerSgl::instance().clearNeighbourShapes((pEdge) pE);
    
      break;
    } 
    case MAd_ECOLLAPSE: {
      // In the edge collapse case, we have to delete the shapes attached to the neighbour elements of the deleted node.
      // clear all shapes in the neighbour of the deleted node
      MeshQualityManagerSgl::instance().clearNeighbourShapes((pVertex) ppp);
  
      break;
    }
    case MAd_FSWAP:{
      // In the face swap case, we have to delete the shapes associated to the old elements
      void * temp = NULL;
      while ( pEntity ent = PList_next(before,&temp) ) {
        MeshQualityManagerSgl::instance().clearShape(ent);
      }
      break;
    } 
    case MAd_ESWAP: {
      // In the edge swap case, we have to delete the shapes associated to the old elements
      void * temp = NULL;
      while ( pEntity ent = PList_next(before,&temp) ) {
        MeshQualityManagerSgl::instance().clearShape(ent);
      }
      break;
    } 
    case MAd_RREMOVE: {
      // delete the shape associated to the deleted region
      MeshQualityManagerSgl::instance().clearShape(ppp);
      break;
    }
    default: {
      printf("Error in MeshQualityManager: Callback function not implemented for mesh modification %d",
             type);  
      throw;
    }
    }
  }

  // -------------------------------------------------------------------
  void MeshQualityManager::initialize(pMesh m, DiscreteSF * sf, evaluationType type)
  {
    mesh = m;
    sizeField = sf;

    dim = M_dim(mesh);
    shapeId = MD_newMeshDataId("");

    switch(type) {
    case MEANRATIO: {
      if(!evaluator) {
        evaluator = new meanRatioEvaluator(sizeField);
      }
      break;
    }
    default: {
      cerr << "Error: unknown evaluation type when initializing the mesh quality manager\n";
      throw;
    }
    }

    CallBackManagerSgl::instance().registerCallBack(QualityDeletingCBFunction,NULL);
    CallBackManagerSgl::instance().registerCallBackMove(QualityDeletingCBFunctionMove,NULL);

    if(!histogram) {
      histogram = new int[10];
    }
    if(!histogramAvg) {
      histogramAvg = new double[10];
    }
  }

  // -------------------------------------------------------------------
  void MeshQualityManager::setMesh(pMesh m)
  {
    mesh = m;
    dim = M_dim(mesh);
  }

  // -------------------------------------------------------------------
  void MeshQualityManager::finalize()
  {
    clearAllShapes();
    MD_deleteMeshDataId(shapeId);
    if (evaluator) {
      delete evaluator;
      evaluator=NULL;
    }
    if (histogram) {
      delete [] histogram;
      histogram=NULL;
    }
    if (histogramAvg) {
      delete [] histogramAvg;
      histogramAvg=NULL;
    }
  }

  // -------------------------------------------------------------------
  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  void MeshQualityManager::evaluateAllShapes() const
  {
    if (!mesh) {
      cout << "Warning: Could not evaluate shapes in the element evaluator: no mesh specified\n";
      return;
    }

    if ( dim == 3 ) {
      RIter rit = M_regionIter(mesh);
      while ( pRegion pr = RIter_next(rit) ) {
        double tmp;
        if ( !getAttachedShape((pEntity)pr,&tmp) ) {
          double shape;
          evaluator->R_shape(pr,&shape);
          attachShape((pEntity)pr, shape);
        }
      }
      RIter_delete(rit);
    }
    else {
      FIter fit = M_faceIter(mesh);
      while ( pFace pf = FIter_next(fit) ) {
        double tmp;
        if ( !getAttachedShape((pEntity)pf,&tmp) ) {
          double shape;
          evaluator->F_shape(pf,0,&shape);
          attachShape((pEntity)pf, shape);
        }
      }
      FIter_delete(fit);
    }
  }

  // -------------------------------------------------------------------
  void MeshQualityManager::evaluateAndReplaceAllShapes() const
  {
    if (!mesh) {
      cout << "Warning: Could not evaluate shapes in the element evaluator: no mesh specified\n";
      return;
    }

    if ( dim == 3 ) {
      RIter rit = M_regionIter(mesh);
      while ( pRegion pr = RIter_next(rit) ) {
        double shape;
        evaluator->R_shape(pr,&shape);
        attachShape((pEntity)pr, shape);
      }
      RIter_delete(rit);
    }
    else {
      FIter fit = M_faceIter(mesh);
      while ( pFace pf = FIter_next(fit) ) {
        double shape;
        evaluator->F_shape(pf,0,&shape);
        attachShape((pEntity)pf, shape);
      }
      FIter_delete(fit);
    }
  }

  // -------------------------------------------------------------------
  int MeshQualityManager::getShape(pFace pf, double normal[3], double * result) const
  {
    return evaluator->F_shape(pf,normal,result);
  }

  // -------------------------------------------------------------------
  int MeshQualityManager::getShapeWithDisp(pFace pf, double normal[3], 
                                           double disp[3][3], double * result) const
  {
    return evaluator->F_shapeWithDisp(pf,normal,disp,result);
  }

  // -------------------------------------------------------------------
  int MeshQualityManager::getShape(pRegion pr, double * result) const
  {
    if ( !getAttachedShape((pEntity)pr, result) ) {
      int flag = evaluator->R_shape(pr,result);
      attachShape((pEntity)pr,*result);
      return flag;
    }

    if ( *result < MAdTOL ) return 0;
    return 1;
  }

  // -------------------------------------------------------------------
  int MeshQualityManager::getShapeWithDisp(pRegion pr, double disp[4][3], 
                                           double * result) const
  {
    evaluator->R_shapeWithDisp(pr,disp,result);
    if ( *result < MAdTOL ) return 0;
    return 1;
  }

  // -------------------------------------------------------------------
  void MeshQualityManager::clearAllShapes() const
  {
    if (!mesh) {
      cout << "Warning: Could not clear shapes in the element evaluator: no mesh specified\n";
      return;
    }

    if ( dim == 3 ) {
      RIter rit = M_regionIter(mesh);
      while ( pRegion pr = RIter_next(rit) ) {
        clearShape((pEntity)pr);
      }
      RIter_delete(rit);
    }
    else {
      FIter fit = M_faceIter(mesh);
      while ( pFace pf = FIter_next(fit) ) {
        clearShape((pEntity)pf);
      }
      FIter_delete(fit);
    }
  }

  // -------------------------------------------------------------------
  void MeshQualityManager::clearShape(pEntity pe) const
  {
    EN_deleteData(pe, shapeId);
  }

  // -------------------------------------------------------------------
  void MeshQualityManager::clearNeighbourShapes(pVertex pv) const
  {
    if ( dim==3 ) {
      pPList regs = V_regions(pv);
      void* temp=0;
      while ( pEntity region = PList_next(regs,&temp) ) {
        clearShape( region );
      }
      PList_delete(regs);
    }
    else {
      pPList faces = V_faces(pv);
      void* temp=0;
      while ( pEntity face = PList_next(faces,&temp) ) {
        clearShape( face );
      }
      PList_delete(faces);
    }
  }

  // -------------------------------------------------------------------
  void MeshQualityManager::clearNeighbourShapes(pEdge pe) const
  {
    if ( dim==3 ) {
      pPList regs = E_regions(pe);
      void* temp=0;
      while ( pEntity region = PList_next(regs,&temp) ) {
        clearShape( region );
      }
      PList_delete(regs);
    }
    else {
      pPList faces = E_faces(pe);
      void* temp=0;
      while ( pEntity face = PList_next(faces,&temp) ) {
        clearShape( face );
      }
      PList_delete(faces);
    }
  }

  // -------------------------------------------------------------------
  int MeshQualityManager::V_worstShape(pVertex vt, double* result) const 
  {
    pRegion region;
    pPList rlist = V_regions(vt);
    double worst = evaluator->bestShapeEver();
    double shape_1;

    void *temp=0;
  
    if( PList_size(rlist)==0 ) {
      // 2D case
      pPList flist = V_faces(vt);
      pFace face;
      while( ( face = (pFace)PList_next(flist,&temp) ) ) {
        getShape(face,0,&shape_1);
        if( shape_1 <= MAdTOL ) {
          PList_delete(flist);
          *result = shape_1;
          return 0;
        }
        if(shape_1 < worst)
          worst = shape_1;
      }
      PList_delete(flist);
    }
    else {
      // 3D case
      while( ( region = (pRegion)PList_next(rlist,&temp) ) ) {
        getShape(region,&shape_1);
        if( shape_1 <= MAdTOL ) {
          PList_delete(rlist);
          *result = shape_1;
          return 0;
        }
        if(shape_1 < worst)
          worst = shape_1;
      }
    }
  
    PList_delete(rlist);
    *result = worst;
    return 1;
  }

  // -------------------------------------------------------------------

  int MeshQualityManager::E_worstShape(pEdge e, double* result) const {
  
    double worst = evaluator->bestShapeEver();
    double shape_1;

    // 2D case
    if( E_numRegions(e)==0 ) {
      pFace face;
      for( int i=0; i<E_numFaces(e); i++ ) {
        face=E_face(e,i);
        if (!getShape(face,0,&shape_1) || shape_1 <= MAdTOL ) {
          *result = 0.0; return 0;
        }
        if( shape_1 < worst)  worst = shape_1;
      }
    }
    // 3D case
    else {
      pPList rlist = E_regions(e);
      pRegion region; void *temp=0;
      while( ( region=(pRegion)PList_next(rlist,&temp) ) ) {
        if ( !getShape(region,&shape_1) || shape_1 <= MAdTOL ) {
          PList_delete(rlist); *result = 0.0; return 0;
        }
        if( shape_1 < worst)  worst = shape_1;
      }
      PList_delete(rlist);
    }

    *result = worst;
    return 1;
  }

  // -------------------------------------------------------------------

  int MeshQualityManager::F_worstShape(pFace f, double* result) const {
  
    double worst = evaluator->bestShapeEver();
    double shape_1;

    pPList rlist = F_regions(f);
    pRegion region; void *temp=0;
    while( ( region=(pRegion)PList_next(rlist,&temp) ) ) {
      if ( !getShape(region,&shape_1) || shape_1 <= MAdTOL ) {
        PList_delete(rlist); *result = 0.0; return 0;
      }
      if( shape_1 < worst)  worst = shape_1;
    }
    PList_delete(rlist);

    *result = worst;
    return 1;
  }

  // -------------------------------------------------------------------
  int MeshQualityManager::FList_worstShape(pPList faces, double * result) const
  {
    double worst = evaluator->bestShapeEver();

    void* tmp=0;
    while( pEntity pe = PList_next(faces,&tmp) ) {
    
      if ( EN_type(pe) != 2 ) {
        cerr << "Error: Received an entity of dimension " << EN_type(pe) << " in FList_worstShape\n";
        throw;
      }
    
      double shape;
      if ( !getShape( (pFace)pe, 0, &shape) ) {
        *result = 0;
        return 0;
      }
      else {
        worst = evaluator->whatsWorst(shape, worst);
      }
    }

    *result = worst;

    return 1;
  }

  // -------------------------------------------------------------------
  int MeshQualityManager::RList_worstShape(pPList regions, double * result) const
  {
    double worst = evaluator->bestShapeEver();

    void* tmp=0;
    while( pEntity pe = PList_next(regions,&tmp) ) {
    
      if ( EN_type(pe) != 3 ) {
        cerr << "Error: Received an entity of dimension " << EN_type(pe) << " in RList_worstShape\n";
        throw;
      }
    
      double shape;
      if ( !getShape( (pRegion)pe, &shape) ) {
        *result = 0.;
        return 0;
      }
      else {
        worst = evaluator->whatsWorst(shape, worst);
      }
    }

    *result = worst;

    return 1;
  }

  // -------------------------------------------------------------------
  void MeshQualityManager::attachShape(pEntity pe, double shape) const
  {
    EN_modifyDataDbl(pe, shapeId, shape);
  }

  // -------------------------------------------------------------------
  bool MeshQualityManager::getAttachedShape(const pEntity pe, double* result) const
  {
    if ( !EN_getDataDbl(pe, shapeId, result) ) return false;
    return true;
  }

  // -------------------------------------------------------------------
  bool MeshQualityManager::checkAttachedShapes() const
  {
    bool ok = true;

    int noData = 0;
    int dataOK = 0;
    int wrongData = 0;

    if (dim == 3) {
      RIter rit = M_regionIter(mesh);
      while ( pRegion pr = RIter_next(rit) ) {
        double attachedShape;
        if ( !getAttachedShape((pEntity)pr,&attachedShape) ) noData++;
        else {
          double shape;
          evaluator->R_shape(pr,&shape);
          if ( fabs(shape - attachedShape) <= 1.e-12 ) dataOK++;
          else { 
            cout<<"Warning: wrong shape found: attached: "<<attachedShape<<", real: "<<shape<<endl;
            wrongData++; 
            ok = false;
          }
        }
      }
      RIter_delete(rit);
    }
    else {
      FIter fit = M_faceIter(mesh);
      while ( pFace pf = FIter_next(fit) ) {
        double attachedShape;
        if ( !getAttachedShape((pEntity)pf,&attachedShape) ) noData++;
        else {
          double shape;
          evaluator->F_shape(pf,0,&shape);
          if ( fabs(shape - attachedShape) <= 1.e-12 ) dataOK++;
          else { 
            cout<<"Warning: wrong shape found: attached: "<<attachedShape<<", real: "<<shape<<endl;
            wrongData++;
            ok = false;
          }
        }
      }
      FIter_delete(fit);
    }

    cout << "\nAttached shapes report:\n\n";
    cout << "Elements checked: "<<noData+dataOK+wrongData<<endl;
    cout << "Elements with no attached shape: "<<noData<<endl;
    cout << "Elements with a correct shape: "<<dataOK<<endl;
    cout << "Elements with a wrong shape: "<<wrongData<<endl;
    cout<<endl;

    return ok;
  }

  // -------------------------------------------------------------------
  void MeshQualityManager::evaluateSizes()
  {
    minAbsoluteSize = MAdBIG;
    maxAbsoluteSize = 0.;
    sizesSum = 0.;

    switch( dim )
      {
      case 3: {
        RIter rit = M_regionIter(mesh);
        while ( pRegion r = RIter_next(rit) ) 
          {
            double vol = R_volume(r);
            minAbsoluteSize = std::min(minAbsoluteSize,vol);
            maxAbsoluteSize = std::max(maxAbsoluteSize,vol);
            sizesSum += vol;
          }
        RIter_delete(rit);
        break;
      }
      case 2: {
        FIter fit = M_faceIter(mesh);
        while ( pFace f = FIter_next(fit) ) 
          {
            double area = F_area(f,0);
            minAbsoluteSize = std::min(minAbsoluteSize,area);
            maxAbsoluteSize = std::max(maxAbsoluteSize,area);
            sizesSum += area;
          }
        FIter_delete(fit);
        break;
      }
      default: {
        cerr << "Dimension " << dim << " not handled by meshEvaluator\n";
        throw;
      }
      }
  }

  // -------------------------------------------------------------------
  void MeshQualityManager::evaluateShapes()
  {
    worstShape = evaluator->bestShapeEver();
    meanShape = 0.;
    for (int i=0; i<10; i++) histogram[i]=0;
    notAcpt = 0;

    switch( dim )
      {
      case 3: {
        nbElem = M_numRegions(mesh);
        RIter rit = M_regionIter(mesh);
        while ( pRegion r = RIter_next(rit) ) 
          {
            double shape;
            if ( !getShape(r,&shape) )  notAcpt++;
            else if (shape < 0.) {cerr << "Error: element with negative quality\n"; throw;}
            else {
              if (shape >= 1.) { shape = 1.-MAdTOL; }
              unsigned int qualityLevel = (unsigned int)floor(shape*10.);
              if ( qualityLevel == 10 ) qualityLevel = 9;
              histogram[qualityLevel]++;
            }
            meanShape += shape;
            worstShape = evaluator->whatsWorst(worstShape,shape);
          }
        RIter_delete(rit);
        break;
      }
      case 2: {
        nbElem = M_numFaces(mesh);
        FIter fit = M_faceIter(mesh);
        while ( pFace f = FIter_next(fit) ) {
          double shape;
          if ( !getShape(f,0,&shape))  histogram[0]++;
          else if (shape < 0.) {cerr << "Error: element with negative quality\n"; throw;}
          else if (shape > 1.) { shape = 1.; }
          else {
            unsigned int qualityLevel = (unsigned int)floor(shape*10);
            histogram[qualityLevel]++;
          }
          meanShape += shape;
          worstShape = evaluator->whatsWorst(worstShape,shape);
        }
        FIter_delete(fit);
        break;
      }
      default: {
        cerr << "Dimension " << dim << " not handled by meshEvaluator\n";
        throw;
      }
      }

    meanShape /= nbElem;
    for (int i=0; i<10; i++) histogramAvg[i] = (double)histogram[i] / nbElem;
  }

  // -------------------------------------------------------------------
  void MeshQualityManager::evaluateStatistics()
  {
    evaluateSizes();
    evaluateShapes();
  }

  // -------------------------------------------------------------------
  void MeshQualityManager::printStatistics(ostream& out) const
  {
    out << "\n ---------- Mesh quality report ----------\n\n";
    out << "Criterion: " << evaluator->getName() << "\n\n";
    out << "  Average element quality\t"<<std::setprecision(4)<<meanShape<<"\n";
    out << "  Worst element quality\t\t"<<std::setprecision(4)<<worstShape<<"\n";
    out << "  Non acceptable elements\t"<< notAcpt << " \n";
    out << std::setprecision(2) << std::setw(0); 
    out << std::fixed;

    out << "\n        --- Histogram ---\n\n";

    for (int i = 0; i < 10; i++) {
      out <<std::setprecision(1)
          << "  " << ((double)i)/10. << "  < Q <  " << ((double)(i+1))/10. << "  :  "
          <<std::setprecision(2)<<std::setw(5)
          << (histogramAvg[i])*100. <<" % "<<std::setw(7)<< histogram[i]<<" elements\n";
    }
    out<<endl;

    out << std::setprecision(6) << std::setw(0); // return to default values
    out.unsetf(std::ios::floatfield);            // return to default values

    out<<"  Smallest ";
    if (dim == 3) out <<"volume ";
    else out<<"area ";
    out <<"(absolute space)\t"<< minAbsoluteSize << "\n";
    out<<"  Biggest ";
    if (dim == 3) out <<"volume ";
    else out<<"area ";
    out <<"(absolute space)\t"<< maxAbsoluteSize << "\n";
    out<<"  Total ";
    if (dim == 3) out <<"volume ";
    else out<<"area ";
    out <<"   (absolute space)\t"<< sizesSum << "\n\n";
  }

  // -------------------------------------------------------------------

}

