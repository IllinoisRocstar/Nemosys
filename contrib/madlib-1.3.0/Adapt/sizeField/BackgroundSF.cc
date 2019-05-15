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

#include "BackgroundSF.h"
#include "IsoMeshSize.h"
#include "MAdDefines.h"

#include <math.h>
#include <cmath>
#include <string.h>
#include <map>

namespace MAd {

  // -------------------------------------------------------------------
  BackgroundSF::BackgroundSF(std::string _name): 
    SizeFieldBase(_name), bgModel(NULL), bgMesh(NULL)
  {
    pMSizeFieldId = MD_newMeshDataId("");
  }

  // -------------------------------------------------------------------
  BackgroundSF::~BackgroundSF()
  {
    MD_deleteMeshDataId(pMSizeFieldId);
    if ( bgMesh ) M_delete(bgMesh);
    if ( bgModel ) GM_delete(bgModel);
  }

  // -------------------------------------------------------------------
  void BackgroundSF::setSize(pEntity ent, double h)
  {
    pMSize pS = new IsoMeshSize(h);
    setSize(ent,pS);
  }

  // -------------------------------------------------------------------
  void BackgroundSF::setSize(pEntity pEnt, pMSize pS)
  {
    void * temp;
    if( EN_getDataPtr(pEnt,pMSizeFieldId,&temp) ) {
      delete (pMSize)temp;
      EN_modifyDataPtr(pEnt,pMSizeFieldId,pS);
    }
    else {
      EN_attachDataPtr(pEnt,pMSizeFieldId,pS);
    }
  }

  // -------------------------------------------------------------------
  // GC warning: very low implementation: should use octrees
  pMSize BackgroundSF::getSize(const pVertex pv) const
  {
    if ( M_dim(bgMesh) < 3 ) {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                  "Background size fields not implemented for dim < 3");
    }

    double xyz[3];
    V_coord(pv,xyz);

    pRegion pr;
    pRegion rBest = NULL;
    {
      RIter rit = M_regionIter(bgMesh);
      while( ( pr = RIter_next(rit) ) )
        {
          if ( !R_inBox(pr,xyz,1.e-6) ) continue;
          if ( !R_contains(pr,xyz,1.e-6) ) continue;
          rBest = pr;
          break;
        }
      RIter_delete(rit);
    }
    if ( !rBest ) {
    
      double best = 1.e12;
      double tmp[4], tmp2;
      RIter rit = M_regionIter(bgMesh);
      while( ( pr = RIter_next(rit) ) )
        {
          R_linearParams(pr,xyz,tmp);
          tmp[3] = 1. - tmp[0] - tmp[1] - tmp[2];
          
          double mi = std::min(tmp[0],std::min(tmp[1],std::min(tmp[2],tmp[3])));
          double ma = std::max(tmp[0],std::max(tmp[1],std::max(tmp[2],tmp[3])));
        
          tmp2 = 0.;
          if ( mi < 0. ) tmp2 = std::abs(mi);
          if ( ma > 1. ) tmp2 = std::max(tmp2,ma-1.);

          if ( tmp2 < best ) {
            best = tmp2;
            rBest = pr;
          }
        }
      RIter_delete(rit);
    }

    void * temp;
    if( EN_getDataPtr((pEntity)rBest,pMSizeFieldId,&temp) ) {
      return MS_copy((pMSize)temp);
    }
    return NULL;
  }

  // -------------------------------------------------------------------
  void BackgroundSF::loadData(std::string mshName)
  {
    GM_create(&bgModel,"bgModel");
    bgMesh = M_new(bgModel);
    M_load(bgMesh,mshName.c_str());

    FILE * fp = fopen(mshName.c_str(),"r");
    if(!fp) {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                  "Unknown file %s",mshName.c_str());
    }

#ifdef PARALLEL
    MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                "Set sizes from file not implemented in parallel"); 
#endif
    
    char str[256] = "XXX";
    double version = 1.0;
    
    std::map<int,pEntity> elId;

    while(1) {

      while(str[0] != '$') {
        if(!fgets(str, sizeof(str), fp) || feof(fp)) {
          MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                      "No data section found in file %s",
                                      mshName.c_str());
        }
      }

      if(!strncmp(&str[1], "MeshFormat", 10)) {
        if(!fgets(str, sizeof(str), fp)) throw;
        int format, size;
        if(sscanf(str, "%lf %d %d", &version, &format, &size) != 3) throw;
      }
      else if(!strncmp(&str[1], "ELM", 3) || !strncmp(&str[1], "Elements", 8)) {
        if(!fgets(str, sizeof(str), fp)) throw;
        int numElements;
        sscanf(str, "%d", &numElements);

        for(int i = 0; i < numElements; i++) {
          int num, type, physical = 0, elementary = 0, partition = 0, numV;
          if(version <= 1.0){
            fscanf(fp, "%d %d %d %d %d", &num, &type, &physical, &elementary, &numV);
          }
          else{
            int numTags;
            fscanf(fp, "%d %d %d", &num, &type, &numTags);
            //	  printf("%d %d %d\n", num, type, numTags);
            //	  getchar();
            for(int j = 0; j < numTags; j++){
              int tag;
              fscanf(fp, "%d", &tag);       
              if(j == 0)      physical = tag;
              else if(j == 1) elementary = tag;
              else if(j == 2) partition = tag;
              // ignore any other tags for now
            }
            if (type == 15) numV = 1;
            else if (type == 1) numV = 2;
            else if (type == 2) numV = 3;
            else if (type == 3) numV = 4;
            else if (type == 4) numV = 4;
            else throw;
          }

          int indices[60];
          for(int j = 0; j < numV; j++) fscanf(fp, "%d", &indices[j]);
	
          if ( M_dim(bgMesh) == 3 && type == 4 ){
            pRegion pr = M_region( bgMesh, indices );
            if ( !pr ) printf("didn't find element with nodes %d %d %d %d\n",indices[0],indices[1],indices[2],indices[3]);
            elId[num] = (pEntity)pr;
          }
          else if ( M_dim(bgMesh) == 2 && type == 2 ){
            pFace pf = M_face(  bgMesh, indices );
            elId[num] = (pEntity)pf;
          }
        }
      }
//       else if(!strncmp(&str[1], "NodeData", 8)) {

//         int numStringTags;
//         fscanf(fp,"%d\n",&numStringTags);
//         char dummy[256];
//         for (int i=0;i<numStringTags;i++){
//           fgets(dummy,256,fp);
//         }

//         int numDoubleTags;
//         fscanf(fp,"%d\n",&numDoubleTags);
//         for (int i=0;i<numDoubleTags;i++){
//           fgets(dummy,256,fp);
//         }

//         int numIntTags;
//         fscanf(fp,"%d\n",&numIntTags);
//         int numDataVertices;
//         for (int i=0;i<numIntTags;i++){	
//           if (i == 2)fscanf(fp,"%d",&numDataVertices);
//           else fgets(dummy,256,fp);
//         }

//         if ( M_numVertices(bgMesh) != numDataVertices ) {
//           MAdMsgSgl::instance().error(__LINE__,__FILE__,
//                                       "Data available for %d vertices in file %s while there are %d vertices in the mesh",
//                                       numDataVertices, mshName.c_str(),
//                                       M_numVertices(bgMesh));
//         }

//         for (int i=0; i<numDataVertices; i++) {
//           double dData;
//           int vId;
//           fscanf(fp,"%d %lf",&vId,&dData);
//           setSize(vId,dData);
//         }

//         break;
//       }
      
      else if(!strncmp(&str[1], "ElementData", 11)) {

        int numStringTags;
        fscanf(fp,"%d\n",&numStringTags);
        char dummy[256];
        for (int i=0;i<numStringTags;i++){
          fgets(dummy,256,fp);
        }

        int numDoubleTags;
        fscanf(fp,"%d\n",&numDoubleTags);
        for (int i=0;i<numDoubleTags;i++){
          fgets(dummy,256,fp);
        }

        int numIntTags;
        fscanf(fp,"%d\n",&numIntTags);
        int numDataEl;
        for (int i=0;i<numIntTags;i++){	
          if (i == 2)fscanf(fp,"%d",&numDataEl);
          else fgets(dummy,256,fp);
        }

        int nElem = M_numRegions(bgMesh);
        if ( !nElem ) nElem = M_numFaces(bgMesh);
        if ( nElem != numDataEl ) {
          MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                      "Data available for %d elements in file %s while there are %d elements in the mesh",
                                      numDataEl, mshName.c_str(), nElem );
        }

        for (int i=0; i<nElem; i++) {
          double dData;
          int eId;
          fscanf(fp,"%d %lf",&eId,&dData);
          pEntity pe = elId[eId];
          setSize(pe,dData);
        }

        break;
      }
    
      do {
        if(!fgets(str, sizeof(str), fp) || feof(fp)) {
          MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                      "No data section found in file %s",
                                      mshName.c_str());
        }
      } while(str[0] != '$');
    }
      
    fclose(fp);
    
  }

  // -------------------------------------------------------------------
  void BackgroundSF::redistributeToVertices()
  {
    VIter vIter = M_vertexIter(bgMesh);
    while( pVertex pV = VIter_next(vIter) )
      {
        pPList vRegs = V_regions(pV);
        void * tmp = NULL;
        pRegion pr;
        void * ms;
        double h = 0.;
        int count = 0;
        while ( ( pr = (pRegion)PList_next(vRegs,&tmp) ) ) {
          EN_getDataPtr( (pEntity)pr, pMSizeFieldId, &ms);
          if ( !ms ) {
            MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                        "No size defined for an entity");
          }
          h += ((pMSize)ms)->getMeanLength();
          count++;
        }
        PList_delete(vRegs);
        h /= (double)count;
        setSize((pEntity)pV,h);
      }
    VIter_delete(vIter);
    
  }

  // -------------------------------------------------------------------
}
