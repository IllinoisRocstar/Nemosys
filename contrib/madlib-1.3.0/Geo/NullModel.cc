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

#ifndef _HAVE_GMSH_

#include "NullModel.h"

#include "MAdMessage.h"
#include "MshTags.h"

#include <stdio.h>
#include <set>
using std::set;
#include <string.h>

namespace MAd {

  // -------------------------------------------------------------------
  NullModel::~NullModel()
  {
    clean();
  }

  // -------------------------------------------------------------------
  void NullModel::clean()
  {
    riter itR = regions.begin();
    for (; itR != regions.end(); itR++) delete (*itR);
    regions.clear();
    fiter itF = faces.begin();
    for (; itF != faces.end(); itF++) delete (*itF);
    faces.clear();
    eiter itE = edges.begin();
    for (; itE != edges.end(); itE++) delete (*itE);
    edges.clear();
    viter itV = vertices.begin();
    for (; itV != vertices.end(); itV++) delete (*itV);
    vertices.clear();
  }

  // -------------------------------------------------------------------
  MAdGEntity * NullModel::getEntityByTag(int dim, int tag)
  {
    switch(dim) {
    case 3: return getRegionByTag(tag);
    case 2: return getFaceByTag(tag);
    case 1: return getEdgeByTag(tag);
    case 0: return getVertexByTag(tag);
    }
    return NULL;
  }

  // -------------------------------------------------------------------
  MAdGRegion *NullModel::getRegionByTag(int tag)
  {
    set<MAdGRegion *,MAdGEntityLessThan>::const_iterator it = regions.begin();
    for (; it != regions.end(); it++) {
      if ( (*it)->tag() == tag ) return *it;
    }
    MAdGRegion * newRegion = new MAdGRegion(tag,this);
    regions.insert(newRegion);
    return newRegion;
  }

  // -------------------------------------------------------------------
  MAdGFace *NullModel::getFaceByTag(int tag)
  {
    set<MAdGFace *,MAdGEntityLessThan>::const_iterator it = faces.begin();
    for (; it != faces.end(); it++) {
      if ( (*it)->tag() == tag ) return *it;
    }
    MAdGFace * newFace = new MAdGFace(tag,this);
    faces.insert(newFace);
    return newFace;
  }

  // -------------------------------------------------------------------
  MAdGEdge *NullModel::getEdgeByTag(int tag)
  {
    set<MAdGEdge *,MAdGEntityLessThan>::const_iterator it = edges.begin();
    for (; it != edges.end(); it++) {
      if ( (*it)->tag() == tag ) return *it;
    }
    MAdGEdge * newEdge = new MAdGEdge(tag,this);
    edges.insert(newEdge);
    return newEdge;
  }

  // -------------------------------------------------------------------
  MAdGVertex *NullModel::getVertexByTag(int tag)
  {
    set<MAdGVertex *,MAdGEntityLessThan>::const_iterator it = vertices.begin();
    for (; it != vertices.end(); it++) {
      if ( (*it)->tag() == tag ) return *it;
    }
    MAdGVertex * newVertex = new MAdGVertex(tag,this);
    vertices.insert(newVertex);
    return newVertex;
  }

  // -------------------------------------------------------------------
  int NullModel::readGEO(const std::string &filename)
  {
    MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                "Reading a model from a .geo file not implemented in NullModel");
  }

  // -------------------------------------------------------------------
  void SwapBytes(char *array, int size, int n)
  {
    char *x = new char[size];
    for(int i = 0; i < n; i++) {
      char *a = &array[i * size];
      memcpy(x, a, size);
      for(int c = 0; c < size; c++)
        a[size - 1 - c] = x[c];
    }
    delete [] x;
  }

  // -------------------------------------------------------------------
  int NullModel::readMSH(const std::string &filename)
  {
    FILE *fp = fopen(filename.c_str(), "rb");
    if(!fp){
      MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                  "Unable to open file '%s'", filename.c_str());
      return 0;
    }

    char str[256] = "XXX";
    double version = 1.0;
    bool binary = false, swap = false;
 
    while(1) {

      while(str[0] != '$'){
        if(!fgets(str, sizeof(str), fp) || feof(fp))
          break;
      }
    
      if(feof(fp))
        break;

      if(!strncmp(&str[1], "MeshFormat", 10)) {

        if(!fgets(str, sizeof(str), fp)) return 0;
        int format, size;
        if(sscanf(str, "%lf %d %d", &version, &format, &size) != 3) return 0;
        if(format){
          binary = true;
          MAdMsgSgl::instance().info(__LINE__,__FILE__,
                                     "Mesh is in binary format");
          int one;
          if(fread(&one, sizeof(int), 1, fp) != 1) return 0;
          if(one != 1){
            swap = true;
            MAdMsgSgl::instance().info(__LINE__,__FILE__,
                                       "Swapping bytes from binary file");
          }
        }

      }
      else if(!strncmp(&str[1], "PhysicalNames", 13)) {
      }
      else if(!strncmp(&str[1], "NO", 2) || !strncmp(&str[1], "Nodes", 5) ||
              !strncmp(&str[1], "ParametricNodes", 15)) {
      }
      else if(!strncmp(&str[1], "ELM", 3) || !strncmp(&str[1], "Elements", 8)) {

        if(!fgets(str, sizeof(str), fp)) return 0;
        int numElements;
        sscanf(str, "%d", &numElements);
        MAdMsgSgl::instance().info(__LINE__,__FILE__,
                                   "%d elements", numElements);
        if(!binary){
          for(int i = 0; i < numElements; i++) {
            int num, type, physical = -1, elementary = -1, partition = 0, numVertices;
            if(version <= 1.0){
              fscanf(fp, "%d %d %d %d %d", &num, &type, &physical, &elementary, &numVertices);
            }
            else{
              int numTags;
              fscanf(fp, "%d %d %d", &num, &type, &numTags);
              for(int j = 0; j < numTags; j++){
                int tag;
                fscanf(fp, "%d", &tag);       
                if(j == 0)      physical = tag;
                else if(j == 1) elementary = tag;
                else if(j == 2) partition = tag;
                // ignore any other tags for now
                numVertices = getNumVerticesForElementTypeMSH(type);
              }
            }
            int tmp;
            for(int j = 0; j < numVertices; j++) fscanf(fp, "%d", &tmp);

            int gTag = elementary;
            if ( physical != -1 ) gTag = physical;
            int gDim = getDimForElementTypeMSH(type);
            getEntityByTag(gDim,gTag);
          }
        }
        else{
          int numElementsPartial = 0;
          while(numElementsPartial < numElements){
            int header[3];
            if(fread(header, sizeof(int), 3, fp) != 3) return 0;
            if(swap) SwapBytes((char*)header, sizeof(int), 3);
            int type = header[0];
            int numElms = header[1];
            int numTags = header[2];
            int numVertices = getNumVerticesForElementTypeMSH(type);
            unsigned int n = 1 + numTags + numVertices;
            int *data = new int[n];
            for(int i = 0; i < numElms; i++) {
              if(fread(data, sizeof(int), n, fp) != n) return 0;
              if(swap) SwapBytes((char*)data, sizeof(int), n);
              int num = data[0];
              int physical = (numTags > 0) ? data[4 - numTags] : -1;
              int elementary = (numTags > 1) ? data[4 - numTags + 1] : -1;
              int partition = (numTags > 2) ? data[4 - numTags + 2] : 0;
              int *indices = &data[numTags + 1];
              int gTag = elementary;
              if ( physical != -1 ) gTag = physical;
              int gDim = getDimForElementTypeMSH(type);
              getEntityByTag(gDim,gTag);
            }
            delete [] data;
            numElementsPartial += numElms;
          }
        }

      }

      do {
        if(!fgets(str, sizeof(str), fp) || feof(fp))
          break;
      } while(str[0] != '$');
    }

    fclose(fp);
    return 1;
  }

  // -------------------------------------------------------------------
  int NullModel::readSTEP(const std::string &filename)
  {
    MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                "Null model cannot read STEP format");
  }

  // -------------------------------------------------------------------
  int NullModel::readBREP(const std::string &filename)
  {
    MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                "Null model cannot read BREP format");
  }

  // -------------------------------------------------------------------
  int NullModel::readIGES(const std::string &filename)
  {
    MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                "Null model cannot read IGES format");
  }

  // -------------------------------------------------------------------

}

#endif

