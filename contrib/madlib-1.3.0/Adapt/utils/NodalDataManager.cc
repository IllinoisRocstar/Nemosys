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

#include "NodalDataManager.h"
#include "CallbackManager.h"
#include "MAdOutput.h"
#include "MeshSizeBase.h"

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
using std::vector;
using std::set;
using std::ostream;
using std::string;

namespace MAd {

  // -------------------------------------------------------------------
  void P1CBFunction (pPList before, pPList after, void *dataNames,
                     operationType type , pEntity ppp) {
  
    set<string>* knownDataNames = static_cast<set<string>*>(dataNames);

    switch (type) {
    case MAd_ESPLIT: {
      // In the edge split case, we have to interpolate the data at the new node

      set<string>::const_iterator tIter = (*knownDataNames).begin();
      set<string>::const_iterator tLast = (*knownDataNames).end();
      for (; tIter != tLast; tIter++ ) {

        pMeshDataId dataId = MD_lookupMeshDataId((*tIter).c_str());
    
        // find the old edge
        void *tmp=0;
        pEntity pE = PList_next(before,&tmp);
        double t = E_linearParams((pEdge)pE,(pVertex)ppp);
     
        // get datas at old nodes
        double data0 = 0.;
        pVertex pV0 = E_vertex((pEdge)pE, 0); 
        int gotit0 = EN_getDataDbl((pEntity)pV0, dataId,  &data0);
      
        double data1 = 0.;
        pVertex pV1 = E_vertex((pEdge)pE, 1); 
        int gotit1 = EN_getDataDbl((pEntity)pV1, dataId,  &data1);
      
        if (!gotit0 || !gotit1) {
          printf("Error: one of the nodes has no data attached to with name %s",
                 (*tIter).c_str());
          throw;
        }

        // Interpolate the data at the new node. 
        double newData = (1.-t) * data0 + t * data1;
      
        // attach it
        EN_attachDataDbl(ppp, dataId, newData);
      }

      break;
    } 
    case MAd_ECOLLAPSE: {
      // In the edge collapse case, we have to delete the datas attached to the deleted node
    
      set<string>::const_iterator tIter = (*knownDataNames).begin();
      set<string>::const_iterator tLast = (*knownDataNames).end();
      for (; tIter != tLast; tIter++ ) {

        pMeshDataId dataId = MD_lookupMeshDataId((*tIter).c_str());
        EN_deleteData(ppp, dataId);
      }
  
      break;
    }
    case MAd_ESWAP:
    case MAd_FSWAP: {
      // nothing to be done (no modification at nodes)
      break;
    }
    case MAd_RREMOVE: {
      void * temp = NULL;
      while ( pEntity pE = PList_next(before,&temp) ) {
        if ( EN_type(pE) == 0 ) {
          set<string>::const_iterator tIter = (*knownDataNames).begin();
          set<string>::const_iterator tLast = (*knownDataNames).end();
          for (; tIter != tLast; tIter++ ) {
            pMeshDataId dataId = MD_lookupMeshDataId((*tIter).c_str());
            EN_deleteData(pE, dataId);
          }
        }
      }
      break;
    }
    case MAd_UNKNOWNOPERATION:
    case MAd_VERTEXMOVE:
    case MAd_FCOLLAPSE:
    case MAd_DESPLTCLPS:
    default: {
      MAdMsgSgl::instance().error(__LINE__,__FILE__,
                                  "Not implemented for mesh modification %d",
                                  type);
    }
    }
  
  }

  // -------------------------------------------------------------------
  void P1VecCBFunction (pPList before, pPList after, void *dataNames,
                        operationType type , pEntity ppp) {
  
    set<string>* knownDataNames = static_cast<set<string>*>(dataNames);

    switch (type) {
    case MAd_ESPLIT: {
      // In the edge split case, we have to interpolate the data at the new node

      set<string>::const_iterator tIter = (*knownDataNames).begin();
      set<string>::const_iterator tLast = (*knownDataNames).end();
      for (; tIter != tLast; tIter++ ) {

        pMeshDataId dataId = MD_lookupMeshDataId((*tIter).c_str());
    
        // find the old edge
        void *tmp=0;
        pEntity pE = PList_next(before,&tmp);
        double t = E_linearParams((pEdge)pE,(pVertex)ppp);
     
        // get coordinates and datas at old nodes
        void * data0 = NULL;
        pVertex pV0 = E_vertex((pEdge)pE, 0);
        int gotit0 = EN_getDataPtr((pEntity)pV0, dataId,  &data0);
        vector<double> * vec0 = (vector<double>*) data0;
      
        void * data1 = NULL;
        pVertex pV1 = E_vertex((pEdge)pE, 1);
        int gotit1 = EN_getDataPtr((pEntity)pV1, dataId,  &data1);
        vector<double> * vec1 = (vector<double>*) data1;
      
        if (!gotit0 || !gotit1) {
          printf("Error: one of the nodes has no data attached to with name %s",
                 (*tIter).c_str());
          throw;
        }

        // Interpolate the data at the new node. 
        vector<double> * newData = new vector<double>;
        for (int i=0; i<3; i++)  (*newData).push_back( (1.-t) * (*vec0)[i] + t * (*vec1)[i] );
      
        // attach it
        EN_attachDataPtr(ppp, dataId, newData);
      }

      break;
    } 
    case MAd_ECOLLAPSE: {
      // In the edge collapse case, we have to delete the datas attached to the deleted node
          
      set<string>::const_iterator tIter = (*knownDataNames).begin();
      set<string>::const_iterator tLast = (*knownDataNames).end();
      for (; tIter != tLast; tIter++ ) {

        pMeshDataId dataId = MD_lookupMeshDataId((*tIter).c_str());
        void * data = NULL;
        int gotit = EN_getDataPtr((pEntity)ppp, dataId,  &data);
        vector<double> * vec = (vector<double>*) data;
        if (gotit && vec) delete vec;
        EN_deleteData(ppp, dataId);
      }
  
      break;
    }
    case MAd_ESWAP:
    case MAd_FSWAP: {
      // nothing to be done (no modification at nodes)
      break;
    }
    case MAd_RREMOVE: {
      void * temp = NULL;
      while ( pEntity pE = PList_next(before,&temp) ) {
        if ( EN_type(pE) == 0 ) {
          set<string>::const_iterator tIter = (*knownDataNames).begin();
          set<string>::const_iterator tLast = (*knownDataNames).end();
          for (; tIter != tLast; tIter++ ) {
            pMeshDataId dataId = MD_lookupMeshDataId((*tIter).c_str());
            void * data = NULL;
            int gotit = EN_getDataPtr((pEntity)ppp, dataId,  &data);
            vector<double> * vec = (vector<double>*) data;
            if (gotit && vec) delete vec;
            EN_deleteData(pE, dataId);
          }
        }
      }
      break;
    }
    case MAd_UNKNOWNOPERATION:
    case MAd_VERTEXMOVE:
    case MAd_FCOLLAPSE:
    case MAd_DESPLTCLPS:
    default: {
      printf("Error in NodalDataManager: Callback function not implemented for mesh modification %d",
             type);  
      throw;
    }
    }
  
  }

  // -------------------------------------------------------------------
  void NodalDataManager::initialize(pMesh m)
  {
    mesh = m;

    prefix    = "NodalDataManagerId_";
    prefixVec = prefix + "Vec_";

    coordNameBase = "StoredCoordinates__";

    if(knownDataNames )
      {
        knownDataNames->clear();
      }
    else
      { 
        knownDataNames= new set<string>();
      }
    if ( knownDataVecNames  )
      {
        knownDataVecNames->clear();
      }
    else
      {
        knownDataVecNames = new set<string>();
      }

    CallBackManagerSgl::instance().registerCallBack(P1CBFunction,   knownDataNames   );
    CallBackManagerSgl::instance().registerCallBack(P1VecCBFunction,knownDataVecNames);
  }

  // -------------------------------------------------------------------
  void NodalDataManager::setMesh(pMesh m)
  {
    mesh = m;
  }

  // -------------------------------------------------------------------
  void NodalDataManager::finalize()
  {
    removeAllData();
    if (knownDataNames) 
      {
        delete knownDataNames;
        knownDataNames = NULL;
      }
    if (knownDataVecNames)
      {
        delete knownDataVecNames;
        knownDataVecNames = NULL;
      }
  }

  // -------------------------------------------------------------------
  void NodalDataManager::diagnostics (ostream& out) const
  {
    out << "\nContent of the nodal data manager:\n";

    out << "\n  - Double scalar : " << knownDataNames->size() << " fields:\n";
    set<string>::iterator iIter;
    set<string>::iterator iLast;
    iIter = (*knownDataNames).begin();
    iLast = (*knownDataNames).end();
    for (; iIter != iLast; iIter++) {
      out << "     - "<< *iIter << "\n";
    }

    out << "\n  - Double vector : " << knownDataVecNames->size() << " vector fields:\n";
    iIter = (*knownDataVecNames).begin();
    iLast = (*knownDataVecNames).end();
    for (; iIter != iLast; iIter++) {
      out << "     - "<< *iIter << "\n";
    }
  
    out << "\n" << endl;
  }

  // -------------------------------------------------------------------
  void NodalDataManager::removeAllData()
  {
    set<string>::iterator iIter;
    set<string>::iterator iLast;

    iIter = (*knownDataNames).begin();
    iLast = (*knownDataNames).end();
    for (; iIter != iLast; iIter++)  removeData(*iIter);

    iIter = (*knownDataVecNames).begin();
    iLast = (*knownDataVecNames).end();
    for (; iIter != iLast; iIter++)  removeVData(*iIter);
  }

  // -------------------------------------------------------------------
  void NodalDataManager::registerData(string name, const vector<double> datas)
  {
    string dataName = prefix + name;

    // check that this name doesn't exists yet
    if ( nameExists(dataName,0) ) {
      cerr << "Error: this scalar data \'"<<name<<"\' is already registered in the NodalDataManager\n";
      throw;
    }

    // check that the number of datas is equal to the number of nodes in the mesh
    int nbDatas = datas.size();
    if ( nbDatas != M_numVertices(mesh) ) {
      cerr<< "Error: wrong number of datas ("<<nbDatas<<") with "<<M_numVertices(mesh)<<" nodes\n";
      throw;
    }

    // add this name to the list
    (*knownDataNames).insert(dataName);

    // attach the data at nodes
    pMeshDataId id = MD_lookupMeshDataId( dataName.c_str() );
    int i = 0;
    VIter vit = M_vertexIter(mesh);
    while (pVertex pv = VIter_next(vit))
      {
        EN_attachDataDbl((pEntity)pv,id,datas[i]);
        i++;
      }
    VIter_delete(vit);
  }

  // -------------------------------------------------------------------
  void NodalDataManager::registerVData(string name, 
                                       const vector<vector<double> > datas)
  {
    string dataName = prefixVec + name;

    // check that this name doesn't exists yet
    if ( nameExists(dataName,1) ) {
      cerr << "Error: this vectorial data \'"<<name<<"\' is already registered in the NodalDataManager\n";
      throw;
    }

    // check that the number of datas is equal to the number of nodes in the mesh
    int nbDatas = datas.size();
    if ( nbDatas != M_numVertices(mesh) ) {
      cerr<< "Error: wrong number of datas ("<<nbDatas<<") with "<<M_numVertices(mesh)<<" nodes\n";
      throw;
    }

    // add this name to the list
    (*knownDataVecNames).insert(dataName);

    // attach the data at nodes
    pMeshDataId id = MD_lookupMeshDataId( dataName.c_str() );
    int i = 0;
    VIter vit = M_vertexIter(mesh);
    while (pVertex pv = VIter_next(vit))
      {
        vector<double> * copy = new vector<double>(datas[i]);
        EN_attachDataPtr((pEntity)pv,id,copy);
        i++;
      }
    VIter_delete(vit);
  }

  // -------------------------------------------------------------------
  void NodalDataManager::getMeshData(string name, vector<double> * result) const
  {
    string dataName = prefix + name;

    // check that this name exists
    if ( !nameExists(dataName,0) ) {
      cerr << "Error: this scalar data \'"<<name<<"\' is not registered in the NodalDataManager\n";
      throw;
    }

    result->clear();

    pMeshDataId id = MD_lookupMeshDataId( dataName.c_str() );
  
    int i = 0;
    VIter vit = M_vertexIter(mesh);
    while (pVertex pv = VIter_next(vit))
      {
        double theDbl;
        EN_getDataDbl((pEntity)pv,id,&theDbl);
        result->push_back(theDbl);
        i++;
      }
    VIter_delete(vit);
  }

  // -------------------------------------------------------------------
  void NodalDataManager::getMeshVData(string name, vector<vector<double> > * result) const
  {
    string dataName = prefixVec + name;

    // check that this name exists
    if ( !nameExists(dataName,1) ) {
      cerr << "Error: this vectorial data \'"<<name<<"\' is not registered in the NodalDataManager\n";
      throw;
    }

    result->clear();

    pMeshDataId id = MD_lookupMeshDataId( dataName.c_str() );
  
    int i = 0;
    VIter vit = M_vertexIter(mesh);
    while (pVertex pv = VIter_next(vit))
      {
        void * tmp = 0;
        EN_getDataPtr((pEntity)pv,id,&tmp);
        vector<double> * vec = (vector<double>*) tmp;
        result->push_back(*vec);
        i++;
      }
    VIter_delete(vit);
  }

  // -------------------------------------------------------------------
  void NodalDataManager::getData(string name, pVertex pv, double* res) const
  {
    string dataName = prefix + name;

    // check that this name exists
    if ( !nameExists(dataName,0) ) {
      cerr << "Error: this scalar data \'"<<name<<"\' is not registered in the NodalDataManager\n";
      throw;
    }

    pMeshDataId id = MD_lookupMeshDataId( dataName.c_str() );
  
    EN_getDataDbl((pEntity)pv,id,res);
  }

  // -------------------------------------------------------------------
  void NodalDataManager::getVData(string name, pVertex pv, vector<double>& res) const
  {
    string dataName = prefixVec + name;

    // check that this name exists
    if ( !nameExists(dataName,1) ) {
      cerr << "Error: this vectorial data \'"<<name<<"\' is not registered in the NodalDataManager\n";
      throw;
    }

    pMeshDataId id = MD_lookupMeshDataId( dataName.c_str() );
  
    void * tmp = 0;
    EN_getDataPtr((pEntity)pv,id,&tmp);
    res.clear();
    for (int i=0; i<3; i++) res.push_back(  (*( (vector<double>*) tmp))[i]  );

  }

  // -------------------------------------------------------------------
  void NodalDataManager::removeData(string name)
  {
    string dataName = prefix + name;

    if ( !nameExists(dataName,0) )  return;

    // remove it from the list of names
    set<string>::iterator iIter = (*knownDataNames).find(dataName);
    (*knownDataNames).erase(iIter);

    // delete the datas at nodes
    pMeshDataId id = MD_lookupMeshDataId( (dataName).c_str() );
    VIter vit = M_vertexIter(mesh);
    while (pVertex pv = VIter_next(vit))
      {
        EN_deleteData((pEntity)pv,id);
      }
    VIter_delete(vit);
  }

  // -------------------------------------------------------------------
  void NodalDataManager::removeVData(string name)
  {
    string dataName = prefixVec + name;

    if ( !nameExists(dataName,1) )  return;

    // remove it from the list of names
    set<string>::iterator iIter = (*knownDataVecNames).find(dataName);
    (*knownDataVecNames).erase(iIter);

    // delete the datas at nodes
    pMeshDataId id = MD_lookupMeshDataId( dataName.c_str() );
    VIter vit = M_vertexIter(mesh);
    while (pVertex pv = VIter_next(vit))
      {
        void * data = NULL;
        int gotit = EN_getDataPtr((pEntity)pv, id,  &data);
        vector<double> * vec = (vector<double>*) data;
        if (gotit && vec) delete vec;
        EN_deleteData((pEntity)pv, id);
      }
    VIter_delete(vit);
  }

  // -------------------------------------------------------------------
  void NodalDataManager::writeData(string name, const char *fn)
  {
    string dataName = prefix + name;

    // check that this name exists
    if ( !nameExists(dataName,0) ) {
      cerr << "Error: this data \'"<<name<<"\' is not registered in the NodalDataManager\n";
      throw;
    }

    pMeshDataId id = MD_lookupMeshDataId( dataName.c_str() );
  
    MAdAttachedNodalDataOutput(mesh, fn, id);
  }

  // MS
  // -------------------------------------------------------------------
  void NodalDataManager::writeDataCSV(string name, const char *fn)
  {
    string dataName = prefix + name;

    // check that this name exists
    if ( !nameExists(dataName,0) ) {
      cerr << "Error: this data \'"<<name<<"\' is not registered in the NodalDataManager\n";
      throw;
    }

    pMeshDataId id = MD_lookupMeshDataId( dataName.c_str() );

    MAdAttachedNodalDataCSVOutput(mesh, fn, id, name);
  }
  // MS END

  // -------------------------------------------------------------------
  void NodalDataManager::writeVData(string name, const char *fn)
  {
//     string dataName = prefixVec + name;

//     // check that this name exists
//     if ( !nameExists(dataName,1) ) {
//       cerr << "Error: this data \'"<<name<<"\' is not registered in the NodalDataManager\n";
//       throw;
//     }

//     pMeshDataId id = MD_lookupMeshDataId( dataName.c_str() );
  
    MAdMsgSgl::instance().warning(__LINE__,__FILE__,
                                  "Output for vectorial attached data is not implemented");
    // it is implemented for double arrays, not vectors. 
    // This class should store double arrays too, to be done (GC).
    //    MAdAttachedNodalDataVecOutput(mesh, fn, id);
  }

  // -------------------------------------------------------------------
  void NodalDataManager::storeCoordinates()
  {
    removeCoordinates();

    // add the coordinates name to the list
    string coordName = prefixVec + coordNameBase;
    (*knownDataVecNames).insert(coordName);

    // attach the data at nodes
    pMeshDataId id = MD_lookupMeshDataId( coordName.c_str() );
    VIter vit = M_vertexIter(mesh);
    while (pVertex pv = VIter_next(vit))
      {
        double xyz[3]; V_coord(pv,xyz);
        vector<double> * vec = new vector<double>;
        for (int j = 0; j<3; j++)  (*vec).push_back(xyz[j]);
        EN_attachDataPtr((pEntity)pv,id,vec);
      }
    VIter_delete(vit);
  }

  // -------------------------------------------------------------------
  void NodalDataManager::getStoredCoordinates(pVertex pv, vector<double>& res)
  {
    getVData(coordNameBase,pv,res);
  }

  // -------------------------------------------------------------------
  void NodalDataManager::removeCoordinates()
  {
    removeVData(coordNameBase);
  }

  // -------------------------------------------------------------------
  bool NodalDataManager::isCoordinates()
  {
    string coordName = prefixVec + coordNameBase;
    return nameExists(coordName,1);
  }

  // -------------------------------------------------------------------

  bool NodalDataManager::nameExists(string dataName, int vec) const
  {
    if ( vec == 0 ) {
      set<string>::iterator iIter = (*knownDataNames).find(dataName);
      if ( iIter == (*knownDataNames).end() ) return false;
      else return true;
    }
    if ( vec == 1 ) {
      set<string>::iterator iIter = (*knownDataVecNames).find(dataName);
      if ( iIter == (*knownDataVecNames).end() ) return false;
      else return true;
    }
    return true;
  }

  // -------------------------------------------------------------------

}
