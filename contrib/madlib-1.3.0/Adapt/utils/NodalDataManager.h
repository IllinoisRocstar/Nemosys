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
// Authors: Gaetan Compere, Jean-Francois Remacle
// -------------------------------------------------------------------

#ifndef _H_NODALDATAMANAGER
#define _H_NODALDATAMANAGER

#include "MAdSingleton.h"
#include "MeshDataBaseInterface.h"

#include <vector>
#include <string>
#include <set>
#include <iostream>

namespace MAd {

  // -------------------------------------------------------------------

  class NodalDataManager {

  public:
  
    NodalDataManager():mesh(NULL),
                       knownDataNames(NULL), knownDataVecNames(NULL){};
    ~NodalDataManager() {};
  
    void initialize(pMesh m);
    void setMesh(pMesh m);
    void finalize();

    void removeAllData();

    void diagnostics (std::ostream& out) const;

    // will attach the datas to the nodes of the mesh
    // order in vector is the node iterator order
    void registerData  (std::string name, const std::vector<double>);
    void registerVData (std::string name, const std::vector<std::vector<double> >);

    // will fill the variable with the datas attached to the nodes
    // this function does not allocate memory
    // order in vector is the node iterator order
    void getMeshData  (std::string name, std::vector<double> *)          const;
    void getMeshVData (std::string name, std::vector<std::vector<double> > *) const;

    // get the data at a node
    void getData  (std::string name, pVertex pv, double* res)         const;
    void getVData (std::string name, pVertex pv, std::vector<double>& res) const;

    // will delete the data attached to the nodes
    void removeData  (std::string name);
    void removeVData (std::string name);

    // write the data in a postprocessing file
    void writeData  (std::string name, const char *fn);
    void writeVData (std::string name, const char *fn);

    // functions to keep track of the initial coordinates
    void storeCoordinates();
    void getStoredCoordinates(pVertex pv, std::vector<double>& res);
    void removeCoordinates();
    bool isCoordinates();
  
  private:

    pMesh mesh;

    // prefix for identification of the datas from this manager
    // ensures there is no clash with other datas
    std::string prefix;
    std::string prefixVec;

    // identifier for the coordinates
    std::string coordNameBase;

    // stored with the prefix
    std::set<std::string>* knownDataNames;
    std::set<std::string>* knownDataVecNames;

  private:

    bool nameExists(std::string name, int vec) const;

  };

  // -------------------------------------------------------------------

  typedef MAdSingleton<NodalDataManager> NodalDataManagerSgl;

}

#endif
