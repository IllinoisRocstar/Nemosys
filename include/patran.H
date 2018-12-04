#ifndef patran_H
#define patran_H

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <vtkCellTypes.h>

class meshBase;

namespace PATRAN
{

class patran
{

    public:
      patran(const std::shared_ptr<meshBase> fullMesh, const std::string inFnameVtk,
        const std::string outFnameNeu, std::map<int, int> faceTypeMap, 
        std::map<int, int> nodeTypeMap,
        std::map<int, bool> nodeStructuralMap, std::map<int, bool> nodeMeshMotionMap,
        std::map<int, bool> nodeThermalMap, std::vector<int> nppItr);

    public:
      void write25(std::ofstream& outputStream);
      void write26(std::ofstream& outputStream);
      void write1(std::ofstream& outputStream);
      void write2(std::ofstream& outputStream);
      void write6(std::ofstream& outputStream);
      void write7(std::ofstream& outputStream);
      void write8(std::ofstream& outputStream);
      void write10(std::ofstream& outputStream);
      void write99(std::ofstream& outputStream);
      bool comparePatch(int i, int j);

    private:
      std::shared_ptr<meshBase> fullMesh;
      std::string inFnameVtk;
      std::string outFnameNeu;
      std::shared_ptr<meshBase> volMeshBase;
      std::shared_ptr<meshBase> surfMeshBase;
      std::map<int,int> faceTypeMap;
      std::map<int,int> nodeTypeMap;
      std::map<int,bool> nodeStructuralMap;
      std::map<int,bool> nodeMeshMotionMap;
      std::map<int,bool> nodeThermalMap;
      std::vector<int> nppVec;

      std::map<int,std::string> face2nodes;
      std::map<int,std::vector<int>> boundaryNodeId2PatchNo;


};

} // end namespace patran
#endif
