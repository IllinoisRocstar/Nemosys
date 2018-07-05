#ifndef ROCPREPDRIVER_H
#define ROCPREPDRIVER_H

#include <NemDriver.H>
#include <vtkProcess.h>
#include <unordered_set>

class RocPrepDriver : public NemDriver
{
  public:
    RocPrepDriver(std::string& fname, int numPartitions);
    ~RocPrepDriver();
    static RocPrepDriver* readJSON(json inputjson);

  private:
    meshBase* mesh;
    std::vector<meshBase*> partitions;

};

// vtkProcess is an abstract class representing a process that can be launched
// by a vtkMultiProcessController. Concrete classes just have to implement
// Execute() method and make sure it set the proper value in ReturnValue.
// The controller has a function consuming vtkProcess objects (SetSingleProcessObject())
//  where the method to execute to is the vtkProcess' Execute() and the data can be 
// encapsulated in the vtkProcess

class GhostGenerator : public vtkProcess
{
  public:
    static GhostGenerator* New();
    // single process function, called by each proc
    virtual void Execute();
    // sets partitioned meshes which are generated outside of mpi world
    void setPartitions(std::vector<meshBase*>& partitions);

    ~GhostGenerator() 
    {
      delete ghostCellMesh;
    }

  protected:
    GhostGenerator();
    std::vector<meshBase*> partitions;          // vector of partitioned meshes
    meshBase* ghostCellMesh;                    // this proc's mesh with ghost info
    std::vector<int> myGlobalNodeIds;           // global node indices of partition
    std::vector<int> myGlobalGhostNodeIds;      // global ghost node indices of partition
    std::vector<int> myGlobalCellIds;           // global cell indices of partition
    std::vector<int> myGlobalGhostCellIds;      // global ghost cell indices of partition
    std::map<int, std::vector<int>> sharedNodes;// <proc, shared nodes>
    std::map<int, std::unordered_set<int>> sentNodes; // <proc, sent nodes> 
    std::map<int, std::unordered_set<int>> sentCells; // <proc, sent cells>
    std::map<int, int> receivedNodesNum;        // <proc, num received nodes>
    std::map<int, int> receivedCellsNum;        // <proc, num received cells>
    std::map<int,int> myGlobToPartNodeMap;      // <global nodeId, local nodeId>
    std::map<int,int> myPartToGlobNodeMap;      // <local nodeId, global nodeId>
    std::vector<int> pConnVec;                  // pConn vector for rocstar
     
 
    // populates the vectors and maps above
    void getPconnInformation(int me, int numProcs);
    // helpers to get pconn
    void getGlobalIds(int me);
    void getGlobalGhostIds(int me); 

    // write pconn vec
    void writeSharedToPconn(const std::string& type);
    void writeSentToPconn(const std::string& type, bool nodeOrCell);
    void writeReceivedToPconn(const std::string& type, int me, bool nodeOrCell);


};

std::vector<std::string> getCgFNames(const std::string& case_dir, 
                                     const std::string& prefix,
                                     const std::string& base_t);


#endif
