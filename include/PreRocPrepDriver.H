#ifndef PREROCPREPDRIVER_H
#define PREROCPREPDRIVER_H
#include <NemDriver.H>
#include <vtkPolyData.h>
#include <unordered_set>

class PreRocPrepDriver : public NemDriver
{
  public:
    PreRocPrepDriver(std::string& fname, int numPartitions, int numPatches);
    ~PreRocPrepDriver();
    static PreRocPrepDriver* readJSON(json inputjson);
    void extractPatches(int numPatches);
    void AddGlobalCellIds(meshBase* mesh);
    // find shared nodes and sent nodes/cells from me to other procs
    void getGhostInformation(int me, bool vol);
    void getGlobalIdsAndMaps(int numPartitions, bool vol);
    void getVirtualCells(int me, int you, bool vol);

    vtkSmartPointer<vtkPolyData>
    deleteInterPartitionSurface(meshBase* fullSurf, vtkSmartPointer<vtkDataSet> partSurf);
  
  private:
    meshBase* mesh;
    meshBase* stitchedSurf;
    std::vector<meshBase*> partitions;
    std::vector<std::map<int, meshBase*>> virtualCellsOfPartitions;
    std::vector<meshBase*> surfacePartitions;
    std::vector<std::map<int, meshBase*>> virtualCellsOfSurfPartitions;
    //std::vector<meshBase*> surfacePartitionPatches;    
    // <proc, <patch, mesh>> 
    std::vector<std::map<int, meshBase*>> patchesOfSurfacePartitions;
   
    // --- the following are cleared before use between vol and surf partitions 
    // global node indices of partition
    std::vector<std::vector<int>> globalNodeIds;             
    // global cell indices of partition  
    std::vector<std::vector<int>> globalCellIds;               
    // <proc, <global nodeId, local nodeId>>
    std::vector<std::map<int,int>> globToPartNodeMap;            
    // <proc, <global cellId, local cellId>>
    std::vector<std::map<int,int>> globToPartCellMap; 
    // <proc, <local nodeId, global nodeId>>
    std::vector<std::map<int,int>> partToGlobNodeMap; 
    // <proc, <local cellId, global cellId>>
    std::vector<std::map<int,int>> partToGlobCellMap;
    // <proc, <global surface nodeId, local surface nodeId>>
    std::vector<std::map<int,int>> globToPartSurfNodeMap;        
    // <proc, <local surface nodeId, global surface nodeId>>
    std::vector<std::map<int,int>> partToGlobSurfNodeMap;        
 
    // --- volume partition ghost information   
    // <proc, <proc, shared nodes>>
    std::map<int, std::map<int, std::vector<int>>> sharedNodes;      
    // <proc, <proc, sent nodes>>
    std::map<int, std::map<int, std::unordered_set<int>>> sentNodes; 
    // <proc, <proc, sent cells>>
    std::map<int, std::map<int, std::unordered_set<int>>> sentCells;
    // <proc, <proc, received nodess>>
    std::map<int, std::map<int, std::unordered_set<int>>> receivedNodes;
    // <proc, <proc, receivedcells>>
    std::map<int, std::map<int, std::unordered_set<int>>> receivedCells;

    // --- surface partition ghost information
    // <proc, <proc, shared nodes>>
    std::map<int, std::map<int, std::vector<int>>> sharedSurfNodes;      
    // <proc, <proc, sent nodes>>
    std::map<int, std::map<int, std::unordered_set<int>>> sentSurfNodes; 
    // <proc, <proc, sent cells>>
    std::map<int, std::map<int, std::unordered_set<int>>> sentSurfCells;
    // <proc, <proc, received nodess>>
    std::map<int, std::map<int, std::unordered_set<int>>> receivedSurfNodes;
    // <proc, <proc, receivedcells>>
    std::map<int, std::map<int, std::unordered_set<int>>> receivedSurfCells;

    // --- patch comm information

    // <proc, <proc, num received nodes>>
    std::map<int, std::map<int, int>> receivedNodesNum;             
    // <proc, < num received cells>>
    std::map<int, std::map<int, int>> receivedCellsNum;              

    // <proc, pConn vector for rocstar>
    std::map<int, std::vector<int>> pConnVec;                        

    // gets shared nodes, sent nodes/cells, received nodes/cells for
    // both partitions me and you
    void getGhostInformation(int me, int you, bool hasShared, bool vol,
                             vtkSmartPointer<vtkIdList> cellIdsList,
                             vtkSmartPointer<vtkGenericCell> genCell);
};


#endif
