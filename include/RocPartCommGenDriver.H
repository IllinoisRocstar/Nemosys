#ifndef ROCPARTCOMMGENDRIVER_H
#define ROCPARTCOMMGENDRIVER_H
#include <NemDriver.H>
#include <vtkPolyData.h>
#include <unordered_set>

class RocPartCommGenDriver : public NemDriver
{
  public:
    RocPartCommGenDriver(std::shared_ptr<meshBase> vol, 
                         std::shared_ptr<meshBase> surf, 
                         std::shared_ptr<meshBase> volWithSol,
                         std::shared_ptr<meshBase> surfWithSol,
                         std::shared_ptr<meshBase> burnSurfWithSol,
                         int numPartitions, const std::string& base_t,
                         bool writeIntermediateFiles, double searchTolerance,
                         const std::string& caseName);
    RocPartCommGenDriver(std::string& volname, std::string& surfname, int numPartitions);
    ~RocPartCommGenDriver();
    static RocPartCommGenDriver* readJSON(json inputjson);

  // --- executor functions
  private:    
    // run the driver (called on construction)
    void execute(int numPartitions);
    // extract patches from each surface partition and get patch virtual cells
    void extractPatches();
    // add global cell ids to provided mesh (used to add to full surface)
    void AddGlobalCellIds(std::shared_ptr<meshBase> mesh);
    // find shared nodes and sent nodes/cells from me to other procs
    void getGhostInformation(int me, bool vol);
    // get global ids and maps which were loaded into vol mesh partitions during partitioning
    void getGlobalIdsAndMaps(int numPartitions, bool vol);
    // get virtual cells for complete (not patch) vol (vol=true) or surf (vol=false) partition
    void getVirtualCells(int me, int you, bool vol);
    // get the shared nodes between patches both intra and inter partition
    void getSharedPatchInformation();

  // --- info for write to cgns and intermediaries
  private:
    // base time step (rocstar convention, eg. 04.210000)
    std::string base_t;
    // trimmed time step (eg. 4.21)
    std::string trimmed_base_t;
    // if true, all intermediate files will be written (i.e. surf partitions, virtuals etc.)
    bool writeAllFiles;
    // Rocstar case name, used to write Rocstar input file names
    std::string caseName;

  // --- write cgns files
  private:
    // write vol cgns (proc is partition num, type is patch num, prefix usually fluid)
    void writeVolCgns(const std::string& prefix, int proc, int type, int gsExists); 
    // write surf cgns (proc is partition num, other info is deduced from maps)
    void writeSurfCgns(const std::string& prefix, int proc); 
    //void writeVolGs(const std::string& prefix, int proc, int type);

  // --- write Rocstar files
  private:
    // write region mapping file
    void mapWriter();
    // write cell mapping file
    void cmpWriter(int proc, int nCells);
    // write communication lists file
    void comWriter(int proc);
    // write dimensions file
    void dimWriter(int proc, std::shared_ptr<meshBase> realMB, std::shared_ptr<meshBase> totalMB);
    void dimSurfWriter(int proc, std::vector<int> cgConnReal, std::vector<int> cgConnVirtual, int patchNo);
    void dimSurfWriter(int proc);
    // write processor mapping file
    void txtWriter();
    // utilities

  // --- real and virtual meshes 
  private:
    // remeshed volume to be partitioned
    std::shared_ptr<meshBase> mesh;
    // stitched surface
    std::shared_ptr<meshBase> remeshedSurf;
    // original volume mesh with solution data
    std::shared_ptr<meshBase> volWithSol;
    // original surface mesh with solution data
    std::shared_ptr<meshBase> surfWithSol;
    // original burning surface mesh with solution data (iburn)
    std::shared_ptr<meshBase> burnSurfWithSol;
    // TODO: add meshBase (mbObj[2) from remeshDriver) for iburn solution mesh (iBurnWithSol)
    //       figure out a way to extract the cells with bflag on each surf/vol partition
    //        transfer data from iBurnWithSol to surf partitions and write them out to cgns
    // partitions <proc, partition>
    std::vector<std::shared_ptr<meshBase>> partitions;
    // <proc, <proc, mesh>>
    std::vector<std::map<int, std::shared_ptr<meshBase>>> virtualCellsOfPartitions;
    // surface partitions <proc, partition>
    std::vector<std::shared_ptr<meshBase>> surfacePartitions;
    // <proc, <proc, mesh>>
    std::vector<std::map<int, std::shared_ptr<meshBase>>> virtualCellsOfSurfPartitions;
    // <proc, <patch, mesh>> 
    std::vector<std::map<int, std::shared_ptr<meshBase>>> patchesOfSurfacePartitions;
    // <proc, <proc, <patch, mesh>>>
    std::vector<std::map<int, std::map<int, std::shared_ptr<meshBase>>>> 
      virtualCellsOfPatchesOfSurfacePartitions;  
  
  // --- tolerance for cell and point search
  private:
    double searchTolerance;

  // --- global to local node and cell index maps
  // --- the following are cleared and reused for vol and surf partitions 
  private: 

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

  // --- volume partition ghost information   
  private:
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
  private:
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
  private:
    /* shared nodes b/w patches and proc [proc][patch][proc][patch]
       eg) partition 0's patch i potentially shares with patch j on proc 0 and 
           possibly the other procss */
    std::vector<std::map<int, std::map<int, std::map<int, std::vector<int>>>>> sharedPatchNodes;

  // --- pconn vectors for vol and surf partitions
  private:
    // pconn for each vol partition 
    // <proc, pconns>
    std::vector<std::vector<int>> volPconns;
    // pconn for each patch of each surf partition
    // <proc, <patch, pconns>>
    std::vector<std::map<int, std::vector<int>>> surfPconns;
    // maximum Pconn values for each proc
    std::map<int,int> pconnProcMax;
    std::map<int,int> pconnProcMin;

  // --- storing number of unique volumetric faces for each partition
    // <proc, number of faces>
    std::map<int, int> nUniqueVolFaces;

  // --- volume to surface node maps for each partition
  private: 
    // <proc, <surface nodeId, volume nodeId>>
    //std::vector<std::map<int,int>> surfToVolNodeMap;            
    // <proc, <volume nodeId, surface nodeId>>
    //std::vector<std::map<int,int>> volToSurfNodeMap; 
    // Stores volume meshes with virtual for each proc
        std::vector<std::shared_ptr<meshBase>> procVolMesh;
    // Stores surface meshes for each proc
    //std::vector<std::vector<std::vector<double>>> procSurfMesh;

  // --- write pconn information into pconn vectors
  private:
    // write shared pconn vec for vol partition proc - returns length
    int writeSharedToPconn(int proc, const std::string& type);

    // write shared pconn vec for all patches in surf pconn proc
    void writeSharedToPconn(int proc);
    // write sent pconn vec for vol partition proc for nodes (nodeOrCell=true) and cells
    void writeSentToPconn(int proc, const std::string& type, bool nodeOrCell);
    // write received pconn vec for vol partition proc for nodes (nodeOrCell=true) and cells
    void writeReceivedToPconn(int proc, const std::string& type, bool nodeOrCell);
    // get shared nodes, sent nodes/cells, received nodes/cells for
    // both partitions me and you
    void getGhostInformation(int me, int you, bool hasShared, bool vol,
                             vtkSmartPointer<vtkIdList> cellIdsList,
                             vtkSmartPointer<vtkGenericCell> genCell);
    // restructure Pconn information so that Rocstar can read it properly
    void restructurePconn(std::vector<int> &pconnVec, int proc, int volOrSurf, std::map<int,int> old2New);
    // Map old to new indices for point data
    void mapOld2NewPointData(std::vector<double> &pointData, const std::map<int,int> &new2Old);
    // function for re-ordering tri indices
    void swapTriOrdering(std::vector<int> &connVec);


  private:
    // --- temporary vectors to assist with writing out pconn info to Rocstar input files
    std::vector<int> neighborProcsTmp;
    std::vector< std::vector<int> > sharedNodesTmp;
    std::vector< std::vector<int> > sentNodesTmp;
    std::vector< std::vector<int> > receivedNodesTmp;
    std::vector< std::vector<int> > sentCellsTmp;
    std::vector< std::vector<int> > receivedCellsTmp;
    std::map<int,int> totalTrisPerPatch;
    // Rocflu
    // <proc<patch,zone>>
    std::vector<std::map<int, int>> surfZoneMap;
    // <proc,zone>
    std::vector<int> volZoneMap;
    // Rocburn
    std::vector<std::map<int, int>> burnSurfZoneMap;

  // --- utility to clear above vectors
  private:
    void clearPconnVectors();
    void clearBorderVector();

  // --- descriptors for cgns writing
  private:  
  std::vector<int> notGhostInPconn;
  // helpers
  private:  
    /* deletes inter partition surfaces which result from extracting surface of vol partition.
       this is done by comparing cells from full surface mesh to extracted surface and
       removing those from the extracted surface which are not found in the full surface */
    vtkSmartPointer<vtkPolyData>
    deleteInterPartitionSurface(std::shared_ptr<meshBase> fullSurf, vtkSmartPointer<vtkDataSet> partSurf);
};


#endif
