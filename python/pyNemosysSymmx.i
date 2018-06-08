%module pyNemosys
%include "std_string.i"
%include "std_vector.i"
%include "std_map.i"
%{
#include "meshBase.H"
#include "vtkMesh.H"
#include "TransferDriver.H"
#include "NemDriver.H"
#include "RefineDriver.H"
#include "MeshGenDriver.H"
#include "MeshQualityDriver.H"
#include "ConversionDriver.H"
#include "RemeshDriver.H"
#include "RocRestartDriver.H"
#include "OrderOfAccuracy.H"
#include "jsoncons/json.hpp"
#include "meshGen.H"
#include "meshingParams.H"
#include "symmxGen.H"
#include "symmxParams.H"
#include "cgnsAnalyzer.H"
#include "meshStitcher.H"
#include "meshPartitioner.H"
%}



%template(vectorString) std::vector<std::string>;
%template(doubleV) std::vector<double>;
%template(intV) std::vector<int>;
%template(doubleVV) std::vector<std::vector<double>>;
%template(cellMap) std::map<int, std::vector<double>>;

class meshBase
{

  // constructors and destructors
  public:

    meshBase():numPoints(0),numCells(0),hasSizeField(0),checkQuality(0);
    ~meshBase();

    static meshBase* Create(std::string fname);
    static meshBase* Create(vtkSmartPointer<vtkDataSet> other, std::string newname);
    static std::shared_ptr<meshBase> CreateShared(std::string fname);
    //static std::unique_ptr<meshBase> CreateUnique(std::string fname);
    static meshBase* generateMesh(std::string fname, std::string meshEngine,
                                  meshingParams* params); 
    virtual std::vector<double> getPoint(int id);
    virtual std::vector<std::vector<double>> getVertCrds() const;
    virtual std::map<int, std::vector<double>> getCell(int id);
    virtual std::vector<std::vector<double>> getCellVec(int id);
    vtkSmartPointer<vtkDataSet> getDataSet();
    virtual vtkSmartPointer<vtkDataSet> extractSurface();
    virtual void inspectEdges(const std::string& ofname);
    virtual void setPointDataArray(const char* name,
                                     const std::vector<std::vector<double>>& data);
    virtual void setCellDataArray(const char* name,
                                    const std::vector<std::vector<double>>& data);
    virtual void setCellDataArray(const char* name,
                                    const std::vector<double>& data);
    virtual void unsetPointDataArray(int arrayID);
    virtual void unsetPointDataArray(const char* name);
    virtual void unsetCellDataArray(int arrayID);
    virtual void unsetCellDataArray(const char* name);
    virtual void unsetFieldDataArray(const char* name);
    virtual std::vector<double> getCellLengths();
    virtual std::vector<double> getCellCenter(int cellID);
    vtkSmartPointer<vtkCellLocator> buildLocator();
    int transfer(meshBase* target, std::string method,
                 const std::vector<int>& arrayIDs);
    virtual int getCellType() const = 0;
    int transfer(meshBase* target, std::string method,
                 const std::vector<std::string>& arrayNames);
    int transfer(meshBase* target, std::string method);
    
    std::vector<std::vector<double>> integrateOverMesh(const std::vector<int>& arrayIDs);
    
    void generateSizeField(std::string method, int arrayID, double dev_mlt, bool maxIsmin, double sizeFactor=1.);

    void setSFBool(bool q);
    bool getSFBool();
    int IsArrayName(std::string name);
    void setOrder(int _order); 
    int getOrder() const; 

    void refineMesh(std::string method, int arrayID,
                    double dev_mult, bool maxIsmin, double edge_scale, std::string ofname, bool transferData,
                    double sizeFactor=1.);
    void refineMesh(std::string method, std::string arrayName,
                    double dev_mult, bool maxIsmin, double edge_scale, std::string ofname, bool transferData, 
                    double sizeFactor=1.);
    void refineMesh(std::string method, double edge_scale, std::string ofname, bool transferData);

    virtual void report();
    int getNumberOfPoints();
    int getNumberOfCells();
    void checkMesh(std::string ofname);

    virtual void write();
    virtual void write(std::string fname);
    void writeMSH(std::ofstream& outputStream);
    void writeMSH(std::string fname);
    void writeMSH(std::ofstream& outputStream, std::string pointOrCell, int arrayID);
    void writeMSH(std::string fname, std::string pointOrCell, int arrayID);
    void writeMSH(std::ofstream& outputStream, std::string pointOrCell, int arrayID,
                    bool onlyVol);
    void writeMSH(std::string fname, std::string pointOrCell, int arrayID,
                    bool onlyVol);

    void setFileName(std::string fname);
    std::string getFileName();
    void setCheckQuality(bool x);
    void setContBool(bool x);
    void setNewArrayNames(const std::vector<std::string>& newnames);
    void unsetNewArrayNames();
};


class vtkMesh : public meshBase
{
  // constructor and destructor
  public:
    vtkMesh();
    vtkMesh(const char* fname);
    vtkMesh(const char* fname1, const char* fname2);
    vtkMesh(vtkSmartPointer<vtkDataSet> other, std::string fname);
    ~vtkMesh();

  // access
  public:
    // get point with id
    std::vector<double> getPoint(int id)  const;
    // get cell with id : returns point indices and respective coordinates
    std::map<int, std::vector<double>> getCell(int id) const;
    std::vector<std::vector<double>> getCellVec(int id) const;
    // get diameter of circumsphere of each cell
    std::vector<double> getCellLengths() const;
    // get center of a cell
    std::vector<double> getCellCenter(int cellID) const;
    // get cell type as an integer
    // assumes all elements are the same type
    int getCellType() const;
 
  // diagnostics
  public:
    void report();
    void write();
    void write(std::string fname); 

  // processing
  public:
    vtkSmartPointer<vtkDataSet> extractSurface();

  // set and get point and cell data
  public:
    // set point data (numComponets per point determined by dim of data[0] 
    void setPointDataArray(const char* name, const std::vector<std::vector<double>>& data);
    // set cell data (numComponents per cell determined by dim of data[0])
    void setCellDataArray(const char* name, const std::vector<std::vector<double>>& data);
    // set scalar cell data
    void setCellDataArray(const char* name, const std::vector<double>& data);
    // remove point data with given id from dataSet if it exists
    void unsetPointDataArray(int arrayID);
    void unsetPointDataArray(const char* name);
    // remove cell data with given id from dataSet if it exists
    void unsetCellDataArray(int arrayID);
    void unsetCellDataArray(const char* name);
    // remove field data with given id from dataSet
    void unsetFieldDataArray(const char* name);
};


int diffMesh(meshBase* mesh1, meshBase* mesh2);

class TransferDriver //: public NemDriver
{
  public:

    TransferDriver(std::string srcmsh, std::string trgmsh, std::string method,
                   std::string ofname, bool checkQuality);

    TransferDriver(std::string srcmsh, std::string trgmsh, std::string method,
                   std::vector<std::string> arrayNames, std::string ofname,
                   bool checkQuality);

    ~TransferDriver();
};

%extend TransferDriver {

    static TransferDriver* py_readJSON(std::string serialized_json, std::string ifname, bool serialized){
      if (serialized) {
        jsoncons::json inputjson = jsoncons::json::parse(serialized_json);
        return TransferDriver::readJSON(inputjson);
      }
      else if (!serialized) {
        return TransferDriver::readJSON(ifname);
      }
      return NULL;
    }

    %pythoncode %{

    @staticmethod
    def readJSON( json_obj):
      import json
      if type(json_obj) is list:
        serialized_json = json.dumps(json_obj[0])
        return TransferDriver.py_readJSON(serialized_json, '', True)
      elif type(json_obj) is dict:
        serialized_json = json.dumps(json_obj)
        return TransferDriver.py_readJSON(serialized_json, '', True)
      elif type(json_obj) is str:
        try:
            json.loads(json_obj)
            serialized_json = json_obj
            return TransferDriver.py_readJSON(serialized_json, '', True)
        except:
            return TransferDriver.py_readJSON('', json_obj, False)

    %}
};

class NemDriver
{
  public:
    NemDriver();
    virtual ~NemDriver();
    static NemDriver* readJSON(json inputjson);
};

%extend NemDriver {

    static NemDriver* py_readJSON(std::string serialized_json, std::string ifname, bool serialized){
      if (serialized) {
        jsoncons::json inputjson = jsoncons::json::parse(serialized_json);
        return NemDriver::readJSON(inputjson);
      }
      else if (!serialized) {
        return NemDriver::readJSON(ifname);
      }
      return NULL;
    }

    %pythoncode %{

    @staticmethod
    def readJSON( json_obj):
      import json
      if type(json_obj) is list:
        serialized_json = json.dumps(json_obj[0])
        return NemDriver.py_readJSON(serialized_json, '', True)
      elif type(json_obj) is dict:
        serialized_json = json.dumps(json_obj)
        return NemDriver.py_readJSON(serialized_json, '', True)
      elif type(json_obj) is str:
        try:
            json.loads(json_obj)
            serialized_json = json_obj
            return NemDriver.py_readJSON(serialized_json, '', True)
        except:
            return NemDriver.py_readJSON('', json_obj, False)

    %}
};

class RefineDriver : public NemDriver
{
  public:
    RefineDriver(std::string _mesh, std::string method, std::string arrayName,
                 double dev_mult, bool maxIsmin, double edgescale, std::string ofname, bool transferData, 
                 double sizeFactor = 1.);
    RefineDriver(std::string _mesh, std::string method,
                 double edgescale, std::string ofname, bool transferData);
    RefineDriver(std::string _mesh, std::string method, std::string arrayName, int order,
                 std::string ofname, bool transferData);

  ~RefineDriver();

};

%extend RefineDriver {

    static RefineDriver* py_readJSON(std::string serialized_json, std::string ifname, bool serialized){
      if (serialized) {
        jsoncons::json inputjson = jsoncons::json::parse(serialized_json);
        return RefineDriver::readJSON(inputjson);
      }
      else if (!serialized) {
        return RefineDriver::readJSON(ifname);
      }
      return NULL;
    }

    %pythoncode %{

    @staticmethod
    def readJSON( json_obj):
      import json
      if type(json_obj) is list:
        serialized_json = json.dumps(json_obj[0])
        return RefineDriver.py_readJSON(serialized_json, '', True)
      elif type(json_obj) is dict:
        serialized_json = json.dumps(json_obj)
        return RefineDriver.py_readJSON(serialized_json, '', True)
      elif type(json_obj) is str:
        try:
            json.loads(json_obj)
            serialized_json = json_obj
            return RefineDriver.py_readJSON(serialized_json, '', True)
        except:
            return RefineDriver.py_readJSON('', json_obj, False)

    %}
};


class MeshGenDriver : public NemDriver
{
  public:

    MeshGenDriver(std::string ifname, std::string meshEngine,
                  meshingParams* params, std::string ofname);

    static MeshGenDriver* readJSON(json inputjson);
};

%extend MeshGenDriver {

    static MeshGenDriver* py_readJSON(std::string serialized_json, std::string ifname, bool serialized){
      if (serialized) {
        jsoncons::json inputjson = jsoncons::json::parse(serialized_json);
        return MeshGenDriver::readJSON(inputjson);
      }
      return NULL;
    }

    %pythoncode %{

    @staticmethod
    def readJSON( json_obj):
      import json
      if type(json_obj) is list:
        serialized_json = json.dumps(json_obj[0])
        return MeshGenDriver.py_readJSON(serialized_json, '', True)
      elif type(json_obj) is dict:
        serialized_json = json.dumps(json_obj)
        return MeshGenDriver.py_readJSON(serialized_json, '', True)
      elif type(json_obj) is str:
        try:
            json.loads(json_obj)
            serialized_json = json_obj
        except:
            with open(json_obj, 'r') as jsonfile:
                serialized_json = json.dumps(json.load(jsonfile))
        return MeshGenDriver.py_readJSON(serialized_json, '', True)

    %}
};


class MeshQualityDriver : public NemDriver
{
  public:

    MeshQualityDriver(std::string _mesh, std::string ofname);
    ~MeshQualityDriver();

    static MeshQualityDriver* readJSON(json inputjson);
};


%extend MeshQualityDriver {

    static MeshQualityDriver* py_readJSON(std::string serialized_json, std::string ifname, bool serialized){
      if (serialized) {
        jsoncons::json inputjson = jsoncons::json::parse(serialized_json);
        return MeshQualityDriver::readJSON(inputjson);
      }
      else if (!serialized) {
        return MeshQualityDriver::readJSON(ifname);
      }
      return NULL;
    }

    %pythoncode %{

    @staticmethod
    def readJSON( json_obj):
      import json
      if type(json_obj) is list:
        serialized_json = json.dumps(json_obj[0])
        return MeshQualityDriver.py_readJSON(serialized_json, '', True)
      elif type(json_obj) is dict:
        serialized_json = json.dumps(json_obj)
        return MeshQualityDriver.py_readJSON(serialized_json, '', True)
      elif type(json_obj) is str:
        try:
            json.loads(json_obj)
            serialized_json = json_obj
            return MeshQualityDriver.py_readJSON(serialized_json, '', True)
        except:
            return MeshQualityDriver.py_readJSON('', json_obj, False)

    %}
};


class ConversionDriver : public NemDriver
{
  public:

    ConversionDriver(std::string srcmsh, std::string trgmsh, std::string method,
                   std::string ofname, json inputjson);

    ~ConversionDriver();
};

%extend ConversionDriver {

    static ConversionDriver* py_readJSON(std::string serialized_json, std::string ifname, bool serialized){
      if (serialized) {
        jsoncons::json inputjson = jsoncons::json::parse(serialized_json);
        return ConversionDriver::readJSON(inputjson);
      }
      else if (!serialized) {
        return ConversionDriver::readJSON(ifname);
      }
      return NULL;
    }

    %pythoncode %{

    @staticmethod
    def readJSON( json_obj):
      import json
      if type(json_obj) is list:
        serialized_json = json.dumps(json_obj[0])
        return ConversionDriver.py_readJSON(serialized_json, '', True)
      elif type(json_obj) is dict:
        serialized_json = json.dumps(json_obj)
        return ConversionDriver.py_readJSON(serialized_json, '', True)
      elif type(json_obj) is str:
        try:
            json.loads(json_obj)
            serialized_json = json_obj
            return ConversionDriver.py_readJSON(serialized_json, '', True)
        except:
            return ConversionDriver.py_readJSON('', json_obj, False)

    %}
};

class RemeshDriver : public NemDriver
{
  public:
    RemeshDriver(const std::vector<std::string>& fluidNames,
                 const std::vector<std::string>& ifluidniNames,
                 const std::vector<std::string>& ifluidnbNames,
                 const std::vector<std::string>& ifluidbNames,
                 const json& remeshjson);
    ~RemeshDriver();
};

%extend RemeshDriver {

    static RemeshDriver* py_readJSON(std::string serialized_json, std::string ifname, bool serialized){
      if (serialized) {
        jsoncons::json inputjson = jsoncons::json::parse(serialized_json);
        return RemeshDriver::readJSON(inputjson);
      }
      else if (!serialized) {
        return RemeshDriver::readJSON(ifname);
      }
      return NULL;
    }

    %pythoncode %{

    @staticmethod
    def readJSON( json_obj):
      import json
      if type(json_obj) is list:
        serialized_json = json.dumps(json_obj[0])
        return RemeshDriver.py_readJSON(serialized_json, '', True)
      elif type(json_obj) is dict:
        serialized_json = json.dumps(json_obj)
        return RemeshDriver.py_readJSON(serialized_json, '', True)
      elif type(json_obj) is str:
        try:
            json.loads(json_obj)
            serialized_json = json_obj
            return RemeshDriver.py_readJSON(serialized_json, '', True)
        except:
            return RemeshDriver.py_readJSON('', json_obj, False)

    %}
};

class RocRestartDriver : public NemDriver
{

  public:
    RocRestartDriver(const std::vector<std::string>& fluidNamesRm,
                     const std::vector<std::string>& ifluidniNamesRm,
                     const std::vector<std::string>& ifluidnbNamesRm,
                     const std::vector<std::string>& ifluidbNamesRm,
                     const std::vector<std::string>& fluidNamesLts,
                     const std::vector<std::string>& ifluidniNamesLts,  
                     const std::vector<std::string>& ifluidnbNamesLts,   
                     const std::vector<std::string>& ifluidbNamesLts); 
 
    ~RocRestartDriver(); 
};

%extend RocRestartDriver {

    static RocRestartDriver* py_readJSON(std::string serialized_json, std::string ifname, bool serialized){
      if (serialized) {
        jsoncons::json inputjson = jsoncons::json::parse(serialized_json);
        return RocRestartDriver::readJSON(inputjson);
      }
      else if (!serialized) {
        return RocRestartDriver::readJSON(ifname);
      }
      return NULL;
    }

    %pythoncode %{

    @staticmethod
    def readJSON( json_obj):
      import json
      if type(json_obj) is list:
        serialized_json = json.dumps(json_obj[0])
        return RocRestartDriver.py_readJSON(serialized_json, '', True)
      elif type(json_obj) is dict:
        serialized_json = json.dumps(json_obj)
        return RocRestartDriver.py_readJSON(serialized_json, '', True)
      elif type(json_obj) is str:
        try:
            json.loads(json_obj)
            serialized_json = json_obj
            return RocRestartDriver.py_readJSON(serialized_json, '', True)
        except:
            return RocRestartDriver.py_readJSON('', json_obj, False)

    %}
};

class OrderOfAccuracy
{

  public: 
    OrderOfAccuracy(meshBase* _f1, meshBase* _f2, meshBase* _f3,
                    const std::vector<int>& _arrayIDs)
      : f1(_f1), f2(_f2), f3(_f3), arrayIDs(_arrayIDs);
    ~OrderOfAccuracy(); 
    std::vector<std::vector<double>> computeOrderOfAccuracy(); 
    std::vector<std::vector<double>> computeGCI_21();
    std::vector<std::vector<double>> computeGCI_32();
    std::vector<std::vector<double>> computeResolution(double gciStar);
    std::vector<std::vector<double>> getOrderOfAccuracy();
    std::vector<std::vector<double>> checkAsymptoticRange();
    std::vector<std::vector<double>> 
    computeDiff(meshBase* mesh, const std::vector<std::string>& newArrNames);
    void computeRichardsonExtrapolation();
    void computeMeshWithResolution(double gciStar, const std::string& ofname);
    std::vector<std::vector<double>> computeDiffF3F1();
    std::vector<std::vector<double>> getDiffF2F1();
  
  private:
    meshBase *f1, *f2, *f3;
    const std::vector<int> arrayIDs;
    std::vector<int> diffIDs;
    std::vector<int> relEIDs;
    std::vector<int> realDiffIDs;
    std::vector<std::string> f3ArrNames, f2ArrNames;    
    double r21, r32;
    std::vector<std::vector<double>> diffF3F2;
    std::vector<std::vector<double>> diffF2F1;  
    std::vector<std::vector<double>> GCI_32;
    std::vector<std::vector<double>> GCI_21;
    std::vector<std::vector<double>> orderOfAccuracy; 
};

class meshingParams
{
  public:
    meshingParams();
    virtual ~meshingParams();    
    
};

class meshGen 
{
  public:
    meshGen();
    virtual ~meshGen();
    
    // creates generator with default parameters
    static meshGen* Create(std::string fname, std::string meshEngine);
    // creates generater with specified parameters
    static meshGen* Create(std::string fname, std::string meshEngine, meshingParams* params);
    virtual int createMeshFromSTL(const char* fname) = 0;
    vtkSmartPointer<vtkDataSet> getDataSet(); 
};

class symmxParams : public meshingParams
{
  
  public:
    symmxParams();
    ~symmxParams();
    std::string logFName;
    std::string features;
    std::string licFName;
    double meshSize;
    double anisoMeshCurv;
    double minCurvSize;
    double glbSizeGradRate;
    double surfMshImprovGradRate;  
    double surfMshImprovMinSize; 
};

class symmxGen : public meshGen
{

  public:
    symmxGen();
    symmxGen(symmxParams* params); 
    ~symmxGen();
    void createMeshFromModel(const char* mdlFName);
    int createModelFromSTL(const char* stlFName);
    int createSurfaceMeshFromSTL(const char* stlFName);
    int createVolumeMeshFromSTL(const char* stlFName);
    int createMeshFromSTL(const char* fname); 
    void convertToVTU(); 
    void saveMesh(const std::string& mshFName);
    void setWriteSurfAndVol(bool b);
};

class meshStitcher
{
  public:
    
    meshStitcher(const std::vector<std::string>& cgFileNames, bool surf);
    ~meshStitcher();
    
    cgnsAnalyzer* getStitchedCGNS();
    meshBase* getStitchedMB();
};

class meshPartitioner
{

public:
  meshPartitioner(int nNde, int nElm, std::vector<int>& elemConn, MeshType_t meshType);
  meshPartitioner(cgnsAnalyzer* inCg);
  meshPartitioner(meshBase* inMB);
  // destructor
  ~meshPartitioner();

  // mesh information
  int partition(int nPartition);
  int partition();
  std::vector<double> getCrds(int iPart, std::vector<double> crds);
  std::vector<int> getConns(int iPart);
  int getNNdePart(int iPart);
  int getNElmPart(int iPart);
};
