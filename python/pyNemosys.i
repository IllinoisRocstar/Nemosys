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
#include "RichardsonExtrapolation.H"
#include "jsoncons/json.hpp"
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
    static std::shared_ptr<meshBase> CreateShared(std::string fname);
    //static std::unique_ptr<meshBase> CreateUnique(std::string fname);
    virtual std::vector<double> getPoint(int id);
    virtual std::map<int, std::vector<double>> getCell(int id);
    virtual std::vector<std::vector<double>> getCellVec(int id);
    vtkSmartPointer<vtkDataSet> getDataSet();
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
    virtual void classifyAndAddBoundaries();
    vtkSmartPointer<vtkCellLocator> buildLocator();
    int transfer(meshBase* target, std::string method,
                 const std::vector<int>& arrayIDs);
    virtual int getCellType() const = 0;
    int transfer(meshBase* target, std::string method,
                 const std::vector<std::string>& arrayNames);
    int transfer(meshBase* target, std::string method);
    std::vector<std::vector<double>> integrateOverMesh(const std::vector<int>& arrayIDs);
    
    void generateSizeField(std::string method, int arrayID, double dev_mlt, bool maxIsmin);

    void setSFBool(bool q);
    bool getSFBool();
    int IsArrayName(std::string name);
    void setOrder(int _order); 
    int getOrder() const; 

    void refineMesh(std::string method, int arrayID,
                    double dev_mult, bool maxIsmin, double edge_scale, std::string ofname, bool transferData);
    void refineMesh(std::string method, std::string arrayName,
                    double dev_mult, bool maxIsmin, double edge_scale, std::string ofname, bool transferData);
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
    void setNewArrayNames(const std::vector<std::string>& newnames);
};


class vtkMesh : public meshBase
{
  // constructor and destructor
  public:
    vtkMesh();
    vtkMesh(const char* fname);
    vtkMesh(const char* fname1, const char* fname2);
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
                 double dev_mult, bool maxIsmin, double edgescale, std::string ofname, bool transferData);
    RefineDriver(std::string _mesh, std::string method,
                 double edgescale, std::string ofname, bool transferData);

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

<<<<<<< HEAD


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
class RichardsonExtrapolation
{

  public:
    RichardsonExtrapolation(meshBase* _fineMesh, meshBase* coarseMesh,
                            double _ref_factor, int _order, 
                            const std::vector<int>& _arrayIDs)
      : fineMesh(_fineMesh), ref_factor(_ref_factor), order(_order),
        arrayIDs(_arrayIDs);

    std::vector<std::vector<double>> computeDiscretizationError();
    std::vector<double> computeObservedOrderOfAccuracy(meshBase* finerMesh);
  private:
    meshBase* fineMesh;
    double ref_factor;
    int order;
    const std::vetor<int> arrayIDs;
    std::vector<std::string> newArrNames;
};
