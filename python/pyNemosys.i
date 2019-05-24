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
#include "OrderOfAccuracy.H"
#include "jsoncons/json.hpp"
#include "meshGen.H"
#include "meshingParams.H"
%}


%template(vectorString) std::vector<std::string>;
%template(doubleV) std::vector<double>;
%template(intV) std::vector<int>;
%template(doubleVV) std::vector<std::vector<double>>;
%template(cellMap) std::map<int, std::vector<double>>;
%template(meshBaseV) std::vector<meshBase *>;


%ignore CreateUnique;
%include "meshBase.H"


%include "vtkMesh.H"


int diffMesh(meshBase *mesh1, meshBase *mesh2);


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
      } else if (!serialized) {
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


class NemDriver {
  public:
    NemDriver();
    virtual ~NemDriver();
    static NemDriver *readJSON(json inputjson);
};

%extend NemDriver {

    static NemDriver* py_readJSON(std::string serialized_json, std::string ifname, bool serialized){
      if (serialized) {
        jsoncons::json inputjson = jsoncons::json::parse(serialized_json);
        return NemDriver::readJSON(inputjson);
      } else if (!serialized) {
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


class RefineDriver : public NemDriver {
  public:
    RefineDriver(std::string _mesh, std::string method, std::string arrayName,
                 double dev_mult, bool maxIsmin, double edgescale,
                 std::string ofname, bool transferData,
                 double sizeFactor = 1.);
    RefineDriver(std::string _mesh, std::string method,
                 double edgescale, std::string ofname, bool transferData);
    RefineDriver(std::string _mesh, std::string method, std::string arrayName,
                 int order,
                 std::string ofname, bool transferData);

    ~RefineDriver();
};

%extend RefineDriver {

    static RefineDriver* py_readJSON(std::string serialized_json, std::string ifname, bool serialized){
      if (serialized) {
        jsoncons::json inputjson = jsoncons::json::parse(serialized_json);
        return RefineDriver::readJSON(inputjson);
      } else if (!serialized) {
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


class MeshGenDriver : public NemDriver {
  public:
    MeshGenDriver(std::string ifname, std::string meshEngine,
                  meshingParams *params, std::string ofname);

    static MeshGenDriver *readJSON(json inputjson);
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


class MeshQualityDriver : public NemDriver {
  public:
    MeshQualityDriver(std::string _mesh, std::string ofname);
    ~MeshQualityDriver();

    static MeshQualityDriver *readJSON(json inputjson);
};


%extend MeshQualityDriver {

    static MeshQualityDriver* py_readJSON(std::string serialized_json, std::string ifname, bool serialized){
      if (serialized) {
        jsoncons::json inputjson = jsoncons::json::parse(serialized_json);
        return MeshQualityDriver::readJSON(inputjson);
      } else if (!serialized) {
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


class ConversionDriver : public NemDriver {
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
      } else if (!serialized) {
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


class OrderOfAccuracy {
  public:
    OrderOfAccuracy(meshBase *_f1, meshBase *_f2, meshBase *_f3,
                    const std::vector<int> &_arrayIDs)
        : f1(_f1), f2(_f2), f3(_f3), arrayIDs(_arrayIDs);
    ~OrderOfAccuracy();
    std::vector<std::vector<double>> computeOrderOfAccuracy();
    std::vector<std::vector<double>> computeGCI_21();
    std::vector<std::vector<double>> computeGCI_32();
    std::vector<std::vector<double>> computeResolution(double gciStar);
    std::vector<std::vector<double>> getOrderOfAccuracy();
    std::vector<std::vector<double>> checkAsymptoticRange();
    std::vector<std::vector<double>>
    computeDiff(meshBase *mesh, const std::vector<std::string> &newArrNames);
    void computeRichardsonExtrapolation();
    void computeMeshWithResolution(double gciStar, const std::string &ofname);
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


class meshingParams {
  public:
    meshingParams();
    virtual ~meshingParams();
};


class meshGen {
  public:
    meshGen();
    virtual ~meshGen();

    // creates generator with default parameters
    static meshGen *Create(std::string fname, std::string meshEngine);
    // creates generater with specified parameters
    static meshGen *
    Create(std::string fname, std::string meshEngine, meshingParams *params);
    virtual int createMeshFromSTL(const char *fname) = 0;
    vtkSmartPointer<vtkDataSet> getDataSet();
};
