%module pyNemosys
%include "std_string.i"
%include "std_vector.i"
%include "std_map.i"
%{
#include "nemosys_export.h"
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


%include "NemDriver.H"

%extend NemDriver {

    static NemDriver *py_readJSON(const std::string &serialized_json,
                                  const std::string &ifname, bool serialized) {
      if (serialized) {
        jsoncons::json inputjson = jsoncons::json::parse(serialized_json);
        return NemDriver::readJSON(inputjson);
      } else {
        return NemDriver::readJSON(ifname);
      }
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


%include "TransferDriver.H"

%extend TransferDriver {

    static TransferDriver *py_readJSON(const std::string &serialized_json,
                                       const std::string &ifname,
                                       bool serialized) {
      if (serialized) {
        jsoncons::json inputjson = jsoncons::json::parse(serialized_json);
        return TransferDriver::readJSON(inputjson);
      } else {
        return TransferDriver::readJSON(ifname);
      }
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


%include "RefineDriver.H"

%extend RefineDriver {

    static RefineDriver *py_readJSON(const std::string &serialized_json,
                                     const std::string &ifname,
                                     bool serialized) {
      if (serialized) {
        jsoncons::json inputjson = jsoncons::json::parse(serialized_json);
        return RefineDriver::readJSON(inputjson);
      } else {
        return RefineDriver::readJSON(ifname);
      }
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


%include "MeshGenDriver.H"

%extend MeshGenDriver {

    static MeshGenDriver *py_readJSON(const std::string &serialized_json,
                                      const std::string &ifname,
                                      bool serialized) {
      if (serialized) {
        jsoncons::json inputjson = jsoncons::json::parse(serialized_json);
        return MeshGenDriver::readJSON(inputjson);
      } else {
        return nullptr;
      }
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


%include "MeshQualityDriver.H"

%extend MeshQualityDriver {

    static MeshQualityDriver *py_readJSON(const std::string &serialized_json,
                                          const std::string &ifname,
                                          bool serialized) {
      if (serialized) {
        jsoncons::json inputjson = jsoncons::json::parse(serialized_json);
        return MeshQualityDriver::readJSON(inputjson);
      } else {
        return nullptr;
      }
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
            except:
                with open(json_obj, 'r') as jsonfile:
                    serialized_json = json.dumps(json.load(jsonfile))
            return MeshQualityDriver.py_readJSON(serialized_json, '', True)

    %}
};


%include "ConversionDriver.H"

%extend ConversionDriver {

    static ConversionDriver *py_readJSON(const std::string &serialized_json,
                                         const std::string &ifname,
                                         bool serialized) {
      if (serialized) {
        jsoncons::json inputjson = jsoncons::json::parse(serialized_json);
        return ConversionDriver::readJSON(inputjson);
      } else {
        return ConversionDriver::readJSON(ifname);
      }
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


%include "OrderOfAccuracy.H"

%include "meshingParams.H"

%include "meshGen.H"
