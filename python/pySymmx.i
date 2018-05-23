%module pySymmx
%include "std_string.i"
%include "std_vector.i"
%include "std_map.i"
%{
#include "symmxGen.H"
%}

%template(vectorString) std::vector<std::string>;
%template(doubleV) std::vector<double>;
%template(intV) std::vector<int>;
%template(doubleVV) std::vector<std::vector<double>>;
%template(cellMap) std::map<int, std::vector<double>>;

class meshSymmx
{

  public:
    meshSymmx(char* logFName="NONE", char* features="geomsim_core,meshsim_surface,meshsim_volume", 
            const char* licFName="simmodsuite.lic");
    ~meshSymmx();

    // symmetrix mesh creation
    // create mesh from symmetrix model file
    void createMeshFromModel(char* mdlFName);
    int createModelFromSTL(char* stlFName);
    int createSurfaceMeshFromSTL(char* stlFName);
    int createVolumeMeshFromSTL(char* stlFName);
    void convertToVTU();
    void saveMesh(char* mshFName);
    vtkSmartPointer<vtkDataSet> getDataSet();
    void setWriteSurfAndVol(bool b);
};
