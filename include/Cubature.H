#ifndef CUBATURE_H
#define CUBATURE_H
#include <meshBase.H>
#include <vtkMesh.H>
#include <Eigen/Dense>

// 3 point symmetric gaussian quadrature for triangle element
extern double TRI3 [];
extern double TRI3W [];
// 4 point symmetric gaussian quadrature for tetrahedral element
extern double TET4 [];
extern double TET4W [];

// pair type for coordinate and data there
typedef std::pair<std::vector<double>,std::vector<std::vector<double>>> pntDataPair;

// holds gauss points and data at these points as pairs
typedef std::vector<pntDataPair> pntDataPairVec;

class GaussCubature
{
  public:
    //GaussCubature(meshBase* _nodeMesh):nodeMesh(_nodeMesh){ buildMap();} 
    GaussCubature(meshBase* _nodeMesh);
    GaussCubature(meshBase* _nodeMesh, const std::vector<int>& arrayIDs);
    ~GaussCubature(){ delete [] dict;}
    // builds a polydata of gauss points
    void constructGaussMesh();

  public:
    // returns coordinates of gauss points at cell
    pntDataPairVec getGaussPointsAndDataAtCell(int cellID);
    std::vector<std::vector<double>> integrateOverAllCells();
    // get interpolated values at gauss points
    void interpolateToGaussPoints();
  // access
  public:
    // returns pointer to mesh object
    meshBase* getNodeMesh() { return nodeMesh; }

    // writes a vtp (polydata) of gauss points 
    // interpolated values also written if interpolateToGaussPoitns(...) has been called
    void writeGaussMesh(const char* name);

  // factory constructors
  public:
    static GaussCubature* Create(meshBase* nodeMesh);
    static GaussCubature* Create(meshBase* nodeMesh, const std::vector<int>& arrayIDs);
    static std::unique_ptr<GaussCubature> CreateUnique(meshBase* nodeMesh);
    static std::unique_ptr<GaussCubature>
      CreateUnique(meshBase* nodeMesh, const std::vector<int>& arrayIDs);
    static std::shared_ptr<GaussCubature>
      CreateShared(meshBase* nodeMesh);
    static std::shared_ptr<GaussCubature>
      CreateShared(meshBase* nodeMesh, const std::vector<int>& arrayIDs);

    void setArrayIDs(const std::vector<int>& _arrayIDs) { arrayIDs = _arrayIDs; }
    double computeScaledJacobian(vtkSmartPointer<vtkGenericCell>genCell, int cellType);

  private:
    // we want gauss points of this mesh
    meshBase* nodeMesh;
    // we put the gauss points and interpolated data into this mesh
    vtkSmartPointer<vtkPolyData> gaussMesh;
    // dictionary relating cell type quadrature scheme info for that type
    vtkQuadratureSchemeDefinition** dict;   
    // array ids of data to be interpolated
    std::vector<int> arrayIDs;
    std::vector<int> numComponents;

    // get offset from nodeMesh for lookup of gauss points in polydata
    int getOffset(int cellID);
    // interpolates provided data (das) to gauss points in cell
    int interpolateToGaussPointsAtCell(const int cellID,
                                       vtkSmartPointer<vtkGenericCell> genCell,
                                       const std::vector<vtkSmartPointer<vtkDataArray>>& das,
                                       std::vector<vtkSmartPointer<vtkDoubleArray>>& daGausses);
    // integrates provided data over cell
    void integrateOverCell(int cellID,
                           vtkSmartPointer<vtkGenericCell> genCell,
                           vtkSmartPointer<vtkPointData> pd,
                           std::vector<vtkSmartPointer<vtkDoubleArray>>& integralData,
                           std::vector<std::vector<double>>& totalIntegralData);



    // disabled
    GaussCubature(const GaussCubature& that) = delete;
    GaussCubature& operator=(const GaussCubature& that) = delete; 

};



#endif
