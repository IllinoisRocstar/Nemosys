#ifndef CUBATURE_H
#define CUBATURE_H
#include <meshBase.H>
#include <vtkMesh.H>

// 3 point symmetric gaussian quadrature for triangle element
extern double TRI3 [];
extern double TRI3W;
// 4 point symmetric gaussian quadrature for tetrahedral element
extern double TET4 [];
extern double TET4W;

class GaussCubature
{
  public:
    GaussCubature(meshBase* _nodeMesh):nodeMesh(_nodeMesh){ buildMap();} 

    ~GaussCubature(){}

  public:
    // returns coordinates of gauss points at cell
    std::vector<std::vector<double>> getGaussPointsAtCell(int cellID);
    // interpolates provided data (das) to gauss points in cell
    int interpolateToGaussPointsAtCell(const int cellID,
                                       vtkSmartPointer<vtkGenericCell> genCell,
                                       const std::vector<vtkSmartPointer<vtkDataArray>>& das,
                                       std::vector<vtkSmartPointer<vtkDoubleArray>>& daGausses,
                                       const std::vector<int>& numComponents, 
                                       const int polyPnt);
    // integrates provided data over cell
    void integrateOverCell(int cellID,const std::vector<int>& arrayIDs);

  public:
    // returns pointer to mesh object
    meshBase* getNodeMesh() { return nodeMesh; }
    // for testing; writes a vtp (polydata) of gauss points and interpolated values
    void constructGaussMesh(const std::vector<int>& arrayIDs);
  

  private:
    // we want gauss points of this mesh
    meshBase* nodeMesh;
    // map relating cell type and number of gauss points for that type
    std::map<int,int> nGaussForCellTMap;
    // used to generate the above map
    int getNumGaussPointsForCellType(int cellType);
    void buildMap();
    // used by constructGaussMesh to get interpolated values at gauss points
    void interpolateToGaussPoints(vtkSmartPointer<vtkPolyData> gaussMesh,
                                  const std::vector<int>& arrayIDs);
    double computeVolume(vtkSmartPointer<vtkGenericCell> genCell, const int cellType); 
    // disabled
    GaussCubature(const GaussCubature& that) = delete;
    GaussCubature& operator=(const GaussCubature& that) = delete; 
};

#endif
