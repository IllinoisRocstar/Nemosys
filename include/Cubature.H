#ifndef CUBATURE_H 
#define CUBATURE_H
#include <meshBase.H>
#include <vtkMesh.H>
#include <vtkPolyData.h>
#include <vtkGenericCell.h>
#include <vtkQuadratureSchemeDefinition.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <Eigen/Dense>

// 3 point symmetric gaussian quadrature for triangle element
extern double TRI3 [];
extern double TRI3W [];
// 4 point symmetric gaussian quadrature for tetrahedral element
extern double TET4 [];
extern double TET4W [];

using Eigen::Matrix3d;
using Eigen::VectorXd;

// pair type for coordinate and data there
//typedef std::pair<std::vector<double>,std::vector<std::vector<double>>> pntDataPair;
typedef std::pair<std::vector<double>,std::vector<double>> pntDataPair;
//typedef std::pair<std::vector<double>, std::vector<VectorXd>> pntDataPair;

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
    // returns coordinates of gauss points and associated data at cell
    pntDataPairVec getGaussPointsAndDataAtCell(int cellID);
    // get interpolated values at gauss points for arrays specified in arrayIDs
    void interpolateToGaussPoints();
    // get interpolated values at gauss points for arrays specified by name
    void interpolateToGaussPoints(const std::vector<std::string>& newArrayNames);
    // integrate arrays specified in arrayIDs over all cells
    // the array returned contains the sum of the integral in each cell
    // taking the square root of its components will give you the L2-norm 
    std::vector<std::vector<double>> integrateOverAllCells();
    // in each cell, integrate arrays specified by name, normalize them by volume 
    // and take the sqrt. the returned array contains the sum over all cells
    // this function is meant to be used by the patch recovery class 
    std::vector<std::vector<double>> 
      integrateOverAllCells(const std::vector<std::string>& newArrayNames, bool computeRMSE);
    void setArrayIDs(const std::vector<int>& _arrayIDs) { arrayIDs = _arrayIDs; }
    double computeJacobian(vtkSmartPointer<vtkGenericCell> genCell, int cellType);
    double computeCellVolume(vtkSmartPointer<vtkGenericCell> genCell, int cellType);

  // access
  public:
    // returns pointer to mesh object
    meshBase* getNodeMesh() { return nodeMesh; }
    vtkSmartPointer<vtkPolyData> getGaussMesh() { return gaussMesh; } 
    vtkQuadratureSchemeDefinition** getDict() { return dict; }
    std::vector<int> getNumComponents(){ return numComponents; }    
    std::vector<int> getArrayIDs() { return arrayIDs; }
    int getTotalComponents() { return totalComponents; }

    // writes a vtp (polydata) of gauss points 
    // interpolated values also written if interpolateToGaussPoints(...) has been called
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

  private:
    // we want gauss points of this mesh
    meshBase* nodeMesh;
    // number of volume elements in nodemesh
    int numVolCells;
    // we put the gauss points and interpolated data into this mesh
    vtkSmartPointer<vtkPolyData> gaussMesh;
    // dictionary relating cell type quadrature scheme info for that type
    vtkQuadratureSchemeDefinition** dict;   
    // array ids of data to be interpolated
    std::vector<int> arrayIDs;
    std::vector<int> numComponents;
    int totalComponents;

    // get offset from nodeMesh for lookup of gauss points in gaussMesh
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
    // integrates data specified by newArrayNames over cell, 
    // normalizes the integral by the cell volume and takes the sqrt
    void integrateOverCell(int cellID,
                           vtkSmartPointer<vtkGenericCell> genCell,
                           vtkSmartPointer<vtkPointData> pd,
                           std::vector<vtkSmartPointer<vtkDoubleArray>>& integralData,
                           std::vector<std::vector<double>>& totalIntegralData,
                           const std::vector<std::string>& newArrayNames, bool normalizeByVol);

    // disabled
    GaussCubature(const GaussCubature& that) = delete;
    GaussCubature& operator=(const GaussCubature& that) = delete; 

};



#endif
