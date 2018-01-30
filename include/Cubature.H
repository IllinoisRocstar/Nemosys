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
    GaussCubature(meshBase* _nodeMesh);  

    ~GaussCubature()
    { 
      delete [] dict;
    }

  public:
    // we don't use the transfer class here because the gauss points in
    // the polydata are the same order of traversal as the cells in nodeMesh
    void interpolateToGaussPoints(const std::vector<int>& arrayIDs);
    void integrateOverCells(const std::vector<int>& arrayIDs);

  public:
    void writeGaussMesh() { writeVTFile<vtkXMLPolyDataWriter> ("gaussTest.vtp",gaussMesh); }  

  private:
    meshBase* nodeMesh;
    vtkSmartPointer<vtkPolyData> gaussMesh;
    vtkQuadratureSchemeDefinition** dict; 
    void constructGaussMesh();
    GaussCubature(const GaussCubature& that) = delete;
    GaussCubature& operator=(const GaussCubature& that) = delete; 

};

#endif
