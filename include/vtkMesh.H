#ifndef vtkMesh_H
#define vtkMesh_H

#include <meshBase.H>

using namespace nemAux;

// read files from any vtk format (except legacy .vtk)
template<class TReader> vtkDataSet* ReadAnXMLFile(const char* fileName)
{
  vtkSmartPointer<TReader> reader =
    vtkSmartPointer<TReader>::New();
  reader->SetFileName(fileName);
  reader->Update();
  reader->GetOutput()->Register(reader);
  return vtkDataSet::SafeDownCast(reader->GetOutput());
}

// read an ASCII legacy vtk file containing an unstructured grid
// the available vtk readers do not suffice for this purpose if the file
// contains more than one point or cell data array. This function has mainlyl been
// tested on vtk output from MFEM
vtkSmartPointer<vtkUnstructuredGrid> ReadALegacyVTKFile(const char* fileName);

// helpers for reading legacy vtk 

// helper that reads/checks header from legacy vtk
bool readLegacyVTKHeader(const std::string& line);

// helper that reads field data from legacy vtk
bool readLegacyVTKFieldData(std::istream& meshStream, std::string& line, 
                            vtkSmartPointer<vtkUnstructuredGrid> dataSet_tmp);
// helper that reads points from legacy vtk
bool readLegacyVTKPoints(std::istream& meshStream, std::string& line, int& numPoints,
                         vtkSmartPointer<vtkPoints> points,
                         vtkSmartPointer<vtkUnstructuredGrid> dataSet_tmp);
// helper that reads cells from legacy vtk
bool readLegacyVTKCells(std::istream& meshStream, std::string& line, int& numCells, 
                        std::vector<vtkSmartPointer<vtkIdList>>& vtkCellIds,
                        vtkSmartPointer<vtkUnstructuredGrid> dataSet_tmp);
// helper that reads point and cell data from legacy vtk
bool readLegacyVTKData(std::ifstream& meshStream, std::string& line, const int numTuple, 
                       const bool pointOrCell, bool& hasPointOrCell,
                       vtkSmartPointer<vtkUnstructuredGrid> dataSet_tmp);
// helper that casts data arrays to type specified in legacy file (not used rn)
void addLegacyVTKData(vtkDataArray* arr, const std::string& type, bool pointOrCell, 
                      vtkSmartPointer<vtkUnstructuredGrid> dataSet_tmp);

template<class TWriter> void writeVTFile(std::string fname, vtkSmartPointer<vtkDataSet> dataSet)
{
  vtkSmartPointer<TWriter> Writer = vtkSmartPointer<TWriter>::New();
  Writer->SetFileName(&fname[0u]);
  Writer->SetInputData(dataSet);
  Writer->Write();
}

class vtkMesh : public meshBase
{
  // constructor and destructor
  public:
    vtkMesh(){};
    vtkMesh(const char* fname);
    ~vtkMesh()
    {
      std::cout << "vtkMesh destroyed" << std::endl;
    } 

  // access
  public:
    // get point with id
    std::vector<double> getPoint(int id);
    // get cell with id : returns point indices and respective coordinates
    std::map<int, std::vector<double>> getCell(int id);
    std::vector<std::vector<double>> getCellVec(int id);
    // get diameter of circumsphere of each cell
    std::vector<double> getCellLengths();
    // get center of a cell
    std::vector<double> getCellCenter(int cellID);
 
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



#endif
