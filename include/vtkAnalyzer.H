#ifndef __VTKANALYZER_H
#define __VTKANALYZER_H

#include "nemosys_export.h"

#include "spheres.H"

// VTK
#include <vtksys/SystemTools.hxx>
#include <vtkSmartPointer.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>

// VTK types
#include <vtkUnstructuredGrid.h>
#include <vtkRectilinearGrid.h>
#include <vtkStructuredGrid.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkFieldData.h>
#include <vtkPolyData.h>
#include <vtkImageData.h>
#include <vtkCell.h>
#include <vtkCellTypes.h>
#include <vtkCellCenters.h>

// VTK Reader
#include <vtkXMLReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkXMLRectilinearGridReader.h>
#include <vtkXMLCompositeDataReader.h>
#include <vtkXMLImageDataReader.h>
#include <vtkDataSetReader.h>

// VTK Writer
#include <vtkXMLWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkXMLRectilinearGridWriter.h>
#include <vtkXMLCompositeDataWriter.h>

// others
#include <vector>
#include <map>


template<class TReader>
TReader *ReadAnXMLFile(const char *fileName)
{
  vtkSmartPointer<TReader> reader =
      vtkSmartPointer<TReader>::New();
  reader->SetFileName(fileName);
  reader->Update();
  reader->GetOutput()->Register(reader);
  return reader;
  //return vtkDataSet::SafeDownCast(reader->GetOutput());
}

template<class TWriter, class TObject>
void WriteAnXMLFile(const char *fileName, TObject Grid)
{
  vtkSmartPointer<TWriter> writer =
      vtkSmartPointer<TWriter>::New();
  writer->SetFileName(fileName);
#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(Grid);
#else
  writer->SetInputData(Grid);
#endif
  writer->Write();
}

class NEMOSYS_EXPORT vtkAnalyzer
{
    typedef std::map<int, int> CellContainer;

  public:
    int numberOfPoints;
    int numberOfCells;
    int numberOfPointData;
    int numberOfCellData;

  public:
    // constructor stuff
    vtkAnalyzer(char *fName) :
        numberOfPoints(0), numberOfCells(0),
        numberOfPointData(0), numberOfCellData(0), xmlFileName(fName),
        dataSet(nullptr), pointData(nullptr), cellData(nullptr) {}

    vtkAnalyzer(vtkDataSet *_dataSet, char *fName) :
        numberOfPointData(0), numberOfCellData(0), xmlFileName(fName),
        dataSet(_dataSet), pointData(nullptr), cellData(nullptr)
  {
      numberOfPoints = dataSet->GetNumberOfPoints();
      numberOfCells = dataSet->GetNumberOfCells();
    };

    ~vtkAnalyzer()
    {
      dataSet->Delete();
    };

    // read file and populate the class
    void read();
    void write(std::string outXMLFileName);

  // inquiry
  public:
    // number of points
    int getNumberOfPoints();
    // number of cells
    int getNumberOfCells();
    // number of point data
    int getNumberOfPointData();
    // number of cell data
    int getNumberOfCellData();

    vtkDataSet *getDataSet() const { return dataSet; }

    std::multimap<int, std::vector<int>> findBoundaryFaces();

    //get given point ID coordinates
    double *getPointCoords(int pntId);
    // returns coordinates of member points in cell with ID
    std::vector<double *> getCellCoords(int pntId, int &numComponent);
    // get number of nonTri
    int getNumberOfNonTri();
    std::vector<std::vector<double *>> getSurfaceTriElements(int &numComponent);
    void writeSurfaceTriElements(std::string fname);
    // returns coordinates of centers of all cells
    std::vector<double> getCellCenters(int &numComponent);
    // returns cells sharing point (topological search)
    std::vector<int> getCellsWithPoint(int pnt);
    // returns min extent
    double getMinExtent(int nDim, const std::vector<double> &pntCrds);
    //return coordinates of all points in mesh
    std::vector<double> getAllPointCoords(int nDim);
    // check for named array in vtk
    int IsArrayName(std::string name);
    // access to data
    // get a point data array
    int getPointDataArray(int id, std::vector<std::vector<double>> &pntData,
                          int &numTuple, int &numComponent);

    // get a cell data array
    int getCellDataArray(int id, std::vector<std::vector<double> > &cllData,
                         int &numTuple, int &numComponent);

    std::vector<std::vector<double> >
    getInterpData(int nDim, int num_neighbors, int numComponent, int numTuple,
                  std::vector<std::vector<double>> &volDataMat,
                  std::vector<double> &PlaneCellCenters,
                  std::vector<double> &VolPointCoords, double tol);

    // consider inclusions in interpolation
    std::vector<std::vector<double>>
    getInterpData(int nDim, int num_neighbors, int numComponent, int numTuple,
                  std::vector<std::vector<double>> &volDataMat,
                  std::vector<double> &PlaneCellCenters,
                  std::vector<double> &VolPointCoords,
                  std::vector<sphere> &spheres,
                  std::vector<double> &maskData, double tol);

    // add data
    // add point data array
    void setPointDataArray(const char *name, int numComponent,
                           std::vector<double> &pntArrData);
    // add a cell data array
    void setCellDataArray(const char *name, int numComponent,
                          std::vector<double> &pntArrData);

    void setCellDataArray(const char *name, int numComponent, int cellType,
                          std::vector<double> &pntArrData);
    // write interpolated data
    // no spheres, don't write coords
    void writeInterpData(const std::vector<std::vector<double>> &interpData,
                         double Mc, double M, double youngs_dom_default,
                         double poisson_dom_default, double T, double R,
                         std::ostream &outputStream, int numNonTri,
                         double conc_conv);
    // no spheres, switch to write coords
    void writeInterpData(const std::vector<std::vector<double>> &interpData,
                         double Mc, double M, double youngs_dom_default,
                         double poison_dom_default, double T, double R,
                         const std::vector<double> &PlaneCellCenters, int nDim,
                         std::ostream &outputStream, bool writeCoord,
                         int numNonTri,
                         double conc_conv);
    void writeInterpData(const std::vector<std::vector<double>> &interpData,
                         double Mc, double M, double youngs_dom_default,
                         double poisson_dom_default, double T, double R,
                         const std::vector<double> &PlaneCellCenters, int nDim,
                         std::string filename, bool writeCoord, int numNonTri,
                         double conc_conv);

    // write interpolated data if spheres
    // has spheres, don't write coords
    void writeInterpData(const std::vector<std::vector<double>> &interpData,
                         double Mc, double M, double youngs_dom_default,
                         double poisson_dom_default, double T, double R,
                         const std::vector<double> &PlaneCellCenters, int nDim,
                         std::vector<sphere> &spheres,
                         std::vector<std::string> &mat_sphere_names,
                         std::vector<std::string> &material_names,
                         std::vector<double> &youngs_inc_default,
                         std::vector<double> &shear_inc_default,
                         std::vector<double> &poisson_inc_default,
                         std::ostream &outputStream, int numNonTri,
                         double conc_conv);
    // has spheres, switch to write coords
    void writeInterpData(const std::vector<std::vector<double>> &interpData,
                         double Mc, double M, double youngs_dom_default,
                         double poisson_dom_default, double T, double R,
                         const std::vector<double> &PlaneCellCenters, int nDim,
                         std::vector<sphere> &spheres,
                         std::vector<std::string> &mat_sphere_names,
                         std::vector<std::string> &material_names,
                         std::vector<double> &youngs_inc_default,
                         std::vector<double> &shear_inc_default,
                         std::vector<double> &poisson_inc_default,
                         std::ostream &outputStream, bool writeCoord,
                         int numNonTri,
                         double conc_conv);

    void writeInterpData(const std::vector<std::vector<double>> &interpData,
                         double Mc, double M, double youngs_dom_default,
                         double poisson_dom_default, double T, double R,
                         const std::vector<double> &PlaneCellCenters, int nDim,
                         std::vector<sphere> &spheres,
                         std::vector<std::string> &mat_sphere_names,
                         std::vector<std::string> &material_names,
                         std::vector<double> &youngs_inc_default,
                         std::vector<double> &shear_inc_default,
                         std::vector<double> &poisson_inc_default,
                         std::string filename, bool writeCoord, int numNonTri,
                         double conc_conv);

    void writeMSH(std::string filename);
    void writeBackgroundMSH(std::string filename, double size);
    // utilities
    // print to output a report about the file
    void report();
    void writeCSV(char *fname, std::vector<double> slnVec);

    vtkSmartPointer<vtkDataSetReader> getReader() { return dataSetReader; }

  protected:
    char *xmlFileName;
    vtkSmartPointer<vtkXMLUnstructuredGridReader> unsGridReader;
    vtkSmartPointer<vtkDataSetReader> dataSetReader;
    // adding support for image file
    vtkSmartPointer<vtkXMLImageDataReader> imageDataSetReader;
    std::string extension;
    vtkDataSet *dataSet;
    vtkPointData *pointData;
    vtkCellData *cellData;
};

#endif
