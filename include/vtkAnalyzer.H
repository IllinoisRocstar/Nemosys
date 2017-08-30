#ifndef __VTKANALYZER_H
#define __VTKANALYZER_H

#include <vtkSmartPointer.h>
#include <vtkXMLReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkXMLRectilinearGridReader.h>
#include <vtkXMLHyperOctreeReader.h>
#include <vtkXMLCompositeDataReader.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkXMLImageDataReader.h>
#include <vtkDataSetReader.h>
#include <vtkDataSet.h>
#include <vtkUnstructuredGrid.h>
#include <vtkRectilinearGrid.h>
#include <vtkHyperOctree.h>
#include <vtkImageData.h>
#include <vtkPolyData.h>
#include <vtkStructuredGrid.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkFieldData.h>
#include <vtkCellTypes.h>
#include <vtksys/SystemTools.hxx>
#include <vtkDoubleArray.h>
#include <vtkXMLWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkXMLRectilinearGridWriter.h>
#include <vtkXMLHyperOctreeWriter.h>
#include <vtkXMLCompositeDataWriter.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkXMLImageDataReader.h>
#include <vtkDataSetReader.h>
#include <vtkCellCenters.h>
#include <vtkCell.h>
#include <baseInterp.H>
#include <spheres.H>
// others
#include <vector>
#include <iostream>
#include <ostream>
#include <fstream>
#include <map>
#include <string>
#include <cstdlib>

template<class TReader> TReader *ReadAnXMLFile(const char*fileName)
{
  vtkSmartPointer<TReader> reader =
    vtkSmartPointer<TReader>::New();
  reader->SetFileName(fileName);
  reader->Update();
  reader->GetOutput()->Register(reader);
  return reader;
  //return vtkDataSet::SafeDownCast(reader->GetOutput());
}

template<class TWriter, class TObject> void WriteAnXMLFile(const char*fileName, TObject Grid)
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

class vtkAnalyzer {

   typedef std::map<int,int> CellContainer;

public:
   int numberOfPoints;
   int numberOfCells;
   int numberOfPointData;
   int numberOfCellData;

public:
   // constructor stuff
   vtkAnalyzer(char* fName): 
     xmlFileName(fName), numberOfPoints(0), numberOfCells(0),
     numberOfPointData(0), numberOfCellData(0), dataSet(NULL),
     pointData(NULL), cellData(NULL)
   {

   };
   ~vtkAnalyzer()
     {};
 
   // read file and populate the class
   void read();
   void write(char* outXMLFileName);

   // inquiry 
   
	 // number of points
   int getNumberOfPoints();
   // number of cells
   int getNumberOfCells();
   // number of point data
   int getNumberOfPointData();
   // number of cell data
   int getNumberOfCellData();
   //get given point ID coordinates
   double* getPointCoords(int pntId);
	 // returns coordinates of member points in cell with ID 
	 std::vector<double* > getCellCoords(int pntId, int& numComponent);
	 // returns coordinates of centers of all cells 
	 std::vector<double > getCellCenters(int& numComponent, int nDim); 
	 //return coordinates of all points in mesh
	 std::vector<double> getAllPointCoords(int nDim); 
	 // check for named array in vtk
	 int IsArrayName(std::string name);
	 // access to data
   // get a point data array
   int getPointDataArray(int id, std::vector<std::vector<double> > &pntData, 
                         int &numTuple, int &numComponent);

   // get a cell data array
   int getCellDataArray(int id, std::vector<std::vector<double> > &cllData, 
                         int &numTuple, int &numComponent);

	 std::vector<std::vector<double> >
	 getInterpData(int nDim, int num_neighbors, int numComponent, int numTuple,     
   		                        std::vector<std::vector<double>> volDataMat,
      		                    std::vector<double> PlaneCellCenters);

	// consider inclusions in interpolation
	std::vector<std::vector<double>>
	getInterpData(int nDim, int num_neighbors, int numComponent, int numTuple,
                          std::vector<std::vector<double>> volDataMat,
                          std::vector<double> PlaneCellCenters,
                          std::vector<sphere> spheres);

 
  // add data
  // add point data array
  void setPointDataArray( const char* name, int numComponent,
                          std::vector<double> &pntArrData);
  // add point data array
  void setCellDataArray( const char* name, int numComponent,
                         std::vector<double> &pntArrData);

	// write interpolated data
	// no spheres, don't write coords
	void writeInterpData(std::vector<std::vector<double>> interpData,
											 double Mc, double M, double E, 
  	                   std::ostream& outputStream);
	// no spheres, switch to write coords
	void writeInterpData(std::vector<std::vector<double>> interpData,
    	                 double Mc, double M, double E,
											 std::vector<double> PlaneCellCenters, int nDim,
                       std::ostream& outputStream, bool writeCoord);
	void writeInterpData(std::vector<std::vector<double>> interpData,
											 double Mc, double M,	double E,
                       std::vector<double> PlaneCellCenters, int nDim,
                       std::string filename, bool writeCoord);

	// write interpolated data if spheres
	// has spheres, don't write coords
	void writeInterpData(std::vector<std::vector<double>> interpData,
											 double Mc, double M,	double E,
                       std::vector<double> PlaneCellCenters, int nDim,
                       std::vector<sphere> spheres,
                       std::ostream& outputStream);
	// has spheres, switch to write coords
	void writeInterpData(std::vector<std::vector<double>> interpData,
											 double Mc, double M,	double E,
                       std::vector<double> PlaneCellCenters, int nDim,
                       std::vector<sphere> spheres,
                       std::ostream& outputStream, bool writeCoord);
	void writeInterpData(std::vector<std::vector<double>> interpData,
											 double Mc, double M,	double E,
                       std::vector<double> PlaneCellCenters, int nDim,
                       std::vector<sphere> spheres,
                       std::string filename, bool writeCoord);

   // utilities
   // print to output a report about the file
   void report();
   void writeCSV(char* fname, std::vector<double> slnVec);

private:
   char* xmlFileName;
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
