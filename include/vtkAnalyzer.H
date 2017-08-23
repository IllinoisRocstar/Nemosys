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
// others
#include <vector>
#include <fstream>
#include <map>

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


   // access to data
   // get a point data array
   int getPointDataArray(int id, std::vector<std::vector<double> > &pntData, 
                         int &numTuple, int &numComponent);

   // get a cell data array
   int getCellDataArray(int id, std::vector<std::vector<double> > &cllData, 
                         int &numTuple, int &numComponent);

 
   // add data
   // add point data array
   void setPointDataArray( const char* name, int numComponent,
                          std::vector<double> &pntArrData);
   // add point data array
   void setCellDataArray( const char* name, int numComponent,
                          std::vector<double> &pntArrData);


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
