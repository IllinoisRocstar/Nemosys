#ifndef vtkMesh_H
#define vtkMesh_H

#include <meshBase.H>

using namespace nemAux;

template<class TReader> vtkDataSet* ReadAnXMLFile(const char* fileName)
{
  vtkSmartPointer<TReader> reader =
    vtkSmartPointer<TReader>::New();
  reader->SetFileName(fileName);
  reader->Update();
  reader->GetOutput()->Register(reader);
  return vtkDataSet::SafeDownCast(reader->GetOutput());
}

template<class TWriter> void writeVTFile(std::string fname, vtkDataSet* dataSet)
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



  // integration
  public:
      
    std::vector<double> getGaussWeightTet4(int I)
    {
      double baryWeights4[2] = {.585410196624968, .138196601125010};
      std::vector<double> sfWeights(4);
      if (I == 1)
      {
        sfWeights[0] = baryWeights4[0];
        sfWeights[1] = baryWeights4[1];
        sfWeights[2] = baryWeights4[1];
        sfWeights[3] = baryWeights4[1];
      }
      else if (I == 2)
      {
        sfWeights[0] = baryWeights4[1];
        sfWeights[1] = baryWeights4[0];
        sfWeights[2] = baryWeights4[1];
        sfWeights[3] = baryWeights4[1];
      }
      else if (I == 3)
      {
        sfWeights[0] = baryWeights4[1];
        sfWeights[1] = baryWeights4[1];
        sfWeights[2] = baryWeights4[0];
        sfWeights[3] = baryWeights4[1];
      }
      else if (I == 4)
      {
        sfWeights[0] = baryWeights4[1];
        sfWeights[1] = baryWeights4[1];
        sfWeights[2] = baryWeights4[1];
        sfWeights[3] = baryWeights4[0];
      }
      return sfWeights;
    }   
 
    std::vector<double> getIntegrationPointAtCell(int I, 
                                                  const std::vector<std::vector<double>>& cellPoints)
        
    {
      std::vector<double> sfWeights = getGaussWeightTet4(I);

      std::vector<double> gaussPoint(3,0);
      for (int j = 0; j < cellPoints.size(); ++j)
      { 
        for (int i = 0; i < 3; ++i)
        {
            gaussPoint[i] += sfWeights[j]*cellPoints[j][i];
        }
      }
      return gaussPoint;
    }    
 
    void getIntegrationPointsAtCell(int cellID)
    {

      if (dataSet->GetCell(cellID)->GetCellType() == VTK_TETRA)
      {
        int numGaussPoints = 4;
        std::vector<std::vector<double>> gaussPoints(numGaussPoints);
        std::vector<std::vector<double>> cellCoords = getCellVec(cellID);
        for (int i = 0; i < numGaussPoints; ++i)
        {
          gaussPoints[i] = getIntegrationPointAtCell(i+1,cellCoords);
          printVec(gaussPoints[i]);
        }

      }
      

      // Get the dictionary key     
      vtkInformationQuadratureSchemeDefinitionVectorKey *key = 
          vtkQuadratureSchemeDefinition::DICTIONARY(); 

      if (dataSet->GetCellData() == NULL)
        std::cout << "NO CELL DATA" << std::endl;

      // Get the cell types used by the data set
      vtkSmartPointer<vtkCellTypes> cellTypes 
        = vtkSmartPointer<vtkCellTypes>::New();
      dataSet->GetCellTypes(cellTypes);
      int nCellTypes = cellTypes->GetNumberOfTypes(); 

      // create offset array and store the dictionary within
      vtkSmartPointer<vtkIdTypeArray> offsets 
        = vtkSmartPointer<vtkIdTypeArray>::New();
      
      std::string basename = "QuadratureOffset";
      vtkSmartPointer<vtkDataArray> data = dataSet->GetCellData()->GetArray(basename.c_str());

      offsets->SetName(basename.c_str());
      dataSet->GetCellData()->AddArray(offsets);
      
      vtkSmartPointer<vtkInformation> info = vtkSmartPointer<vtkInformation>::New();
      info = offsets->GetInformation();

      for (int typeId = 0; typeId < nCellTypes; ++typeId)
      {
        int cellType = cellTypes->GetCellType(typeId);
        // Initialize quadrature scheme definition for given cell type
        vtkSmartPointer<vtkQuadratureSchemeDefinition> def
          = vtkSmartPointer<vtkQuadratureSchemeDefinition>::New();
        switch(cellType)
        {
          case VTK_TRIANGLE:
            def->Initialize(VTK_TRIANGLE,3,3,TRI3);
            break;
          case VTK_TETRA:
            def->Initialize(VTK_TETRA, 4, 4, TET4);
            break;
          default:
            std::cout << "Error: Cell type: " << cellType << "found "
                      << "with no quadrature definition provided" << std::endl;
            exit(1);
        }
        // the definition must apear in the dictionary associated with 
        // the offset array
        key->Set(info, def, cellType); 
      }

      // get dictionary size 
      int dictSize = key->Size(info);
      vtkQuadratureSchemeDefinition** dict 
       = new vtkQuadratureSchemeDefinition * [dictSize];
      key->GetRange(info, dict, 0, 0, dictSize);
      
      offsets->SetNumberOfTuples(dataSet->GetNumberOfCells());
      vtkIdType offset = 0;
      
      for (int cellid = 0; cellid < dataSet->GetNumberOfCells() ; ++cellid)
      {
        offsets->SetValue(cellid, offset);
        vtkCell* cell = dataSet->GetCell(cellid);
        int cellType = cell->GetCellType();
        vtkQuadratureSchemeDefinition* celldef = dict[cellType];
        offset += celldef->GetNumberOfQuadraturePoints();
      }  
    


      vtkSmartPointer<vtkQuadraturePointsGenerator> pointGen = 
        vtkSmartPointer<vtkQuadraturePointsGenerator>::New();

      pointGen->SetInputArrayToProcess
                (0, 0, 0, 
                 vtkDataObject::FIELD_ASSOCIATION_CELLS, 
                 "QuadratureOffset");
      pointGen->SetInputData(dataSet);

      vtkSmartPointer<vtkPolyData> output = vtkSmartPointer<vtkPolyData>::New();
      output = vtkPolyData::SafeDownCast(pointGen->GetOutput());
      pointGen->Update();

      // get desired point data array to be interpolated to gauss points
      vtkSmartPointer<vtkDataArray> da = dataSet->GetPointData()->GetArray(7);
      // get tuple length of given data
      int numComponent = da->GetNumberOfComponents();
      // declare data array to be populated for polydata
      vtkSmartPointer<vtkDoubleArray> daGauss = vtkSmartPointer<vtkDoubleArray>::New();
      // names and sizing
      daGauss->SetName(dataSet->GetPointData()->GetArrayName(7));
      daGauss->SetNumberOfComponents(numComponent);
      // generic cell to store given cell in dataSet
      vtkSmartPointer<vtkGenericCell> genCell = vtkSmartPointer<vtkGenericCell>::New();       
      // number of points in poly data to which data has been interpolated
      int polyPnt = 0;
      for (int i = 0; i < numCells; ++i)
      {
        // putting current cell into genCell
        dataSet->GetCell(i, genCell);
        // getting cellType information for lookup in dictionary
        int cellType = dataSet->GetCellType(i);
        // number of gauss points in this cell
        int numGaussPoints = dict[cellType]->GetNumberOfQuadraturePoints();
          
        for (int j = 0; j < numGaussPoints; ++j)
        {
          double x[3];
          output->GetPoint(j+polyPnt,x); 
          // parameters for interpolation
          int subId; // not used
          double minDist2; // not used
          double pcoords[3];
          double comps[numComponent];
          double weights[genCell->GetNumberOfPoints()];
          genCell->EvaluatePosition(x,NULL,subId,pcoords,minDist2,weights); 
          double interps[numComponent];
          for (int m = 0; m < genCell->GetNumberOfPoints(); ++m)
          {
            int pntId = genCell->GetPointId(m);
            da->GetTuple(pntId, comps);
            for (int h = 0; h < numComponent; ++h)
            {
              interps[h] += comps[h]*weights[m]; 
            }
          }
          daGauss->InsertNextTuple(interps);
        } 
        polyPnt += numGaussPoints;
      }
      std::cout << "NUM POLY POINTS: " << polyPnt << std::endl; 
      output->GetPointData()->SetActiveScalars(dataSet->GetPointData()->GetArrayName(7));
      output->GetPointData()->SetScalars(daGauss);        

      writeVTFile<vtkXMLPolyDataWriter> ("refined_beamGaussPoints.vtp",output);
      
      delete [] dict;
      



      //vtkDataSet* dataSet_tmp = fieldInterp->GetOutput();
      //fieldInterp->Update();

      //dataSet->CopyAttributes(dataSet_tmp);    
      //       
     //// vtkSmartPointer<vtkUnstructuredGrid> dataSet_tmp
     ////   = vtkSmartPointer<vtkUnstructuredGrid>::New();
     //// dataSet_tmp = vtkUnstructuredGrid::SafeDownCast(fieldInterp->GetOutput());
     //// fieldInterp->Update();
      //
      //writeVTFile<vtkXMLUnstructuredGridWriter> ("refined_beamGaussFieldInterp.vtu",dataSet);

      //vtkSmartPointer<vtkQuadraturePointsGenerator> pointGen = 
      //  vtkSmartPointer<vtkQuadraturePointsGenerator>::New();
      //pointGen->SetInputArrayToProcess
      //          (0, 0, 0, 
      //           vtkDataObject::FIELD_ASSOCIATION_CELLS, 
      //           "QuadratureOffset");
      //pointGen->SetInputData(dataSet);
      ////vtkSmartPointer<vtkPolyData> output = vtkSmartPointer<vtkPolyData>::New();
      ////pointGen->Update();
      //vtkPolyData* output=vtkPolyData::SafeDownCast(pointGen->GetOutput());

      //writeVTFile<vtkXMLPolyDataWriter> ("refined_beamGaussPoints.vtp",output);
      
    
    }
 
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
