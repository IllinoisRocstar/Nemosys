#include <vtkAnalyzer.H>

class meshPhys: public vtkAnalyzer
{
  public:
    // constructors
    meshPhys(char* fname):vtkAnalyzer(fname)
    {
      read();
      
      if (dataSet)
      {
        numArr = dataSet->GetPointData()->GetNumberOfArrays();
        for (int i = 0; i < numArr; ++i) 
        {
          int numTuple, numComponent;
          std::vector<std::vector<double>> pntData_tmp;
          getPointDataArray(i, pntData_tmp, numTuple, numComponent);
          std::string array_name = dataSet->GetPointData()->GetArrayName(i);
          pntData[array_name] = pntData_tmp;
        }
      } 
    };

    ~meshPhys(){};

  public:
    // relevant computations (gradients etc.)


  public:

    writeBackgroundMSH(string filename, std::vector<double> sizes)
    {
      
      std::ofstream outputStream(filename.c_str());
      if(!outputStream.good()) 
      {
        std::cout << "Output file stream is bad" << std::endl;
        exit(1);
      }
     
      if (!dataSet) 
      {
        std::cout << "No data to write" << std::endl;
        exit(2);
      }
    
      // ---- get number of points and number of elements ---- //
      getNumberOfCells();
      getNumberOfPoints();
      
      if (sizes.size() != numberOfCells)
      {
        std::cout << "number of sizes is unequal to number of cells!" << std::endl; 
        exit(3);
      }
      

      // ---------  writing gmsh header ----------- //
      outputStream << "$MeshFormat" << std::endl
                   << "2.2 0 8" << std::endl
                   << "$EndMeshFormat" << std::endl; 

      // -------- ensure all cell types are tri/tet or below -------------- //
      int num_bad = 0;
      for (int i = 0; i < numberOfCells; i++)
      {
        int type_id = dataSet->GetCellType(i);
        if (!(type_id == 3 || type_id == 5 || type_id == 10))
        {
          std::cout << "Error: Only tetrahedral and triangular" 
                    << " meshes can be written to gmsh format" << std::endl;
          exit(3);
        }
        if (!(type_id == 5 || type_id == 10))
          num_bad+=1;
      }

      // ------------------------ write point coords -------------------------- //
      outputStream << "$Nodes" << std::endl << numberOfPoints << std::endl;
      for (int i = 0; i < numberOfPoints; ++i)
      {
        double* pntcrds = getPointCoords(i);
        outputStream << i + 1 << " "
                     << pntcrds[0] << " "
                     << pntcrds[1] << " "
                     << pntcrds[2] << " " << std::endl;
      }
      outputStream << "$EndNodes" << std::endl;

      // ------------- write element type and connectivity --------------------- //
      outputStream << "$Elements" << std::endl << numberOfCells-num_bad << std::endl;
      int k = 0;
      for (int i = 0; i < numberOfCells; ++i)
      {

        vtkIdList* point_ids = dataSet->GetCell(i)->GetPointIds();
        int numComponent = point_ids->GetNumberOfIds();
        int type_id = dataSet->GetCellType(i);
        if (type_id==5 || type_id == 10)
        {
          outputStream << k + 1 << " ";
          switch(numComponent)
          {
            case 2:
            {
              break;
            }
            case 3:
            {
              outputStream << 2 << " " << 2 << " " << 1 << " " << 1 << " ";
              break;
            }
            case 4:
            {
              outputStream << 4 << " " << 2 << " " << 1 << " " << 1 << " ";
              break;
            }
          
            default: 
            {  
              std::cerr << "Components in cell should be less than 4"<< std::endl;
              exit(1);
            }
          }
          for (int j = 0; j < numComponent; ++j)
             outputStream << point_ids->GetId(j) + 1 << " ";
          outputStream << std::endl;
          k+=1;
        }
      }
      outputStream << "$EndElements" << std::endl;
      
      // -------------------------- write cell data ---------------------------- // 

      std::string tmpname = "BackgroundSF";
      outputStream << "$ElementData" << std::endl
                   << 1 << std::endl // 1 string tag
                   << "\"" << (tmpname) // name of view
                   << "\"" << std::endl 
                   << 0 << std::endl // 0 real tag
                   << 3 << std::endl // 3 int tags (dt index, dim of field, number of fields)
                   << 0 << std::endl // dt index
                   << 1 << std::endl // dim of field
                   << numberOfCells-num_bad << std::endl; // number of fields
      int i = 0;
      for (int j = 0; j < numberOfCells; ++j)
      {
        int type_id = dataSet->GetCellType(j);
        if (type_id == 5 || type_id == 10) {
          outputStream << i+1 << " ";
          outputStream << sizes[j] << " ";
        outputStream << std::endl;
          i+=1;
        }
      }
      outputStream << "$EndElementData" << std::endl;    
    }


    private:
      // point data
      int numArr;
      std::map<std::string, std::vector<std::vector<double>>> pntData ;
      
};

 
  /*  meshPhys(char* fname):vtkAnalyzer(fname)
    {
      vtkAnalyzer::read();
      if (vtkAnalyzer::dataSet)
      {
        numArr = vtkAnalyzer::dataSet->GetPointData()->GetNumberOfArrays();
        pntData.resize(numArr);
        for (int i = 0; i < numArr; ++i) 
        {
          int numTuple, numComponent;
          std::vector<std::vector<double>> pntData_tmp;
          vtkAnalyzer::getPointDataArray(i, pntData_tmp, numTuple, numComponent);
          pntData[i].resize(numTuple);
          for (int j = 0; j < numTuple; ++j)
          {
            pntData[i][j].resize(numComponent);
            for (int k = 0; k < numComponent; ++k)
              pntData[i][j][k] = pntData_tmp[j][k];
          }
        }
        std::cout << "number of data arrays: " << numArr << std::endl;
        std::cout << "size of array0: " << pntData[0].size() << std::endl;
        std::cout << "random elem from array0: " << pntData[0][1000][1] << std::endl;
      } 
    };

        std::map<std::string, std::vector<std::vector<double>>>::iterator it = pntData.begin();
        while( it != pntData.end())
        {
          std::cout << it->first << ": " << (it->second)[1000][2] << std::endl;
          it++;
        }
    */
