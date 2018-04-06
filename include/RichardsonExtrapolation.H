#ifndef RICHARDSON_EXTRAPOLATION
#define RICHARDSON_EXTRAPOLATION

#include <meshBase.H>

class RichardsonExtrapolation
{

  public:
    RichardsonExtrapolation(meshBase* _fineMesh, meshBase* _coarseMesh,
                            double _ref_factor, int _order, 
                            const std::vector<int>& _arrayIDs)
      : fineMesh(_fineMesh), coarseMesh(_coarseMesh), ref_factor(_ref_factor), order(_order),
        arrayIDs(_arrayIDs)
    {
      // set names for array from coarse mesh to be transferred to fine
      std::vector<std::string> newArrNames(arrayIDs.size());
      for (int i = 0; i < arrayIDs.size(); ++i)
      {
        std::string name(fineMesh->getDataSet()->GetPointData()->GetArrayName(arrayIDs[i]));
        std::string ext("New");
        name += ext;
        newArrNames[i] = name;
      } 
      coarseMesh->setNewArrayNames(newArrNames);
      // transfer data from coarse mesh to fine
      coarseMesh->transfer(fineMesh, "Finite Element", arrayIDs);
    }

    ~RichardsonExtrapolation(){}
    
    std::vector<std::vector<double>> computeDiscretizationError()
    {
      int numArr = arrayIDs.size();
      std::vector<vtkSmartPointer<vtkDoubleArray>> fineDatas(numArr);
      std::vector<vtkSmartPointer<vtkDoubleArray>> coarseDatas(numArr);
      std::vector<vtkSmartPointer<vtkDoubleArray>> diffDatas(numArr);
  
      vtkSmartPointer<vtkPointData> finePD 
        = fineMesh->getDataSet()->GetPointData();
    
      for (int id = 0; id < numArr; ++id)
      {
        fineDatas[id] 
          = vtkDoubleArray::SafeDownCast(finePD->GetArray(arrayIDs[id]));
        coarseDatas[id] 
          = vtkDoubleArray::SafeDownCast(finePD->GetArray(arrayIDs[id] + numArr));
        vtkSmartPointer<vtkDoubleArray> diffData = vtkSmartPointer<vtkDoubleArray>::New();
        diffData->SetNumberOfComponents(fineDatas[id]->GetNumberOfComponents());
        diffData->SetNumberOfTuples(fineMesh->getNumberOfPoints());
        std::string name(finePD->GetArrayName(arrayIDs[id]));
        name += "Diff";
        diffData->SetName(&name[0u]);
        diffDatas[id] = diffData;
      }
      
      for (int i = 0; i < fineMesh->getNumberOfPoints(); ++i)
      {
        for (int id = 0; id < numArr; ++id)
        {
          int numComponent = fineDatas[id]->GetNumberOfComponents();
          double fine_comps[numComponent];
          double coarse_comps[numComponent];
          fineDatas[id]->GetTuple(i,fine_comps);
          coarseDatas[id]->GetTuple(i,coarse_comps);
          double diff[numComponent];
          for (int j = 0; j < numComponent; ++j)
          {
            diff[j] = pow((fine_comps[j] - coarse_comps[j])/(pow(ref_factor,order)-1),2);  
          } 
          diffDatas[id]->SetTuple(i,diff);
        }  
      }
      
      std::vector<int> diffIDs(numArr);
      for (int id = 0; id < numArr; ++id)
      {
        finePD->AddArray(diffDatas[id]);
        diffIDs[id] = arrayIDs[id] + 2*numArr; 
      }
         
      std::vector<std::vector<double>> discretization_error(fineMesh->integrateOverMesh(diffIDs)); 
      for (int i = 0; i < discretization_error.size(); ++i)
      {
        for (int j = 0; j < discretization_error[i].size(); ++j)
        {
          discretization_error[i][j] = log(discretization_error[i][j]);
          std::cout << discretization_error[i][j] << std::endl;
        }
      } 
      return discretization_error; 
    }
      

  private:
    meshBase* fineMesh;
    meshBase* coarseMesh;
    double ref_factor;
    int order;
    const std::vector<int> arrayIDs;

};



#endif
