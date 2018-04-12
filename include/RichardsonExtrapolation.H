#ifndef RICHARDSON_EXTRAPOLATION
#define RICHARDSON_EXTRAPOLATION

#include <meshBase.H>

class RichardsonExtrapolation
{

  public:
    RichardsonExtrapolation(meshBase* _fineMesh, meshBase* coarseMesh,
                            double _ref_factor, int _order, 
                            const std::vector<int>& _arrayIDs)
      : fineMesh(_fineMesh), ref_factor(_ref_factor), order(_order),
        arrayIDs(_arrayIDs)
    {
      // set names for array from coarse mesh to be transferred to fine
      newArrNames.resize((arrayIDs.size()));
      for (int i = 0; i < arrayIDs.size(); ++i)
      {
        std::string name(fineMesh->getDataSet()->GetPointData()->GetArrayName(arrayIDs[i]));
        std::string ext("Coarse");
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
   
      std::vector<std::string> names(numArr);
      for (int id = 0; id < numArr; ++id)
      {
        fineDatas[id] 
          = vtkDoubleArray::SafeDownCast(finePD->GetArray(arrayIDs[id]));
        coarseDatas[id] 
          = vtkDoubleArray::SafeDownCast(finePD->GetArray(&(newArrNames[id])[0u]));
        vtkSmartPointer<vtkDoubleArray> diffData = vtkSmartPointer<vtkDoubleArray>::New();
        diffData->SetNumberOfComponents(fineDatas[id]->GetNumberOfComponents());
        diffData->SetNumberOfTuples(fineMesh->getNumberOfPoints());
        std::string name(finePD->GetArrayName(arrayIDs[id]));
        name += "Diff";
        names[id] = name;
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
            double error = (coarse_comps[j] -fine_comps[j]);
            diff[j] = error*error; 
          } 
          diffDatas[id]->SetTuple(i,diff);
        }  
      }
      
      std::vector<int> diffIDs(numArr); 
      for (int id = 0; id < numArr; ++id)
      {
        finePD->AddArray(diffDatas[id]);
        finePD->GetArray(&(names[id])[0u], diffIDs[id]);
      }
         
      std::vector<std::vector<double>> discretization_error(fineMesh->integrateOverMesh(diffIDs)); 
      std::vector<std::vector<double>> discretization_errorAll(discretization_error.size()*2);
      int k = 0;
      for (int i = 0; i < discretization_error.size(); ++i)
      {
        discretization_errorAll[k].resize(discretization_error[i].size());
        discretization_errorAll[k+1].resize(discretization_error[i].size());
        for (int j = 0; j < discretization_error[i].size(); ++j)
        {
          discretization_errorAll[k][j] = 0.5*log(discretization_error[i][j])
                                          - log(pow(ref_factor,order) - 1); 
          discretization_errorAll[k+1][j] = discretization_errorAll[k][j] + order*log(ref_factor);
        }
        k += 2;
      }
 
      for (int i = 0; i < discretization_errorAll.size(); ++i)
      {
        for (int j = 0; j < discretization_errorAll[i].size(); ++j)
        {
          std::cout << discretization_errorAll[i][j] << " ";
        }
        std::cout << std::endl;
      } 

      for (int id = 0; id < names.size(); ++id)
      {
        std::string intgrl_name(names[id]);
        intgrl_name += "Integral";
        fineMesh->unsetPointDataArray(&(names[id])[0u]);
        fineMesh->unsetCellDataArray(&intgrl_name[0u]);
      }
      for (int id = 0; id < numArr; ++id)
      {
        fineMesh->unsetPointDataArray(&(newArrNames[id])[0u]);
      }
      return discretization_errorAll; 
    }

    std::vector<double> computeObservedOrderOfAccuracy(meshBase* finerMesh)
    {
      std::vector<std::vector<double>> discretization_error 
        = computeDiscretizationError();
      std::unique_ptr<RichardsonExtrapolation> finerExtrapolation
        = std::unique_ptr<RichardsonExtrapolation>
          (new RichardsonExtrapolation(finerMesh, fineMesh, ref_factor,order,arrayIDs));
      std::vector<std::vector<double>> discretization_error1 
        = finerExtrapolation->computeDiscretizationError();
      fineMesh->report();
      finerMesh->report();      

      std::vector<double> observed_order;
      for (int i = 0; i < discretization_error.size(); i+=2)
      {
        for (int j = 0; j < discretization_error[i].size(); ++j)
        {
          observed_order.push_back((discretization_error[i][j] 
                                    - discretization_error1[i][j])/log(ref_factor)); 
        }
      }
      printVec(observed_order); 
      return observed_order;
    }
      

  private:
    meshBase* fineMesh;
    double ref_factor;
    int order;
    const std::vector<int> arrayIDs;
    std::vector<std::string> newArrNames;
};



#endif
