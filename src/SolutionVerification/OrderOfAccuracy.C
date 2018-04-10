#include <OrderOfAccuracy.H>

// TODO: RUN THE SIMULATIONS TO ACHIEVE ASYMPTOTIC GRID CONVERGENCE REGIME

OrderOfAccuracy::OrderOfAccuracy(meshBase* _f1, meshBase* _f2, meshBase* _f3,
                const std::vector<int>& _arrayIDs)
  : f1(_f1), f2(_f2), f3(_f3), arrayIDs(_arrayIDs)
{
  // set names for array from coarse mesh to be transferred to fine
  f3ArrNames.resize((arrayIDs.size()));
  f2ArrNames.resize((arrayIDs.size()));
  diffIDs.resize(arrayIDs.size());
  relEIDs.resize(arrayIDs.size());
  for (int i = 0; i < arrayIDs.size(); ++i)
  {
    std::string name3(f1->getDataSet()->GetPointData()->GetArrayName(arrayIDs[i]));
    std::string name2(name3);
    std::string coarse("f3");
    std::string fine("f2");
    name3 += coarse;
    name2 += fine;
    f3ArrNames[i] = name3;
    f2ArrNames[i] = name2;
  }
  f3->setNewArrayNames(f3ArrNames);
  f2->setNewArrayNames(f2ArrNames);
  f3->transfer(f2, "Finite Element", arrayIDs);
  f2->transfer(f1, "Finite Element", arrayIDs); 

  std::vector<double> f3CellLengths = f3->getCellLengths();
  std::vector<double> f2CellLengths = f2->getCellLengths();
  std::vector<double> f1CellLengths = f1->getCellLengths();
  double h3 = std::accumulate(f3CellLengths.begin(),f3CellLengths.end(),0.0)
              /f3->getNumberOfCells();
  double h2 = std::accumulate(f2CellLengths.begin(),f2CellLengths.end(),0.0)
              /f2->getNumberOfCells();
  double h1 = std::accumulate(f1CellLengths.begin(),f1CellLengths.end(),0.0)
              /f1->getNumberOfCells();
  r21 = h2/h1;
  r32 = h3/h2;
  diffF3F2 = computeDiffAndRelativeError(f2,f3ArrNames);
  diffF2F1 = computeDiffAndRelativeError(f1,f2ArrNames);

  f1->report();
  f2->report();
  f3->report();
}

std::vector<std::vector<double>> OrderOfAccuracy::computeOrderOfAccuracy()
{
  std::vector<std::vector<double>> orderOfAccuracy(diffF3F2.size());
  for (int i = 0; i < diffF3F2.size(); ++i)
  {
    orderOfAccuracy[i].resize(diffF3F2[i].size());
    for (int j = 0; j < diffF3F2[i].size(); ++j)
    {
      // picard-type iteration to solve transendental equation for p
      double q_p = 0;
      double f32_f21 = diffF3F2[i][j]/diffF2F1[i][j];
      int s = (f32_f21 > 0) - (f32_f21 < 0);
      double p = -1;
      double old_p = 1;
      for (int k = 0; k < 100; ++k)
      {
        p = std::fabs(log(f32_f21) + q_p)/log(r21);
        q_p = log((pow(r21,p) - s)/(pow(r32,p)-s));
        if (std::fabs(old_p - p) < 1e-16)
          break;
        else
          old_p = p;
      }
      orderOfAccuracy[i][j] = p;
    }
  }
  return orderOfAccuracy;
}


std::vector<std::vector<double>> OrderOfAccuracy::computeGridConvergenceIndex()
{
  std::vector<std::vector<double>> orderOfAccuracy(computeOrderOfAccuracy());
  std::vector<std::vector<double>> relativeError(f1->integrateOverMesh(relEIDs));
  std::vector<std::vector<double>> GCI(orderOfAccuracy.size());
  for (int i = 0; i < orderOfAccuracy.size(); ++i)
  {
    GCI[i].resize(orderOfAccuracy[i].size());
    for (int j = 0; j < orderOfAccuracy[i].size(); ++j)
    {
      GCI[i][j] = 1.25*std::sqrt(relativeError[i][j])/(pow(r21,orderOfAccuracy[i][j])-1);  
    }
  }
  return GCI;
}

std::vector<std::vector<double>> 
OrderOfAccuracy::computeDiffAndRelativeError
  (meshBase* mesh, const std::vector<std::string>& newArrNames)
{
  int numArr = arrayIDs.size();
  std::vector<vtkSmartPointer<vtkDoubleArray>> fineDatas(numArr);
  std::vector<vtkSmartPointer<vtkDoubleArray>> coarseDatas(numArr);
  std::vector<vtkSmartPointer<vtkDoubleArray>> diffDatas(numArr);
  std::vector<vtkSmartPointer<vtkDoubleArray>> relativeErrors(numArr);

  vtkSmartPointer<vtkPointData> finePD 
    = mesh->getDataSet()->GetPointData();

  std::vector<std::string> names(numArr);
  std::vector<std::string> namesRE(numArr);
  for (int id = 0; id < numArr; ++id)
  {
    fineDatas[id] 
      = vtkDoubleArray::SafeDownCast(finePD->GetArray(arrayIDs[id]));
    coarseDatas[id] 
      = vtkDoubleArray::SafeDownCast(finePD->GetArray(&(newArrNames[id])[0u]));
    vtkSmartPointer<vtkDoubleArray> diffData = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> relativeError = vtkSmartPointer<vtkDoubleArray>::New();
    diffData->SetNumberOfComponents(fineDatas[id]->GetNumberOfComponents());
    diffData->SetNumberOfTuples(mesh->getNumberOfPoints());
    std::string name(finePD->GetArrayName(arrayIDs[id]));
    name += "Diff";
    names[id] = name;
    diffData->SetName(&name[0u]);
    diffDatas[id] = diffData;
    relativeError->SetNumberOfComponents(fineDatas[id]->GetNumberOfComponents());
    relativeError->SetNumberOfTuples(mesh->getNumberOfPoints());
    std::string nameRE(finePD->GetArrayName(arrayIDs[id]));
    nameRE += "RelativeError";
    namesRE[id] = nameRE;
    relativeError->SetName(&nameRE[0u]);
    relativeErrors[id] = relativeError;
  }
  for (int i = 0; i < mesh->getNumberOfPoints(); ++i)
  {
    for (int id = 0; id < numArr; ++id)
    {
      int numComponent = fineDatas[id]->GetNumberOfComponents();
      double fine_comps[numComponent];
      double coarse_comps[numComponent];
      fineDatas[id]->GetTuple(i,fine_comps);
      coarseDatas[id]->GetTuple(i,coarse_comps);
      double diff[numComponent];
      double diffRel[numComponent];
      for (int j = 0; j < numComponent; ++j)
      { 
        double error = (coarse_comps[j] - fine_comps[j]);
        diff[j] = error*error;
        diffRel[j] = (error < 1e-8 ? 0 : pow(error/fine_comps[j],2)); 
      } 
      diffDatas[id]->SetTuple(i,diff);
      relativeErrors[id]->SetTuple(i,diffRel); 
    }  
  }
 
   
  for (int id = 0; id < numArr; ++id)
  {
    finePD->AddArray(diffDatas[id]);
    finePD->GetArray(&(names[id])[0u], diffIDs[id]);
    finePD->AddArray(relativeErrors[id]);
    finePD->GetArray(&(namesRE[id])[0u],relEIDs[id]);
  }
     
  std::vector<std::vector<double>> diff_integral(mesh->integrateOverMesh(diffIDs)); 
  for (int i = 0; i < diff_integral.size(); ++i)
  {
    for (int j = 0; j < diff_integral[i].size(); ++j)
    {
      diff_integral[i][j] = std::sqrt(diff_integral[i][j]);
      std::cout << diff_integral[i][j] << " ";
    }
    std::cout << std::endl;
  }
  return diff_integral; 
}  
  

//for (int id = 0; id < names.size(); ++id)
  //{
  //  std::string intgrl_name(names[id]);
  //  intgrl_name += "Integral";
  //  mesh->unsetPointDataArray(&(names[id])[0u]);
  //  mesh->unsetCellDataArray(&intgrl_name[0u]);
  //}
  //for (int id = 0; id < numArr; ++id)
  //{
  //  mesh->unsetPointDataArray(&(newArrNames[id])[0u]);
  //}
