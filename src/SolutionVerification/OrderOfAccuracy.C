#include <OrderOfAccuracy.H>

OrderOfAccuracy::OrderOfAccuracy(meshBase* _f1, meshBase* _f2, meshBase* _f3,
                const std::vector<int>& _arrayIDs)
  : f1(_f1), f2(_f2), f3(_f3), arrayIDs(_arrayIDs)
{
  // set names for array from coarse mesh to be transferred to fine
  f3ArrNames.resize((arrayIDs.size()));
  f2ArrNames.resize((arrayIDs.size()));
  diffIDs.resize(arrayIDs.size());
  relEIDs.resize(arrayIDs.size());
  realDiffIDs.resize(arrayIDs.size());
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

  r21 = pow(f1->getNumberOfPoints()/f2->getNumberOfPoints(),1./3.);
  r32 = pow(f2->getNumberOfPoints()/f3->getNumberOfPoints(),1./3.);

  std::cout << r21 << " " << r32 << std::endl;

  diffF3F2 = computeDiff(f2,f3ArrNames);
  diffF2F1 = computeDiff(f1,f2ArrNames);
  f1->report();
}

OrderOfAccuracy::~OrderOfAccuracy()
{

}

std::vector<std::vector<double>> OrderOfAccuracy::computeOrderOfAccuracy()
{
  orderOfAccuracy.resize(diffF3F2.size());
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


std::vector<std::vector<double>> OrderOfAccuracy::computeGCI_21()
{
  if (GCI_21.empty())
  {
    if (orderOfAccuracy.empty())
    {
      orderOfAccuracy = computeOrderOfAccuracy();
    }
    std::vector<std::vector<double>> f1L2(f1->integrateOverMesh(relEIDs));
    GCI_21.resize(orderOfAccuracy.size());
    for (int i = 0; i < orderOfAccuracy.size(); ++i)
    {
      GCI_21[i].resize(orderOfAccuracy[i].size());
      for (int j = 0; j < orderOfAccuracy[i].size(); ++j)
      {
        double relativeError = diffF2F1[i][j]/std::sqrt(f1L2[i][j]);
        GCI_21[i][j] = 1.25*relativeError/(pow(r21,orderOfAccuracy[i][j])-1); 
      }
    }
  }
  return GCI_21;
}

std::vector<std::vector<double>> OrderOfAccuracy::computeGCI_32()
{
  if (GCI_32.empty())
  {
    if (orderOfAccuracy.empty())
    {
      orderOfAccuracy = computeOrderOfAccuracy();
    }
    std::vector<std::vector<double>> f2L2(f2->integrateOverMesh(relEIDs));
    GCI_32.resize(orderOfAccuracy.size());
    for (int i = 0; i < orderOfAccuracy.size(); ++i)
    {
      GCI_32[i].resize(orderOfAccuracy[i].size());
      for (int j = 0; j < orderOfAccuracy[i].size(); ++j)
      {
        double relativeError = diffF3F2[i][j]/std::sqrt(f2L2[i][j]);
        GCI_32[i][j] = 1.25*relativeError/(pow(r32,orderOfAccuracy[i][j])-1);
      }
    }
  }
  return GCI_32;
}

std::vector<std::vector<double>> OrderOfAccuracy::computeResolution(double gciStar)
{
  std::vector<std::vector<double>> GCI(computeGCI_32());
  std::vector<std::vector<double>> res(GCI.size());
  for (int i = 0; i < GCI.size(); ++i)
  {
    res[i].resize(GCI[i].size());
    for (int j = 0; j < GCI[i].size(); ++j)
    {
      res[i][j] = pow(gciStar/GCI[i][j],1./orderOfAccuracy[i][j]);
    }
  }
  return res;
}

std::vector<std::vector<double>> OrderOfAccuracy::checkAsymptoticRange()
{
  if (GCI_21.empty())
    computeGCI_21();
  if (GCI_32.empty())
    computeGCI_32();
  std::vector<std::vector<double>> ratios(GCI_21.size());
  for (int i = 0; i < GCI_21.size(); ++i)
  {
    ratios[i].resize(GCI_21[i].size());
    for (int j = 0; j < GCI_21[i].size(); ++j)
    {
      ratios[i][j] = GCI_32[i][j]/(pow(r21,orderOfAccuracy[i][j])*GCI_21[i][j]);
    }
  }
  return ratios;
}

//OrderOfAccuracy::computeRichardsonExtrapolation()
//{
//  
//}

std::vector<std::vector<double>> 
OrderOfAccuracy::computeDiff
  (meshBase* mesh, const std::vector<std::string>& newArrNames)
{
  int numArr = arrayIDs.size();
  std::vector<vtkSmartPointer<vtkDoubleArray>> fineDatas(numArr);
  std::vector<vtkSmartPointer<vtkDoubleArray>> coarseDatas(numArr);
  std::vector<vtkSmartPointer<vtkDoubleArray>> diffDatas(numArr);
  std::vector<vtkSmartPointer<vtkDoubleArray>> fineDatasSqr(numArr);
  std::vector<vtkSmartPointer<vtkDoubleArray>> realDiffDatas(numArr);
  
  vtkSmartPointer<vtkPointData> finePD 
    = mesh->getDataSet()->GetPointData();

  std::vector<std::string> names(numArr);
  std::vector<std::string> names2(numArr);
  std::vector<std::string> names3(numArr);
  for (int id = 0; id < numArr; ++id)
  {
    fineDatas[id] 
      = vtkDoubleArray::SafeDownCast(finePD->GetArray(arrayIDs[id]));
    coarseDatas[id] 
      = vtkDoubleArray::SafeDownCast(finePD->GetArray(&(newArrNames[id])[0u]));
    vtkSmartPointer<vtkDoubleArray> diffData = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> fineDataSqr = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> realDiffData = vtkSmartPointer<vtkDoubleArray>::New();
    diffData->SetNumberOfComponents(fineDatas[id]->GetNumberOfComponents());
    diffData->SetNumberOfTuples(mesh->getNumberOfPoints());
    std::string name(finePD->GetArrayName(arrayIDs[id]));
    name += "DiffSqr";
    names[id] = name;
    diffData->SetName(&name[0u]);
    diffDatas[id] = diffData;
    fineDataSqr->SetNumberOfComponents(fineDatas[id]->GetNumberOfComponents());
    fineDataSqr->SetNumberOfTuples(mesh->getNumberOfPoints());
    std::string name2(finePD->GetArrayName(arrayIDs[id]));
    name2 += "Sqr";
    names2[id] = name2;
    fineDataSqr->SetName(&name2[0u]);
    fineDatasSqr[id] = fineDataSqr;
    realDiffData->SetNumberOfComponents(fineDatas[id]->GetNumberOfComponents());
    realDiffData->SetNumberOfTuples(mesh->getNumberOfPoints());
    std::string name3(finePD->GetArrayName(arrayIDs[id]));
    name3 += "Diff";
    names3[id] = name3;
    realDiffData->SetName(&name3[0u]);
    realDiffDatas[id] = realDiffData;
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
      double fsqr[numComponent];
      double realdiff[numComponent];
      for (int j = 0; j < numComponent; ++j)
      { 
        double error = (coarse_comps[j] - fine_comps[j]);
        diff[j] = error*error;
        fsqr[j] = fine_comps[j]*fine_comps[j];
        realdiff[j] = fine_comps[j] = coarse_comps[j];
      } 
      diffDatas[id]->SetTuple(i,diff);
      fineDatasSqr[id]->SetTuple(i,fsqr);
      realDiffDatas[id]->SetTuple(i,realdiff); 
    }  
  }
 
   
  for (int id = 0; id < numArr; ++id)
  {
    finePD->AddArray(diffDatas[id]);
    finePD->GetArray(&(names[id])[0u], diffIDs[id]);
    finePD->AddArray(fineDatasSqr[id]);
    finePD->GetArray(&(names2[id])[0u],relEIDs[id]);
    finePD->AddArray(realDiffDatas[id]);
    finePD->GetArray(&(names3[id])[0u],realDiffIDs[id]);
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
