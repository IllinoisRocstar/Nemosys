#include "AuxiliaryFunctions.H"
#include "SolutionVerification/OrderOfAccuracy.H"
#include "Drivers/TransferDriver.H"

#include <vtkDoubleArray.h>
#include <vtkPointData.h>

OrderOfAccuracy::OrderOfAccuracy(meshBase *_f3, meshBase *_f2, meshBase *_f1,
                                 std::vector<int> _arrayIDs,
                                 std::string transferType, double targetGCI)
    : f3(_f3),
      f2(_f2),
      f1(_f1),
      arrayIDs(std::move(_arrayIDs)),
      targetGCI(targetGCI) {
  // set names for array from coarse mesh to be transferred to fine
  int numArr = arrayIDs.size();

  f3ArrNames.resize(numArr);
  f2ArrNames.resize(numArr);
  diffIDs.resize(numArr);
  relEIDs.resize(numArr);
  realDiffIDs.resize(numArr);

  for (int i = 0; i < numArr; ++i) {
    // get names from coarsest mesh
    std::string name =
        f3->getDataSet()->GetPointData()->GetArrayName(arrayIDs[i]);
    std::string coarseName = name + "f3";
    std::string fineName = name + "f2";
    f3ArrNames[i] = coarseName;
    f2ArrNames[i] = fineName;
  }

  auto f3f2Transfer =
      NEM::DRV::TransferDriver::CreateTransferObject(f3, f2, transferType);
  f3f2Transfer->transferPointData(arrayIDs, f3ArrNames);
  auto f2f1Transfer =
      NEM::DRV::TransferDriver::CreateTransferObject(f2, f1, transferType);
  f2f1Transfer->transferPointData(arrayIDs, f2ArrNames);

  diffF3F2 = computeDiff(f2, f3ArrNames);
  diffF2F1 = computeDiff(f1, f2ArrNames);

  // TODO: Double-check the integer division below.
  r21 = pow(f1->getNumberOfPoints() / f2->getNumberOfPoints(), 1. / 3.);
  r32 = pow(f2->getNumberOfPoints() / f3->getNumberOfPoints(), 1. / 3.);
  std::cout << "Refinement ratio from 2-to-1 is " << r21 << " "
            << " and from 3-to-2 is " << r32 << std::endl;
}

std::vector<std::vector<double>> OrderOfAccuracy::computeOrderOfAccuracy() {
  orderOfAccuracy.resize(diffF3F2.size());
  for (int i = 0; i < diffF3F2.size(); ++i) {
    orderOfAccuracy[i].resize(diffF3F2[i].size());
    for (int j = 0; j < diffF3F2[i].size(); ++j) {
      // Picard-type iteration to solve transcendental equation for p
      double q_p = 0;
      double f32_f21 = diffF3F2[i][j] / diffF2F1[i][j];
      int s = (f32_f21 > 0) - (f32_f21 < 0);
      double p = -1;
      double old_p = 1;
      for (int k = 0; k < 1000; ++k) {
        p = std::fabs(log(std::fabs(f32_f21)) + q_p) / log(r21);
        q_p = log((pow(r21, p) - s) / (pow(r32, p) - s));
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

std::vector<std::vector<double>> OrderOfAccuracy::computeGCI_21() {
  if (GCI_21.empty()) {
    if (orderOfAccuracy.empty()) {
      computeOrderOfAccuracy();
    }
    std::vector<std::vector<double>> f1L2(f1->integrateOverMesh(relEIDs));
    GCI_21.resize(orderOfAccuracy.size());
    for (int i = 0; i < orderOfAccuracy.size(); ++i) {
      GCI_21[i].resize(orderOfAccuracy[i].size());
      for (int j = 0; j < orderOfAccuracy[i].size(); ++j) {
        double relativeError = diffF2F1[i][j] / std::sqrt(f1L2[i][j]);
        GCI_21[i][j] =
            1.25 * relativeError / (pow(r21, orderOfAccuracy[i][j]) - 1);
      }
    }
  }
  return GCI_21;
}

std::vector<std::vector<double>> OrderOfAccuracy::computeGCI_32() {
  if (GCI_32.empty()) {
    if (orderOfAccuracy.empty()) {
      computeOrderOfAccuracy();
    }
    std::vector<std::vector<double>> f2L2(f2->integrateOverMesh(relEIDs));
    GCI_32.resize(orderOfAccuracy.size());
    for (int i = 0; i < orderOfAccuracy.size(); ++i) {
      GCI_32[i].resize(orderOfAccuracy[i].size());
      for (int j = 0; j < orderOfAccuracy[i].size(); ++j) {
        double relativeError = diffF3F2[i][j] / std::sqrt(f2L2[i][j]);
        GCI_32[i][j] =
            1.25 * relativeError / (pow(r32, orderOfAccuracy[i][j]) - 1);
      }
    }
  }
  return GCI_32;
}

std::vector<std::vector<double>> OrderOfAccuracy::computeResolution(
    double gciStar) {
  if (GCI_32.empty()) computeGCI_32();
  std::vector<std::vector<double>> res(GCI_32.size());
  for (int i = 0; i < GCI_32.size(); ++i) {
    res[i].resize(GCI_32[i].size());
    for (int j = 0; j < GCI_32[i].size(); ++j) {
      res[i][j] = pow(gciStar / GCI_32[i][j], 1. / orderOfAccuracy[i][j]);
    }
  }
  return res;
}

void OrderOfAccuracy::computeMeshWithResolution(double gciStar,
                                                const std::string &ofname) {
  std::vector<std::vector<double>> ratios(checkAsymptoticRange());
  for (const auto &ratio : ratios) {
    for (double j : ratio) {
      if (std::abs(j - 1) > 1) {
        std::cerr
            << "WARNING: Solutions are not in the asymptotic convergence range"
            << std::endl;
      }
    }
  }

  std::vector<double> resolution = nemAux::flatten(computeResolution(gciStar));
  auto minmax = std::minmax_element(resolution.begin(), resolution.end());
  double ave = (*minmax.first + *minmax.second) / 2;
  f3->refineMesh("uniform", ave, ofname, false);
  meshBase *refined = meshBase::Create(ofname);
  // f3->unsetNewArrayNames();
  // f3->transfer(refined, "Consistent Interpolation", arrayIDs);
  auto f3refinedTransfer = NEM::DRV::TransferDriver::CreateTransferObject(
      f3, refined, "Consistent Interpolation");
  f3refinedTransfer->transferPointData(arrayIDs);
  delete f3;
  f3 = refined;
  computeRichardsonExtrapolation();
  f3->write(ofname);
}

std::vector<std::vector<double>> OrderOfAccuracy::checkAsymptoticRange() {
  if (GCI_21.empty()) computeGCI_21();
  if (GCI_32.empty()) computeGCI_32();
  std::vector<std::vector<double>> ratios(GCI_21.size());
  for (int i = 0; i < GCI_21.size(); ++i) {
    ratios[i].resize(GCI_21[i].size());
    for (int j = 0; j < GCI_21[i].size(); ++j) {
      ratios[i][j] =
          GCI_32[i][j] / (pow(r21, orderOfAccuracy[i][j]) * GCI_21[i][j]);
    }
  }
  return ratios;
}

void OrderOfAccuracy::computeRichardsonExtrapolation() {
  if (orderOfAccuracy.empty()) computeOrderOfAccuracy();

  int numArr = arrayIDs.size();

  std::vector<vtkSmartPointer<vtkDoubleArray>> fineDatas(numArr);
  std::vector<vtkSmartPointer<vtkDoubleArray>> diffDatas(numArr);
  std::vector<vtkSmartPointer<vtkDoubleArray>> richardsonDatas(numArr);

  vtkSmartPointer<vtkPointData> finePD = f1->getDataSet()->GetPointData();

  std::vector<std::string> names(numArr);
  for (int id = 0; id < numArr; ++id) {
    diffDatas[id] =
        vtkDoubleArray::SafeDownCast(finePD->GetArray(realDiffIDs[id]));
    fineDatas[id] =
        vtkDoubleArray::SafeDownCast(finePD->GetArray(arrayIDs[id]));

    auto richardsonData = vtkSmartPointer<vtkDoubleArray>::New();
    richardsonData->SetNumberOfComponents(
        diffDatas[id]->GetNumberOfComponents());
    richardsonData->SetNumberOfTuples(f1->getNumberOfPoints());
    std::string name(finePD->GetArrayName(arrayIDs[id]));
    name += "_richExtrap";
    names[id] = name;
    richardsonData->SetName(&name[0u]);
    richardsonDatas[id] = richardsonData;
  }

  for (int i = 0; i < f1->getNumberOfPoints(); ++i) {
    for (int id = 0; id < numArr; ++id) {
      int numComponent = diffDatas[id]->GetNumberOfComponents();

      auto *fine_comps = new double[numComponent];
      fineDatas[id]->GetTuple(i, fine_comps);

      auto *diff_comps = new double[numComponent];
      diffDatas[id]->GetTuple(i, diff_comps);

      auto *val = new double[numComponent];
      for (int j = 0; j < numComponent; ++j) {
        double comp = fine_comps[j] +
                      diff_comps[j] / (pow(r21, orderOfAccuracy[id][j]) - 1);
        val[j] = comp;
      }
      richardsonDatas[id]->SetTuple(i, val);

      delete[] fine_comps;
      delete[] diff_comps;
      delete[] val;
    }
  }

  std::vector<int> richExtrapIDs(numArr);
  for (int id = 0; id < numArr; ++id) {
    finePD->AddArray(richardsonDatas[id]);
    finePD->GetArray(names[id].c_str(), richExtrapIDs[id]);
  }
  // f1->transfer(f3, "Consistent Interpolation", richExtrapIDs);
  auto f1f3Transfer = NEM::DRV::TransferDriver::CreateTransferObject(
      f1, f3, "Consistent Interpolation");
  f1f3Transfer->transferPointData(arrayIDs, f1->getNewArrayNames());
}

std::vector<std::vector<double>> OrderOfAccuracy::computeDiff(
    meshBase *mesh, const std::vector<std::string> &newArrNames) {
  // diff is computed across multiple data items, specified by arrayIDs
  int numArr = arrayIDs.size();
  std::vector<vtkSmartPointer<vtkDoubleArray>> fineDatas(numArr);
  std::vector<vtkSmartPointer<vtkDoubleArray>> coarseDatas(numArr);
  std::vector<vtkSmartPointer<vtkDoubleArray>> diffDatas(numArr);
  std::vector<vtkSmartPointer<vtkDoubleArray>> fineDatasSqr(numArr);
  std::vector<vtkSmartPointer<vtkDoubleArray>> realDiffDatas(numArr);

  vtkSmartPointer<vtkPointData> finePD = mesh->getDataSet()->GetPointData();

  std::vector<std::string> names(numArr);
  std::vector<std::string> names2(numArr);
  std::vector<std::string> names3(numArr);

  std::string diffSqrName, sqrName, diffName;
  for (int id = 0; id < numArr; ++id) {
    fineDatas[id] =
        vtkDoubleArray::SafeDownCast(finePD->GetArray(arrayIDs[id]));
    coarseDatas[id] =
        vtkDoubleArray::SafeDownCast(finePD->GetArray(newArrNames[id].c_str()));

    // initialize
    vtkSmartPointer<vtkDoubleArray> diffData =
        vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> fineDataSqr =
        vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> realDiffData =
        vtkSmartPointer<vtkDoubleArray>::New();

    diffData->SetNumberOfComponents(fineDatas[id]->GetNumberOfComponents());
    diffData->SetNumberOfTuples(mesh->getNumberOfPoints());
    std::string name(finePD->GetArrayName(arrayIDs[id]));
    name += "DiffSqr";
    names[id] = name;
    diffData->SetName(name.c_str());
    diffDatas[id] = diffData;

    fineDataSqr->SetNumberOfComponents(fineDatas[id]->GetNumberOfComponents());
    fineDataSqr->SetNumberOfTuples(mesh->getNumberOfPoints());
    std::string name2(finePD->GetArrayName(arrayIDs[id]));
    name2 += "Sqr";
    names2[id] = name2;
    fineDataSqr->SetName(name2.c_str());
    fineDatasSqr[id] = fineDataSqr;

    realDiffData->SetNumberOfComponents(fineDatas[id]->GetNumberOfComponents());
    realDiffData->SetNumberOfTuples(mesh->getNumberOfPoints());
    std::string name3(finePD->GetArrayName(arrayIDs[id]));
    name3 += "Diff";
    names3[id] = name3;
    realDiffData->SetName(name3.c_str());
    realDiffDatas[id] = realDiffData;
  }

  int numPts = mesh->getNumberOfPoints();
  for (int id = 0; id < numArr; ++id) {
    auto fineData = fineDatas[id];
    auto coarseData = coarseDatas[id];
    auto diffData = diffDatas[id];
    auto fineDataSqr = fineDatasSqr[id];
    auto realDiffData = realDiffDatas[id];

    int numComps = fineDatas[id]->GetNumberOfComponents();

    double *fine_comps = new double[numComps];
    double *coarse_comps = new double[numComps];
    double *diff = new double[numComps];
    double *fsqr = new double[numComps];
    double *realdiff = new double[numComps];

    for (int i = 0; i < numPts; ++i) {
      fineData->GetTuple(i, fine_comps);
      coarseData->GetTuple(i, coarse_comps);

      for (int j = 0; j < numComps; ++j) {
        double error = coarse_comps[j] - fine_comps[j];
        diff[j] = error * error;
        fsqr[j] = fine_comps[j] * fine_comps[j];
        realdiff[j] = fine_comps[j] - coarse_comps[j];
      }

      diffData->SetTuple(i, diff);
      fineDataSqr->SetTuple(i, fsqr);
      realDiffData->SetTuple(i, realdiff);
    }
    delete[] fine_comps;
    delete[] coarse_comps;
    delete[] diff;
    delete[] fsqr;
    delete[] realdiff;
  }

  for (int id = 0; id < numArr; ++id) {
    finePD->AddArray(diffDatas[id]);
    finePD->GetArray(names[id].c_str(), diffIDs[id]);
    finePD->AddArray(fineDatasSqr[id]);
    finePD->GetArray(names2[id].c_str(), relEIDs[id]);
    finePD->AddArray(realDiffDatas[id]);
    finePD->GetArray(names3[id].c_str(), realDiffIDs[id]);
  }

  std::vector<std::vector<double>> diff_integral(
      mesh->integrateOverMesh(diffIDs));
  for (auto &&i : diff_integral) {
    for (double &j : i) {
      j = std::sqrt(j);
      std::cout << j << " ";
    }
    std::cout << std::endl;
  }
  return diff_integral;
}
