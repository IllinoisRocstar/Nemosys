#include "SizeFieldGen.H"

#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkPointData.h>

#include "AuxiliaryFunctions.H"
#include "GradSizeField.H"
#include "ValSizeField.H"
#include "Z2ErrorSizeField.H"

namespace NEM {
namespace ADP {

using nemAux::operator*;  // for vector multiplication.

SizeFieldBase::SizeFieldBase(vtkDataSet *_ds, int arrayID, double _dev_mult,
                             bool _maxIsmin, const std::string &arrName)
    : ds(_ds), dev_mult(_dev_mult), maxIsmin(_maxIsmin), sizeFactor(1.) {
  // checking for point data
  int numArr = ds->GetPointData()->GetNumberOfArrays();
  if (arrayID >= numArr) {
    std::cerr << "ERROR: arrayID is out of bounds\n";
    std::cerr << "There are " << numArr << " point data arrays" << std::endl;
    exit(1);
  } else if (numArr < 1) {
    std::cerr << "no point data found" << std::endl;
    exit(1);
  }
  // setting data array member
  if (arrName != "Z2ErrorSF") da = ds->GetPointData()->GetArray(arrayID);
  // setting name of size field
  std::string array_name = ds->GetPointData()->GetArrayName(arrayID);
  sfname = array_name.append(arrName);

  {  // checking for name conflicts and removing SF with same name if it exists
    vtkSmartPointer<vtkCellData> cd = ds->GetCellData();
    if (cd->GetNumberOfArrays()) {
      for (int i = 0; i < cd->GetNumberOfArrays(); ++i) {
        std::string currname = cd->GetArrayName(i);
        if (sfname == currname) {
          std::cout << "Found size field identifier in cell data: " << currname
                    << std::endl;
          std::cout << "Removing " << currname << " from dataSet" << std::endl;
          ds->GetCellData()->RemoveArray(i);
          break;
        }
      }
    }
  }
  std::cout << "SizeFieldBase constructed" << std::endl;
}

SizeFieldBase *SizeFieldBase::Create(vtkDataSet *_dataSet,
                                     const std::string &method, int arrayID,
                                     double _dev_mult, bool _maxIsmin,
                                     double sizeFactor, int _order) {
  SizeFieldBase *sf;
  if (method == "value") {
    sf = new ValSizeField(_dataSet, arrayID, _dev_mult, _maxIsmin);
  } else if (method == "gradient") {
    sf = new GradSizeField(_dataSet, arrayID, _dev_mult, _maxIsmin);
  } else if (method == "Z2 Error Estimator") {
    sf = new Z2ErrorSizeField(_dataSet, arrayID, _order);
  } else {
    std::cerr << "Specified method " << method << " is not supported\n";
    std::cerr << "Available methods are gradient, value and error" << std::endl;
    exit(1);
  }
  sf->setSizeFactor(sizeFactor);
  sf->computeSizeField(_dataSet->GetPointData()->GetArray(arrayID));
  return sf;
}

std::unique_ptr<SizeFieldBase> SizeFieldBase::CreateUnique(
    vtkDataSet *_dataSet, const std::string &method, int arrayID,
    double _dev_mult, bool _maxIsmin, double _sizeFactor, int _order) {
  return std::unique_ptr<SizeFieldBase>(SizeFieldBase::Create(
      _dataSet, method, arrayID, _dev_mult, _maxIsmin, _sizeFactor, _order));
}

// identifies cells to refine and mutates current size values
// into a compatible size field for the mesh
void SizeFieldBase::mutateValues(std::vector<double> &values) const {
  std::cout << "Size Factor = " << sizeFactor << std::endl;
  // get circumsphere diameter of all cells
  std::vector<double> lengths(ds->GetNumberOfCells());
  for (vtkIdType i = 0; i < ds->GetNumberOfCells(); ++i)
    lengths[i] = std::sqrt(ds->GetCell(i)->GetLength2());
  // find minmax of diameters
  std::vector<double> lengthminmax = nemAux::getMinMax(lengths);
  // redefine minmax values for appropriate size definition reference
  if (maxIsmin)
    lengthminmax[1] = lengthminmax[0];
  else
    lengthminmax[1] *= 0.65;
  lengthminmax[0] -= lengthminmax[0] / 2.;

  // min/max length
  std::cout << "Min Elm Length Scale : " << lengthminmax[0] << "\n"
            << "Max Elm Length Scale : " << lengthminmax[1] << std::endl;
  // get mean and stdev of values
  std::vector<double> meanStdev = nemAux::getMeanStdev(values);
  // get bool array of which cells to refine based on multiplier of stdev
  // std::vector<bool> cells2Refine
  //    = nemAux::cellsToRefine(values, meanStdev[0] + meanStdev[1] * dev_mult);
  // std::vector<bool> cells2Refine
  //    = nemAux::cellsToRefineStdev(values, meanStdev[0],
  //                                 meanStdev[1] * dev_mult);
  std::vector<bool> cells2Refine =
      nemAux::cellsToRefineMaxdev(values, dev_mult);
  // normalize values by mean
  std::vector<double> values_norm = (1. / meanStdev[0]) * values;
  // take the reciprocal of values for size definition (high value -> smaller
  // size)
  if (!nemAux::hasZero(values)) nemAux::reciprocal_vec(values);

  // scale values to min max circumsphere diam of cells
  // now, values represents a size field
  std::vector<double> valuesMinMax = nemAux::getMinMax(values);
  nemAux::scale_vec_to_range(values, valuesMinMax, lengthminmax);

  // setting sizes (values) to f*max element diam based on return of
  // cellsToRefine function
  std::ofstream elmLst;
  elmLst.open("refineCellList.csv");
  if (!elmLst.good()) {
    std::cerr << "Error opening the stream for refinement list." << std::endl;
    exit(1);
  }
  bool isFirstElmIdx = true;
  for (int i = 0; i < values.size(); ++i) {
    if (!cells2Refine[i])
      values[i] = lengthminmax[1];  // if cell shouldn't be refined, size set to
                                    // min of diams
    else {
      if (isFirstElmIdx) {
        elmLst << i;
        isFirstElmIdx = false;
      } else
        elmLst << "," << i;
      values[i] = sizeFactor * values[i];
    }
  }
  elmLst.close();
}

}  // namespace ADP
}  // namespace NEM
