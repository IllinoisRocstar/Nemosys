#include <SizeFieldBase.H>
#include <GradSizeField.H>
#include <ValSizeField.H>
#include <Z2ErrorSizeField.H>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <AuxiliaryFunctions.H>

using namespace nemAux;

SizeFieldBase* SizeFieldBase::Create(meshBase* _mesh, std::string method, int arrayID,
                                   double _dev_mult, bool _maxIsmin, double sizeFactor)
{
  if (!method.compare("value"))
  {
    ValSizeField* valsf = new ValSizeField(_mesh, arrayID, _dev_mult, _maxIsmin);
    valsf->setSizeFactor(sizeFactor);
    valsf->computeSizeField(arrayID);
    return valsf;
  }
  else if (!method.compare("gradient"))
  {
    GradSizeField* gradsf = new GradSizeField(_mesh, arrayID, _dev_mult, _maxIsmin);
    gradsf->setSizeFactor(sizeFactor);
    gradsf->computeSizeField(arrayID);
    return gradsf;
  }
  else if (!method.compare("Z2 Error Estimator"))
  {
    Z2ErrorSizeField* z2errorsf = new Z2ErrorSizeField(_mesh, arrayID);
    z2errorsf->setSizeFactor(sizeFactor);
    z2errorsf->computeSizeField(arrayID);
    return z2errorsf;
  }
  else
  {
    std::cout << "Specified method " << method << " is not supported" << std::endl;
    std::cout << "Available methods are gradient, value and error" << std::endl;
    exit(1);
  }
  
}

std::unique_ptr<SizeFieldBase> 
SizeFieldBase::CreateUnique(meshBase* _mesh, std::string method, int arrayID, double _dev_mult, 
                            bool _maxIsmin, double _sizeFactor)
{
  std::cout << __FILE__ << __LINE__ << std::endl;
  std::cout << "Size Factor = " << _sizeFactor << std::endl;
  return std::unique_ptr<SizeFieldBase>(
          SizeFieldBase::Create(_mesh,method,arrayID,_dev_mult,_maxIsmin, _sizeFactor));
}
// initializes derived class, only difference is tmp arrName 
void SizeFieldBase::initialize(meshBase* _mesh, int arrayID, double _dev_mult, 
                               bool _maxIsmin, const std::string& arrName)
{
  // setting private vars
  mesh = _mesh;
  dev_mult = _dev_mult;
  maxIsmin = _maxIsmin;
  // checking for point data
  int numArr = mesh->getDataSet()->GetPointData()->GetNumberOfArrays();
  if (arrayID >= numArr)
  {
    std::cout << "ERROR: arrayID is out of bounds" << std::endl;
    std::cout << "There are " << numArr << " point data arrays" << std::endl;
    exit(1);
  }
  else if (numArr < 1)
  {
    std::cout << "no point data found" << std::endl;
    exit(1);
  }
  // setting data array member
  if (arrName.compare("Z2ErrorSF"))
    da = mesh->getDataSet()->GetPointData()->GetArray(arrayID); 
  // setting name of size field
  std::string array_name = mesh->getDataSet()->GetPointData()->GetArrayName(arrayID);
  sfname = array_name.append(arrName);
  
  { // checking for name conflicts and removing SF with same name if it exists
    vtkSmartPointer<vtkCellData> cd = mesh->getDataSet()->GetCellData();
    if (cd->GetNumberOfArrays())
    { 
      for (int i = 0; i < cd->GetNumberOfArrays(); ++i)
      {
        std::string currname = cd->GetArrayName(i);
        if (!sfname.compare(currname))
        {
          std::cout << "Found size field identifier in cell data: " << currname << std::endl;
          std::cout << "Removing " << currname << " from dataSet" << std::endl;
          mesh->unsetCellDataArray(i);
          break;
        }
      }
    }
  }
}

// identifies cells to refine and mutates current size values
// into a compatible size field for the mesh
void SizeFieldBase::mutateValues(std::vector<double>& values)
{
  std::cout << "Size Factor = " << sizeFactor << std::endl;
  // get circumsphere diameter of all cells 
  std::vector<double> lengths = mesh->getCellLengths();
  // find minmax of diameters
  std::vector<double> lengthminmax = getMinMax(lengths);
  // redefine minmax values for appropriate size definition reference
  if (maxIsmin)
    lengthminmax[1] = lengthminmax[0];
  else
    lengthminmax[1] *= 0.65;
  lengthminmax[0] -= lengthminmax[0]/2.; 

  // min/max length
  std::cout << "Min Elm Lenght Scale : " << lengthminmax[0]
            << "\nMax Elm Lenght Scale : " << lengthminmax[1]
            << std::endl;
  // get mean and stdev of values 
  std::vector<double> meanStdev = getMeanStdev(values);
  // get bool array of which cells to refine based on multiplier of stdev
  //std::vector<int> cells2Refine = cellsToRefine(values, meanStdev[0]+meanStdev[1]*dev_mult);
  //std::vector<int> cells2Refine = cellsToRefineStdev(values, meanStdev[0], meanStdev[1]*dev_mult);
  std::vector<int> cells2Refine = cellsToRefineMaxdev(values, dev_mult);
  // normalize values by mean
  std::vector<double> values_norm = (1./meanStdev[0])*values;  
  // take the reciprocal of values for size definition (high value -> smaller size)
  if (!hasZero(values)){
    reciprocal_vec(values);
  }
  // scale values to min max circumsphere diam of cells 
  // now, values represents a size field
  std::vector<double> valuesMinMax = getMinMax(values);
  scale_vec_to_range(values, valuesMinMax, lengthminmax);

  // setting sizes (values) to f*max element diam based on return of cellsToRefine function
  std::ofstream elmLst;
  elmLst.open("refineCellList.csv");
  if (!elmLst.good())
  {
    std::cerr << "Error opening the stream for refinement list.\n";
    exit(1);
  }
  bool isFirstElmIdx = true;
  for (int i = 0; i < values.size(); ++i)
  {
    if (!cells2Refine[i])
      values[i] = lengthminmax[1]; // if cell shouldn't be refined, size set to min of diams
    else
    {
      if (isFirstElmIdx)
      {
        elmLst << i;
        isFirstElmIdx = false;
      }
      else
        elmLst << "," << i;
      values[i] = sizeFactor*values[i];
    }
  }
  elmLst.close();
}
