#include <SizeFieldBase.H>
#include <GradSizeField.H>
#include <ValSizeField.H>

SizeFieldBase* SizeFieldBase::Create(meshBase* _mesh, std::string method, int arrayID,
                                   double _dev_mult, bool _maxIsmin)
{
  if (!method.compare("value"))
  {
    ValSizeField* valsf = new ValSizeField(_mesh, arrayID, _dev_mult, _maxIsmin);
    valsf->computeSizeField(arrayID);
    return valsf;
  }
  else if (!method.compare("gradient"))
  {
    GradSizeField* gradsf = new GradSizeField(_mesh, arrayID, _dev_mult, _maxIsmin);
    gradsf->computeSizeField(arrayID);
    return gradsf;
  }
  else if (!method.compare("error")){}
  else
  {
    std::cout << "Specified method " << method << " is not supported" << std::endl;
    std::cout << "Available methods are gradient, value and error" << std::endl;
    exit(1);
  }
  
}

// identifies cells to refine and mutates current size values
// into a compatible size field for the mesh
void SizeFieldBase::mutateValues(std::vector<double>& values)
{
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

  // get mean and stdev of values 
  std::vector<double> meanStdev = getMeanStdev(values);
  // get bool array of which cells to refine based on multiplier of stdev
  std::vector<int> cells2Refine = cellsToRefine(values, meanStdev[0]+meanStdev[1]*dev_mult);
  // normalize values by mean
  std::vector<double> values_norm = (1./meanStdev[0])*values;  
  // take the reciprocal of values for size definition (high value -> smaller size)
  reciprocal_vec(values);
  // scale values to min max circumsphere diam of cells 
  // now, values represents a size field
  std::vector<double> valuesMinMax = getMinMax(values);
  scale_vec_to_range(values, valuesMinMax, lengthminmax);
  // setting sizes (values) to f*max element diam based on return of cellsToRefine function
  for (int i = 0; i < values.size(); ++i)
  {
    if (!cells2Refine[i])
      values[i] = lengthminmax[1]; // if cell shouldn't be refined, size set to min of diams
  }
}
