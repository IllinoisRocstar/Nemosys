//standard
#include <iostream>
#include <fstream>
#include<cstring>
#include <sstream>
// third party
#include <spheres.H>
#include <vtkAnalyzer.H>
#include <ModelInterface.h>
#include <MAdLib.h>
#include <NodalDataManager.h>
#include <MeshDataBaseIO.h>
#include <CheckOrientation.h>
#include <GmshEntities.h>

/********************************************************
 * Usage: vol2planeTransfer file.inp                    *
 * G = rho*R*T(1 - Mc/M) is calculated                  *
 * Given E by user, V = EG/(3*(3*G - E)                 *
 * Given V by user, E = 2G(1 + V)                       *
 *******************************************************/

/************************************************************************/
/* TODO: Check if plane is contained in volume and act accordingly if not
   TODO: Generalize for n Dimensions. nDim is useless atm  
   TODO: Unit conversion
   TODO: Refactor by storing opened objects in inputs class and redefining function 
         prototypes accordingly
   TODO: Since the determination of a whether a plane point is in the sphere is made during the interpolation,
         we can modify the return values of interpolate to include a bbmask type array for the plane.
         This can be used in the writeInterpData functions instead of the determinations currently made there. 
   RECENT DEVELOPMENT:
    - added support to consider radius in nn search of k-d tree
      as an overload to interpolate function
    - improved parameter passing (const ref)
    - implemented extent discovery to set tol in radius for nn search 
    - added user input field for fraction of min extent to consider
      for radius in nn search 
    - reading material names from geo file 
    - added support for specified/unspecified material and domain params 
    - added exception handling for when neighbors of points on plane outside of
      inclusion are inside inclusion 
    - no longer considers inclusion points in weight calculation
    - moved sphere checks to interpolate and removed during write*/
/************************************************************************
 *Auxilliary classes and functions
************************************************************************/

/* map between element type and node number 
   as given in GMSH ASCII documentation */
int nodes_per_type(int type);

/* get array of physical names and numeric identifiers 
   from gmsh .msh file*/
std::vector<std::string> get_phys_names(std::string mshFileName);

/* get array of numeric identifiers for physical names
   from gmsh .msh file */
std::vector<int> get_phys_nums(std::string mshFileName);

/* pair struct for physical constants */
struct physical_constants
{
  long double R;
  double T;  
};

//unit conversion consts to convert RocLB unit to SI kg/m/s  
struct unit_convt_consts {
    //const double convt_conc_rocLB = 10.0;  //10nmol/mm^3 to mol/m^3
    //const double convt_time_rocLB = 2.592e6; //month to second
    //const double convt_len_rocLB = 1.e-3;   //mm to m
    double convt_conc_rocLB;   
    double convt_time_rocLB;  
    double convt_len_rocLB;  
};

// class to read input file
class inputs {
public:
  inputs(){};
  ~inputs(){};
  inputs(std::string input_file);
  void print();
  void validate();
public:
  std::string vol_file;
  std::string geo_file;
  std::string plane_file;
  std::string mask_file;
  std::string out_file;
  vector<std::string> material_names;
  std::string cross_link_name;
  int write_coords, has_spheres;
  vector<double> youngs_inc_default;
  vector<double> shear_inc_default;
  vector<double> poisson_inc_default;
  double youngs_dom_default;
  double poisson_dom_default;
  double M_weight, Mc_weight, NN_TOL;
  double T;
  bool writePlaneMesh;
  double len_convt;
  double conc_convt;
  double t_convt;
};

/*************************************************************************/

int main(int argc, char* argv[])
{
  if (argc < 2) {
    std::cerr << "3D-2D mesh transfer utility" << std::endl
              << " Usage: " << argv[0] << " inputFile" << std::endl << std::endl
              << "inputFile must specify the following, in order: " << std::endl
              << "Vol_vti_File" << std::endl 
              << "geo_File" << std:: endl
              << "Plane_vtk_File" << std::endl
              << "bbmask_vti_File" << std::endl
              << "outputFile" << std::endl
              << "material_names" << std::endl
              << "cross_link_name" << std::endl
              << "write_coords" << std::endl
              << "has_spheres" << std::endl
              << "youngs_inc_default" << std::endl
              << "shear_inc_default" << std::endl
              << "poisson_inc_default" << std::endl
              << "youngs_dom_default" << std::endl
              << "poisson_dom_default" << std::endl
              << "M_weight" << std::endl
              << "Mc_weight" << std::endl
              << "NN_TOL" << std::endl
              << "Temperature" << std::endl
              << "writePlaneMesh" << std::endl
              << "len_convt" << std::endl
              << "conc_convt" << std::endl
              << "t_convt" << std::endl;
    exit(1);
  }

  // parsing input file and reading/loading stuff 
  inputs inp(argv[1]);
  inp.validate();
  // read geo file for sphere locations
  sphere_string spherestring = readSpheres(inp.geo_file);
  vector<sphere> spheres = spherestring.spheres;
  vector<std::string> mat_sphere_names = spherestring.strings;
  // read mesh data
  vtkAnalyzer* VolMesh;
  vtkAnalyzer* PlaneMesh;
  vtkAnalyzer* maskMesh;
  VolMesh = new vtkAnalyzer((char*) &(inp.vol_file)[0u]);
  PlaneMesh = new vtkAnalyzer((char*) &(inp.plane_file)[0u]);
  maskMesh = new vtkAnalyzer((char*) &(inp.mask_file)[0u]);
  VolMesh->read();
  PlaneMesh->read();
  maskMesh->read();

  int numComponent; 
  int num_plane_points = PlaneMesh->getNumberOfPoints();
  int num_vol_points = VolMesh->getNumberOfPoints();
  int nDim = 3;

  // TEST STUFF

  /*std::vector<std::vector<double>> maskDatas_tmp;
  int tmpTuples, tmpComponents;
  maskMesh->getPointDataArray(0, maskDatas_tmp, tmpTuples, tmpComponents);
  std::vector<double> maskData_tmp(tmpTuples);
  for (int i = 0; i < tmpTuples; ++i)
    maskData_tmp[i] = maskDatas_tmp[i][0];
  
  vtkDataSet* data = maskMesh->getdataSet();
  
  for (int i = 0; i < maskMesh->getNumberOfCells(); ++i) {
    vtkIdList* point_ids = data->GetCell(i)->GetPointIds();
    for (int j = 0; j < point_ids->GetNumberOfIds() ; ++j) {
      if (maskData_tmp[point_ids->GetId(j)] == 1.0) {
        double* pntcoords = maskMesh->getPointCoords(point_ids->GetId(j));
        std::cerr << pntcoords[0] << " " 
                  << pntcoords[1] << " " 
                  << pntcoords[2] << std::endl;
      }
    }
  }

   // Reading gmsh file and info
    std::string mshFileName = argv[6];
    std::ifstream mshStream(mshFileName.c_str());

    MAd::pGModel model = NULL;
    MAd::GM_create(&model,"initial");
    MAd::pMesh mesh = MAd::M_new(model);
    MAd::M_load(mesh, mshFileName.c_str());
    std::cout << "Using tags?: " << MAd::GM_physical(model) << std::endl;
  
  //  std::cout << model.getNumVertices() << std::endl;

    std::vector<std::vector<double>> curds;
    MAd::VIter vit = M_vertexIter(mesh);
    MAd::pVertex pv;
    while (pv = VIter_next(vit)) {
      MAd::pPoint pp = MAd::V_point(pv);
      std::cout << pp->X << " " << pp->Y << " " << pp->Z << std::endl;
    }
  
    //MAd::checkMesh(mesh);
    //MAd::M_info(mesh);

*/


  // get centers of cells in plane mesh
  std::vector<double> PlaneCellCenters = 
  PlaneMesh->getCellCenters(numComponent);
  
  // get all coordinates in vol mesh
  std::vector<double> VolPointCoords = 
  VolMesh->getAllPointCoords(nDim);

  // get min extent
  double minExtent = VolMesh->getMinExtent(nDim, VolPointCoords);

  // setting tol for knn search in k-d tree
  double tol = inp.NN_TOL*minExtent;

  // defining physical constants
  physical_constants phys_const;
  //phys_const.R = .0000000083144598; // in m^3*GPa/(K*mol)
  phys_const.R = 8.3144598; // in J/K/mol)
  phys_const.T = inp.T;
  
  // unit conversion consts
  unit_convt_consts unit_consts;
  unit_consts.convt_conc_rocLB = inp.conc_convt;
  
  // creating interpolation for cross_link 
  int species_id = VolMesh->IsArrayName(inp.cross_link_name);
  if (species_id >= 0) {
    std::vector<std::vector<double>> volDataMat;
    std::vector<std::vector<double>> maskDatas;
    int numTuple, numComponent, tmpTuple, tmpComponent;
    VolMesh->getPointDataArray(species_id, volDataMat, numTuple, numComponent);
    maskMesh->getPointDataArray(0, maskDatas, tmpTuple, tmpComponent);
    std::vector<double> maskData(tmpTuple);
    for (int i = 0; i < tmpTuple; ++i)
      maskData[i] = maskDatas[i][0];
  

    std::vector<std::vector<double> > interpData;
    // check if spheres are present and do accordingly
    switch(inp.has_spheres) {
      case 0: { interpData=
                  VolMesh->getInterpData(nDim, 10, numComponent, numTuple,
                                         volDataMat, PlaneCellCenters,
                                         VolPointCoords, tol);
                
                VolMesh->writeInterpData(interpData, inp.Mc_weight, 
                                         inp.M_weight, inp.youngs_dom_default,
                                         inp.poisson_dom_default,
                                         phys_const.T, phys_const.R,
                                         PlaneCellCenters, nDim, 
                                         inp.out_file, (bool) inp.write_coords, PlaneMesh->getNumberOfNonTri(),
                                         unit_consts.convt_conc_rocLB);
                // writing plane cell data to vtk
                if (inp.writePlaneMesh) 
                {
                  PlaneMesh->setCellDataArray("crosslink", 1, interpData[0]); 
                  PlaneMesh->write("crosslink_on_plane.vtu");
                }
                break;
              }

      case 1: { interpData=
                  VolMesh->getInterpData(nDim, 10, numComponent, numTuple,
                                         volDataMat, PlaneCellCenters,
                                         VolPointCoords, spheres, maskData, tol);
                VolMesh->writeInterpData(interpData, inp.Mc_weight,
                                         inp.M_weight, inp.youngs_dom_default,
                                         inp.poisson_dom_default,
                                         phys_const.T, phys_const.R,
                                         PlaneCellCenters, nDim, spheres, 
                                         mat_sphere_names, inp.material_names,
                                         inp.youngs_inc_default, inp.shear_inc_default,
                                         inp.poisson_inc_default,
                                         inp.out_file, (bool) inp.write_coords, PlaneMesh->getNumberOfNonTri(),
                                         unit_consts.convt_conc_rocLB);

                if (inp.writePlaneMesh)
                {
                  // writing plane cell data to vtk
                  PlaneMesh->setCellDataArray("crosslink",1,interpData[0]); 
                  PlaneMesh->write("crosslink_on_plane.vtu");
                }
                break;
              }
      default: {  std::cerr << "has_spheres must be 0|1" << std::endl;
                  std::cerr << "correct in input file " << std::endl;
                  exit(1);
               }
      

    }
  }
  else {
    std::cerr << "No point data arrays found in " << inp.vol_file << std::endl;
    exit(1);
  }



  delete VolMesh;
  delete PlaneMesh;
  delete maskMesh;
  return EXIT_SUCCESS;

}

/******************************************************************************/

// Auxillary Functions

// map between element type and node number 
// as given in GMSH ASCII documentation
int nodes_per_type(int type)
{
  std::vector<int> nodes = {2,3,4,4,8,6,5,3,6,9,10,
                            27,18,14,1,8,20,15,13,9,
                            10,12,15,15,21,4,5,6,20,
                            35,56,64,125};
  if (type < 92) 
    return nodes[type-1];
  else if (type == 92) 
    return nodes[nodes.size() - 2]; 
  else if (type == 93) 
    return nodes[nodes.size() - 1]; 
  else {
    std::cerr << "Invalid type identifier" << std::endl;
    std::cerr << "See GMSH ASCII Docs for Valid Types" << std::endl;
    exit(1);
  }
}


// get array of physical names and numeric identifiers 
// from gmsh .msh file
std::vector<std::string> get_phys_names(std::string mshFileName)
{
  std::ifstream mshStream(mshFileName.c_str());
  if (!mshStream.good()) {
    cout << "Error: " << mshFileName<< " not found" << endl;
    exit(1);
  }
  std::string line;
  std::vector<std::string> phys_names;
  while(getline(mshStream, line)) {
    if (line.find("$PhysicalNames") != -1) {
      getline(mshStream, line);
      std::stringstream currline(line);
      int num_names;
      currline >> num_names;
      for (int i = 0; i < num_names; ++i) {
        getline(mshStream, line);
        std::stringstream nextline(line);
        std::string name;
        int phys_dim;
        int phys_num;

        nextline >> phys_dim >> phys_num >> name;
        phys_names.push_back(name);
      }   
      break;
    }   
  }
    mshStream.close();
    return phys_names;
}


// get array of numeric identifiers for physical names
// from gmsh .msh file
std::vector<int> get_phys_nums(std::string mshFileName)
{   
  std::ifstream mshStream(mshFileName.c_str());
  if (!mshStream.good()) {
    cout << "Error: " << mshFileName<< " not found" << endl;
    exit(1);
  }
  std::string line;
  std::vector<int> phys_nums;
  while(getline(mshStream, line)) {
    if (line.find("$PhysicalNames") != -1) {
      getline(mshStream, line);
      std::stringstream currline(line);
      int num_names;
      currline >> num_names;
      for (int i = 0; i < num_names; ++i) {
        getline(mshStream, line);
        std::stringstream nextline(line);
        std::string name;
        int phys_dim;
        int phys_num;

        nextline >> phys_dim >> phys_num >> name;
        phys_nums.push_back(phys_num);
      }
      break;
    }
  }
  mshStream.close();
  return phys_nums;
}

// input reader constructor with input file
using std::string; using std::vector; using std::size_t;
inputs::inputs(string input_file)
{
  std::ifstream inputStream(input_file.c_str());
  if (!inputStream.good()) {
    std::cout << "Error reading " << input_file << std::endl;
    exit(1);
  }
  string line;
  int i = 0;
  while (getline(inputStream, line)) {
    if (line.find("=") != -1) {
      size_t beg = line.find("{");
      size_t end = line.find("}");
      if (beg != -1 && end != -1) {
        string data = line.substr(beg+1, end-beg-1);
        switch(i) {
          case 0: vol_file=data;
                  break;
          case 1: geo_file=data;
                  break;
          case 2: plane_file=data;
                  break;
          case 3: mask_file=data;
                  break;
          case 4: out_file=data;
                  break;
          case 5: { std::stringstream ss(data);
                    while (ss.good()) {
                      string tmp;
                      getline(ss, tmp, ',');
                      material_names.push_back(tmp);
                    }
                    break;
                  } 
          case 6: cross_link_name=data;
                  break;
          case 7: { std::stringstream ss(data);
                    ss >> write_coords;
                    break;
                  } 
          case 8: { std::stringstream ss(data);
                    ss >> has_spheres;
                    break;
                  } 
          case 9: { std::stringstream ss(data);
                    double tmp; 
                    while(ss >> tmp) {
                      youngs_inc_default.push_back(tmp);
                      if (ss.peek() == ',')
                        ss.ignore();
                    }
                    break;
                  } 
          case 10:  { std::stringstream ss(data);
                      double tmp; 
                      while(ss >> tmp) {
                      shear_inc_default.push_back(tmp);
                        if (ss.peek() == ',')
                          ss.ignore();
                      }
                    break;
                    } 
          case 11:  { std::stringstream ss(data);
                      double tmp; 
                      while(ss >> tmp) {
                      poisson_inc_default.push_back(tmp);
                        if (ss.peek() == ',')
                          ss.ignore();
                      }
                    break;
                    } 
          case 12:  { std::stringstream ss(data);
                      ss >> youngs_dom_default;
                      break;
                    } 
          case 13:  { std::stringstream ss(data);
                      ss >> poisson_dom_default;
                      break;
                    } 
          case 14:  { std::stringstream ss(data);
                      ss >> M_weight;
                      break;
                    } 
          case 15:  { std::stringstream ss(data);
                      ss >> Mc_weight;
                      break;
                    }
          case 16:  { std::stringstream ss(data);
                      ss >> NN_TOL;
                      break;
                    }
          case 17:  { std::stringstream ss(data);
                      ss >> T;
                      break;
                    }
          case 18:  { std::stringstream ss(data);
                      ss >> writePlaneMesh; 
                      break;
                    }
          case 19:  { std::stringstream ss(data);
                      ss >> conc_convt;
                      break;
                    }
          case 20:  { std::stringstream ss(data);
                      ss >> len_convt;
                      break;
                    }
          case 21:  { std::stringstream ss(data);
                      ss >> t_convt;
                      break;
                    }
        }
      }
      i+=1;
    }
  }
}

void inputs::validate()
{
  if (material_names.size() != youngs_inc_default.size() ||
      material_names.size() != shear_inc_default.size() ||
      material_names.size() != poisson_inc_default.size()) {
    std::cerr << "Inconsistency with default specification" << std::endl;
    exit(3);
  }

  for (int i = 0; i < material_names.size(); ++i) {
    if (youngs_inc_default[i] == -1 && poisson_inc_default[i] == -1) {
      std::cerr << "Either Youngs or Poisson defaults must be specified" << std::endl; 
      exit(3);
    }
  }
    
  if (youngs_dom_default == -1 && poisson_dom_default == -1) {
    std::cerr << "Either Youngs or Poisson defaults must be specified" << std::endl; 
    exit(3);
    }
}


void inputs::print()
{
  std::cout << "Inputs: " << std::endl
            << "Vol_file: " << vol_file << std::endl
            << "geo_file: " << geo_file << std::endl
            << "plane_file: " << plane_file << std::endl
            << "mask_file: " << mask_file << std::endl
            << "out_file: " << out_file << std::endl;
  std::cout << "Material names: ";
  for (int i = 0; i < material_names.size() ; ++i) 
    std::cout << material_names[i] << " ";
  std::cout << std::endl
            << "cross_link name: " << cross_link_name << std::endl
            << "write coords: " << write_coords << std::endl
            << "has spheres: " << has_spheres << std::endl;
  std::cout << "youngs_inc_default: ";
  for (int i = 0; i < youngs_inc_default.size(); ++i)
    std::cout << youngs_inc_default[i] << " ";
  std::cout << std::endl;
  std::cout << "shear_inc_default: ";
  for (int i = 0; i < shear_inc_default.size(); ++i)
    std::cout << shear_inc_default[i] << " ";
  std::cout << std::endl;

  std::cout << "poisson_inc_default: ";
  for (int i = 0; i < poisson_inc_default.size(); ++i)
    std::cout << poisson_inc_default[i] << " ";
  std::cout << std::endl;

  std::cout << "youngs default: " << youngs_dom_default << std::endl
            << "poisson default: " << poisson_dom_default << std::endl
            << "M_weight: " << M_weight << std::endl
            << "Mc_weight: " << Mc_weight << std::endl
            << "NN_TOL: " << NN_TOL << std::endl
            << "Temperature: " << T << std::endl
            << "write Plane Mesh" << writePlaneMesh << std::endl; 
}
