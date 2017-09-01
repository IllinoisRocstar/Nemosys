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
 * Usage: vt_mesh_map file.inp                          *
 * G = rho*R*T(1 - Mc/M) is calculated                  *
 * Given E by user, V = EG/(3*(3*G - E)                 *
 * Given V by user, E = 2G(1 + V)                       *
 * CURRENTLY ASSUMING USER PROVIDES YOUNG'S MODULUS (E) * 
 *******************************************************/

/************************************************************************/
/* TODO: Check if plane is contained in volume and act accordingly if not
	 TODO: input file validation		
	 TODO: If point is in sphere, return index of sphere and map that to
				 the material name- pull material names from geo file 
	 TODO: Generalize for n Dimensions. nDim is useless atm  
	 TODO: Add support for checking plane points in 
				 sphere shell using extents as ref
	
	 RECENT DEVELOPMENT:
		- added support to consider radius in nn search of k-d tree
		  as an overload to interpolate function
		- improved parameter passing (const ref)
		- implemented extent discovery to set tol in radius for nn search 
		- added user input field for fraction of min extent to consider
			for radius in nn search */
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

// class to read input file
class inputs {
public:
	inputs(){};
	~inputs(){};
	inputs(std::string input_file);
public:
  std::string vol_file;
  std::string msh_file;
  std::string geo_file;
  std::string plane_file;
  std::string out_file;
  vector<std::string> material_names;
  std::string cross_link_name;
  int write_coords, has_spheres;
  vector<double> youngs_default;
  vector<double> shear_default;
  double M_weight, Mc_weight, NN_TOL;
};

/*************************************************************************/

int main(int argc, char* argv[])
{
	if (argc < 2) {
		std::cerr << "3D-2D mesh transfer utility" << std::endl
							<< " Usage: " << argv[0] << " inputFile" << std::endl << std::endl
							<< "inputFile must specify the following, in order: " << std::endl
							<< "Vol_vti_File" << std::endl 
							<< "Vol_msh_File" << std::endl
							<< "Vol_geo_File" << std:: endl
							<< "Plane_vtk_File" << std::endl
							<< "outputFile" << std::endl
							<< "material_names" << std::endl
							<< "cross_link_name" << std::endl
							<< "write_coords" << std::endl
							<< "has_spheres" << std::endl
							<< "youngs_default" << std::endl
							<< "shear_default" << std::endl
							<< "M_weight" << std::endl
							<< "Mc_weight" << std::endl
							<< "NN_TOL" << std::endl;
		exit(1);
	}

	// parsing input file and reading/loading stuff	
	inputs inp(argv[1]);
	vector<sphere> spheres = readSpheres(inp.geo_file);
	vtkAnalyzer* VolMesh;
	vtkAnalyzer* PlaneMesh;
	VolMesh = new vtkAnalyzer((char*) &(inp.vol_file)[0u]);
	PlaneMesh = new vtkAnalyzer((char*) &(inp.plane_file)[0u]);
	VolMesh->read();
	PlaneMesh->read();

	int numComponent;	
	int num_plane_points = PlaneMesh->getNumberOfPoints();
	int num_vol_points = VolMesh->getNumberOfPoints();
	int nDim = 3;

	// get centers of cells in plane mesh
	std::vector<double> PlaneCellCenters = 
	PlaneMesh->getCellCenters(numComponent, nDim);
	
	// get all coordinates in vol mesh
	std::vector<double> VolPointCoords = 
	VolMesh->getAllPointCoords(nDim);

	// get min extent
	double minExtent = VolMesh->getMinExtent(nDim, VolPointCoords);

	// setting tol for knn search in k-d tree
	double tol = inp.NN_TOL*minExtent;

	// creating interpolation for cross_link 
	// cross_link 'species' name should be passed by cmdline
	//int species_id = VolMesh->IsArrayName(std::string(argv[3]));
	int species_id = VolMesh->IsArrayName(inp.cross_link_name);
	if (species_id >= 0) {
		std::vector<std::vector<double>> volDataMat;
		int numTuple, numComponent;
		VolMesh->getPointDataArray(species_id, volDataMat, numTuple, numComponent);

		std::vector<std::vector<double> > interpData;
		// check if spheres are present and do accordingly
		switch(inp.has_spheres) {
			case 0: {	interpData=
									VolMesh->getInterpData(nDim, 10, numComponent, numTuple,
																				 volDataMat, PlaneCellCenters,
																				 VolPointCoords, tol);
								
								VolMesh->writeInterpData(interpData, inp.Mc_weight, 
								  											 inp.M_weight, inp.youngs_default[0],
																				 PlaneCellCenters, nDim, 
																				 inp.out_file, (bool) inp.write_coords);
								break;
							}

			case 1: {	interpData=
									VolMesh->getInterpData(nDim, 10, numComponent, numTuple,
																				 volDataMat, PlaneCellCenters,
																				 VolPointCoords, spheres, tol);
								VolMesh->writeInterpData(interpData, inp.Mc_weight,
																				 inp.M_weight, inp.youngs_default[0],
																				 PlaneCellCenters, nDim, spheres, 
																				 inp.out_file, (bool) inp.write_coords);
								break;
							}
			default: {	std::cerr << "has_spheres must be 0|1" << std::endl;
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
	return EXIT_SUCCESS;

		//writeSpheres(spheres, std::cout);
	
	// TODO: Testing gmsh/madlib stuff
	/*if (argc == 7) { 
		std::string mshFileName = argv[6];
		std::vector<std::string> phys_names = get_phys_names(mshFileName);
		std::vector<int> phys_nums = get_phys_nums(mshFileName);		

		std::ifstream mshStream(mshFileName.c_str());
    if (!mshStream.good()) {
      cout << "Error: " << mshFileName<< " not found" << endl;
      exit(1);
    }
	
		std::string line;
		std::vector<int> node_list;
		std::vector<int> geom_list;
		while (getline(mshStream, line)) {
			if (line.find("$Elements") != -1) {
				getline(mshStream, line);
				std::stringstream currline(line);
				int num_elem;
				currline >> num_elem;
				for (int i = 0; i < num_elem; ++i) {
					getline(mshStream, line);
					std::stringstream nextline(line);
					int elm_num, elm_type, num_tags, phys_tag, geom_tag;
					nextline >> elm_num >> elm_type >> num_tags >> phys_tag >> geom_tag;
					// skipping other tags
					for (int j = 0; j < num_tags-2; ++j) { 
						int tmp_tag;
						nextline >> tmp_tag;
					}
					if (phys_tag == 2) { Determine RHS based on user input
						geom_list.push_back(geom_tag);
						for (int k = 0; k < nodes_per_type(elm_type); ++k) { 
							int node;
							nextline >> node;	
							node_list.push_back(node);
						}
					}
				}
				break;
			}
		}

		std::cout << node_list.size() << std::endl;
		for (int i = 0; i < phys_names.size(); ++i)
			std::cout << phys_names[i] << std::endl;
		for (int i = 0; i < node_list.size()/3; ++i) {
			std::cout << node_list[i*3] << " " << node_list[i*3 + 1] << " " << node_list[i*3 + 2] <<  " " << geom_list[i] << std::endl;		
		}*/
		// TODO: Should the inclusions have interpolated values?
		
	//}
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
          case 1: msh_file=data;
                  break;
          case 2: geo_file=data;
                  break;
          case 3: plane_file=data;
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
                      youngs_default.push_back(tmp);
                      if (ss.peek() == ',')
                        ss.ignore();
                    }
                    break;
                  } 
          case 10:  { std::stringstream ss(data);
                      double tmp; 
                      while(ss >> tmp) {
                      shear_default.push_back(tmp);
                        if (ss.peek() == ',')
                          ss.ignore();
                      }
                    break;
                    } 
          case 11:  { std::stringstream ss(data);
                      ss >> M_weight;
                      break;
                    } 
          case 12:  { std::stringstream ss(data);
                      ss >> Mc_weight;
                      break;
                    }
					case 13: 	{	std::stringstream ss(data);
											ss >> NN_TOL;
											break;
									 	}
        }
      }
    }
    i+=1;
  }
}

