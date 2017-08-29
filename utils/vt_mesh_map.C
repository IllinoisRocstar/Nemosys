//standard
#include <iostream>
#include <fstream>
#include<cstring>
#include <sstream>
// third party
#include <vtkAnalyzer.H>
#include <ModelInterface.h>
#include <MAdLib.h>
#include <NodalDataManager.h>
#include <MeshDataBaseIO.h>
#include <CheckOrientation.h>
#include <GmshEntities.h>

// Auxilliary functions
// map between element type and node number 
// as given in GMSH ASCII documentation
int nodes_per_type(int type);
int main(int argc, char* argv[])
{
  if (argc < 4)
	{
		std::cerr << " VTK 3D-2D mesh transfer utility" << std::endl
							<< " Usage: " << argv[0] << 
									" VolVTKFile PlaneVTKFile species_name [write_coords=0|1] [outfile] [material_names]" 
							<< std::endl;
		exit(1);
	}

	// TODO: Testing gmsh/madlib stuff
	if (argc == 7) { 
		std::string mshFileName = argv[6];
		std::ifstream mshStream(mshFileName.c_str());
  	if (!mshStream.good()) {
    	cout << "Error: " << mshFileName<< " not found" << endl;
    	exit(1);
  	}
		std::string line;
		std::vector<std::string> phys_names;
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
					phys_names.push_back(name);
					phys_nums.push_back(phys_num);
				}
				break;
			}
		}
	
		std::vector<int> node_list;
		while (getline(mshStream, line)) {
			if (line.find("$Elements") != -1) {
				getline(mshStream, line);
				std::stringstream currline(line);
				int num_elem;
				currline >> num_elem;
				for (int i = 0; i < num_elem; ++i) {
					getline(mshStream, line);
					std::stringstream nextline(line);
					int elm_num, elm_type, num_tags, phys_tag;
					nextline >> elm_num >> elm_type >> num_tags >> phys_tag;
					// skipping other tags
					for (int j = 0; j < num_tags-1; ++j) { 
						int tmp_tag;
						nextline >> tmp_tag;
					}
					if (phys_tag == 2) { //TODO: Determine RHS based on user input
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
		for (int i = 0; i < node_list.size()/3; ++i)
			std::cout << node_list[i*3] << " " << node_list[i*3 + 1] << " " << node_list[i*3 + 2] << std::endl;		
		
						
							
		
	// Reading gmsh file and info
  	/*MAd::pGModel model = NULL;
		MAd::GM_create(&model,"initial");
  	MAd::pMesh mesh = MAd::M_new(model);
  	MAd::M_load(mesh, mshFileName.c_str());
		std::cout << "Using tags?: " << MAd::GM_physical(model) << std::endl;
	*/
	//	std::cout << model.getNumVertices() << std::endl;

		/*std::vector<std::vector<double>> curds;
		MAd::VIter vit = M_vertexIter(mesh);
		MAd::pVertex pv;
		while (pv = VIter_next(vit)) {
			MAd::pPoint pp = MAd::V_point(pv);
			std::cout << pp->X << " " << pp->Y << " " << pp->Z << std::endl;
		}*/
	
		//MAd::checkMesh(mesh);
  	//MAd::M_info(mesh);
	}


	// Parsing args and Reading files from command line
	vtkAnalyzer* VolMesh;
	vtkAnalyzer* PlaneMesh;
	VolMesh = new vtkAnalyzer(argv[1]);
	PlaneMesh = new vtkAnalyzer(argv[2]);
	VolMesh->read();
	PlaneMesh->read();

	int numComponent;	
	int num_plane_points = PlaneMesh->getNumberOfPoints();
	int num_vol_points = VolMesh->getNumberOfPoints();
	int nDim = 3;

	// get centers of cells in plane mesh
	std::vector<double> PlaneCellCenters = 
	PlaneMesh->getCellCenters(numComponent, nDim);

	// TODO: HOW MANY NEIGHBORS?

	// creating interpolation for cross_link 
	// cross_link 'species' name should be passed by cmdline
	//int species_id = VolMesh->IsArrayName(std::string(argv[3]));
	int species_id = VolMesh->IsArrayName(argv[3]);
	if (species_id >= 0) {
		std::vector<std::vector<double>> volDataMat;
		int numTuple, numComponent;
		VolMesh->getPointDataArray(species_id, volDataMat, numTuple, numComponent);
		std::vector<std::vector<double> > interpData = 
			VolMesh->getInterpData(nDim, 10, numComponent, numTuple, 
											         volDataMat, PlaneCellCenters);

		// check if coordinate writing on and write accordingly
		if (argc >= 5) { 
			int write_coord = atoi(argv[4]);
			switch(write_coord) {
				case 0: { if (argc >= 6)
										VolMesh->writeInterpData(interpData, argv[5]);
									else
										VolMesh->writeInterpData(interpData, std::cout);
									break;
								}
				case 1: { if (argc >= 6)
										VolMesh->writeInterpData(interpData, PlaneCellCenters, nDim, argv[5]);
									else
										VolMesh->writeInterpData(interpData, PlaneCellCenters, nDim, std::cout);
									break;
								}
				default: { std::cerr << "write_coord must be 0|1" << std::endl;
									 std::cerr << "defaulting to 0\n" << std::endl;
									 VolMesh->writeInterpData(interpData,std::cout);
									 break;
								 }
			}
		}
		else
			VolMesh->writeInterpData(interpData, std::cout);
	}
	else {
		std::cerr << "No point data arrays found in " << argv[1] << std::endl;
		exit(1);
	}


	delete VolMesh;
	delete PlaneMesh;
	return EXIT_SUCCESS;
}

// Auxillary Functions
int nodes_per_type(int type)
{
  std::vector<int> nodes = {2,3,4,4,8,6,5,3,6,9,10,27,18,14,1,8,20,15,13,9,10,12,15,15,21,4,5,6,20,35,56,64,125};
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
