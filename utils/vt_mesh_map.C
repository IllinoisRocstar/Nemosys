//standard
#include <iostream>
// third party
#include <vtkAnalyzer.H>

int main(int argc, char* argv[])
{
  if (argc < 4)
	{
		std::cerr << " VTK 3D-2D mesh transfer utility" << std::endl
							<< " Usage: " << argv[0] << 
									" VolMeshFile PlaneMeshFile species_name [write_coords=0|1] [outfile]" 
							<< std::endl;
		exit(1);
	}

	// Parsing args and Reading files from command line
	vtkAnalyzer* VolMesh;
	vtkAnalyzer* PlaneMesh;
	VolMesh = new vtkAnalyzer(argv[1]);
	PlaneMesh = new vtkAnalyzer(argv[2]);
	VolMesh->read();
	PlaneMesh->read();

	
	
	// Report mesh statistics
//	std::cout << "Reporting on VolMesh:" << std::endl;
//	VolMesh->report();
//	std::cout << std::endl;
//	std::cout << "Reporting on PlaneMesh:" << std::endl;
//	PlaneMesh->report();
//	std::cout << std::endl;
	
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
	int species_id = VolMesh->IsArrayName(std::string(argv[3]));
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
				case 0: { if (argc == 6)
										VolMesh->writeInterpData(interpData, argv[5]);
									else
										VolMesh->writeInterpData(interpData, std::cout);
									break;
								}
				case 1: { if (argc==6)
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

