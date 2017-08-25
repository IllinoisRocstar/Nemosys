//standard
#include <iostream>
// third party
#include <vtkAnalyzer.H>

int main(int argc, char* argv[])
{
  if (argc < 4)
	{
		std::cerr << " VTK 3D-2D mesh transfer utility" << std::endl
							<< " Usage: " << argv[0] << " VolMeshFile PlaneMeshFile species_name" 
							<< std::endl;
		exit(1);
	}

	// Reading files from command line
	vtkAnalyzer* VolMesh;
	vtkAnalyzer* PlaneMesh;
	VolMesh = new vtkAnalyzer(argv[1]);
	PlaneMesh = new vtkAnalyzer(argv[2]);
	VolMesh->read();
	PlaneMesh->read();
	
	// Report mesh statistics
	VolMesh->report();
	std::cout << std::endl;
	PlaneMesh->report();
	
	int numComponent;	
	//int num_plane_cells = PlaneMesh->getNumberOfCells();
	int num_plane_points = PlaneMesh->getNumberOfPoints();
	int num_vol_points = VolMesh->getNumberOfPoints();
	int nDim = 3;

	// get centers of cells in plane mesh
	std::vector<double> PlaneCellCenters = 
	PlaneMesh->getCellCenters(numComponent, nDim);

	//get coordinates of points of cells in vol mesh
	std::vector<double> VolPointCoords = 
	VolMesh->getAllPointCoords(nDim);		

	// TODO: HOW MANY NEIGHBORS?

	// creating interpolation for cross_link 
	// cross_link 'species' name should be passed by cmdline
	int species_id = VolMesh->IsArrayName(std::string(argv[3]));
	if (species_id >= 0) {
		basicInterpolant* VolPointInterp = 
			new basicInterpolant(nDim, num_vol_points, 10, VolPointCoords);
		std::vector<std::vector<double>> volDataMat;
		int numTuple, numComponent;
		VolMesh->getPointDataArray(species_id, volDataMat, numTuple, numComponent);
		for (int j = 0; j < numComponent; ++j) {
			std::vector<double > volData(numTuple);
			for (int i = 0; i < numTuple; ++i) {
					volData[i] = volDataMat[i][j];
			}
			std::vector<double> interpData;
			VolPointInterp->interpolate(PlaneCellCenters.size()/nDim, 
																	PlaneCellCenters, volData, interpData);
			
			for (int i = 0; i < PlaneCellCenters.size()/nDim; ++i)
				std::cout << interpData[i] << std::endl;			

		}
		delete VolPointInterp;
	}



	delete VolMesh;
	delete PlaneMesh;
	return EXIT_SUCCESS;
}

