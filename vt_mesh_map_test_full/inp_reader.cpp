#include<iostream>
#include<fstream>
#include<ostream>
#include<stdlib.h>
#include<string>
#include<vector>
#include<cstring>
#include <sstream>

using std::string; using std::vector; using std::size_t;

class inputs {
public:
	// constructors
	inputs(){}
	~inputs(){};
	inputs(string input_file)
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
						case 5: {	std::stringstream ss(data);
											while (ss.good()) {
												string tmp;
												getline(ss, tmp, ',');
												material_names.push_back(tmp);
												// if comma encountered, move on
												//if (ss.peek() == ',')
												//	ss.ignore();

											}
											break;
										}	
						case 6: cross_link_name=data;
										break;
						case 7: {	std::stringstream ss(data);
											ss >> write_coords;
											break;
										}
						case 8: {	std::stringstream ss(data);
											ss >> has_spheres;
											break;
										}
						case 9:	{ std::stringstream ss(data);
											double tmp;
											while(ss >> tmp) {
												youngs_default.push_back(tmp);
												if (ss.peek() == ',')
													ss.ignore();
											}
											break;
										}
						case 10:	{	std::stringstream ss(data);
												double tmp;
												while(ss >> tmp) {
												shear_default.push_back(tmp);
													if (ss.peek() == ',')
														ss.ignore();
												}
											break;
											}
						case 11:	{	std::stringstream ss(data);
										 		ss >> M_weight;
												break;
											}
						case 12: 	{	std::stringstream ss(data);
												ss >> Mc_weight;
												break;
											}
					}
				}
			}
			i+=1;
		}
	}	

public:
	string vol_file;
	string msh_file;
	string geo_file;
	string plane_file;	
	string out_file;
	vector<string> material_names;
	string cross_link_name;
	int write_coords, has_spheres;
	vector<double> youngs_default;
	vector<double> shear_default;
	double M_weight, Mc_weight; 
};

int main (int argc, char* argv[])
{

	inputs my_inp(argv[1]);
	std::cout << my_inp.vol_file << std::endl;
	
	std::cout << my_inp.msh_file << std::endl;
	std::cout << my_inp.geo_file << std::endl;
	std::cout << my_inp.plane_file << std::endl;
	std::cout << my_inp.out_file << std::endl;
	for (int i = 0; i < my_inp.material_names.size(); ++i)
		std::cout << my_inp.material_names[i] << std::endl;
	std::cout << my_inp.cross_link_name << std::endl;
	std::cout << my_inp.write_coords << std::endl;
	std::cout << my_inp.has_spheres << std::endl;
	for (int i = 0; i < my_inp.youngs_default.size(); ++i) {
		std::cout << my_inp.youngs_default[i] << std::endl;
		std::cout << my_inp.shear_default[i] << std::endl;
	}
	std::cout << my_inp.M_weight << std::endl;
	std::cout << my_inp.Mc_weight << std::endl;
	return 0;
}
