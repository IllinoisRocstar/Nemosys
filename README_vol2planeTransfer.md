vol2planeTransfer
-----------------
vol2planeTransfer is utility to transfer data on a 3D mesh (output from RocLB)
to an arbitrary planar cross-section of the 3D mesh. For each cell in the
cross-section, a nearest neighbor search filtered by a user defined radius is
performed in the 3D mesh. The data at these neighbors is then averaged, weighted
by the inverse of their distance from the query point on the cross-section, and 
the result is assigned to the query point.

The driver for the utility is utils/vol2planeTransfer.C

The following classes were modified or added during development:

- vtkAnalyzer
- baseInterp
- spheres (added)

I didn't modify the original interfaces, only overloaded
and added members.
 
The executable 'vol2planeTransfer' is placed in Nemosys/build/bin. 
Append this path to $PATH so the executable name is within scope of the shell.


### Building ###

All utilities, including vol2planeTransfer will be built if cmake is invoked as
follows from within the build directory:

```
$ CMAKE_PREFIX_PATH=../install/madlib:../install/gmsh:../install/cgns cmake ..
$ make -j8 (change number of cores if needed)
```
To turn off (or turn back on) utility building, pass the 'BUILD_UTILS' option:

```
$ CMAKE_PREFIX_PATH=../install/madlib:../install/gmsh:../install/cgns cmake -DBUILD_UTILS=OFF ..
$ make -j8
```

To test the utility, run:
```
$	make test
```
### Usage ###

```
$ vol2planeTransfer file.inp
```
### Input File Specification ###

The name of the input file is the only argument passed to the utility.
Below is an example of the required specifications in the file with
descriptions of each option:

``` 

	Vol_vti_File          =   {vtk_diff1_002000.vti}    // 3D mesh file
 	Vol_geo_File          =   {TEST.geo}								// 2D geo file (with inclusion definitions)
	Plane_vtk_File        =   {TEST.vtk}              	// 2D mesh file
	mask_file             =   {vtk_bbmask_diff1.vti}		// 3D bbmask outupt with RocLB
	outputFile            =   {has_sphere_has_coord1}   // output file name
	material_names        =   {sphere_1.00,sphere_2.00} // names of physical surfaces
	cross_link_name       =   {species3_conc}           // name of crosslink array in 3D mesh file
	write_coords          =   {1}                       // option to write coordinates
	has_spheres           =   {1}                       // option to consider spherical inclusions
	youngs_inc_default    =   {70,-1}                   // youngs moduli defaults for each physical surface 
	shear_inc_default     =   {26,36}										// shear moduli defaults for each physical surface
	poisson_inc_default   =   {-1,100}                  // poisson ratio defaults for each physical surcace
	youngs_dom_default    =   {-1}											// youngs modulus default for domain (non-inclusions)
	poisson_dom_default   =   {105}                     // poisson ratio default for domain (non-inclusions)
	M_weight              =   {1}                       // Molecular weight of primary polymer before cross linking
	Mc_weight             =   {0}                       // Molecular weight of polymer chain between cross links
	NN_TOL                =   {.01}											// Multiplier of min-extent to define radius in k-NN search
```

NOTE: The order of the options matters and must be as given above.

NOTE: DO NOT put spaces between elements in a list.

NOTE: The order of default values MUST correspond to the order of
			the physical names listed in material_names.

NOTE: A default value of "-1" indicates an unknown parameter that 
			will be calculated by the utility and requires non-negative
			values for the other two defaults.

############################ Output Specification ################################

The output file will have at most 7 columns, as in:

```
id    X    Y    Z    rho    E    V
 |    |    |    |     |     |    |
```

id 	- index of cell on cross-section
X		- X coordinate of cell center
Y		- Y coordinate of cell center
Z		- Z coordinate of cell center
rho	- cross-link density at cell center
E		- Young's modulus at cell center
V		- Poisson ratio at cell center

If write_coords is set to 0 in the input file, the output will 
exclude the X, Y and Z columns. 

An example case can be found in testing/test_data/vol2planeTransfer_test