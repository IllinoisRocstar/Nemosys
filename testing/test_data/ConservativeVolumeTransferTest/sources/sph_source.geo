SetFactory("OpenCASCADE");
Sphere(1) = {0,0,0,1.0};
Mesh.CharacteristicLengthMin = 0.1;
Mesh.CharacteristicLengthMax = 0.1;
Mesh 3;
Save "sph_source.vtk";
