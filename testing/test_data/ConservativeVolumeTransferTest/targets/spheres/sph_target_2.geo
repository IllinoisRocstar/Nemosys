SetFactory("OpenCASCADE");
Sphere(1) = {0.5,0.,0.,0.1};
Physical Volume(1) = { 1 };
Mesh.CharacteristicLengthMin = 0.05;
Mesh.CharacteristicLengthMax = 0.1;
Mesh.SaveAll = 0;
Mesh 3;
Save "sph_target_2.vtk";
