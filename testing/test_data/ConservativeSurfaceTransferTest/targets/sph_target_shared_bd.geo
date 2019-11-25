SetFactory("OpenCASCADE");
Sphere(1) = {0,0,0,1.0};
Mesh.CharacteristicLengthMin = 0.03;
Mesh.CharacteristicLengthMax = 0.07;
Mesh 2;
Save "sph_target_shared_bd.vtk";
