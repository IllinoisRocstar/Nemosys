SetFactory("OpenCASCADE");
Sphere(1) = {0,0,0,.995};
Mesh.CharacteristicLengthMin = 0.03;
Mesh.CharacteristicLengthMax = 0.07;
Mesh 2;
Save "sph_target.vtk";
