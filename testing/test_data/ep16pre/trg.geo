SetFactory("OpenCASCADE");

Box(0) = {-1,-1,-1, 2,2,2};  // 0
v[] = Volume "*";
size = #v[];
// Translate to (0,0,0)
//Translate {1,1,1}{Volume{v[]};}

//v[] = BooleanFragments{Volume{v[]}; Delete;}{};
p[] = PointsOf{Volume{v[]};};
Characteristic Length{p[]} =0.5;
 
Mesh.ElementOrder = 1;
Mesh.Algorithm = 6;   // Frontal
Mesh.Algorithm3D = 1;   // Delaunay
Mesh.CharacteristicLengthMin = 0.001;
Mesh.CharacteristicLengthMax = 1;
Mesh.OptimizeThreshold = 1.0;
Mesh.SaveAll = 1;
