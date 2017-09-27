//Plane Definition
SetFactory("OpenCASCADE");
z = 0.800000;
lc = 0.01;
Point(1) = {0.000000, 0.000000, 0.800000, lc};
Point(2) = {0.000000, 2.000000, 0.800000, lc};
Point(3) = {2.000000, 2.000000, 0.800000, lc};
Point(4) = {2.000000, 0.000000, 0.800000, lc};
Line(1) = {1,2} ;
Line(2) = {2,3} ;
Line(3) = {3,4} ;
Line(4) = {4,1} ;
Line Loop(1) = {1,2,3,4} ;
Plane Surface(1) = {1};
Plane Surface(2) = {1};
SetFactory("OpenCASCADE");
Sphere(5) = {1.723565,0.958345,0.437759,0.568406};
Sphere(6) = {-0.276435,0.958345,0.437759,0.568406};
Sphere(7) = {1.187667,1.958345,0.365962,0.568406};
Sphere(8) = {1.187667,-0.041655,0.365962,0.568406};
Sphere(9) = {1.651768,0.422446,1.437759,0.568406};
Sphere(10) = {-0.348232,0.422446,1.437759,0.568406};
Sphere(11) = {-0.348232,0.422446,-0.562241,0.568406};
Sphere(12) = {1.651768,0.422446,-0.562241,0.568406};
Sphere(13) = {0.723565,0.886548,1.901861,0.568406};
Sphere(14) = {0.723565,0.886548,-0.098139,0.568406};
Sphere(15) = {0.651768,1.886548,1.365962,0.568406};
Sphere(16) = {0.651768,-0.113452,1.365962,0.568406};
Physical Volume("sphere_1.00") = {5:16};
d = BooleanDifference{ Surface{1}; Delete;}{Volume{5:16}; };
Physical Surface("Domain") = d();
s0 = BooleanIntersection{ Surface{2}; }{Volume{5:16}; Delete;};
Physical Surface("sphere_1.00") = s0();
ltmp = Line In BoundingBox {-0.001000,-0.001000,0.799000,0.001000,2.001000,0.801000};
Physical Line ("Left") = {ltmp()};
ltmp = Line In BoundingBox {1.999000,-0.001000,0.799000,2.001000,2.001000,0.801000};
Physical Line ("Right") = {ltmp()};
ltmp = Line In BoundingBox {-0.001000,1.999000,0.799000,2.001000,2.001000,0.801000};
Physical Line ("Top") = {ltmp()};
ltmp = Line In BoundingBox {-0.001000,-0.001000,0.799000,2.001000,0.001000,0.801000};
Physical Line ("Bottom") = {ltmp()};
Recursive Delete{Surface {2};} 
//Meshhhh
Mesh.Algorithm = 6; //Frontal
Mesh.Optimize = 1; // Gmsh smoother, works with boundary layers (netgen version does not).
Mesh.CharacteristicLengthFromPoints = 0;
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.CharacteristicLengthMax = 0.05;
Mesh.RemeshAlgorithm = 1; //(0) nosplit (1) automatic (2) split only with metis
Mesh 2;
Coherence Mesh;
Save StrCat(StrPrefix(General.FileName), ".msh");
Save StrCat(StrPrefix(General.FileName), ".vtk");
