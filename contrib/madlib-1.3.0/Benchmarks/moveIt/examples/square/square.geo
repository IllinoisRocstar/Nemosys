// -*- C++ -*-

lc = 0.1;
C = 2.;

nx = 5;
ny = 10;

Point(1) = {-C/2,-C/2,0,lc};
Point(2) = {C/2,-C/2,0,lc};
Point(3) = {C/2,C/2,0,lc};
Point(4) = {-C/2,C/2,0,lc};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(5) = {2,3,4,1};
Plane Surface(6) = {5};

//Transfinite Line {1,3} = nx; // permet de couper la ligne en segments reguliers (avec 'nx' noeuds)
//Transfinite Line {2,4} = ny;
//Transfinite Surface {6} = {1,2,3,4}; // la surface sera egalement maillee de maniere reguliere

//Recombine Surface (6); // combine les triangles pour former des rectangles

Physical Point(1000) = {1,2,3,4};
Physical Line(10) = {1,2,3,4};
Physical Surface(100) = {6};
