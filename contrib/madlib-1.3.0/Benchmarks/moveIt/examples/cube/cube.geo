lc = 0.3;
//nz = 23;//mettre n+1 points
C = 1.;

Point(1) = {0.0,0.0,0.0,lc};
Point(2) = {0.0,  C,0.0,lc};
Point(3) = {0.0,  C,C  ,lc};
Point(4) = {0.0,0.0,C  ,lc};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(4) = {1,2,3,4};
Plane Surface(5) = {4};
//Transfinite Line {1,2,3,4} = nz; // Using Bump 3.5;
//Transfinite Surface {5} = {1,2,3,4};
//Recombine Surface {5};



//Extrude Surface {5, {1,0,0.0}}
//{
//  // Recombine;
//   Layers { nz-1, 9000, 1 } ;
//};

//Physical Line(1) = {1,2,3,4};
//Physical Surface(99) = {5};


Extrude {C,0,0} {
  Surface{5};
}


//Recombine Surface (1);



Physical Surface(1) = {5,18,14};
Physical Surface(2) = {27,26,22};
Physical Volume(99) = {1};