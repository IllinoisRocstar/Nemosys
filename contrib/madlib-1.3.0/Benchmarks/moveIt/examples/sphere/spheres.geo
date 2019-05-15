lc = 0.3;//0.3
lc2 = 0.3;//0.3
R = 1.0;
R2 = 3.;

Point(1) = {0,0,0,lc};
Point(2) = {R,0,0,lc};
Point(3) = {0,R,0,lc};
Point(4) = {0,-R,0,lc};
Point(5) = {-R,0,0,lc};
Point(6) = {0,0,R,lc};
Point(7) = {0,0,-R,lc};

Circle(1) = {2,1,3};
Circle(2) = {3,1,5};
Circle(3) = {5,1,4};
Circle(4) = {4,1,2};
Circle(5) = {2,1,6};
Circle(6) = {6,1,3};
Circle(7) = {3,1,7};
Circle(8) = {7,1,4};
Circle(9) = {4,1,6};
Circle(10) = {6,1,5};
Circle(11) = {5,1,7};
Circle(12) = {7,1,2};

Point(8) = {R2,0,0,lc2};
Point(9) = {0,R2,0,lc2};
Point(10) = {0,-R2,0,lc2};
Point(11) = {-R2,0,0,lc2};
Point(12) = {0,0,R2,lc2};
Point(13) = {0,0,-R2,lc2};

Delete {
  Point{1};
}
Circle(13) = {8, 1, 9};
Circle(14) = {8, 1, 13};
Circle(15) = {13, 1, 11};
Circle(16) = {11, 1, 12};
Circle(17) = {12, 1, 8};
Circle(18) = {8, 1, 10};
Circle(19) = {10, 1, 11};
Circle(20) = {11, 1, 9};
Circle(21) = {10, 1, 12};
Circle(22) = {12, 1, 9};
Circle(23) = {9, 1, 13};
Circle(24) = {13, 1, 10};
Line Loop(25) = {12, 1, 7};
Ruled Surface(26) = {25};
Line Loop(27) = {7, -11, -2};
Ruled Surface(28) = {27};
Line Loop(29) = {2, -10, 6};
Ruled Surface(30) = {29};
Line Loop(31) = {6, -1, 5};
Ruled Surface(32) = {31};
Line Loop(33) = {12, -4, -8};
Ruled Surface(34) = {33};
Line Loop(35) = {8, -3, 11};
Ruled Surface(36) = {35};
Line Loop(37) = {3, 9, 10};
Ruled Surface(38) = {37};
Line Loop(39) = {5, -9, 4};
Ruled Surface(40) = {39};
Line Loop(41) = {13, 23, -14};
Ruled Surface(42) = {41};
Line Loop(43) = {15, 20, 23};
Ruled Surface(44) = {43};
Line Loop(45) = {20, -22, -16};
Ruled Surface(46) = {45};
Line Loop(47) = {22, -13, -17};
Ruled Surface(48) = {47};
Line Loop(49) = {17, 18, 21};
Ruled Surface(50) = {49};
Line Loop(51) = {14, 24, -18};
Ruled Surface(52) = {51};
Line Loop(53) = {24, 19, -15};
Ruled Surface(54) = {53};
Line Loop(55) = {16, -21, 19};
Ruled Surface(56) = {55};
Surface Loop(57) = {50, 48, 46, 44, 54, 52, 42, 56};
Surface Loop(58) = {38, 36, 34, 26, 32, 30, 28, 40};
Volume(59) = {57, 58};

Physical Surface(1) = {26, 28, 30, 32, 34, 36, 38, 40};

Physical Volume(2) = {59};