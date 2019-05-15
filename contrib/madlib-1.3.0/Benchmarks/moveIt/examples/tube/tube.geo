// -*- C++ -*-

lc1 = 0.1;//0.1
lc2 = 0.2;//0.4
lcBox = 0.6;//0.6

// cylinders parameters
r = 0.9;
R = 1.;
dR = 0.2;
L1 = 3.;
L2 = 3.;
dL = 1.;

// external box parameters
B1 = 2.;
B2 = 2.;
B3 = 3.;

// first cylinder (small one)
Point(1) = {0,0,0,lc1};
Point(2) = {r,0,0,lc1};
Point(3) = {0,r,0,lc1};
Point(4) = {0,-r,0,lc1};
Point(5) = {-r,0,0,lc1};
Circle(1) = {2,1,3};
Circle(2) = {3,1,5};
Circle(3) = {5,1,4};
Circle(4) = {4,1,2};
Extrude {0,0,L1} {
  Line{3,2,1,4};
}
Line Loop(21) = {13,9,5,17};
Plane Surface(22) = {21};

// tube
Point(20) = {R,0,L1+dL,lc2};
Point(21) = {-R,0,L1+dL,lc2};
Point(22) = {0,R,L1+dL,lc2};
Point(23) = {0,-R,L1+dL,lc2};
Point(24) = {0,-R-dR,L1+dL,lc2};
Point(25) = {0,+R+dR,L1+dL,lc2};
Point(26) = {R+dR,0,L1+dL,lc2};
Point(27) = {-R-dR,0,L1+dL,lc2};
Point(28) = {0,0,L1+dL,lc2};
Circle(29) = {23,28,20};
Circle(30) = {20,28,22};
Circle(31) = {22,28,21};
Circle(32) = {21,28,23};
Circle(33) = {24,28,26};
Circle(34) = {26,28,25};
Circle(35) = {25,28,27};
Circle(36) = {27,28,24};
Extrude {0,0,L2} {
  Line{36,33,34,35,31,30,29,32};
}
Line Loop(69) = {45,49,37,41};
Line Loop(70) = {53,65,61,57};
Plane Surface(71) = {69,70};
Line Loop(104) = {35,36,33,34};
Line Loop(105) = {32,29,30,31};
Plane Surface(106) = {104,105};

// // external box
Point(62) = {B3,B3,-B1,lcBox};
Point(63) = {-B3,B3,-B1,lcBox};
Point(64) = {-B3,-B3,-B1,lcBox};
Point(65) = {B3,-B3,-B1,lcBox};
Line(82) = {64,63};
Line(83) = {63,62};
Line(84) = {62,65};
Line(85) = {65,64};
Extrude {0,0,B1+L1+dL+L2+B2} {
  Line{82,83,84,85};
}
Line Loop(123) = {107,111,115,119};
Plane Surface(124) = {123};
Line Loop(125) = {2,3,4,1};
Plane Surface(126) = {125};
Line Loop(130) = {84,85,82,83};
Plane Surface(131) = {130};
Surface Loop(132) = {124,110,131,118,122,114};
Surface Loop(133) = {52,106,40,44,48,71,60,56,68,64};
Surface Loop(134) = {12,126,8,20,16,22};
Volume(135) = {132,133,134};
