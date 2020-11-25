SetFactory("OpenCASCADE");

//------------ Parameters --------------//
//--------------------------------------//
x = 1;
y = 1;

num_x = 2;
num_y = 2;

Rectangle(1) = {0,0,0,x,y};

Transfinite Line {4,2} = num_y;
Transfinite Line {1,3} = num_x;

Transfinite Surface {1};

