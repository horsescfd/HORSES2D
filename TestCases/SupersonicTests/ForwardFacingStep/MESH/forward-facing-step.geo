// Inputs

height = 1.0;
length = 3.0;
step_height = 0.2;
step_depth = 0.6;
gridsize = 0.2/10;

// Point coordinates
Point(1) = {0,0,0,gridsize};
Point(2) = {step_depth,0,0,gridsize};
Point(3) = {step_depth,step_height,0,gridsize};
Point(4) = {length,step_height,0,gridsize};
Point(5) = {length,height,0,gridsize};
Point(6) = {0,height,0,gridsize};

Line(7) = {1,2};
Line(8) = {2,3};
Line(9) = {3,4};
Line(10) = {4,5};
Line(11) = {5,6};
Line(12) = {6,1};

Line Loop(13) = {7,8,9,10,11,12};

Plane Surface(14) = 13;

Transfinite Line{7,8,10,12} = 30;
Transfinite Line{9,11} = 70;
Recombine Surface{14};

Physical Line("inlet") = {12};
Physical Line("step") = {7,8,9};
Physical Line("outlet") = {10};
Physical Line("top") = {11};
Physical Surface("interior") = {14};
