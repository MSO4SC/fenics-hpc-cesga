
H = 0.5;
h = H/75;
h_a = H/15;


//+ Big cube parameters
XMIN = -4.0;
XMAX =  4.0;
YMIN = -4.0;
YMAX =  4.0;
ZMIN = -4.0;
ZMAX =  4.0;


//+ Small cube parameters
xmin = -0.25;
xmax =  0.25;
ymin = -0.25;
ymax =  0.25;
zmin = -0.25;
zmax =  0.25;


Point(1) = {XMIN, YMIN, ZMIN, H};
Point(2) = {XMAX, YMIN, ZMIN, H};
Line(3) = {1, 2};

Point(4) = {xmax,  ymin, zmax, h};
Point(5) = {xmax,  ymax, zmax, h};
Line(4) = {4, 5};


e1[] = Extrude{0, YMAX-YMIN, 0}{Line{3};};
e2[] = Extrude{0, 0, ZMAX-ZMIN}{Surface{e1[1]};};

e3[] = Extrude{-(xmax-xmin),0, 0}{Line{4};};
e4[] = Extrude{0, 0, -(zmax-zmin)}{Surface{e3[1]};};
//Delete extruded volumes
Delete {Volume{e2[1]};} 
Delete {Volume{e4[1]};}


Surface Loop(57) = {25, 8, 17, 21, 30, 29};
Surface Loop(58) = {47, 34, 43, 56, 51, 55};
//Build volume with hole from extruded surfaces
Volume(59) = {57, 58};

//Points of plane in Y = 0
rx_Y=0.5;
rz_Y=0.5;

Point(3) = {  rx_Y,0.0, - rz_Y, 1.0};
Point(30) = { rx_Y,0.0,   rz_Y, 1.0};
Point(31) = {-rx_Y,0.0,   rz_Y, 1.0};
Point(32) = {-rx_Y,0.0, - rz_Y, 1.0};


//Points of plane in X = 0
ry_X=0.35;
rz_X=0.35;
Point(33) = { 0.0, ry_X, - rz_X, 1.0};
Point(34) = { 0.0, ry_X,   rz_X, 1.0};
Point(35) = { 0.0, -ry_X,   rz_X, 1.0};
Point(36) = { 0.0, -ry_X, - rz_X, 1.0};

//Points of plane in Z = 0
rx_Z=0.35;
ry_Z=0.35;
Point(37) = {  rx_Z,  ry_Z,  0.0, 1.0};
Point(38) = { -rx_Z,  ry_Z,  0.0, 1.0};
Point(39) = { -rx_Z, -ry_Z, 0.0, 1.0};
Point(40) = {  rx_Z, -ry_Z, 0.0, 1.0};

//Plan Y=0 Lines
Line(60) = {3, 30};
Line(61) = {30, 31};
Line(62) = {31, 32};
Line(63) = {32, 3};
//Plan X=0 Lines
Line(64) = {33, 34};
Line(65) = {34, 35};
Line(66) = {33, 36};
Line(67) = {36, 35};

//Plan Z=0 Lines
Line(68) = {38, 39};
Line(69) = {39, 40};
Line(70) = {40, 37};
Line(71) = {37, 38};

Field[1] = Attractor;
Field[1].EdgesList = {60, 61, 62, 63 ,64, 65, 66, 67,68,69,70,71};
Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = h_a;
Field[2].LcMax = H;
Field[2].DistMin = 0.25;
Field[2].DistMax = 0.7;

// Use minimum of all the fields as the background field
Field[3] = Min;
Field[3].FieldsList = {2};
Background Field = 3;
