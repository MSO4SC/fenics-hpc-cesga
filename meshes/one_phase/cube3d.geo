// Gmsh project created on Mon Dec 19 12:10:33 2016

//+ Element length near the big cube
H = 0.8;

//+ Element length near the big cube
h = 0.01;

//+ Big cube parameters
XMIN = -4.0;
XMAX =  4.0;
YMIN = -2.0;
YMAX =  2.0;
ZMIN = -2.0;
ZMAX =  2.0;

//+ Small cube parameters
xmin = -0.25;
xmax =  0.25;
ymin = -0.25;
ymax =  0.25;
zmin = -0.25;
zmax =  0.25;


//+ BEGIN: BIG CUBE DEFINITION
Point(1) = {XMIN, YMIN, ZMAX, H};
Point(2) = {XMIN, YMAX, ZMAX, H};
Point(3) = {XMAX, YMAX, ZMAX, H};
Point(4) = {XMAX, YMIN, ZMAX, H};
Point(5) = {XMIN, YMIN, ZMIN, H};
Point(6) = {XMIN, YMAX, ZMIN, H};
Point(7) = {XMAX, YMAX, ZMIN, H};
Point(8) = {XMAX, YMIN, ZMIN, H};
//+
Line(1) = {4, 3};
Line(2) = {3, 2};
Line(3) = {2, 1};
Line(4) = {1, 4};
//+
Line Loop(5) = {3, 4, 1, 2};
Plane Surface(6) = {5};
//+
//+
Line(7) = {7, 8};
Line(8) = {8, 5};
Line(9) = {5, 6};
Line(10) = {6, 7};
Line(11) = {3, 7};
Line(12) = {2, 6};
Line(13) = {4, 8};
Line(14) = {1, 5};
//+
Line Loop(15) = {11, -10, -12, -2};
Plane Surface(16) = {15};
//+
Line Loop(17) = {3, 14, 9, -12};
Plane Surface(18) = {17};
//+
Line Loop(19) = {11, 7, -13, 1};
Plane Surface(20) = {19};
//+
Line Loop(21) = {4, 13, 8, -14};
Plane Surface(22) = {21};
//+
Line Loop(23) = {10, 7, 8, 9};
Plane Surface(24) = {23};
//+ END: BIG CUBE DEFINITION

//+ BEGIN: SMALL CUBE DEFINITION
Point(9) = {xmin, ymin, zmax, h};
Point(10) = {xmin, ymax, zmax, h};
Point(11) = {xmax, ymax, zmax, h};
Point(12) = {xmax, ymin, zmax, h};
Point(13) = {xmin, ymin, zmin, h};
Point(14) = {xmin, ymax, zmin, h};
Point(15) = {xmax, ymax, zmin, h};
Point(16) = {xmax, ymin, zmin, h};
//+
Line(25) = {13, 14};
Line(26) = {14, 15};
Line(27) = {15, 16};
Line(28) = {16, 13};
Line(29) = {13, 9};
Line(30) = {9, 12};
Line(31) = {12, 16};
Line(32) = {9, 10};
Line(33) = {10, 11};
Line(34) = {11, 12};
Line(35) = {10, 14};
Line(36) = {11, 15};
//+
Line Loop(37) = {32, 35, -25, 29};
Plane Surface(38) = {37};
//+
Line Loop(39) = {28, 25, 26, 27};
Plane Surface(40) = {39};
//+
Line Loop(41) = {29, 30, 31, 28};
Plane Surface(42) = {41};
//+
Line Loop(43) = {35, 26, -36, -33};
Plane Surface(44) = {43};
//+
Line Loop(45) = {31, -27, -36, 34};
Plane Surface(46) = {45};
//+
Line Loop(47) = {30, -34, -33, -32};
Plane Surface(48) = {47};
//+ END: SMALL CUBE DEFINITION

//+ Define the whole volume
Surface Loop(49) = {18, 6, 22, 20, 16, 24};
Surface Loop(50) = {40, 42, 38, 48, 46, 44};
Volume(51) = {49, 50};
