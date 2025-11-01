//---------------------------------------------------------------------------
// In this section we create the points that make up the outer domain

// Compute the length of the line closest to the edge of the slit:
b = .2 * Tan((60 * Pi)/180);

// Points on the outer domain
Point(1) = {-3.5, -3.36, 0, 1.0};
Point(2) = {3.5, -3.36, 0, 1.0};
Point(3) = {3.5, 3.36, 0, 1.0};
Point(4) = {-3.5, 3.36, 0, 1.0};
Point(5) = {-3.5, 0.2, 0, 1.0};
Point(6) = {-0.98 - b, 0.2, 0, 1.0};
Point(7) = {-0.98, 0, 0, 1.0};
Point(8) = {-0.98 - b, -0.2, 0, 1.0};
Point(9) = {-3.5, -0.2, 0, 1.0};
Point(10) = {3.5, 0, 0, 1.0}; // Tip of slit on left hand side

// The 2 points below are for separating the left and right side of the mesh
Point(36) = {-0.98, 3.36, 0, 1.0};
Point(37) = {-0.98, -3.36, 0, 1.0};

// Points for the right hand side slit
Point(107) = {3.5, -1.88, 0, 1.0};
Point(108) = {3.5, -1.48, 0, 1.0};
Point(109) = {1.2, -1.88, 0, 1.0};
Point(110) = {1.8, -1.48, 0, 1.0};

// Points for the triangle on the right hand side
Point(111) = {3.5, 1.18, 0, 1.0};
Point(112) = {3.5, 2.18, 0, 1.0};
Point(113) = {3, 1.68, 0, 1.0};

//---------------------------------------------------------------------------
// This section contains the lines that make up the outer domain

// Lines for the left hand side of the domain
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 1};
Line(10) = {7, 10}; // Line that bisects domain from top/bottom

// Lines on the bottom/top domain
Line(38) = {1, 37};
Line(39) = {37, 2};
Line(40) = {3, 36};
Line(41) = {36, 4};

// Lines that bisect the domain from left/right
Line(42) = {37, 7};
Line(43) = {7, 36};

// Lines on the right hand side of the domain
Line(52) = {2, 107};
Line(53) = {107, 109};
Line(54) = {109, 110};
Line(55) = {110, 108};
Line(56) = {108, 10};
Line(57) = {10, 111};
Line(58) = {112, 3};
Line(59) = {111, 113};
Line(60) = {113, 112};

//---------------------------------------------------------------------------
// Here are the points and lines that create the two holes in the domain

// Relevant points for bottom left hole
Point(31) = {-2.1, -1.54, 0, 1.0}; // origin
Point(32) = {-2.1, -0.84, 0, 1.0}; // top point
Point(33) = {-2.1, -2.24, 0, 1.0}; // bottom point 
Point(103) = {-1.4, -1.54, 0, 1.0}; // right point
Point(104) = {-2.8, -1.54, 0, 1.0}; // left point

// Bottom left hole arcs
Circle(48) = {103, 31, 32};
Circle(49) = {32, 31, 104};
Circle(50) = {104, 31, 33};
Circle(51) = {33, 31, 103};

// Similarly, the upper right circle
radius1 = .5;
Point(19) = {0.42, 2.0, 0, 1.0}; // origin 
Point(20) = {0.42, 2.0 - radius1, 0, 1.0};
Point(21) = {0.42 + radius1, 2.0, 0, 1.0};
Point(22) = {0.42, 2.0 + radius1, 0, 1.0};
Point(23) = {0.42 - radius1, 2.0, 0, 1.0};
Circle(12) = {20, 19, 21};
Circle(13) = {21, 19, 22};
Circle(14) = {22, 19, 23};
Circle(15) = {23, 19, 20};

//---------------------------------------------------------------------------
// This section describes the "Plane Surfaces", i.e., the 2D surfaces for meshing

// The surface of the top right partition
Curve Loop(1) = {40, -43, 10, 57, 59, 60, 58};
Curve Loop(2) = {12, 13, 14, 15};
Plane Surface(1) = {1, 2};

// Surface of the bottom right partition
Curve Loop(3) = {39, 52, 53, 54, 55, 56, -10, -42};
Plane Surface(2) = {3};

// Surface of bottom left mesh
Curve Loop(4) = {42, 7, 8, 9, 38};
Curve Loop(5) = {49, 50, 51, 48};
Plane Surface(3) = {4, 5};

// Surface of top left mesh
Curve Loop(6) = {6, 43, 41, 4, 5};
Plane Surface(4) = {6};

// Creates a physical surface.  The expression list on the right hand side is the list
// of elementary surfaces created above.  This is what makes our mesh 2D.
Physical Surface(2) = {1, 2, 3, 4};

//---------------------------------------------------------------------------
// This section describes the physical IDs of certain objects.

// Assign boundary ID of 0 to outer boundary
Physical Curve(0) = {38, 39, 52, 53, 54, 55, 56, 57, 59, 60, 58, 40, 41, 4, 5, 6, 7, 8, 9};

// Assign boundary ID of 1 top right circle
Physical Curve(1) = {12};

// Assign boundary IDs of 2 & 3 to the top & bottom half of the bottom left
Physical Curve(2) = {48, 49};
Physical Curve(3) = {50, 51};

//---------------------------------------------------------------------------
// Parameters for the meshing

Mesh.Algorithm = 3;
Mesh.RecombineAll = 1;
Mesh.CharacteristicLengthFactor = .6;
Mesh.SubdivisionAlgorithm = 1;
Mesh.Smoothing = 20;
Show "*";
