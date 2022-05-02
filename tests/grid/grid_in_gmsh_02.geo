meshsize = 2;

Point(1) = {-1, 0, 0, meshsize};
Point(2) = {1, 0, 0, meshsize};
Point(3) = {1, 2, 0, meshsize};
Point(4) = {-1, 2, 0, meshsize};
Point(5) = {-1, 4, 0, meshsize};
Point(6) = {1, 4, 0, meshsize};
Line(1) = {5, 6};
Line(2) = {4, 3};
Line(3) = {1, 2};
Line(4) = {5, 4};
Line(5) = {4, 1};
Line(6) = {6, 3};
Line(7) = {3, 2};

Line Loop(9) = {1, 6, -2, -4};
Plane Surface(9) = {9};
Line Loop(11) = {3, -7, -2, 5};
Plane Surface(11) = {11};

Physical Curve("1") = {1};
Physical Curve("2") = {3};
Physical Surface(1) = {9, 11};

Mesh.SubdivisionAlgorithm = 1;
