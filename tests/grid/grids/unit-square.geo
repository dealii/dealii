// square_partitioned.geo

// Parameters
L = 1.0;

// Points (with unique IDs)
Point(1) = {0, 0, 0, 1.0};
Point(2) = {L, 0, 0, 1.0};
Point(3) = {L, L, 0, 1.0};
Point(4) = {0, L, 0, 1.0};

// Lines (edges of square)
Line(1) = {1, 2};  // bottom
Line(2) = {2, 3};  // right
Line(3) = {3, 4};  // top
Line(4) = {4, 1};  // left

// Curve loop and surface
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};

// Physical groups for deal.II
Physical Line(11) = {1}; // bottom
Physical Line(12) = {2}; // right
Physical Line(13) = {3}; // top
Physical Line(14) = {4}; // left

Physical Surface(100) = {6}; // domain tag
