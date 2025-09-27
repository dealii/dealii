// L-shaped domain

Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {1, 0.5, 0, 1.0};
Point(4) = {0.5, 0.5, 0, 1.0};
Point(5) = {0.5, 1, 0, 1.0};
Point(6) = {0, 1, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};

Line Loop(1) = {1, 2, 3, 4, 5, 6};
Plane Surface(1) = {1};

// --- Physical groups ---
// Whole domain as one material
Physical Surface(1) = {1};

// Boundaries
Physical Line(1) = {1, 6, 5}; // Dirichlet
Physical Line(2) = {2, 3, 4}; // Neumann

// Mesh size (triangles only)
Mesh.CharacteristicLengthMin = 0.1;
Mesh.CharacteristicLengthMax = 0.2;
