// square-quad.geo
// Create points of the square
Point(1) = {0, 0, 0, 0.1};
Point(2) = {1, 0, 0, 0.1};
Point(3) = {1, 1, 0, 0.1};
Point(4) = {0, 1, 0, 0.1};

// Create lines connecting the points
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Create the line loop and surface
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Define a structured grid (transfinite)
Transfinite Line {1, 3} = 10; // 10 divisions horizontally
Transfinite Line {2, 4} = 10; // 10 divisions vertically
Transfinite Surface {1};
Recombine Surface {1}; // Make quads instead of triangles

// Physical groups for deal.II
Physical Curve(11) = {1, 2, 3, 4};

// Material id for the domain
Physical Surface(111) = {1};

// Optional: force quad meshing globally
Mesh.RecombineAll = 1;
