Rectangle(1) = {0.0, 0, 0, 0.5, 1, 0};
Rectangle(2) = {0.5, 0, 0, 0.5, 1, 0};
Recombine Surface{1};
Physical Surface("All") = {1, 2};
Mesh 2;
Coherence Mesh;
Save "box_2D_mixed.msh";
