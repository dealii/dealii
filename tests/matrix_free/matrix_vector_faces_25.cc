// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



// tests matrix-free face evaluation, matrix-vector products as compared to
// the same implementation with MeshWorker. Similar to matrix_vector_faces_04
// but using a mesh with all sorts of different orientations, using a mesh
// somewhat related to the one in tests/grid/grid_tools_06.cc but with the
// cubes patched together appropriately

#include <deal.II/base/function.h>

#include <deal.II/fe/fe_dgq.h>

#include "../tests.h"

#include "matrix_vector_faces_common.h"



void generate_grid(Triangulation<3> &triangulation, int orientation)
{
  Point<3>              vertices_1[] = {Point<3>(-1., -1., -3.),
                           Point<3>(+1., -1., -3.),
                           Point<3>(-1., +1., -3.),
                           Point<3>(+1., +1., -3.),
                           Point<3>(-1., -1., -1.),
                           Point<3>(+1., -1., -1.),
                           Point<3>(-1., +1., -1.),
                           Point<3>(+1., +1., -1.),
                           Point<3>(-1., -1., +1.),
                           Point<3>(+1., -1., +1.),
                           Point<3>(-1., +1., +1.),
                           Point<3>(+1., +1., +1.)};
  std::vector<Point<3>> vertices(&vertices_1[0], &vertices_1[12]);

  std::vector<CellData<3>> cells(2, CellData<3>());

  /* cell 0 */
  int cell_vertices_0[GeometryInfo<3>::vertices_per_cell] = {
    0, 1, 2, 3, 4, 5, 6, 7};

  /* cell 1 */
  int cell_vertices_1[8][GeometryInfo<3>::vertices_per_cell] = {
    {4, 5, 6, 7, 8, 9, 10, 11},
    {5, 7, 4, 6, 9, 11, 8, 10},
    {7, 6, 5, 4, 11, 10, 9, 8},
    {6, 4, 7, 5, 10, 8, 11, 9},
    {9, 8, 11, 10, 5, 4, 7, 6},
    {8, 10, 9, 11, 4, 6, 5, 7},
    {10, 11, 8, 9, 6, 7, 4, 5},
    {11, 9, 10, 8, 7, 5, 6, 4}};

  for (const unsigned int j : GeometryInfo<3>::vertex_indices())
    {
      cells[0].vertices[j] = cell_vertices_0[j];
      cells[1].vertices[j] = cell_vertices_1[orientation][j];
    }
  cells[0].material_id = 0;
  cells[1].material_id = 0;


  triangulation.create_triangulation(vertices, cells, SubCellData());
}



template <int dim, int fe_degree>
void
test()
{
  if (dim == 2)
    return;

  Triangulation<3> tria;
  for (unsigned int orientation = 0; orientation < 8; ++orientation)
    {
      deallog << "Testing orientation case " << orientation << std::endl;
      tria.clear();
      generate_grid(tria, orientation);

      FE_DGQ<3>     fe(fe_degree);
      DoFHandler<3> dof(tria);
      dof.distribute_dofs(fe);
      AffineConstraints<double> constraints;
      constraints.close();

      do_test<3, fe_degree, fe_degree + 1, double>(dof, constraints, false);

      tria.refine_global(1);
      dof.distribute_dofs(fe);
      do_test<3, fe_degree, fe_degree + 1, double>(dof, constraints, false);
    }
}
