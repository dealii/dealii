// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Distribute FE_WedgeP on a DoFHandler.

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_wedge_p.h>

#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


void
test_3()
{
  const int dim      = 3;
  const int spacedim = 3;

  std::vector<Point<spacedim>> vertices;
  std::vector<CellData<dim>>   cells;

#if true
  vertices.emplace_back(0.0, 0.0, 0.0);
  vertices.emplace_back(1.0, 0.0, 0.0);
  vertices.emplace_back(0.0, 1.0, 0.0);
  vertices.emplace_back(0.0, 0.0, 1.0);
  vertices.emplace_back(1.0, 0.0, 1.0);
  vertices.emplace_back(0.0, 1.0, 1.0);

  {
    CellData<dim> cell;
    cell.vertices = {0, 1, 2, 3, 4, 5};
    cells.push_back(cell);
  }
#else
  vertices.emplace_back(0.0, 0.0, 0.0);
  vertices.emplace_back(1.0, 0.0, 0.0);
  vertices.emplace_back(0.0, 1.0, 0.0);
  vertices.emplace_back(1.0, 1.0, 0.0);
  vertices.emplace_back(0.0, 0.0, 1.0);
  vertices.emplace_back(1.0, 0.0, 1.0);
  vertices.emplace_back(0.0, 1.0, 1.0);
  vertices.emplace_back(1.0, 1.0, 1.0);

  {
    CellData<dim> cell;
    cell.vertices = {0, 1, 2, 4, 5, 6};
    cells.push_back(cell);
  }
  {
    CellData<dim> cell;
    cell.vertices = {2, 1, 3, 6, 5, 7};
    cells.push_back(cell);
  }
#endif

  Triangulation<dim, spacedim> tria;
  tria.create_triangulation(vertices, cells, SubCellData());

  GridOut       grid_out;
  std::ofstream out("mesh.vtk");
  grid_out.write_vtk(tria, out);

  FE_WedgeP<dim, spacedim> fe(2);

  DoFHandler<dim, spacedim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  deallog << dof_handler.n_dofs() << std::endl;

  std::vector<types::global_dof_index> dof_indices;
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      for (const auto face : cell->face_indices())
        {
#if false
          dof_indices.resize(face <= 1 ? 3 : 4);
#else
          dof_indices.resize(fe.n_dofs_per_face(face));
#endif
          cell->face(face)->get_dof_indices(dof_indices);

          for (const auto i : dof_indices)
            deallog << i << ' ';
          deallog << std::endl;
        }
    }
}

int
main()
{
  initlog();

  test_3();
}
