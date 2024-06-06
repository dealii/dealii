// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check that GridGenerator::merge_triangulations does not apply any manifold
// unless it is told to do so. Test case adapted from merge_triangulations_03


#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <iostream>

#include "../tests.h"


template <int dim>
void
manifold_info(const Triangulation<dim> &tria, const std::string &case_name)
{
  deallog << "Manifold info for " << case_name << ':' << std::endl
          << " dimension: " << dim << std::endl
          << " no. of cells: " << tria.n_active_cells() << std::endl;

  deallog << "Cell manifolds: ";
  for (const auto &cell : tria.active_cell_iterators())
    {
      deallog << cell->manifold_id() << " ";
      // also check that the manifold can be queried
      (void)cell->get_manifold();
    }
  deallog << std::endl;

  deallog << "Face manifolds: ";
  for (const auto &cell : tria.active_cell_iterators())
    for (const auto face : cell->face_iterators())
      deallog << face->manifold_id() << " ";
  deallog << std::endl << std::endl;
}


template <int dim>
void
test()
{
  Point<dim> center;

  Triangulation<dim> ball;
  Triangulation<dim> shell;
  Triangulation<dim> tria_out;

  // create the two meshes
  GridGenerator::hyper_ball(ball, center, 0.5);
  GridGenerator::hyper_shell(
    shell,
    center,
    0.5,
    1,
    2 * dim,
    true); // colorize flag set to true so outer bnd is 1 and inner is 0

  // set boundaries
  static const SphericalManifold<dim> boundary(center);

  ball.set_manifold(0, boundary);
  shell.set_manifold(0, boundary);
  shell.set_manifold(1, boundary);

  manifold_info(shell, "shell");
  manifold_info(ball, "ball");

  GridGenerator::merge_triangulations(ball, shell, tria_out);
  manifold_info(tria_out, "merged_without_manifold");

  tria_out.clear();
  GridGenerator::merge_triangulations(ball, shell, tria_out, 1e-12, true);
  for (const auto tria : {&ball, &shell})
    for (const auto &cell : tria->active_cell_iterators())
      if (cell->manifold_id() != numbers::flat_manifold_id)
        tria_out.set_manifold(cell->manifold_id(), cell->get_manifold());
  manifold_info(tria_out, "merged_with_manifold");
}


int
main()
{
  initlog();
  test<2>();
  test<3>();

  return 0;
}
