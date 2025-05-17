// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test neighbor_child_on_subface for triangles

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>

#include "../tests.h"


template <int dim>
void
test()
{
  Triangulation<dim> triangulation;
  GridGenerator::subdivided_hyper_cube_with_simplices(triangulation, 1);
  triangulation.begin_active()->set_refine_flag(
    RefinementCase<dim>::isotropic_refinement);
  triangulation.execute_coarsening_and_refinement();


  auto cell = triangulation.begin(0);
  ++cell;

  auto child_1 = cell->face(1)->child(0);
  auto child_2 = cell->face(1)->child(1);

  auto child_cell_1 = cell->neighbor_child_on_subface(1, 0);
  auto child_cell_2 = cell->neighbor_child_on_subface(1, 1);

  deallog << child_cell_1->vertex(1) - child_1->vertex(0) << " "
          << child_cell_1->vertex(2) - child_1->vertex(1) << std::endl;

  deallog << child_cell_2->vertex(1) - child_2->vertex(0) << " "
          << child_cell_2->vertex(2) - child_2->vertex(1) << std::endl;
}

int
main()
{
  using namespace dealii;
  initlog();
  test<2>();

  return 0;
}
