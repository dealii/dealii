// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check GridTools::get_vertex_connectivity_of_cells for two different meshes

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include "../tests.h"


template <int dim>
void
test()
{
  // First check connectivity of plain mesh
  Triangulation<dim> tria;
  GridGenerator::subdivided_hyper_cube(tria, 3);
  DynamicSparsityPattern dsp;
  GridTools::get_vertex_connectivity_of_cells(tria, dsp);
  dsp.print(deallog.get_file_stream());

  // Refine a few cells and check again
  tria.begin_active()->set_refine_flag();
  tria.last()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  GridTools::get_vertex_connectivity_of_cells(tria, dsp);
  dsp.print(deallog.get_file_stream());
}


int
main()
{
  initlog();

  test<1>();
  test<2>();
  test<3>();

  return 0;
}
