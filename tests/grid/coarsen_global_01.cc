// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test Triangulation::coarsen_global

#include <deal.II/base/geometry_info.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include "../tests.h"

template <int dim>
void
check(int number)
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  deallog << "dim = " << dim << ", number = " << number << std::endl;

  deallog << "before n_active_cells: " << tria.n_active_cells() << std::endl;

  tria.coarsen_global(number);

  deallog << "after n_active_cells: " << tria.n_active_cells() << std::endl;
}


int
main()
{
  initlog();

  check<1>(1);
  check<1>(2);
  check<1>(3);

  check<2>(1);
  check<2>(2);
  check<2>(3);

  check<3>(1);
  check<3>(2);
  check<3>(3);
}
