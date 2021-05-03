// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2018 by the deal.II authors
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
