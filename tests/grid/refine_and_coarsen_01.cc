// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check that refine_and_coarsen_fixed_fraction behaves correctly
// if all the indicators are the same. This was fixed in 28858.

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include "../tests.h"



template <int dim>
void
check()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  Vector<float> estimated_error_per_cell;
  estimated_error_per_cell.reinit(tria.n_active_cells());
  for (unsigned int j = 0; j < estimated_error_per_cell.size(); ++j)
    estimated_error_per_cell(j) = 1.;

  deallog << "n_active_cells: " << tria.n_active_cells() << std::endl;


  GridRefinement::refine_and_coarsen_fixed_fraction(tria,
                                                    estimated_error_per_cell,
                                                    0.25,
                                                    0);
  tria.execute_coarsening_and_refinement();

  deallog << "n_active_cells: " << tria.n_active_cells() << std::endl;

  deallog << "OK for " << dim << 'd' << std::endl;
}


int
main()
{
  initlog();

  check<1>();
  check<2>();
  check<3>();
}
