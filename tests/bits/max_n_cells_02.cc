// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2007 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test the max_n_cells argument to
// GridRefinement::refine_and_coarsen_fixed_fraction


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"



template <int dim>
void
test()
{
  deallog << dim << "d:" << std::endl;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  for (unsigned int cycle = 0; cycle < 7 * (4 - dim) * (4 - dim); ++cycle)
    {
      deallog << "cycle=" << cycle << ", n_cells=" << tria.n_active_cells()
              << std::endl;

      Vector<float> criteria(tria.n_active_cells());
      for (unsigned int i = 0; i < tria.n_active_cells(); ++i)
        criteria(i) = i;

      GridRefinement::refine_and_coarsen_fixed_fraction(
        tria, criteria, 0.8, 0.03, 10000);
      tria.execute_coarsening_and_refinement();
    }
}



int
main()
{
  initlog();

  test<1>();
  test<2>();
  test<3>();
}
