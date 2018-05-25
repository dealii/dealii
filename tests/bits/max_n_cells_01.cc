// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// test the max_n_cells argument to
// GridRefinement::refine_and_coarsen_fixed_number


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

      GridRefinement::refine_and_coarsen_fixed_number(
        tria, criteria, 1. / 4, 1. / 64, 10000);
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
