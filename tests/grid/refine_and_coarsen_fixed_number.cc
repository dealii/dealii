// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2017 by the deal.II Authors
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

// This test ensures that all cells get refined (resp. coarsened) when
// the value 1.0 is passed as top_fraction (resp. bottom_fraction) to
// the function GridRefinement::refine_and_coarsen_fixed_number.

#include "../tests.h"
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/lac/vector.h>

using namespace dealii;

int
main(int argc, const char *argv[])
{
  initlog();

  Triangulation<2> tria;
  GridGenerator::hyper_cube (tria);
  tria.refine_global(4);

  Vector<float> indicator(tria.n_active_cells());
  for (int i = 0; i != indicator.size(); ++i)
    {
      indicator[i] = i;
    }

  deallog << "n_active_cells: " << tria.n_active_cells() << std::endl;

  GridRefinement::refine_and_coarsen_fixed_number(tria, indicator, 1.0, 0.0);
  tria.execute_coarsening_and_refinement();

  deallog << "n_active_cells: " << tria.n_active_cells() << std::endl;

  indicator.reinit(tria.n_active_cells());
  for (int i = 0; i != indicator.size(); ++i)
    {
      indicator[i] = i;
    }

  GridRefinement::refine_and_coarsen_fixed_number(tria, indicator, 0.0, 1.0);
  tria.execute_coarsening_and_refinement();

  deallog << "n_active_cells: " << tria.n_active_cells() << std::endl;

  return 0;
}
