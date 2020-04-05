// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2018 by the deal.II authors
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

// When setting refine_fraction=0.6 and coarsen_fraction=0.4, these
// numbers may add up to slightly larger than a total of 1.0 in
// floating point arithmetic. Make sure we allow this nevertheless.

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

#include <limits>

#include "../tests.h"


int
main(int argc, const char *argv[])
{
  initlog();

  Triangulation<2> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(4);

  Vector<float> indicator(tria.n_active_cells());
  for (unsigned int i = 0; i != indicator.size(); ++i)
    {
      indicator[i] = i;
    }

  deallog << "n_active_cells: " << tria.n_active_cells() << std::endl;

  GridRefinement::refine_and_coarsen_fixed_number(
    tria,
    indicator,
    0.6,
    1. - 0.6 + std::numeric_limits<double>::epsilon() * 5);
  tria.execute_coarsening_and_refinement();

  deallog << "n_active_cells: " << tria.n_active_cells() << std::endl;

  indicator.reinit(tria.n_active_cells());
  for (unsigned int i = 0; i != indicator.size(); ++i)
    {
      indicator[i] = i;
    }

  GridRefinement::refine_and_coarsen_fixed_number(
    tria,
    indicator,
    0.4,
    1. - 0.4 + std::numeric_limits<double>::epsilon() * 5);
  tria.execute_coarsening_and_refinement();

  deallog << "n_active_cells: " << tria.n_active_cells() << std::endl;

  return 0;
}
