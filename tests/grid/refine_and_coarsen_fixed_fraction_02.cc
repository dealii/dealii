// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

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

  GridRefinement::refine_and_coarsen_fixed_fraction(
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

  GridRefinement::refine_and_coarsen_fixed_fraction(
    tria,
    indicator,
    0.4,
    1. - 0.4 + std::numeric_limits<double>::epsilon() * 5);
  tria.execute_coarsening_and_refinement();

  deallog << "n_active_cells: " << tria.n_active_cells() << std::endl;

  return 0;
}
