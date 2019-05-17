/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2003 - 2018 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------
 *
 * Authors: Fabian Castelli, KIT
 *          Timo Heister, University of Utah
 */

// @note: This is work in progress and will be an example for a non-linear
// problem solved in parallel with matrix-free geometric multigrid. For now,
// this is just step-1.

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <iostream>
#include <fstream>
#include <cmath>

using namespace dealii;

void first_grid()
{
  Triangulation<2> triangulation;

  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(4);

  std::ofstream out("grid-1.eps");
  GridOut       grid_out;
  grid_out.write_eps(triangulation, out);
  std::cout << "Grid written to grid-1.eps" << std::endl;
}


int main()
{
  first_grid();
}
