/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2018 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------
 */

/**
 * An example of in-place smoothing of a distorted grid using Mesquite
 */


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <Mesquite_all_headers.hpp>

#include <fstream>

#include "../tests.h"

#include "mesquite_utilities.h"


int
main()
{
  initlog();

  const unsigned int dim                 = 2;
  const unsigned int n_elements_per_side = 4;

  // Make a distorted grid
  Triangulation<dim>        tria;
  std::vector<unsigned int> repetitions(2, n_elements_per_side);
  GridGenerator::subdivided_hyper_rectangle(tria,
                                            repetitions,
                                            Point<2>(0.0, 0.0),
                                            Point<2>(1.0, 1.0));
  GridTools::distort_random(0.3, tria, true);

  // Output the distorted grid
  deallog << "Distorted mesh" << std::endl;
  output_mesh(tria, "grid.01");

  // Smooth the mesh in-place using Mesquite
  smooth_mesquite_in_place_arraymesh(tria);

  // Output the smoothed grid
  deallog << "Smoothed mesh" << std::endl;
  output_mesh(tria, "grid.02");

  deallog << "OK" << std::endl;
}
