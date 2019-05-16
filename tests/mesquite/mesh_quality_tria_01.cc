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
 * A 2d example of measuring the quality of a mesh after smoothing a distorted
 * grid using the Mesquite interface. If the grid is output, then the results
 * can be compared to those computed in Paraview.
 *
 * This test is based off of mesh_interface_tria_01.cc
 */


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_mesquite.h>
#include <deal.II/grid/tria.h>

#include <fstream>

#include "../tests.h"
#include "mesquite_utilities.h"


int
main()
{
  initlog();

  const unsigned int dim                 = 2;
  const unsigned int n_elements_per_side = 4;

  // Make a distored grid
  Triangulation<dim>        tria;
  std::vector<unsigned int> repetitions(dim, n_elements_per_side);
  GridGenerator::subdivided_hyper_rectangle(tria,
                                            repetitions,
                                            Point<2>(0.0, 0.0),
                                            Point<2>(1.0, 1.0));
  GridTools::distort_random(0.3, tria, true);

  // Output the distorted grid
  output_mesh(tria, "grid.01");

  // Initialise mesh smoother
  GridTools::MesquiteMeshInterface<dim> mesquite_interface;
  mesquite_interface.initialize(tria, true /*fix_all_boundary_vertices*/);
  deallog << "Distorted mesh: Minimum cell quality: "
          << mesquite_interface.get_minimum_cell_quality() << std::endl;

  // Perform mesh smoothing
  mesquite_interface.execute(GridTools::MesquiteWrapperTypes::laplace);
  deallog << "Smoothed mesh: Minimum cell quality: "
          << mesquite_interface.get_minimum_cell_quality() << std::endl;

  // Output smoothed mesh
  output_mesh(tria, "grid.02");

  deallog << "OK" << std::endl;
}
