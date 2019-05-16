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
 * A 2d example of smoothing of a distorted grid using the Mesquite interface.
 * The grid is locally refined, and so it has hanging nodes.
 * This test is similar to mesh_interface_tria_06, except that now one
 * of the smoother wrappers are used.
 *
 * - The vertices of the triangulation are moved to a distorted position
 * - Mesh adaption is performed
 *     - A predefined smoother is used
 *     - The boundary vertices are considered fixed
 *     - The global mesh optmizer is used
 *     - We rely on the automatic enforcement of hanging vertex constraints at
 *       the end of MesquiteMeshInterface::execute()
 * - The vertices of the triangulation are updated to their optimial locations
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
  const unsigned int n_refinement_cycles = 3;

  // Make a distored grid
  Triangulation<dim>        tria(Triangulation<dim>::maximum_smoothing);
  std::vector<unsigned int> repetitions(dim, n_elements_per_side);
  GridGenerator::subdivided_hyper_rectangle(tria,
                                            repetitions,
                                            Point<2>(0.0, 0.0),
                                            Point<2>(1.0, 1.0));
  GridTools::distort_random(0.3, tria, true);
  for (unsigned int c = 0; c < n_refinement_cycles; ++c)
    {
      tria.begin_active(c)->set_refine_flag();
      tria.execute_coarsening_and_refinement();
    }

  // Output the distorted grid
  std::cout << "Distorted mesh" << std::endl;
  output_mesh(tria, "grid.01");
  output_mesh_vertex_numbering(tria, "grid.01");

  // Perform mesh smoothing
  GridTools::MesquiteMeshInterface<dim> mesquite_interface(
    tria, true /*fix_all_boundary_vertices*/);
  // Note: We choose a Laplace smoother specifically because it doesn't produce
  //       the optimial result (unlike MesquiteWrapperTypes::shape_improvement,
  //       which works very well). This allows us to better see that hanging
  //       vertex constraints are indeed enforced correctly as the resulting
  //       grid is not Cartesian aligned.
  GridTools::MesquiteMeshInterface<dim>::Settings settings;
  settings.vertex_level_algorithm = GridTools::MesquiteMeshInterface<
    dim>::Settings::VertexConsideration::global;
  mesquite_interface.execute(GridTools::MesquiteWrapperTypes::laplace,
                             settings);

  // Move the triangulation's vertices
  mesquite_interface.move_triangulation_vertices(tria);

  // Output the smoothed grid
  std::cout << "Smoothed mesh" << std::endl;
  output_mesh(tria, "grid.02", true /*write_to_deallog*/);

  std::cout << "OK" << std::endl;
}
