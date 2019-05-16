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
 * The grid is in its coarsest state, with no hanging nodes.
 *
 * - The vertices of the triangulation are moved to a distorted position
 * - Mesh adaption is performed
 *     - A predefined smoother is used
 *     - The boundary vertices are considered fixed
 * - A vertex displacement map is constructed
 * - A view of the triangulation in the current configuration is output
 *   without moving the triangulation vertices
 */


#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_mesquite.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <fstream>

#include "../tests.h"
#include "mesquite_utilities.h"


int
main()
{
  initlog();

  const unsigned int dim = 2;

  // Make a standard grid
  Triangulation<dim>        tria;
  std::vector<unsigned int> repetitions(dim);
  repetitions[0] = 14;
  repetitions[1] = 4;
  GridGenerator::subdivided_hyper_rectangle(tria,
                                            repetitions,
                                            Point<2>(0.0, 0.0),
                                            Point<2>(14.0, 2.0));

  std::cout << "Reference mesh" << std::endl;
  output_mesh(tria, "grid.01");

  // Prescribe the mesh motion
  const unsigned int   first_u_dof = 0;
  const Vector<double> eulerian_vertex_map_exact =
    GridTools::get_eulerian_vertex_positions(tria,
                                             TransformationSinusoid<dim>(false),
                                             first_u_dof);
  const Vector<double> eulerian_vertex_map_ptrb =
    GridTools::get_eulerian_vertex_positions(tria,
                                             TransformationSinusoid<dim>(true),
                                             first_u_dof);

  // Undo the perturbations of the boundary vertices
  Vector<double> eulerian_vertex_map(dim * tria.n_vertices());
  fix_eulerian_boundary_vertex_map(eulerian_vertex_map,
                                   tria,
                                   eulerian_vertex_map_exact,
                                   eulerian_vertex_map_ptrb);

  // Output the exact correct grid
  std::cout << "Exact mesh" << std::endl;
  output_mesh(tria, eulerian_vertex_map_exact, "grid.02");

  // Output the distorted grid
  std::cout << "Distorted mesh" << std::endl;
  output_mesh(tria, eulerian_vertex_map, "grid.03");

  // Perform mesh smoothing
  GridTools::MesquiteMeshInterface<dim> mesquite_interface(
    tria, eulerian_vertex_map, true /*fix_all_boundary_vertices*/);
  mesquite_interface.execute(GridTools::MesquiteWrapperTypes::untangler);
  mesquite_interface.execute(
    GridTools::MesquiteWrapperTypes::shape_improvement);

  // Create the solution field that describes the mesh motion
  Vector<double> vertex_displacement_map(dim * tria.n_vertices());
  mesquite_interface.update_vertex_displacement_map(vertex_displacement_map,
                                                    tria);

  // Output the smoothed grid
  std::cout << "Smoothed mesh" << std::endl;
  output_mesh(tria,
              vertex_displacement_map,
              "grid.04",
              true /*write_to_deallog*/);

  std::cout << "OK" << std::endl;
}
