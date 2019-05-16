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
 * This test is the same as mesh_interface_dofhandler_01, except that now
 * higher order polynomials are chosen for the displacement field.
 *
 * - A regular, structured mesh is given to the triangulation
 * - A DoFHandler that represents the mesh motion is constructed
 *   - A distorted deformation map is defined
 * - Mesh adaption is performed
 *   - The view of the underlying mesh is constructed from the DoFHandler and
 *     the prescribed deformation map
 *   - A predefined smoother is used
 *   - The boundary vertices are considered fixed
 * - The vertices of the triangulation are updated to their optimal locations
 *    - The view of the underlying mesh is constructed from the reference
 *      mesh vertex positions and the optimized form of the input deformation
 * map
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
main(int argc, char **argv)
{
  initlog();
  // Prevent threading so that VTK output is deterministic
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

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

  // A FE that describes the mesh motion field
  const unsigned int  degree = 2;
  const FESystem<dim> fe(FE_Q<dim>(degree), dim);
  const unsigned int  first_u_component = 0;

  // Distribute some mesh motion DoFs
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  // The Eulerian vertex position field
  Vector<double> solution_x(dof_handler.n_dofs());
  output_mesh(dof_handler, solution_x, first_u_component, "grid.01");

  // Prescribe the mesh motion
  // First create the completely distorted Eulerian solution
  // but fix the boundary mesh to the "undistorted" solution.
  ConstraintMatrix                 constraints_x;
  const FEValuesExtractors::Vector displacement(first_u_component);
  const ComponentMask component_mask_u = fe.component_mask(displacement);
  VectorTools::interpolate_boundary_values(dof_handler,
                                           0,
                                           TransformationSinusoid<dim>(false),
                                           constraints_x,
                                           component_mask_u);
  constraints_x.close();
  VectorTools::project(dof_handler,
                       constraints_x,
                       QGauss<dim>(degree + 2),
                       TransformationSinusoid<dim>(true, 0.75),
                       solution_x);

  // Next we get the reference mesh vertex postions
  Vector<double> solution_X(dof_handler.n_dofs());
  VectorTools::interpolate(dof_handler,
                           GridTools::internal::ReferencePosition<dim>(),
                           solution_X,
                           component_mask_u);

  // Next we compute the mesh motion solution
  Vector<double> solution_u(solution_x);
  solution_u -= solution_X;

  // Output the distorted grid
  std::cout << "Distorted mesh" << std::endl;
  output_mesh(dof_handler, solution_u, first_u_component, "grid.02");

  // Perform mesh smoothing
  GridTools::MesquiteMeshInterface<dim> mesquite_interface(
    dof_handler,
    solution_u,
    first_u_component,
    true /*fix_all_boundary_vertices*/);
  mesquite_interface.execute(GridTools::MesquiteWrapperTypes::untangler);
  mesquite_interface.execute(
    GridTools::MesquiteWrapperTypes::shape_improvement);

  // Move the triangulation's vertices
  mesquite_interface.move_triangulation_vertices(tria);

  // Output the smoothed grid
  std::cout << "Smoothed mesh" << std::endl;
  output_mesh(tria, "grid.03", true /*write_to_deallog*/);

  std::cout << "OK" << std::endl;
}
