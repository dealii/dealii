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
 * This test is the a combination of mesh_interface_dofhandler_03b and
 * mesh_interface_dofhandler_06b, and tests the use of a BlockVector instead
 * of a normal Vector. There's also an extra dummy block just to mix things up.
 *
 * - A regular, structured mesh is given to the triangulation
 * - A DoFHandler that represents the mesh motion is constructed
 *   - A distorted deformation map is defined
 *   - An extra dummy field is introduced for the FESystem
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
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_mesquite.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/block_vector.h>
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
  const unsigned int  degree_u         = 2;
  const unsigned int  n_dummy_fields_0 = 2;
  const unsigned int  n_dummy_fields_1 = 1;
  const unsigned int  n_u_fields       = dim;
  const unsigned int  n_dummy_fields_2 = 1;
  const FESystem<dim> fe(FE_Q<dim>(1),
                         n_dummy_fields_0, // A dummy field
                         FE_DGQ<dim>(3),
                         n_dummy_fields_1, // Another dummy field
                         FE_Q<dim>(degree_u),
                         n_u_fields, // *** The displacement field ***
                         FE_DGQ<dim>(2),
                         n_dummy_fields_2); // Another dummy field


  // Setup DoF components
  const unsigned int first_dummy_0_component = 0;
  const unsigned int first_dummy_1_component =
    first_dummy_0_component + n_dummy_fields_0;
  const unsigned int first_u_component =
    first_dummy_1_component + n_dummy_fields_1;
  const unsigned int first_dummy_2_component = first_u_component + n_u_fields;
  const unsigned int n_dof_components =
    first_dummy_2_component + n_dummy_fields_2;
  const unsigned int n_dof_components_total = n_dof_components;

  // Setup block components
  enum
  {
    dummy_0_block = 0,
    dummy_1_block,
    u_block,
    dummy_3_block,
    n_blocks
  };
  std::vector<unsigned int> block_component(n_dof_components);
  for (unsigned int i = 0; i < n_dof_components; ++i)
    {
      if (i < first_dummy_1_component)
        block_component[i] = dummy_0_block;
      else if (i < first_u_component)
        block_component[i] = dummy_1_block;
      else if (i < first_dummy_2_component)
        block_component[i] = n_u_fields;
      else
        {
          AssertThrow(i < n_dof_components,
                      ExcIndexRange(i,
                                    first_dummy_2_component,
                                    n_dof_components));
          block_component[i] = dummy_3_block;
        }
    }


  // Distribute some mesh motion DoFs
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);
  DoFRenumbering::Cuthill_McKee(dof_handler);
  DoFRenumbering::component_wise(dof_handler, block_component);
  //  Assert(n_dof_components == dof_handler.get_fe_collection().n_components(),
  //         ExcDimensionMismatch(n_dof_components,dof_handler.get_fe_collection().n_components()));

  std::vector<types::global_dof_index> dofs_per_block(n_blocks);
  DoFTools::count_dofs_per_block(dof_handler, dofs_per_block, block_component);

  // The Eulerian vertex position field
  BlockVector<double> solution_x(dofs_per_block);
  output_mesh(dof_handler, solution_x, u_block, first_u_component, "grid.01");

  // Prescribe the mesh motion
  // First create the completely distorted Eulerian solution
  // but fix the boundary mesh to the "undistorted" solution.
  ConstraintMatrix                 constraints_x;
  const FEValuesExtractors::Vector displacement(first_u_component);
  const ComponentMask component_mask_u = fe.component_mask(displacement);
  VectorTools::interpolate_boundary_values(
    dof_handler,
    0,
    TransformationSinusoid<dim>(
      false, 0.0, n_dof_components_total, first_u_component),
    constraints_x,
    component_mask_u);
  constraints_x.close();
  VectorTools::project(dof_handler,
                       constraints_x,
                       QGauss<dim>(degree_u + 2),
                       TransformationSinusoid<dim>(
                         true, 0.75, n_dof_components_total, first_u_component),
                       solution_x);

  // Next we get the reference mesh vertex postions
  BlockVector<double> solution_X(dofs_per_block);
  VectorTools::interpolate(dof_handler,
                           GridTools::internal::ReferencePosition<dim>(
                             n_dof_components_total, first_u_component),
                           solution_X,
                           component_mask_u);

  // Next we compute the mesh motion solution
  BlockVector<double> solution_u(solution_x);
  solution_u -= solution_X;

  // Output the distorted grid
  std::cout << "Distorted mesh" << std::endl;
  output_mesh(dof_handler, solution_u, u_block, first_u_component, "grid.02");

  // Perform mesh smoothing
  GridTools::MesquiteMeshInterface<dim> mesquite_interface(
    dof_handler,
    solution_u,
    first_u_component,
    true /*fix_all_boundary_vertices*/);
  mesquite_interface.execute(GridTools::MesquiteWrapperTypes::untangler);
  mesquite_interface.execute(
    GridTools::MesquiteWrapperTypes::shape_improvement);

  // Adjust the solution field that describes the mesh motion
  mesquite_interface.update_mesh_displacement_solution(solution_u,
                                                       dof_handler,
                                                       first_u_component);

  // Output the smoothed grid
  std::cout << "Smoothed mesh (DoFHandler)" << std::endl;
  output_mesh(dof_handler,
              solution_u,
              u_block,
              first_u_component,
              "grid.03",
              true /*write_to_deallog*/);

  // Move the triangulation's vertices
  mesquite_interface.move_triangulation_vertices(tria);

  // Output the smoothed grid
  std::cout << "Smoothed mesh (Moved triangulation)" << std::endl;
  output_mesh(tria, "grid.04");

  std::cout << "OK" << std::endl;
}
