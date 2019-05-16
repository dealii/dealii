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
 * This test is the same as mesh_interface_dofhandler_07, except that it uses
 * a parallel::shared::Triangulation and has some local refinement too.
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


#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_mesquite.h>

#include <deal.II/lac/trilinos_parallel_block_vector.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"
#include "mesquite_utilities.h"
//#include <deal.II/lac/trilinos_block_vector.h> // Note: cannot have this and
// t_p_b_v.h in same compilation unit!

#include <algorithm>
#include <fstream>


int
main(int argc, char **argv)
{
  typedef dealii::TrilinosWrappers::MPI::BlockVector BlockVector;

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  const MPI_Comm &   mpi_communicator = MPI_COMM_WORLD;
  const unsigned int this_mpi_process =
    Utilities::MPI::this_mpi_process(mpi_communicator);
  const unsigned int n_mpi_processes =
    Utilities::MPI::n_mpi_processes(mpi_communicator);
  ConditionalOStream pcout(std::cout, this_mpi_process == 0);

  const unsigned int dim                 = 2;
  const unsigned int n_refinement_cycles = 2;

  // Make a standard grid
  parallel::shared::Triangulation<dim> tria(
    mpi_communicator, Triangulation<dim>::maximum_smoothing);
  std::vector<unsigned int> repetitions(dim);
  repetitions[0] = 14;
  repetitions[1] = 4;
  GridGenerator::subdivided_hyper_rectangle(tria,
                                            repetitions,
                                            Point<2>(0.0, 0.0),
                                            Point<2>(14.0, 2.0));

  const double threshold = 0.1;
  for (unsigned int c = 0; c < n_refinement_cycles; ++c)
    {
      for (auto cell : tria.active_cell_iterators_on_level(c))
        {
          if (random_value() > (1.0 - threshold))
            cell->set_refine_flag();
        }
      tria.execute_coarsening_and_refinement();
    }

  pcout << "Reference mesh" << std::endl;
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


  // Count DoFs in each block
  std::vector<types::global_dof_index> dofs_per_block(n_blocks);
  DoFTools::count_dofs_per_block(dof_handler, dofs_per_block, block_component);

  // Compute locally owned degrees-of-freedom for all processes
  std::vector<IndexSet> all_locally_owned_dofs =
    DoFTools::locally_owned_dofs_per_subdomain(dof_handler);
  std::vector<IndexSet> all_locally_relevant_dofs =
    DoFTools::locally_relevant_dofs_per_subdomain(dof_handler);

  // Add entries if not all processors are assigned DoFs.
  if (all_locally_owned_dofs.size() < n_mpi_processes)
    all_locally_owned_dofs.resize(n_mpi_processes, IndexSet());
  if (all_locally_relevant_dofs.size() < n_mpi_processes)
    all_locally_relevant_dofs.resize(n_mpi_processes, IndexSet());

  // Extract locally owned degrees-of-freedom
  IndexSet              locally_owned_dofs;
  std::vector<IndexSet> locally_owned_partitioning;
  Assert(all_locally_owned_dofs.size() > this_mpi_process, ExcInternalError());
  locally_owned_dofs = all_locally_owned_dofs[this_mpi_process];

  // Extract locally relevant degrees-of-freedom
  IndexSet              locally_relevant_dofs;
  std::vector<IndexSet> locally_relevant_partitioning;
  Assert(all_locally_relevant_dofs.size() > this_mpi_process,
         ExcInternalError());
  locally_relevant_dofs = all_locally_relevant_dofs[this_mpi_process];

  // Compute the block partitioning
  locally_owned_partitioning.reserve(n_blocks);
  locally_relevant_partitioning.reserve(n_blocks);
  for (unsigned int b = 0; b < n_blocks; ++b)
    {
      const types::global_dof_index idx_begin =
        std::accumulate(dofs_per_block.begin(),
                        std::next(dofs_per_block.begin(), b),
                        0);
      const types::global_dof_index idx_end =
        std::accumulate(dofs_per_block.begin(),
                        std::next(dofs_per_block.begin(), b + 1),
                        0);
      locally_owned_partitioning.push_back(
        locally_owned_dofs.get_view(idx_begin, idx_end));
      locally_relevant_partitioning.push_back(
        locally_relevant_dofs.get_view(idx_begin, idx_end));
    }

  // The mesh motion field
  BlockVector solution_u(locally_owned_partitioning, mpi_communicator);
  output_mesh(dof_handler,
              solution_u,
              locally_relevant_partitioning,
              u_block,
              first_u_component,
              "grid.01",
              mpi_communicator);

  // Prescribe the mesh motion
  // Since VectorTools::project is not impemented for parallel::Triangulation,
  // we will have to be creative an do this another way...
  // Note: Because a random number generator is used to do the perturbation, the
  //       resulting deformation field will look different to of
  //       mesh_interface_dofhandler_07
  {
    const Vector<double> eulerian_vertex_map_exact =
      GridTools::get_eulerian_vertex_positions(
        tria,
        TransformationSinusoid<dim>(
          false, 0.0, n_dof_components_total, first_u_component),
        first_u_component);
    const Vector<double> eulerian_vertex_map_ptrb =
      GridTools::get_eulerian_vertex_positions(
        tria,
        TransformationSinusoid<dim>(
          true, 0.75, n_dof_components_total, first_u_component),
        first_u_component);

    // Undo the perturbations of the boundary vertices
    Vector<double> vertex_displacement_map(dim * tria.n_vertices());
    fix_eulerian_boundary_vertex_map(vertex_displacement_map,
                                     tria,
                                     eulerian_vertex_map_exact,
                                     eulerian_vertex_map_ptrb);

    // Substract the original vertex positions to get the final displacement
    // vector
    vertex_displacement_map -= GridTools::get_eulerian_vertex_positions(
      tria,
      GridTools::internal::ReferencePosition<dim>(n_dof_components_total,
                                                  first_u_component),
      first_u_component);

    GridTools::internal::MeshMotionData<dim> mm_data(dof_handler);
    mm_data.insert_vertex_displacement_map(solution_u,
                                           dof_handler,
                                           vertex_displacement_map,
                                           first_u_component);
  }

  // Output the distorted grid
  pcout << "Distorted mesh" << std::endl;
  output_mesh(dof_handler,
              solution_u,
              locally_relevant_partitioning,
              u_block,
              first_u_component,
              "grid.02",
              mpi_communicator);

  // Perform mesh smoothing
  GridTools::MesquiteMeshInterface<dim> mesquite_interface(
    dof_handler,
    solution_u,
    first_u_component,
    true /*fix_all_boundary_vertices*/);
  mesquite_interface.execute(GridTools::MesquiteWrapperTypes::untangler);
  mesquite_interface.execute(GridTools::MesquiteWrapperTypes::laplace);
  const unsigned int n_cycles = 3;
  for (unsigned int c = 0; c < n_cycles;
       ++c) // Distortion is so severe that it needs multiple cycles!
    mesquite_interface.execute(
      GridTools::MesquiteWrapperTypes::shape_improvement);

  // Adjust the solution field that describes the mesh motion
  mesquite_interface.update_mesh_displacement_solution(solution_u,
                                                       dof_handler,
                                                       first_u_component);

  // Output the smoothed grid
  pcout << "Smoothed mesh (DoFHandler)" << std::endl;
  output_mesh(dof_handler,
              solution_u,
              locally_relevant_partitioning,
              u_block,
              first_u_component,
              "grid.03",
              mpi_communicator,
              true /*write_to_deallog*/,
              OutputFormat::vtk);
  //  output_mesh (dof_handler, solution_u, locally_relevant_partitioning,
  //               u_block, first_u_component, "grid.03", mpi_communicator,
  //               false /*write_to_deallog*/, OutputFormat::vtu);

  // Move the triangulation's vertices
  mesquite_interface.move_triangulation_vertices(tria);

  // Output the smoothed grid
  pcout << "Smoothed mesh (Moved triangulation)" << std::endl;
  output_mesh(tria,
              "grid.04." + Utilities::int_to_string(this_mpi_process),
              true /*write_to_deallog*/);

  pcout << "OK" << std::endl;
}
