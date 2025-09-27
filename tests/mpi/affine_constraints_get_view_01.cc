// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Check AffineConstraints<double>.get_view()
//
// The way this test works is by setting up a Stokes system with a
// Q2xQ1 element on an adaptively refined mesh for which we compute
// the corresponding constraints (which will contain constraints for
// both the velocity and pressure degrees of freedom). Then we create
// a random vector to which we apply distribute(). This is all
// existing functionality.
//
// Next, we create views into the constraints object for both the
// velocity variable only and the pressure variable, and apply these
// views to the individual blocks of the (original) random
// vector. This must result in the same output vector as above.

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



template <int dim>
void
test()
{
  deallog << "dim=" << dim << std::endl;

  const unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const unsigned int n_processes =
    Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  // Set up of triangulation and dof_handler

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  const double R0 = 6371000. - 2890000.;
  const double R1 = 6371000. - 35000.;
  GridGenerator::hyper_shell(tr, Point<dim>(), R0, R1, 12, true);

  tr.refine_global(1);
  for (unsigned int step = 0; step < 5; ++step)
    {
      typename Triangulation<dim>::active_cell_iterator cell =
                                                          tr.begin_active(),
                                                        endc = tr.end();

      for (; cell != endc; ++cell)
        if (Testing::rand() % 42 == 1)
          cell->set_refine_flag();

      tr.execute_coarsening_and_refinement();
    }
  deallog << "Cells on process #" << myid << ": " << tr.n_active_cells()
          << std::endl;

  FESystem<dim>   fe(FE_Q<dim>(2), dim, FE_Q<dim>(1), 1);
  DoFHandler<dim> dof_handler(tr);

  dof_handler.distribute_dofs(fe);

  std::vector<unsigned int> stokes_sub_blocks(dim + 1, 0);
  stokes_sub_blocks[dim] = 1;
  DoFRenumbering::component_wise(dof_handler, stokes_sub_blocks);

  deallog << "Total number of DoFs: " << dof_handler.n_dofs() << std::endl;

  const std::vector<types::global_dof_index> stokes_dofs_per_block =
    DoFTools::count_dofs_per_fe_block(dof_handler, stokes_sub_blocks);

  const types::global_dof_index n_u = stokes_dofs_per_block[0],
                                n_p = stokes_dofs_per_block[1];

  const IndexSet &stokes_locally_owned_index_set =
    dof_handler.locally_owned_dofs();
  const IndexSet stokes_locally_relevant_set =
    DoFTools::extract_locally_relevant_dofs(dof_handler);

  std::vector<IndexSet> stokes_partitioning;
  stokes_partitioning.push_back(
    stokes_locally_owned_index_set.get_view(0, n_u));
  stokes_partitioning.push_back(
    stokes_locally_owned_index_set.get_view(n_u, n_u + n_p));

  std::vector<IndexSet> stokes_relevant_partitioning;
  stokes_relevant_partitioning.push_back(
    stokes_locally_relevant_set.get_view(0, n_u));
  stokes_relevant_partitioning.push_back(
    stokes_locally_relevant_set.get_view(n_u, n_u + n_p));

  // Set up a constraints object from hanging node constraints and
  // boundary conditions.
  AffineConstraints<double> stokes_constraints(stokes_locally_owned_index_set,
                                               stokes_locally_relevant_set);

  DoFTools::make_hanging_node_constraints(dof_handler, stokes_constraints);

  const FEValuesExtractors::Vector velocity_components(0);
  VectorTools::interpolate_boundary_values(
    dof_handler,
    0,
    Functions::ZeroFunction<dim>(dim + 1),
    stokes_constraints,
    fe.component_mask(velocity_components));

  std::set<types::boundary_id> no_normal_flux_boundaries;
  no_normal_flux_boundaries.insert(1);
  VectorTools::compute_no_normal_flux_constraints(dof_handler,
                                                  0,
                                                  no_normal_flux_boundaries,
                                                  stokes_constraints);
  stokes_constraints.close();
  deallog << "Number of constrains: " << stokes_constraints.n_constraints()
          << std::endl;

  // Now set up a block vector for this situation. Fill it with random numbers.
  TrilinosWrappers::MPI::BlockVector random_vector;
  random_vector.reinit(stokes_partitioning, MPI_COMM_WORLD);
  for (const unsigned int i : stokes_locally_owned_index_set)
    random_vector(i) = random_value<double>(-100, 100);
  random_vector.compress(VectorOperation::insert);


  // Check 1: Call distribute() on a vector that contains
  // everything, with an AffineConstraint object that contains all
  // constraints.
  TrilinosWrappers::MPI::BlockVector check_vector_1;
  {
    check_vector_1.reinit(stokes_partitioning, MPI_COMM_WORLD);
    check_vector_1 = random_vector;

    stokes_constraints.distribute(check_vector_1);
  }

  // Check 2: Call distribute() on the two blocks individually, with
  // AffineConstraint objects that only contain views.
  TrilinosWrappers::MPI::BlockVector check_vector_2;
  {
    check_vector_2.reinit(stokes_partitioning, MPI_COMM_WORLD);
    check_vector_2 = random_vector;


    IndexSet velocity_dofs(dof_handler.n_dofs());
    velocity_dofs.add_range(0, n_u);
    velocity_dofs.compress();

    const AffineConstraints<double> velocity_constraints =
      stokes_constraints.get_view(velocity_dofs);


    IndexSet pressure_dofs(dof_handler.n_dofs());
    pressure_dofs.add_range(n_u, n_u + n_p);
    pressure_dofs.compress();

    const AffineConstraints<double> pressure_constraints =
      stokes_constraints.get_view(pressure_dofs);

    velocity_constraints.distribute(check_vector_2.block(0));
    pressure_constraints.distribute(check_vector_2.block(1));
  }

  // Finally check that these two vectors are the same. Note that
  // when we compile with native optimizations we might have a slight
  // difference in results:
  {
    TrilinosWrappers::MPI::BlockVector result = check_vector_1;
    result -= check_vector_2;
    Assert(result.l2_norm() / check_vector_1.l2_norm() < 1e-8,
           ExcInternalError());
  }

  deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  mpi_initlog();

  test<2>();
  test<3>();
}
