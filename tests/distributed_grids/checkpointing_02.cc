/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2021 - 2023 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 */

// Test (de)serialization of Trilinos vectors with checkpointing files >4GB. The
// test here is of course simplified to run quickly.

// Set this to true to run a test that generates an 8GB file. Run with
// 5 MPI ranks to check correctness with a total file size above 4GB.
// Run with 2 MPI ranks to test a single process above 2GB. Both cases were
// broken before this test was made. Warning, you probably need in the
// order of 32GB of RAM to run this.
const bool run_big = false;

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

#include "../tests.h"


template <int dim>
class LaplaceProblem
{
public:
  LaplaceProblem();

  void
  run(unsigned int n_cycles_global, unsigned int n_cycles_adaptive);

private:
  void
  setup_system();
  void
  refine_grid();
  void
  output_results(const unsigned int cycle) const;

  MPI_Comm mpi_communicator;

  parallel::distributed::Triangulation<dim> triangulation;

  FE_Q<dim>       fe;
  DoFHandler<dim> dof_handler;

  IndexSet locally_owned_dofs;
  IndexSet locally_relevant_dofs;
};

template <int dim>
LaplaceProblem<dim>::LaplaceProblem()
  : mpi_communicator(MPI_COMM_WORLD)
  , triangulation(
      mpi_communicator,
      typename Triangulation<dim>::MeshSmoothing(
        Triangulation<dim>::smoothing_on_refinement |
        Triangulation<dim>::smoothing_on_coarsening),
      parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy)
  , fe(2)
  , dof_handler(triangulation)
{}


template <int dim>
void
LaplaceProblem<dim>::setup_system()
{
  dof_handler.distribute_dofs(fe);
  locally_owned_dofs    = dof_handler.locally_owned_dofs();
  locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);
}

template <int dim>
void
LaplaceProblem<dim>::refine_grid()
{
  // refine into a corner
  Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
  for (auto &cell : triangulation.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        estimated_error_per_cell(cell->active_cell_index()) =
          cell->center().norm() * std::pow(cell->diameter(), 0.5);
    }

  parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number(
    triangulation, estimated_error_per_cell, 0.20, 0.0);

  triangulation.execute_coarsening_and_refinement();
}



template <int dim>
void
LaplaceProblem<dim>::run(unsigned int n_cycles_global,
                         unsigned int n_cycles_adaptive)
{
  using VectorType = TrilinosWrappers::MPI::Vector;

  for (unsigned int cycle = 0; cycle < n_cycles_adaptive; ++cycle)
    {
      deallog << "Cycle " << 1 + cycle << " / " << n_cycles_adaptive << ':'
              << std::endl;

      if (cycle == 0)
        {
          GridGenerator::subdivided_hyper_cube(triangulation, 10);
          triangulation.refine_global(n_cycles_global);
        }
      else
        refine_grid();

      deallog << "n_global_active_cells: "
              << triangulation.n_global_active_cells()
              << " n_global_levels: " << triangulation.n_global_levels()
              << " ghost_owners.size: " << triangulation.ghost_owners().size()
              << " level_ghost_owners.size: "
              << triangulation.level_ghost_owners().size() << std::endl;

      setup_system();

      const unsigned int n_vectors = (run_big) ? 70 : 2;
      {
        deallog << "checkpointing..." << std::endl;
        std::vector<VectorType> vectors(n_vectors);
        VectorType              x(locally_owned_dofs, mpi_communicator);
        for (unsigned int i = 0; i < n_vectors; ++i)
          {
            vectors[i].reinit(locally_owned_dofs,
                              locally_relevant_dofs,
                              mpi_communicator);
            vectors[i] = x;
            x.add(1.0);
          }

        std::vector<const VectorType *> x_system(n_vectors);
        int                             i = 0;
        for (auto &v : x_system)
          {
            v = &vectors[i];
            ++i;
          }

        VectorType y(locally_owned_dofs,
                     locally_relevant_dofs,
                     mpi_communicator);
        y = x;

        // To be sure, use two SolutionTransfer objects, because the second one
        // will get a large offset
        SolutionTransfer<dim, VectorType> system_trans(dof_handler);
        SolutionTransfer<dim, VectorType> trans2(dof_handler);

        system_trans.prepare_for_serialization(x_system);
        trans2.prepare_for_serialization(y);
        triangulation.save("restart.mesh");
      }

      {
        deallog << "resume..." << std::endl;
        std::vector<VectorType> vectors(n_vectors);
        for (unsigned int i = 0; i < n_vectors; ++i)
          vectors[i].reinit(locally_owned_dofs, mpi_communicator);

        std::vector<VectorType *> x_system(n_vectors);
        int                       i = 0;
        for (auto &v : x_system)
          {
            v = &vectors[i];
            ++i;
          }

        VectorType y(locally_owned_dofs, mpi_communicator);

        triangulation.coarsen_global(99);
        triangulation.load("restart.mesh");

        SolutionTransfer<dim, VectorType> system_trans(dof_handler);
        SolutionTransfer<dim, VectorType> trans2(dof_handler);
        system_trans.deserialize(x_system);
        trans2.deserialize(y);

        for (unsigned int i = 0; i < n_vectors; ++i)
          deallog << "vec " << i << ": " << vectors[i].linfty_norm()
                  << std::endl;
        deallog << "vec y: " << y.linfty_norm() << std::endl;
      }

      deallog << std::endl;
    }
}


int
main(int argc, char *argv[])
{
  unsigned int n_cycles_global   = (run_big) ? 3 : 1;
  unsigned int n_cycles_adaptive = 1;

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  LaplaceProblem<3> laplace_problem;
  laplace_problem.run(n_cycles_global, n_cycles_adaptive);
}
