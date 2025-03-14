/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2021 by the deal.II authors
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
 */

// Test (de)serialization of LinearAlgebra::distributed::Vectors

// clang-format off
/* This currently fails with

An error occurred in line <624> of file <../include/deal.II/base/partitioner.templates.h> in function
    void dealii::Utilities::MPI::Partitioner::import_from_ghosted_array_finish(dealii::VectorOperation::values, const dealii::ArrayView<const ElementType, MemorySpaceType>&, const dealii::ArrayView<ElementType, MemorySpace>&, const dealii::ArrayView<ElementType, MemorySpace>&, std::vector<int, std::allocator<int> >&) const [with Number = double; MemorySpaceType = dealii::MemorySpace::Host]
The violated condition was:
    *read_position == Number() || internal::get_abs(locally_owned_array[j] - *read_position) <= internal::get_abs(locally_owned_array[j] + *read_position) * 100000. * std::numeric_limits<typename numbers::NumberTraits< Number>::real_type>::epsilon()
Additional information:
    Called compress(VectorOperation::insert), but the element received
    from a remote processor, value 1, does not match with the value 0 on
    the owner processor 1
*/
// clang-format on

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>
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
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
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

  ConditionalOStream pcout;
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
  , pcout(std::cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
{}


template <int dim>
void
LaplaceProblem<dim>::setup_system()
{
  dof_handler.distribute_dofs(fe);
  locally_owned_dofs = dof_handler.locally_owned_dofs();
  DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
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
  for (unsigned int cycle = 0; cycle < n_cycles_adaptive; ++cycle)
    {
      pcout << "Cycle " << 1 + cycle << " / " << n_cycles_adaptive << ':'
            << std::endl;

      if (cycle == 0)
        {
          GridGenerator::subdivided_hyper_cube(triangulation, 10);
          triangulation.refine_global(n_cycles_global);
        }
      else
        refine_grid();

      pcout << "n_global_active_cells: "
            << triangulation.n_global_active_cells()
            << " n_global_levels: " << triangulation.n_global_levels()
            << " ghost_owners.size: " << triangulation.ghost_owners().size()
            << " level_ghost_owners.size: "
            << triangulation.level_ghost_owners().size() << std::endl;

      setup_system();

      const unsigned int n_vectors = 2;
      {
        pcout << "checkpointing..." << std::endl;
        std::vector<LinearAlgebra::distributed::Vector<double>> vectors(
          n_vectors);
        for (unsigned int i = 0; i < n_vectors; ++i)
          {
            vectors[i].reinit(locally_owned_dofs,
                              locally_relevant_dofs,
                              mpi_communicator);
            vectors[i].add((double)i);
            vectors[i].update_ghost_values();
          }

        std::vector<const LinearAlgebra::distributed::Vector<double> *>
            x_system(n_vectors);
        int i = 0;
        for (auto &v : x_system)
          {
            v = &vectors[i];
            ++i;
          }

        parallel::distributed::
          SolutionTransfer<dim, LinearAlgebra::distributed::Vector<double>>
            system_trans(dof_handler);

        system_trans.prepare_for_serialization(x_system);
        triangulation.save("restart.mesh");
      }

      {
        pcout << "resume..." << std::endl;
        std::vector<LinearAlgebra::distributed::Vector<double>> vectors(
          n_vectors);
        for (unsigned int i = 0; i < n_vectors; ++i)
          {
            vectors[i].reinit(locally_owned_dofs,
                              locally_relevant_dofs,
                              mpi_communicator);
            vectors[i].zero_out_ghost_values();
          }

        std::vector<LinearAlgebra::distributed::Vector<double> *> x_system(
          n_vectors);
        int i = 0;
        for (auto &v : x_system)
          {
            v = &vectors[i];
            ++i;
          }

        triangulation.coarsen_global(99);
        triangulation.load("restart.mesh");

        parallel::distributed::
          SolutionTransfer<dim, LinearAlgebra::distributed::Vector<double>>
            system_trans(dof_handler);
        system_trans.deserialize(x_system);


        for (unsigned int i = 0; i < n_vectors; ++i)
          {
            pcout << "vec " << i << " " << vectors[i].linfty_norm()
                  << std::endl;
            vectors[i].update_ghost_values();
          }

        // data out
        if (false)
          {
            pcout << "data out..." << std::endl;
            DataOut<dim> data_out;
            data_out.attach_dof_handler(dof_handler);
            for (unsigned int i = 0; i < n_vectors; ++i)
              data_out.add_data_vector(vectors[i],
                                       "u" + Utilities::to_string(i));

            data_out.build_patches();

            pcout << "writing..." << std::endl;

            data_out.write_vtu_with_pvtu_record(
              "./", "solution", cycle, mpi_communicator, 8);
            pcout << "done." << std::endl;
          }
      }

      pcout << std::endl;
    }
}


int
main(int argc, char *argv[])
{
  unsigned int n_cycles_global   = 2;
  unsigned int n_cycles_adaptive = 2;

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  LaplaceProblem<2> laplace_problem_2d;
  laplace_problem_2d.run(n_cycles_global, n_cycles_adaptive);
}
