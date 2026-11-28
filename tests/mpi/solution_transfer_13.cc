// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// Like solution_transfer_12, but with hp capabilities: cells use FE_Q(1) or
// FE_Q(2) and change their finite element at random. SolutionTransfer then
// sends variable size data to the new owners of the cells.


#include <deal.II/base/function.h>

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/hp/fe_collection.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

#include <random>

#include "../tests.h"


template <int dim>
void
test(const bool average_values)
{
  parallel::shared::Triangulation<dim> tria(
    MPI_COMM_WORLD,
    Triangulation<dim>::limit_level_difference_at_vertices,
    true,
    parallel::shared::Triangulation<dim>::partition_metis);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(dim == 2 ? 3 : 2);

  // The same sequence of random numbers on all processes:
  std::mt19937 random(42);

  const hp::FECollection<dim> fe(FE_Q<dim>(1), FE_Q<dim>(2));
  DoFHandler<dim>             dof_handler(tria);
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      const auto r = random() % 2;
      if (cell->is_locally_owned())
        cell->set_active_fe_index(r);
    }
  dof_handler.distribute_dofs(fe);

  const ScalarFunctionFromFunctionObject<dim> function([](const Point<dim> &p) {
    double value = 1.;
    for (unsigned int d = 0; d < dim; ++d)
      value += (d + 1.) * p[d];
    return value;
  });

  using VectorType = LinearAlgebra::distributed::Vector<double>;
  VectorType solution(dof_handler.locally_owned_dofs(),
                      DoFTools::extract_locally_relevant_dofs(dof_handler),
                      MPI_COMM_WORLD);
  VectorTools::interpolate(dof_handler, function, solution);

  for (unsigned int cycle = 0; cycle < 4; ++cycle)
    {
      for (const auto &cell : dof_handler.active_cell_iterators())
        {
          const auto r = random() % 5;
          if (!cell->is_locally_owned())
            continue;
          if (r == 0)
            cell->set_refine_flag();
          else if (r == 1 && cell->level() > 1)
            cell->set_coarsen_flag();
          else if (r == 2)
            cell->set_future_fe_index(1 - cell->active_fe_index());
        }

      solution.update_ghost_values();
      SolutionTransfer<dim, VectorType> solution_transfer(dof_handler,
                                                          average_values);
      solution_transfer.prepare_for_coarsening_and_refinement(solution);
      tria.execute_coarsening_and_refinement();

      dof_handler.distribute_dofs(fe);
      solution.reinit(dof_handler.locally_owned_dofs(),
                      DoFTools::extract_locally_relevant_dofs(dof_handler),
                      MPI_COMM_WORLD);
      solution_transfer.interpolate(solution);

      VectorType reference(solution);
      VectorTools::interpolate(dof_handler, function, reference);
      for (const auto i : solution.locally_owned_elements())
        AssertThrow(std::abs(solution[i] - reference[i]) < 1e-12,
                    ExcInternalError());
    }

  deallog << "OK" << std::endl;
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  deallog.push("2d");
  test<2>(false);
  test<2>(true);
  deallog.pop();
  deallog.push("3d");
  test<3>(false);
  test<3>(true);
  deallog.pop();
}
