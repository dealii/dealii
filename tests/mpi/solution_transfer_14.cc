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

// Verify that Triangulation::save() and load() work with SolutionTransfer on
// a parallel::shared::Triangulation with artificial cells, with and without
// hp capabilities. METIS partitions the mesh, so the locally owned cells of a
// process are not contiguous, and every process only stores the values of its
// locally relevant degrees of freedom. A linear function is restored exactly.


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

#include "../tests.h"


template <int dim>
void
test(const bool hp)
{
  using VectorType = LinearAlgebra::distributed::Vector<double>;

  const std::string filename =
    "save_load_" + std::to_string(dim) + "d" + (hp ? "_hp" : "");

  const ScalarFunctionFromFunctionObject<dim> function([](const Point<dim> &p) {
    double value = 1.;
    for (unsigned int d = 0; d < dim; ++d)
      value += (d + 1.) * p[d];
    return value;
  });

  const hp::FECollection<dim> fe_collection(FE_Q<dim>(1), FE_Q<dim>(2));

  // The same finite elements before saving and after loading:
  const auto distribute_dofs = [&](DoFHandler<dim> &dof_handler) {
    if (hp)
      {
        for (const auto &cell : dof_handler.active_cell_iterators())
          if (cell->is_locally_owned())
            cell->set_active_fe_index(cell->center()[0] < 0.5 ? 0 : 1);
        dof_handler.distribute_dofs(fe_collection);
      }
    else
      dof_handler.distribute_dofs(FE_Q<dim>(2));
  };

  const auto create_vector = [](const DoFHandler<dim> &dof_handler) {
    return VectorType(dof_handler.locally_owned_dofs(),
                      DoFTools::extract_locally_relevant_dofs(dof_handler),
                      MPI_COMM_WORLD);
  };

  {
    parallel::shared::Triangulation<dim> tria(
      MPI_COMM_WORLD,
      Triangulation<dim>::limit_level_difference_at_vertices,
      true,
      parallel::shared::Triangulation<dim>::partition_metis);
    GridGenerator::hyper_cube(tria);
    tria.refine_global(3);
    for (const auto &cell : tria.active_cell_iterators())
      if (cell->is_locally_owned() && cell->center().norm() < 0.5)
        cell->set_refine_flag();
    tria.execute_coarsening_and_refinement();

    DoFHandler<dim> dof_handler(tria);
    distribute_dofs(dof_handler);

    VectorType solution = create_vector(dof_handler);
    VectorTools::interpolate(dof_handler, function, solution);
    solution.update_ghost_values();

    SolutionTransfer<dim, VectorType> solution_transfer(dof_handler);
    solution_transfer.prepare_for_serialization(solution);
    tria.save(filename);

    deallog << "cells: " << tria.n_global_active_cells() << std::endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);

  {
    parallel::shared::Triangulation<dim> tria(
      MPI_COMM_WORLD,
      Triangulation<dim>::limit_level_difference_at_vertices,
      true,
      parallel::shared::Triangulation<dim>::partition_metis);
    GridGenerator::hyper_cube(tria);
    tria.load(filename);

    DoFHandler<dim> dof_handler(tria);
    distribute_dofs(dof_handler);

    VectorType                        solution = create_vector(dof_handler);
    SolutionTransfer<dim, VectorType> solution_transfer(dof_handler);
    solution_transfer.deserialize(solution);

    VectorType reference = create_vector(dof_handler);
    VectorTools::interpolate(dof_handler, function, reference);
    for (const auto i : solution.locally_owned_elements())
      AssertThrow(std::abs(solution[i] - reference[i]) < 1e-12,
                  ExcInternalError());

    deallog << "OK" << std::endl;
  }
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
