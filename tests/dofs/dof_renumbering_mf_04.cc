// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check DoFRenumbering::matrix_free_data_locality on a hypercube mesh in
// parallel on the multigrid levels

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

template <int dim>
void
test(const unsigned int degree)
{
  deallog << "Test in " << dim << "D with degree " << degree << std::endl;
  parallel::distributed::Triangulation<dim> tria(
    MPI_COMM_WORLD,
    Triangulation<dim>::limit_level_difference_at_vertices,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);
  GridGenerator::hyper_cube(tria, -1., 1.);
  tria.refine_global(5 - dim);
  FE_Q<dim>       fe(degree);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  dof.distribute_mg_dofs();

  using MatrixFreeType = MatrixFree<dim, float, VectorizedArray<float, 1>>;
  typename MatrixFreeType::AdditionalData mf_data;
  mf_data.tasks_parallel_scheme = MatrixFreeType::AdditionalData::none;

  AffineConstraints<double> constraints;
  for (unsigned int level = 0; level < tria.n_global_levels(); ++level)
    {
      mf_data.mg_level = level;

      {
        const auto renumber =
          DoFRenumbering::compute_matrix_free_data_locality(dof,
                                                            constraints,
                                                            mf_data);

        deallog << "Level " << level << std::endl;
        deallog << "Renumbering no constraints: " << std::endl;
        for (unsigned int i = 0; i < renumber.size(); ++i)
          {
            deallog << renumber[i] << " ";
            if (i % 16 == 15)
              deallog << std::endl;
          }
        deallog << std::endl;
      }

      DoFRenumbering::matrix_free_data_locality(dof, constraints, mf_data);
      std::vector<types::global_dof_index> dof_indices(fe.dofs_per_cell);
      deallog << "New dof indices on cells: " << std::endl;
      for (const auto &cell : dof.mg_cell_iterators_on_level(level))
        if (cell->is_locally_owned_on_level())
          {
            cell->get_active_or_mg_dof_indices(dof_indices);
            for (const auto i : dof_indices)
              deallog << i << " ";
            deallog << std::endl;
          }
      deallog << std::endl;
    }
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  MPILogInitAll                    log;

  test<2>(2);
  test<2>(4);
  test<3>(2);
  test<3>(4);
}
