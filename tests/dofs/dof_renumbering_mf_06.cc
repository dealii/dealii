// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
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
// parallel for FE_SimplexP instead of FE_Q

#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_simplex_p.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

template <int dim>
void
test(const unsigned int degree)
{
  deallog << "Test in " << dim << "D with degree " << degree << std::endl;

  parallel::fullydistributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  {
    const auto serial_grid_generator =
      [](dealii::Triangulation<dim, dim> &tria_serial) {
        dealii::Triangulation<dim, dim> temp;
        GridGenerator::subdivided_hyper_cube(temp, 2);
        GridGenerator::convert_hypercube_to_simplex_mesh(temp, tria_serial);
        tria_serial.refine_global(4 - dim);
      };
    const auto serial_grid_partitioner =
      [&](dealii::Triangulation<dim, dim> &tria_serial,
          const MPI_Comm                   comm,
          const unsigned int) {
        dealii::GridTools::partition_triangulation_zorder(
          dealii::Utilities::MPI::n_mpi_processes(comm), tria_serial);
      };

    const unsigned int group_size = 8;

    typename dealii::TriangulationDescription::Settings
      triangulation_description_setting =
        dealii::TriangulationDescription::default_setting;
    const auto description = dealii::TriangulationDescription::Utilities::
      create_description_from_triangulation_in_groups<dim, dim>(
        serial_grid_generator,
        serial_grid_partitioner,
        tria.get_communicator(),
        group_size,
        dealii::Triangulation<dim>::none,
        triangulation_description_setting);

    tria.create_triangulation(description);
  }

  FE_SimplexP<dim> fe(degree);
  DoFHandler<dim>  dof(tria);
  dof.distribute_dofs(fe);

  using MatrixFreeType = MatrixFree<dim, double, VectorizedArray<double, 1>>;
  typename MatrixFreeType::AdditionalData mf_data;
  mf_data.tasks_parallel_scheme = MatrixFreeType::AdditionalData::none;

  AffineConstraints<double> constraints;

  {
    const auto renumber =
      DoFRenumbering::compute_matrix_free_data_locality(dof,
                                                        constraints,
                                                        mf_data);

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
  for (const auto &cell : dof.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        cell->get_dof_indices(dof_indices);
        for (const auto i : dof_indices)
          deallog << i << " ";
        deallog << std::endl;
      }
  deallog << std::endl;

  dof.distribute_dofs(fe);
  VectorTools::interpolate_boundary_values(dof,
                                           0,
                                           Functions::ZeroFunction<dim>(),
                                           constraints);
  constraints.close();

  {
    const auto renumber =
      DoFRenumbering::compute_matrix_free_data_locality(dof,
                                                        constraints,
                                                        mf_data);

    deallog << "Renumbering Dirichlet constraints: " << std::endl;
    for (unsigned int i = 0; i < renumber.size(); ++i)
      {
        deallog << renumber[i] << " ";
        if (i % 16 == 15)
          deallog << std::endl;
      }
    deallog << std::endl;
  }

  DoFRenumbering::matrix_free_data_locality(dof, constraints, mf_data);
  deallog << "New dof indices on cells: " << std::endl;
  for (const auto &cell : dof.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        cell->get_dof_indices(dof_indices);
        for (const auto i : dof_indices)
          deallog << i << " ";
        deallog << std::endl;
      }
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  MPILogInitAll                    log;

  test<2>(2);
  test<3>(1);
  test<3>(3);
}
