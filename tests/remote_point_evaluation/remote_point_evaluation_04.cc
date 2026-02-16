// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test Utilities::MPI::RemotePointEvaluation with n_components>1.

#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi_remote_point_evaluation.h>

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>

#include "../tests.h"



template <int dim, int n_components>
void
test(const MPI_Comm comm)
{
  deallog << "* n_components = " << n_components << std::endl;
  const unsigned int spacedim      = dim;
  const unsigned int n_refinements = 1;

  parallel::shared::Triangulation<dim, spacedim> tria(
    comm,
    Triangulation<dim, spacedim>::none,
    true,
    parallel::shared::Triangulation<dim, spacedim>::Settings::partition_zorder);

  GridGenerator::hyper_cube(tria);

  tria.refine_global(n_refinements);

  MappingQ1<spacedim> mapping;

  Utilities::MPI::RemotePointEvaluation<spacedim, spacedim> eval;

  std::vector<Point<spacedim>> points;
  std::vector<unsigned int>    data;
  if (Utilities::MPI::this_mpi_process(comm) == 0)
    {
      points.resize(2);
      points[1](0) = 0.75;

      data.resize(n_components * points.size());
      std::iota(data.begin(), data.end(), 1);
    }

  eval.reinit(points, tria, mapping);

  const auto fu = [&](const auto &values, const auto &cell_data) {
    for (const auto cell_index : cell_data.cell_indices())
      {
        const auto cell        = cell_data.get_active_cell_iterator(cell_index);
        const auto unit_points = cell_data.get_unit_points(cell_index);
        const ArrayView<const unsigned int> data =
          cell_data.get_data_view(cell_index, values, n_components);

        deallog << "cell index = " << cell_index
                << ", center = " << cell->center() << ":\n";
        for (unsigned int q = 0; q < unit_points.size(); ++q)
          {
            deallog << "\tpoint " << q << ", coordinates = " << unit_points[q]
                    << ", data:\n\t";
            for (unsigned int c = 0; c < n_components; ++c)
              deallog << data[n_components * q + c] << ' ';
            deallog << '\n';
          }
        deallog << std::endl;
      }
  };

  eval.template process_and_evaluate<unsigned int, n_components>(data, fu);

  deallog << "OK." << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  MPILogInitAll                    all;

  const int n_ranks = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  const int my_rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  test<2, 1>(MPI_COMM_WORLD);
  test<2, 2>(MPI_COMM_WORLD);
}
