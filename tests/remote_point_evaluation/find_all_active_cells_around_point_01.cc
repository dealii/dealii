// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test GridTools::find_active_cell_around_point() and
// GridTools::find_all_active_cells_around_point for simplices. These
// functions are used in find_all_locally_owned_active_cells_around_point(),
// which is used by distributed_compute_point_locations().

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

template <int dim>
void
test()
{
  const auto  reference_cell = ReferenceCells::get_simplex<dim>();
  const auto &mapping =
    reference_cell.template get_default_linear_mapping<dim>();

  Triangulation<dim> tria;
  if (false)
    GridGenerator::reference_cell(tria, reference_cell);
  else
    GridGenerator::subdivided_hyper_cube_with_simplices(tria, 2);

  GridTools::Cache<dim> cache(tria, mapping);

  const unsigned n_subdivisions = 8;

  for (unsigned int i = 0; i <= n_subdivisions; ++i)
    for (unsigned int j = 0; j <= n_subdivisions; ++j)
      {
        Point<dim> p(1.0 * i / n_subdivisions, 1.0 * j / n_subdivisions);

        const auto first_point =
          GridTools::find_active_cell_around_point(mapping, tria, p, {}, 1e-6);

        const auto result = GridTools::find_all_active_cells_around_point(
          mapping, tria, p, 1e-6, {first_point.first, first_point.second});

        deallog << result.size() << std::endl; // should be 2
        for (const auto &cell : result)
          deallog << cell.first->id() << " " << cell.second << std::endl;
        deallog << std::endl;
      }
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  deallog.precision(8);

  test<2>();
}
