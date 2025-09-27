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



// test that coarsening/refinement flags are communicated correctly

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/grid/grid_generator.h>

#include "../tests.h"


template <int dim>
void
test()
{
  using Number              = double;
  using VectorizedArrayType = VectorizedArray<Number>;

  parallel::shared::Triangulation<dim> tria(
    MPI_COMM_WORLD,
    ::Triangulation<dim>::none,
    true,
    parallel::shared::Triangulation<dim>::partition_custom_signal);

  tria.signals.create.connect([&]() {
    for (const auto &cell : tria.active_cell_iterators())
      if (cell->center()[1] < 0.5)
        cell->set_subdomain_id(0);
      else
        cell->set_subdomain_id(1);
  });

  const auto print = [&]() {
    for (unsigned int l = 0; l < tria.n_levels(); ++l)
      deallog << tria.n_cells(l) << " ";
    deallog << std::endl;
  };

  GridGenerator::hyper_cube(tria);

  print();

  tria.refine_global();

  print();

  for (const auto &cell : tria.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        const auto center = cell->center();

        if (center[0] < 0.5 && center[1] < 0.5)
          cell->set_refine_flag();
        else if (center[0] > 0.5 && center[1] > 0.5)
          cell->set_coarsen_flag();
      }

  tria.execute_coarsening_and_refinement();

  print();
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll log;

  AssertDimension(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD), 2);

  test<2>();
}
