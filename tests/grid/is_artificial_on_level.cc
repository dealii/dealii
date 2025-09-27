// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test CellAccessor::is_artificial_on_level() on a subdivided hypercube.

#include <deal.II/base/mpi.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>

#include "../tests.h"


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;


  parallel::distributed::Triangulation<2> tria(
    MPI_COMM_WORLD,
    dealii::Triangulation<2>::none,
    parallel::distributed::Triangulation<2>::construct_multigrid_hierarchy);
  GridGenerator::subdivided_hyper_cube(tria, 4);

  AssertDimension(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD), 2);

  const auto rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  for (const auto &cell : tria.active_cell_iterators())
    {
      if (rank == 0)
        {
          if (cell->center()[1] > 0.75)
            {
              Assert(cell->is_artificial_on_level(), ExcInternalError());
            }
          else
            {
              Assert(!cell->is_artificial_on_level(), ExcInternalError());
            }
        }
      else if (rank == 1)
        {
          if (cell->center()[1] < 0.25)
            {
              Assert(cell->is_artificial_on_level(), ExcInternalError());
            }
          else
            {
              Assert(!cell->is_artificial_on_level(), ExcInternalError());
            }
        }
      else
        {
          Assert(false, ExcInternalError());
        }
    }

  deallog << "OK!" << std::endl;
}
