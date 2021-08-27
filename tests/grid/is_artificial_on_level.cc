// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// Test CellAccessor::is_artificial_on_level() on a subdivided hypercube.

#include <deal.II/base/mpi.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>

#include "../tests.h"

using namespace dealii;

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
