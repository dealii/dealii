// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/mpi.h>

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/grid/grid_generator.h>

#include <iostream>

#include "../tests.h"


// Test to check whether Zoltan PHG throws warning because
// PHG_EDGE_SIZE_THRESHOLD value is low.

template <int dim>
void
test(const MPI_Comm mpi_communicator)
{
  parallel::shared::Triangulation<dim> triangulation(
    mpi_communicator, Triangulation<dim>::limit_level_difference_at_vertices);

  GridGenerator::subdivided_hyper_cube(triangulation, 2, 0, 1);

  // Some dummy output.
  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    std::cout << "Number of vertices: " << triangulation.n_vertices()
              << std::endl;
}

int
main(int argc, char *argv[])
{
  try
    {
      dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(
        argc, argv, dealii::numbers::invalid_unsigned_int);

      test<1>(MPI_COMM_WORLD);
      test<2>(MPI_COMM_WORLD);
      test<3>(MPI_COMM_WORLD);
    }
  catch (const std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      throw;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      throw;
    }

  return 0;
}
