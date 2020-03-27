// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2018 by the deal.II authors
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

#include <deal.II/base/mpi.h>

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/grid/grid_generator.h>

#include <iostream>

#include "../tests.h"


// Test to check whether Zoltan PHG throws warning because
// PHG_EDGE_SIZE_THRESHOLD value is low.

template <int dim>
void
test(const MPI_Comm &mpi_communicator)
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
  catch (std::exception &exc)
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
