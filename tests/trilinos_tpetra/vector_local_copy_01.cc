// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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



// Test that the creation of a local copy of a vector does not lead
// to a deadlock if not performed on all MPI processes

#include <deal.II/lac/trilinos_tpetra_vector.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"


int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());
  MPI_Comm           mpi_communicator = MPI_COMM_WORLD;
  const unsigned int this_mpi_process =
    Utilities::MPI::this_mpi_process(mpi_communicator);
  const unsigned int n_mpi_processes =
    Utilities::MPI::n_mpi_processes(mpi_communicator);

  MPI_Comm           serial_communicator;
  const unsigned int colour = this_mpi_process;
  const unsigned int key    = this_mpi_process;
  if (n_mpi_processes > 1)
    MPI_Comm_split(mpi_communicator, colour, key, &serial_communicator);
  else
    serial_communicator = mpi_communicator;

  if (this_mpi_process == 0)
    {
      const LinearAlgebra::TpetraWrappers::Vector<double> tril_vec(
        complete_index_set(10), serial_communicator);

      // Check copy constructor
      const Vector<double> local_vector_1(tril_vec);

      // Check equality operator
      Vector<double> local_vector_2;
      local_vector_2 = tril_vec;
    }

  deallog << "OK" << std::endl;
}
