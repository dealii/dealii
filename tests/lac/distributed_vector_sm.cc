// ---------------------------------------------------------------------
//
// Copyright (C) 2021 - 2022 by the deal.II authors
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

// Test LinearAlgebra::distributed::Vector's shared-memory capability on a
// Cartesian virtual topology.

#include <deal.II/base/mpi.h>

#include <deal.II/lac/la_parallel_vector.h>

#include "../tests.h"

using namespace dealii;



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  MPILogInitAll                    all;

  AssertDimension(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD), 4);

  const auto my_rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  MPI_Comm sm_comm;
  MPI_Comm_split_type(
    MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, my_rank, MPI_INFO_NULL, &sm_comm);
  AssertDimension(Utilities::MPI::n_mpi_processes(sm_comm), 4);

  MPI_Comm row_comm;
  MPI_Comm_split(MPI_COMM_WORLD, my_rank / 2, my_rank, &row_comm);

  const auto my_row_rank = Utilities::MPI::this_mpi_process(row_comm);

  MPI_Comm column_comm;
  MPI_Comm_split(MPI_COMM_WORLD, my_rank % 2, my_rank, &column_comm);

  IndexSet is_local(20);
  IndexSet is_ghost(20);

  if (my_row_rank == 0)
    {
      is_local.add_range(0, 10);
      is_ghost.add_range(10, 15);
    }
  else
    {
      is_local.add_range(10, 20);
      is_ghost.add_range(5, 10);
    }

  const auto test = [&](const auto sm_comm) {
    const auto partitioner =
      std::make_shared<Utilities::MPI::Partitioner>(is_local,
                                                    is_ghost,
                                                    row_comm);

    LinearAlgebra::distributed::Vector<double> vector;
    vector.reinit(partitioner, sm_comm);

    for (unsigned int i = 0; i < partitioner->locally_owned_size(); ++i)
      vector.local_element(i) = 10 * my_rank + i;

    vector.update_ghost_values();

    const auto local_size =
      partitioner->locally_owned_size() + partitioner->n_ghost_indices();
    const auto sizes = Utilities::MPI::all_gather(sm_comm, local_size);

    for (unsigned int i = 0; i < sizes.size(); ++i)
      {
        for (unsigned int j = 0; j < sizes[i]; ++j)
          deallog << vector.shared_vector_data()[i][j] << " ";

        deallog << std::endl;
      }
    deallog << std::endl;
  };

  test(MPI_COMM_WORLD);
  test(row_comm);
  test(column_comm);

  MPI_Comm_free(&column_comm);
  MPI_Comm_free(&row_comm);
  MPI_Comm_free(&sm_comm);
}
