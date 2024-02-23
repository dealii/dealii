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


// Test Utilities::MPI::create_mpi_data_type_n_bytes

// To not require a lot of memory, we run this test with smaller
// size. Change this bool to test the real 64bit size communication:
const bool run_big = false;

#include <deal.II/base/mpi.h>

#include "../tests.h"


void
test_data_type(const std::uint64_t n_bytes)
{
  const auto bigtype = Utilities::MPI::create_mpi_data_type_n_bytes(n_bytes);

  deallog << "checking size " << n_bytes << ':';

  int size32;
  int ierr = MPI_Type_size(*bigtype, &size32);
  AssertThrowMPI(ierr);

  if (size32 == MPI_UNDEFINED)
    deallog << " size32="
            << "UNDEFINED (too big)";
  else
    deallog << " size32=" << size32;

  MPI_Count size64;
  ierr = MPI_Type_size_x(*bigtype, &size64);
  AssertThrowMPI(ierr);

  deallog << " size64=" << size64;


  deallog << std::endl;
}



void
test_send_recv(MPI_Comm comm)
{
  unsigned int        myid    = Utilities::MPI::this_mpi_process(comm);
  const std::uint64_t n_bytes = (run_big) ? ((1ULL << 31) + 37) : (37ULL);

  if (myid == 0)
    {
      std::vector<char> buffer(n_bytes, 'A');
      buffer[n_bytes - 1] = 'B';
      const auto bigtype =
        Utilities::MPI::create_mpi_data_type_n_bytes(buffer.size());
      int ierr =
        MPI_Send(buffer.data(), 1, *bigtype, 1 /* dest */, 0 /* tag */, comm);
      AssertThrowMPI(ierr);
    }
  else if (myid == 1)
    {
      std::vector<char> buffer(n_bytes, '?');
      const auto        bigtype =
        Utilities::MPI::create_mpi_data_type_n_bytes(buffer.size());
      int ierr = MPI_Recv(buffer.data(),
                          1,
                          *bigtype,
                          0 /* src */,
                          0 /* tag */,
                          comm,
                          MPI_STATUS_IGNORE);
      AssertThrowMPI(ierr);

      AssertThrow(buffer[0] == 'A', ExcInternalError());
      AssertThrow(buffer[n_bytes - 1] == 'B', ExcInternalError());
    }
  deallog << "send_recv: OK" << std::endl;
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  test_data_type(0);
  test_data_type(1);
  test_data_type(1ULL << 30);
  test_data_type(1ULL << 31);
  test_data_type(1ULL << 32);
  test_data_type(1ULL << 33);
  test_data_type(1ULL << 55);
  test_data_type(58493729485ULL);

  test_send_recv(MPI_COMM_WORLD);
}
