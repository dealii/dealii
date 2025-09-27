// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/utilities.h>

#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/trilinos_epetra_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"

// Check LinearAlgebra::EpetraWrappers::Vector all_zero.

void
test()
{
  IndexSet     parallel_partitioner(10);
  unsigned int rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  if (rank == 0)
    parallel_partitioner.add_range(0, 5);
  else
    parallel_partitioner.add_range(5, 10);
  parallel_partitioner.compress();
  LinearAlgebra::EpetraWrappers::Vector a(parallel_partitioner, MPI_COMM_WORLD);

  AssertThrow(a.all_zero() == true, ExcInternalError());

  IndexSet                               read_write_index_set(10);
  LinearAlgebra::ReadWriteVector<double> read_write(parallel_partitioner);
  if (rank == 0)
    read_write[0] = 1.;
  a.import_elements(read_write, VectorOperation::insert);
  AssertThrow(a.all_zero() == false, ExcInternalError());
}


int
main(int argc, char **argv)
{
  initlog();
  deallog.depth_console(0);

  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);

  test();

  deallog << "OK" << std::endl;

  return 0;
}
