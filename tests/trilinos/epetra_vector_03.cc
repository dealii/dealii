// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/trilinos_epetra_vector.h>
#include <iostream>
#include <vector>

// Check LinearAlgebra::EpetraWrappers::Vector all_zero.

void
test()
{
  IndexSet     parallel_partitioner(10);
  unsigned int rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  if(rank == 0)
    parallel_partitioner.add_range(0, 5);
  else
    parallel_partitioner.add_range(5, 10);
  parallel_partitioner.compress();
  LinearAlgebra::EpetraWrappers::Vector a(parallel_partitioner, MPI_COMM_WORLD);

  AssertThrow(a.all_zero() == true, ExcInternalError());

  IndexSet                               read_write_index_set(10);
  LinearAlgebra::ReadWriteVector<double> read_write(parallel_partitioner);
  if(rank == 0)
    read_write[0] = 1.;
  a.import(read_write, VectorOperation::insert);
  AssertThrow(a.all_zero() == false, ExcInternalError());
}

int
main(int argc, char** argv)
{
  initlog();
  deallog.depth_console(0);

  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);

  test();

  deallog << "OK" << std::endl;

  return 0;
}
