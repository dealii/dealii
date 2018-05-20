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

// Check import function

#include "../tests.h"
#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/read_write_vector.h>
#include <iostream>
#include <vector>

void
test()
{
  unsigned int my_id   = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int n_procs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  IndexSet locally_owned(n_procs * 2);
  locally_owned.add_range(my_id * 2, my_id * 2 + 2);
  IndexSet read_write_owned(4);
  read_write_owned.add_range(my_id * 2, my_id * 2 + 2);

  LinearAlgebra::distributed::Vector<double> v(locally_owned, MPI_COMM_WORLD);
  LinearAlgebra::ReadWriteVector<double> read_write_vector(read_write_owned);
  read_write_vector.local_element(0) = 1.;
  read_write_vector.local_element(1) = 2.;

  v.import(read_write_vector, VectorOperation::insert);

  AssertThrow(v.local_element(0) == 1., ExcInternalError());
  AssertThrow(v.local_element(1) == 2., ExcInternalError());

  read_write_owned.clear();
  read_write_owned.add_index(1);
  read_write_owned.add_index(2);
  read_write_vector.reinit(read_write_owned);
  read_write_vector.local_element(0) = 1.;
  read_write_vector.local_element(1) = 2.;

  v.import(read_write_vector, VectorOperation::insert);

  AssertThrow(v.local_element(0) == my_id + 1, ExcInternalError());
  AssertThrow(v.local_element(1) == my_id + 1, ExcInternalError());

  if(my_id == 0)
    deallog << "OK" << std::endl;
}

int
main(int argc, char** argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  deallog.push(Utilities::int_to_string(myid));

  if(myid == 0)
    {
      initlog();
      deallog << std::setprecision(4);

      test();
    }
  else
    test();
}
