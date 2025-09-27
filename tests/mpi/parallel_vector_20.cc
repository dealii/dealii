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

// Check import function

#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/read_write_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"

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

  v.import_elements(read_write_vector, VectorOperation::insert);

  AssertThrow(v.local_element(0) == 1., ExcInternalError());
  AssertThrow(v.local_element(1) == 2., ExcInternalError());

  read_write_owned.clear();
  read_write_owned.add_index(1);
  read_write_owned.add_index(2);
  read_write_vector.reinit(read_write_owned);
  read_write_vector.local_element(0) = 1.;
  read_write_vector.local_element(1) = 2.;

  v.import_elements(read_write_vector, VectorOperation::insert);

  AssertThrow(v.local_element(0) == my_id + 1, ExcInternalError());
  AssertThrow(v.local_element(1) == my_id + 1, ExcInternalError());

  if (my_id == 0)
    deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      initlog();
      deallog << std::setprecision(4);

      test();
    }
  else
    test();
}
