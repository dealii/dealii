// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check LinearAlgebra::TpetraWrappers::Vector::import
// for VectorOperation::add and check LinearAlgebra::ReadWriteVector
// for VectorOperation::add/min/max.

#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/trilinos_tpetra_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"

template <typename Number>
void
test()
{
  unsigned int my_id   = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int n_procs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  IndexSet locally_owned(n_procs * 2);
  locally_owned.add_range(my_id * 2, my_id * 2 + 2);
  locally_owned.compress();

  LinearAlgebra::TpetraWrappers::Vector<Number, MemorySpace::Default> v(
    locally_owned, MPI_COMM_WORLD);

  IndexSet     workaround_set(n_procs * 2);
  unsigned int my_first_index  = (my_id * 2 + 2) % (n_procs * 2);
  unsigned int my_second_index = my_id * 2 + 1;
  workaround_set.add_index(my_first_index);
  workaround_set.add_index(my_second_index);
  workaround_set.add_index(0);
  workaround_set.compress();
  LinearAlgebra::ReadWriteVector<Number> rw_vector(workaround_set);

  rw_vector(my_first_index)  = my_id + 10;
  rw_vector(my_second_index) = my_id + 100;
  rw_vector(0)               = 1.;
  // rw_vector(2) = 1.;
  v.import_elements(rw_vector, VectorOperation::add);
  deallog << "Tpetra first import add:" << std::endl;
  v.print(deallog.get_file_stream());
  rw_vector.print(deallog.get_file_stream());

  rw_vector(my_first_index)  = my_id + 20;
  rw_vector(my_second_index) = my_id + 200;
  rw_vector(0)               = 2.;
  // rw_vector(2) = 3.;
  v.import_elements(rw_vector, VectorOperation::add);
  deallog << "Tpetra second import add:" << std::endl;
  v.print(deallog.get_file_stream());
  rw_vector.print(deallog.get_file_stream());

  rw_vector.import_elements(v, VectorOperation::add);
  deallog << "ReadWrite import add:" << std::endl;
  rw_vector.print(deallog.get_file_stream());

  rw_vector(my_first_index)  = my_id + 100;
  rw_vector(my_second_index) = 1;
  rw_vector(0)               = 4.;
  // rw_vector(2) = 3.;
  rw_vector.import_elements(v, VectorOperation::min);
  deallog << "ReadWrite import min:" << std::endl;
  rw_vector.print(deallog.get_file_stream());

  rw_vector(my_first_index)  = my_id + 100;
  rw_vector(my_second_index) = 1;
  rw_vector(0)               = 4.;
  // rw_vector(2) = 3.;
  rw_vector.import_elements(v, VectorOperation::max);
  deallog << "ReadWrite import max:" << std::endl;
  rw_vector.print(deallog.get_file_stream());
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  MPILogInitAll log;
  test<double>();
}
