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

// test: RWV::import_elements() from TrilinosWrappers::MPI::Vector

#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/trilinos_vector.h>

#include <vector>

#include "../tests.h"

void
test()
{
  IndexSet     is(8);
  unsigned int rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  if (rank == 0)
    is.add_range(0, 4);
  else if (rank == 1)
    is.add_range(4, 8);
  else
    AssertThrow(false, ExcNotImplemented());

  is.compress();

  deallog << "IS: ";
  is.print(deallog);

  TrilinosWrappers::MPI::Vector          tril_vector(is);
  LinearAlgebra::ReadWriteVector<double> readwrite(is);

  for (auto idx : is)
    readwrite[idx] = idx;

  deallog << "RWVector contents:" << std::endl;
  readwrite.print(deallog.get_file_stream());

  // import RWV->Trilinos
  tril_vector.import_elements(readwrite, VectorOperation::insert);
  deallog << "trilinos vec:" << std::endl;
  tril_vector.print(deallog.get_file_stream());

  // test that ::add also works
  tril_vector.import_elements(readwrite, VectorOperation::add);
  deallog << "trilinos vec (2x):" << std::endl;
  tril_vector.print(deallog.get_file_stream());

  // import again overwriting the contents
  tril_vector.import_elements(readwrite, VectorOperation::insert);
  deallog << "trilinos vec (1x):" << std::endl;
  tril_vector.print(deallog.get_file_stream());

  deallog << "OK" << std::endl;
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  MPILogInitAll log;

  test();
}
