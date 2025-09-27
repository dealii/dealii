// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check correct behavior of checking all Trilinos vector entries are
// non-negative works for distributed vectors

#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/trilinos_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"


void
test()
{
  const unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const unsigned int n_processes =
    Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  // create a vector that consists of elements indexed from 0 to n
  TrilinosWrappers::MPI::Vector vec;
  {
    IndexSet is(100 * n_processes);
    is.add_range(100 * myid, 100 * myid + 100);
    vec.reinit(is, MPI_COMM_WORLD);
  }
  AssertThrow(vec.locally_owned_size() == 100, ExcInternalError());
  AssertThrow(vec.local_range().first == 100 * myid, ExcInternalError());
  AssertThrow(vec.local_range().second == 100 * myid + 100, ExcInternalError());
  for (unsigned int i = vec.local_range().first; i < vec.local_range().second;
       ++i)
    vec(i) = i;
  vec.compress(VectorOperation::insert);

  // verify correctness so far
  {
    bool exact_non_negative = true;
    AssertThrow(vec.is_non_negative() == exact_non_negative,
                ExcInternalError());
  }

  if (vec.in_local_range(vec.size() / 2))
    vec[vec.size() / 2] = -1;
  {
    bool exact_non_negative = false;
    AssertThrow(vec.is_non_negative() == exact_non_negative,
                ExcInternalError());
  }

  deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      initlog();

      test();
    }
  else
    test();
}
