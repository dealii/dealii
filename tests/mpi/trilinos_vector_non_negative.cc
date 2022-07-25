// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2021 by the deal.II authors
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
  AssertThrow(vec.local_size() == 100, ExcInternalError());
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
