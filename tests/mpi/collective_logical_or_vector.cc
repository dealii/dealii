// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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



// check Utilities::MPI::logical_or() for vectors
//
// std::vector<bool> does not necessarily store elements in contiguous array,
// use std::vector<char> instead

#include <deal.II/base/mpi.templates.h>
#include <deal.II/base/utilities.h>

#include "../tests.h"

void
test()
{
  const unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  std::vector<char> values = {0, (myid % 2) == 0};
  std::vector<char> results(2);
  Utilities::MPI::logical_or(values, MPI_COMM_WORLD, results);
  Assert(results[0] == false, ExcInternalError());
  Assert(results[1] == true, ExcInternalError());

  if (myid == 0)
    deallog << std::boolalpha << static_cast<bool>(results[0]) << ' '
            << static_cast<bool>(results[1]) << std::endl;
}


int
main(int argc, char *argv[])
{
#ifdef DEAL_II_WITH_MPI
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());
#else
  (void)argc;
  (void)argv;
  compile_time_error;

#endif

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      initlog();

      deallog.push("mpi");
      test();
      deallog.pop();
    }
  else
    test();
}
