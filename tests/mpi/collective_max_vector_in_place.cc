// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2017 by the deal.II authors
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

// check Utilities::MPI::max() for vectors, but with input=output

#include "../tests.h"
#include <deal.II/base/utilities.h>

void
test()
{
  unsigned int       myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const unsigned int numprocs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  unsigned int              values_[2] = {myid, numprocs + myid};
  std::vector<unsigned int> maxima(&values_[0], &values_[2]);
  Utilities::MPI::max(maxima, MPI_COMM_WORLD, maxima);
  Assert(maxima[0] == numprocs - 1, ExcInternalError());
  Assert(maxima[1] == 2 * numprocs - 1, ExcInternalError());

  if(myid == 0)
    deallog << maxima[0] << ' ' << maxima[1] << std::endl;
}

int
main(int argc, char* argv[])
{
#ifdef DEAL_II_WITH_MPI
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());
#else
  (void) argc;
  (void) argv;
  compile_time_error;

#endif

  if(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      initlog();

      deallog.push("mpi");
      test();
      deallog.pop();
    }
  else
    test();
}
