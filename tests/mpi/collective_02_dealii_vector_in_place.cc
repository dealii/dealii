// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2018 by the deal.II authors
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



// check Utilities::MPI::sum() for dealii::Vector, but with input=output

#include <deal.II/base/utilities.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"

void
test()
{
  unsigned int       myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const unsigned int numprocs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  {
    Vector<int> sums(2);
    sums[0] = 1;
    sums[1] = 2;
    Utilities::MPI::sum(sums, MPI_COMM_WORLD, sums);
    Assert((unsigned int)sums[0] == numprocs, ExcInternalError());
    Assert((unsigned int)sums[1] == 2 * numprocs, ExcInternalError());

    if (myid == 0)
      deallog << sums[0] << ' ' << sums[1] << std::endl;
  }

  {
    Vector<double> sums(2);
    sums[0] = 1.5;
    sums[1] = 2.5;
    Utilities::MPI::sum(sums, MPI_COMM_WORLD, sums);
    Assert(sums[0] == 1.5 * numprocs, ExcInternalError());
    Assert(sums[1] == 2.5 * numprocs, ExcInternalError());

    if (myid == 0)
      deallog << sums[0] << ' ' << sums[1] << std::endl;
  }
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
