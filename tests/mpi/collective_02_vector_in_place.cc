// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2015 by the deal.II authors
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



// check Utilities::MPI::sum() for vectors, but with input=output

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <fstream>

void test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  const unsigned int numprocs = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  unsigned int values_[2] = { 1, 2 };
  std::vector<unsigned int> sums(&values_[0], &values_[2]);
  Utilities::MPI::sum (sums,
                       MPI_COMM_WORLD,
                       sums);
  Assert (sums[0] == numprocs, ExcInternalError());
  Assert (sums[1] == 2*numprocs, ExcInternalError());

  if (myid==0)
    deallog << sums[0] << ' ' << sums[1] << std::endl;
}


int main(int argc, char *argv[])
{
#ifdef DEAL_II_WITH_MPI
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, numbers::invalid_unsigned_int);
#else
  (void)argc;
  (void)argv;
  compile_time_error;

#endif

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    {
      std::ofstream logfile("output");
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      deallog.push("mpi");
      test();
      deallog.pop();
    }
  else
    test();
}
