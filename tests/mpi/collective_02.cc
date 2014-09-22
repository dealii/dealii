// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
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



// check Utilities::MPI::sum()

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <fstream>

void test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  const unsigned int numprocs = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  int int_sum;
  unsigned int uint_sum;
  double double_sum;
  float float_sum;

  int_sum
    =
      Utilities::MPI::sum<int>(myid+1,
                               MPI_COMM_WORLD);
  uint_sum
    =
      Utilities::MPI::sum<unsigned int>(myid+1,
                                        MPI_COMM_WORLD);
  float_sum
    =
      Utilities::MPI::sum<float>(myid+1,
                                 MPI_COMM_WORLD);
  double_sum
    =
      Utilities::MPI::sum<double>(myid+1,
                                  MPI_COMM_WORLD);

  if (myid==0)
    deallog << int_sum << ' ' << uint_sum << ' ' << double_sum << ' ' << float_sum << std::endl;
}


int main(int argc, char *argv[])
{
#ifdef DEAL_II_WITH_MPI
  Utilities::MPI::MPI_InitFinalize mpi (argc, argv);
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
