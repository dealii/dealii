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



// check Utilities::MPI::min_max_avg()

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <fstream>

void print_it(Utilities::MPI::MinMaxAvg &result)
{
  deallog << "sum: " << result.sum
          << " avg: " << result.avg
          << " min: " << result.min << " @" << result.min_index
          << " max: " << result.max << " @" << result.max_index
          << std::endl;
}

void test()
{
  Assert( Utilities::System::job_supports_mpi(), ExcInternalError());

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  const unsigned int numprocs = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  if (myid==0)
    deallog << "Running on " << numprocs << " CPU(s)." << std::endl;

  Utilities::System::MinMaxAvg result;

  double value = 0.0;
  result = Utilities::MPI::min_max_avg(value, MPI_COMM_WORLD);
  if (myid==0)
    print_it(result);
  Assert(result.sum == 0.0, ExcInternalError());
  Assert(result.min == 0.0, ExcInternalError());
  Assert(result.max == 0.0, ExcInternalError());
  Assert(result.avg == 0.0, ExcInternalError());

  value = 1.0;
  result = Utilities::MPI::min_max_avg(value, MPI_COMM_WORLD);
  if (myid==0)
    print_it(result);
  Assert(result.sum == numprocs, ExcInternalError());
  Assert(result.min == 1.0, ExcInternalError());
  Assert(result.max == 1.0, ExcInternalError());
  Assert(result.avg == 1.0, ExcInternalError());

  value = 0.0;
  if (myid==0)
    value = 1.0;

  result = Utilities::MPI::min_max_avg(value, MPI_COMM_WORLD);
  if (myid==0)
    print_it(result);
  Assert(result.sum == 1.0, ExcInternalError());
  Assert(result.min == (numprocs>1)?0.0:1.0, ExcInternalError());
  Assert(result.max == 1.0, ExcInternalError());
  Assert(result.max_index == 0, ExcInternalError());

  if (myid==0)
    deallog << "done" << std::endl;
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
