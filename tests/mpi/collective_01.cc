//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


// check Utilities::System::calculate_collective_mpi_min_max_avg()

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <fstream>

void print_it(Utilities::System::MinMaxAvg & result)
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
  Utilities::System::calculate_collective_mpi_min_max_avg(MPI_COMM_WORLD,
							  value,
							  result);
  if (myid==0)
    print_it(result);  
  Assert(result.sum == 0.0, ExcInternalError());
  Assert(result.min == 0.0, ExcInternalError());
  Assert(result.max == 0.0, ExcInternalError());
  Assert(result.avg == 0.0, ExcInternalError());

  value = 1.0;
  Utilities::System::calculate_collective_mpi_min_max_avg(MPI_COMM_WORLD,
							  value,
							  result);
  if (myid==0)
    print_it(result);
  Assert(result.sum == numprocs, ExcInternalError());
  Assert(result.min == 1.0, ExcInternalError());
  Assert(result.max == 1.0, ExcInternalError());
  Assert(result.avg == 1.0, ExcInternalError());

  value = 0.0;
  if (myid==0)
    value = 1.0;
  
  Utilities::System::calculate_collective_mpi_min_max_avg(MPI_COMM_WORLD,
							  value,
							  result);
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
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
  Utilities::System::MPI_InitFinalize mpi (argc, argv);
#else
  (void)argc;
  (void)argv;
  compile_time_error;

#endif

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    {
      std::ofstream logfile(output_file_for_mpi("collective_01").c_str());
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
