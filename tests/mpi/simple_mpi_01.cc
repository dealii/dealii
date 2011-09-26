//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


// check if mpi is working

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>

#include <fstream>
//#include <mpi.h>


void test_mpi()
{
  Assert( Utilities::System::job_supports_mpi(), ExcInternalError());


  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  const unsigned int numprocs = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  if (myid==0)
    deallog << "Running on " << numprocs << " CPU(s)." << std::endl;

  for (unsigned int i=1;i<numprocs;++i)
    {
      MPI_Barrier(MPI_COMM_WORLD);
//      system("sleep 1");

      if (myid==0)
	{
	  unsigned int buf=numbers::invalid_unsigned_int;
	  MPI_Status status;
	  MPI_Recv(&buf, 1, MPI_UNSIGNED, i, 1, MPI_COMM_WORLD, &status);
	  deallog << "got message '" << buf << "' from CPU " << i+1 << "!" << std::endl;
	  Assert(buf == i, ExcInternalError());

	}
      else if (myid==i )
	{
	  MPI_Send(&myid, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
	}


    }
  if (myid==0)
    deallog << "done" << std::endl;

}


int main(int argc, char *argv[])
{
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
  Utilities::MPI::MPI_InitFinalize mpi (argc, argv);
#else
  (void)argc;
  (void)argv;
  compile_time_error;

#endif

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    {
      std::ofstream logfile(output_file_for_mpi("simple_mpi_01").c_str());
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      deallog.push("mpi");
      test_mpi();
      deallog.pop();
    }
  else
    test_mpi();
}
