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

  for (unsigned int i=1; i<numprocs; ++i)
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
  Utilities::MPI::MPI_InitFinalize mpi (argc, argv);

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    {
      std::ofstream logfile("output");
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
