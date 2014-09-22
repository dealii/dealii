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

  // select a few destinations
  std::vector<unsigned int> destinations;
  for (unsigned int i=0; i<3+myid/3; ++i)
    if ((myid + 17*i) % numprocs != myid)
      destinations.push_back((myid + 17*i) % numprocs);

  if (myid == 0)
    {
      deallog << "Processor 0 wants to send to ";
      for (unsigned int i=0; i<destinations.size(); ++i)
        deallog << destinations[i] << ' ';
      deallog << std::endl;

      for (unsigned int p=1; p<numprocs; ++p)
        {
          MPI_Status status;
          unsigned int size = 0;
          MPI_Recv (&size, 1, MPI_UNSIGNED, p, 0, MPI_COMM_WORLD, &status);

          std::vector<unsigned int> dest (size);
          MPI_Recv (&dest[0], size, MPI_UNSIGNED, p, 0, MPI_COMM_WORLD, &status);

          deallog << "Processor " << p << " wants to send to ";
          for (unsigned int i=0; i<size; ++i)
            deallog << dest[i] << ' ';
          deallog << std::endl;
        }
    }
  else
    {
      unsigned int size = destinations.size();
      MPI_Send (&size, 1, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD);
      MPI_Send (&destinations[0], size, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD);
    }


  if (myid == 0)
    deallog << "Exchanging data..." << std::endl;

  std::vector<unsigned int> origins
    = Utilities::MPI::compute_point_to_point_communication_pattern (MPI_COMM_WORLD,
        destinations);

  if (myid == 0)
    {
      deallog << "Processor 0 will receive from ";
      for (unsigned int i=0; i<origins.size(); ++i)
        deallog << origins[i] << ' ';
      deallog << std::endl;

      for (unsigned int p=1; p<numprocs; ++p)
        {
          MPI_Status status;
          unsigned int size = 0;
          MPI_Recv (&size, 1, MPI_UNSIGNED, p, 0, MPI_COMM_WORLD, &status);

          std::vector<unsigned int> orig (size);
          MPI_Recv (&orig[0], size, MPI_UNSIGNED, p, 0, MPI_COMM_WORLD, &status);

          deallog << "Processor " << p << " will receive from ";
          for (unsigned int i=0; i<size; ++i)
            deallog << orig[i] << ' ';
          deallog << std::endl;
        }
    }
  else
    {
      unsigned int size = origins.size();
      MPI_Send (&size, 1, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD);
      MPI_Send (&origins[0], size, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD);
    }
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
