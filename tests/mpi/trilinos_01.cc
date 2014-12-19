// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2013 by the deal.II authors
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



// check Trilinos has_ghost_elements() if run on one CPU only

#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/trilinos_vector.h>
#include <fstream>
#include <iostream>
#include <vector>


void test ()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  if (myid==0) deallog << "numproc=" << numproc << std::endl;


  // each processor owns 2 indices and all
  // are ghosting element 1 (the second)
  IndexSet local_active(numproc*2);
  local_active.add_range(myid*2,myid*2+2);
  IndexSet local_relevant(numproc*2);
  local_relevant = local_active;
  local_relevant.add_range(1,2);

  TrilinosWrappers::MPI::Vector v(local_active, MPI_COMM_WORLD);
  TrilinosWrappers::MPI::Vector v_tmp(local_relevant, MPI_COMM_WORLD);

  // only one CPU checks has_ghost_elements()
  MPI_Barrier(MPI_COMM_WORLD);
  if (myid==0)
    {
      deallog << v_tmp.has_ghost_elements() << std::endl;
      int dummy=0;
      MPI_Send(&dummy, 1, MPI_INT, 1, 12345, MPI_COMM_WORLD);
    }

  if (myid==1)
    {
      MPI_Status status;
      MPI_Request request;
      int flag=0;
      int tests=0;

      while (!flag && tests<10)
        {
          tests++;
          MPI_Iprobe(0, 12345, MPI_COMM_WORLD, &flag, &status);
          std::cout << flag << std::endl;

          sleep(1);
        }

      Assert(flag!=0, ExcMessage("hang in has_ghost_elements()"));
    }


  MPI_Barrier(MPI_COMM_WORLD);

  deallog << v_tmp.has_ghost_elements() << std::endl;

  if (myid==0)
    deallog << "OK" << std::endl;
}



int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile("output");
      deallog.attach(logfile);
      deallog << std::setprecision(4);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      test();
    }
  else
    test();

}
