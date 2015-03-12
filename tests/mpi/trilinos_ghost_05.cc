// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2015 by the deal.II authors
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



// check correct behaviour of Trilinos ghosted vectors
// check if assignment from a normal vector works correctly and updates the ghost values

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

  if (myid==0)
    deallog << "numproc=" << numproc << std::endl;

  unsigned int ghostel=(numproc>1)?2:1;

  // each processor owns 2 indices and all
  // are ghosting one element
  IndexSet local_active(numproc*2);
  local_active.add_range(myid*2,myid*2+2);
  IndexSet local_relevant(numproc*2);
  local_relevant = local_active;
  local_relevant.add_range(ghostel,ghostel+1);

  TrilinosWrappers::MPI::Vector x(local_active, MPI_COMM_WORLD);
  TrilinosWrappers::MPI::Vector v(local_relevant, MPI_COMM_WORLD);

  // set local values
  x(myid*2)=myid*2.0;
  x(myid*2+1)=myid*2.0+1.0;

  // transfer to ghosted vector v and check
  x.compress(VectorOperation::insert);
  v=x;

  Assert(v(myid*2) == myid*2.0, ExcInternalError());
  Assert(v(myid*2+1) == myid*2.0+1.0, ExcInternalError());
  Assert(v(ghostel) == ghostel, ExcInternalError());

  // change x, transfer, and check again
  x*=2.0;
  v=x;

  Assert(v(myid*2) == myid*4.0, ExcInternalError());
  Assert(v(myid*2+1) == myid*4.0+2.0, ExcInternalError());
  Assert(v(ghostel) == 2.0*ghostel, ExcInternalError());

  if (myid == 0)
    {
      deallog << myid*2 << ":" << v(myid*2) << std::endl;
      deallog << myid*2+1 << ":" << v(myid*2+1) << std::endl;
    }

  if (myid == 0)
    deallog << "OK" << std::endl;
}



int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, numbers::invalid_unsigned_int);

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
