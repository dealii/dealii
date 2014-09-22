// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2013 by the deal.II authors
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


// check that operator= resets ghosts, both if they have been set and if they
// have not been set

#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/parallel_vector.h>
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
  IndexSet local_owned(numproc*2);
  local_owned.add_range(myid*2,myid*2+2);
  IndexSet local_relevant(numproc*2);
  local_relevant = local_owned;
  local_relevant.add_range(1,2);

  parallel::distributed::Vector<double> v(local_owned, local_relevant, MPI_COMM_WORLD);

  // set local values and check them
  v(myid*2)=myid*2.0;
  v(myid*2+1)=myid*2.0+1.0;

  v.compress(VectorOperation::insert);
  v*=2.0;

  Assert(v(myid*2) == myid*4.0, ExcInternalError());
  Assert(v(myid*2+1) == myid*4.0+2.0, ExcInternalError());

  // set ghost dof on remote process, no
  // compress called
  if (myid > 0)
    v(1) = 7;

  Assert(v(myid*2) == myid*4.0, ExcInternalError());
  Assert(v(myid*2+1) == myid*4.0+2.0, ExcInternalError());

  if (myid > 0)
    Assert (v(1) == 7.0, ExcInternalError());

  // reset to zero
  v = 0;

  Assert(v(myid*2) == 0., ExcInternalError());
  Assert(v(myid*2+1) == 0., ExcInternalError());

  // check that everything remains zero also
  // after compress
  v.compress(VectorOperation::add);

  Assert(v(myid*2) == 0., ExcInternalError());
  Assert(v(myid*2+1) == 0., ExcInternalError());

  // set element 1 on owning process to
  // something nonzero
  if (myid == 0)
    v(1) = 2.;
  if (myid > 0)
    Assert (v(1) == 0., ExcInternalError());

  // check that all processors get the correct
  // value again, and that it is erased by
  // operator=
  v.update_ghost_values();

  Assert (v(1) == 2.0, ExcInternalError());

  v = 0;
  Assert (v(1) == 0.0, ExcInternalError());

  if (myid == 0)
    deallog << "OK" << std::endl;
}



int main (int argc, char **argv)
{
  Utilities::System::MPI_InitFinalize mpi_initialization(argc, argv);

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
