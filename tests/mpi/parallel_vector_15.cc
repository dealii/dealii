// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2015 by the deal.II authors
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


// check that handling of ghost elements in parallel distributed vectors works
// appropriately when creating a vector from a non-ghosted source vector using
// the assignment operator

#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <fstream>
#include <iostream>
#include <vector>


void test ()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  if (myid==0) deallog << "numproc=" << numproc << std::endl;

  // processor 0 and 1 own 2 indices each, higher processors nothing, all are
  // ghosting global elements 1 and 3
  IndexSet local_owned(std::min(numproc*2, 4U));
  if (myid < 2)
    local_owned.add_range(myid*2,myid*2+2);
  IndexSet local_relevant(local_owned.size());
  local_relevant = local_owned;
  local_relevant.add_range(1,2);
  if (numproc > 1)
    local_relevant.add_range(3,4);

  LinearAlgebra::distributed::Vector<double> v(local_owned, local_relevant, MPI_COMM_WORLD);

  // set local values
  if (myid < 2)
    {
      v(myid*2)=myid*2.0;
      v(myid*2+1)=myid*2.0+1.0;
    }

  v.compress(VectorOperation::insert);

  if (myid==0)
    deallog << "v has ghost elements: " << v.has_ghost_elements() << std::endl;

  LinearAlgebra::distributed::Vector<double> w, x;
  w = v;
  if (myid==0)
    deallog << "w has ghost elements: " << w.has_ghost_elements() << std::endl;

  v.update_ghost_values();
  w = v;
  if (myid==0)
    deallog << "w has ghost elements: " << w.has_ghost_elements() << std::endl;

  v.zero_out_ghosts();
  w = v;
  if (myid==0)
    deallog << "w has ghost elements: " << w.has_ghost_elements() << std::endl;

  w.zero_out_ghosts();
  w = v;
  if (myid==0)
    deallog << "w has ghost elements: " << w.has_ghost_elements() << std::endl;

  v.update_ghost_values();
  x = v;
  if (myid==0)
    deallog << "x has ghost elements: " << x.has_ghost_elements() << std::endl;

  x.zero_out_ghosts();
  if (myid==0)
    deallog << "x has ghost elements: " << x.has_ghost_elements() << std::endl;

  if (myid==0)
    deallog << "OK" << std::endl;
}



int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, testing_max_num_threads());

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile("output");
      deallog.attach(logfile);
      deallog << std::setprecision(4);
      deallog.threshold_double(1.e-10);

      test();
    }
  else
    test();

}
