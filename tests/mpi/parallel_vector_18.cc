// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
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


// use same partition as in prarallel_vector_06 to check add() with pointers.

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


  // each processor from processor 1 to 8
  // owns 2 indices (the other processors do
  // not own any dof), and all processors are
  // ghosting element 1 (the second)
  IndexSet local_owned(std::min(16U,numproc*2));
  if (myid < 8)
    local_owned.add_range(myid*2,myid*2+2);
  IndexSet local_relevant(numproc*2);
  local_relevant = local_owned;
  local_relevant.add_range(1,2);

  LinearAlgebra::distributed::Vector<double> v(local_owned, local_owned, MPI_COMM_WORLD);

  // set local values
  if (myid < 8)
    {
      types::global_dof_index n_elements = 2;
      types::global_dof_index indices[n_elements];
      indices[0] = myid*2;
      indices[1] = myid*2+1;
      float values[2];
      values[0] = myid*2.0;
      values[1] = myid*2.0+1.0;
      v.add (n_elements, &indices[0], &values[0]);
    }
  v.compress(VectorOperation::insert);
  v*=2.0;
  if (myid < 8)
    {
      AssertThrow(v(myid*2) == myid*4.0, ExcInternalError());
      AssertThrow(v(myid*2+1) == myid*4.0+2.0, ExcInternalError());
    }

  if (myid == 0)
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
