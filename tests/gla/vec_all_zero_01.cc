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



// check that all_zero is working correctly (as a collective operation)

#include "../tests.h"
#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/base/index_set.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>

#include "gla.h"

template <class LA>
void test ()
{
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  if (myid==0)
    deallog << "numproc=" << numproc << std::endl;

  // each processor owns 2 indices and all
  // are ghosting Element 1 (the second)

  IndexSet local_active(numproc*2);
  local_active.add_range(myid*2,myid*2+2);
  IndexSet local_relevant(numproc*2);
  local_relevant.add_range(1,2);


  typename LA::MPI::Vector x;  
  x.reinit(local_active, MPI_COMM_WORLD);
  x=0;
  typename LA::MPI::Vector g(local_active, local_relevant, MPI_COMM_WORLD);
  g=x;
  deallog << "all_zero? " << g.all_zero() << " (should be true)" << std::endl;
  if (myid==0)
    x(0)=1.0;
  x.compress(VectorOperation::insert);
  g=x;
  deallog << "all_zero? " << g.all_zero() << " (should be false)" << std::endl;
  
  // done
  if (myid==0)
    deallog << "OK" << std::endl;
}



int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
  MPILogInitAll log;
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);

  {
    deallog.push("PETSc");
    test<LA_PETSc>();
    deallog.pop();
    deallog.push("Trilinos");
    test<LA_Trilinos>();
    deallog.pop();
  }

  if (myid==9999)
    test<LA_Dummy>(); // don't execute

}
