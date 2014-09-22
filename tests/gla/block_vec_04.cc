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



// allow operator= for BlockVectors if dest is same size or empty

#include "../tests.h"
#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/constraint_matrix.h>
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

  IndexSet block1(10);
  if (myid==0)
    block1.add_range(0,7);
  if (myid==1)
    block1.add_range(7,10);

  IndexSet block2(numproc);
  block2.add_index(myid);

  std::vector<IndexSet> partitioning;
  partitioning.push_back(block1);
  typename LA::MPI::BlockVector v_1(partitioning, MPI_COMM_WORLD);
  partitioning.push_back(block2);

  std::vector<IndexSet> relevant = partitioning;
  relevant[0].add_index(0);
  relevant[1].add_range(0,numproc);
  
  typename LA::MPI::BlockVector v_2(partitioning, MPI_COMM_WORLD);

  {
    typename LA::MPI::BlockVector x(partitioning, MPI_COMM_WORLD);
    x=v_2;
  }
  {
    typename LA::MPI::BlockVector x;
    x=v_2;
  }
  {
    deal_II_exceptions::disable_abort_on_exception();
    try
      {	
	typename LA::MPI::BlockVector x=v_1;
	x=v_2; // error
      }
    catch (const ExceptionBase &e)
      {
	deallog << "Exception: " << e.get_exc_name() << std::endl;
      }
  }
  
  // done
  if (myid==0)
    deallog << "OK" << std::endl;
}



int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll log;
  {
    deallog.push("PETSc");
    test<LA_PETSc>();
    deallog.pop();
    deallog.push("Trilinos");
    test<LA_Trilinos>();
    deallog.pop();
  }

}
