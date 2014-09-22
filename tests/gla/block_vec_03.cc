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



// copying of ghosted vectors

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
  partitioning.push_back(block2);

  std::vector<IndexSet> relevant = partitioning;
  relevant[0].add_index(0);
  relevant[1].add_range(0,numproc);
  
  typename LA::MPI::BlockVector v(partitioning, MPI_COMM_WORLD);
  typename LA::MPI::BlockVector v2(partitioning, relevant, MPI_COMM_WORLD);

  Assert(!v.has_ghost_elements(), ExcInternalError());
  Assert(v2.has_ghost_elements(), ExcInternalError());
  Assert(!v.block(0).has_ghost_elements(), ExcInternalError());
  Assert(!v.block(1).has_ghost_elements(), ExcInternalError());
  Assert(v2.block(0).has_ghost_elements(), ExcInternalError());
  Assert(v2.block(1).has_ghost_elements(), ExcInternalError());

  v.reinit(partitioning, relevant, MPI_COMM_WORLD);
  Assert(v.has_ghost_elements(), ExcInternalError());
  Assert(v.block(0).has_ghost_elements(), ExcInternalError());
  Assert(v.block(1).has_ghost_elements(), ExcInternalError());
  v.reinit(partitioning, MPI_COMM_WORLD);
  Assert(!v.has_ghost_elements(), ExcInternalError());
  

  typename LA::MPI::BlockVector v3 = v2;
  Assert(v3.has_ghost_elements(), ExcInternalError());

  v3 = v; // just copy data, keep ghosts
  Assert(v3.has_ghost_elements(), ExcInternalError());

  typename LA::MPI::Vector x = v2.block(0);
  Assert(x.has_ghost_elements(), ExcInternalError());
  
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
