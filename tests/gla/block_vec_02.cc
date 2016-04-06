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



// block vectors: test ConstraintMatrix::distribute (fails currently)

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
  if (numproc==1)
    block1.add_range(0,10);

  if (myid==0)
    block1.add_range(0,7);
  if (myid==1)
    block1.add_range(7,10);

  IndexSet block2(6);
  if (numproc==1)
    block2.add_range(0,6);

  if (myid==0)
    block2.add_range(0,2);
  if (myid==1)
    block2.add_range(2,6);

  std::vector<IndexSet> partitioning;
  partitioning.push_back(block1);
  partitioning.push_back(block2);

  typename LA::MPI::BlockVector v(partitioning, MPI_COMM_WORLD);
  v = 0.1;

  v(myid) = myid;
  v.compress(VectorOperation::insert);

  deallog << "size: " << v.size() << std::endl;
  deallog << "block(0).size: " << v.block(0).size() << std::endl;
  deallog << "block(1).size: " << v.block(1).size() << std::endl;
  if (block1.n_elements()>0)
    deallog << "my first entry: " << get_real_assert_zero_imag(v(block1.nth_index_in_set(0))) << std::endl;


  ConstraintMatrix cm;
  cm.add_line(0);
  cm.add_entry(0, 1, 3.0);
  cm.close();

  if (myid==0)
    deallog << "before: " << get_real_assert_zero_imag(v(0)) << std::endl;
  cm.distribute(v); // this should set x(0)= 3.0 * x(1) = 3.0
  if (myid==0)
    deallog << "after: " << get_real_assert_zero_imag(v(0)) << std::endl;

  // done
  if (myid==0)
    deallog << "OK" << std::endl;
}



int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);
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
