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



// create LA::MPI::SparseMatrix where one CPU has no DoFs

// this bug got fixed in PETScWrappers in r29776

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

  IndexSet local_active(10);
  if (myid==0)
    local_active.add_range(0,10);
  IndexSet local_relevant= local_active;
  if (myid==1)
    local_relevant.add_range(5,10);
  local_relevant.add_range(0,1);

  CompressedSimpleSparsityPattern csp (local_relevant);

  for (unsigned int i=0; i<10; ++i)
    if (local_relevant.is_element(i))
      csp.add(i,i);

  csp.add(0,1);

  typename LA::MPI::SparseMatrix mat;
  mat.reinit (local_active, local_active, csp, MPI_COMM_WORLD);

  Assert(mat.n()==10, ExcInternalError());
  Assert(mat.m()==10, ExcInternalError());

  mat.set(0,0,0.1);
  mat.set(1,1,0.1);

  MPI_Barrier(MPI_COMM_WORLD);

  mat.compress(VectorOperation::insert);

  mat.add(0,1,1.0);

  mat.compress(VectorOperation::add);

  // check local values
  if (myid==0)
    {
      deallog << "1,1 : " << get_real_assert_zero_imag(mat(1,1)) << std::endl;
      deallog << "0,1 : " << get_real_assert_zero_imag(mat(0,1)) << std::endl;
    }

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

  // compile, don't run
  //if (myid==9999)
  //  test<LA_Dummy>();


}
