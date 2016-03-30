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



// reinit and copy_from for LA::MPI::SparseMatrix

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
  IndexSet local_relevant= local_active;
  local_relevant.add_range(0,1);

  CompressedSimpleSparsityPattern csp (local_relevant);

  for (unsigned int i=0; i<2*numproc; ++i)
    if (local_relevant.is_element(i))
      csp.add(i,i);

  csp.add(0,1);


  typename LA::MPI::SparseMatrix mat;
  mat.reinit (local_active, local_active, csp, MPI_COMM_WORLD);

  Assert(mat.n()==numproc*2, ExcInternalError());
  Assert(mat.m()==numproc*2, ExcInternalError());

  // set local values
  mat.set(myid*2,myid*2, myid*2.0);
  mat.set(myid*2+1,myid*2+1, myid*2.0+1.0);

  mat.compress(VectorOperation::insert);

  mat.add(0,1,1.0);

  mat.compress(VectorOperation::add);

  deallog << "l1-norm: " << mat.l1_norm() << std::endl;
  if (myid==0)
    deallog << "mat(0,1): " << get_real_assert_zero_imag(mat(0,1)) << std::endl;

  {
    // disabled:
    //typename LA::MPI::SparseMatrix mat2(mat);
  }
  {
    // disabled:
    //typename LA::MPI::SparseMatrix mat2;
    //mat2 = mat;
  }

  {
    deallog << "* reinit(other):" << std::endl;
    typename LA::MPI::SparseMatrix mat2;
    mat2.reinit(mat);
    deallog << "l1-norm: " << mat2.l1_norm() << std::endl;
    if (myid==0)
      deallog << "mat(0,1): " << get_real_assert_zero_imag(mat2(0,1)) << std::endl;
    mat2.add(0,1,0.3);
    mat2.compress(VectorOperation::add);
    if (myid==0)
      deallog << "after adding mat(0,1): " << get_real_assert_zero_imag(mat2(0,1)) << std::endl;
  }
  {
    deallog << "* copy_from():" << std::endl;
    typename LA::MPI::SparseMatrix mat2;
    // PETSc needs a reinit() befor copy_from()
    mat2.reinit (local_active, local_active, csp, MPI_COMM_WORLD);
    mat2.copy_from(mat);
    deallog << "l1-norm: " << mat2.l1_norm() << std::endl;
    if (myid==0)
      deallog << "mat(0,1): " << get_real_assert_zero_imag(mat2(0,1)) << std::endl;
    mat2.add(0,1,0.3);
    mat2.compress(VectorOperation::add);
    if (myid==0)
      deallog << "after adding mat(0,1): " << get_real_assert_zero_imag(mat2(0,1)) << std::endl;
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
