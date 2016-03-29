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



// assignment of ghost vectors

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

  // each processor owns 2 indices and all
  // are ghosting Element 1 (the second)

  IndexSet local_active(numproc*2);
  local_active.add_range(myid*2,myid*2+2);
  IndexSet local_relevant(numproc*2);
  local_relevant.add_range(1,2);

  typename LA::MPI::Vector vb(local_active, MPI_COMM_WORLD);
  typename LA::MPI::Vector v(local_active, local_relevant, MPI_COMM_WORLD);

  vb = 1.0;

  // set local values
  vb(myid*2)=myid*2.0;
  vb(myid*2+1)=myid*2.0+1.0;

  vb.compress(VectorOperation::insert);
  vb*=2.0;
  v=vb;

  Assert(vb.size() == numproc*2, ExcInternalError());
  Assert(v.size() == numproc*2, ExcInternalError());

  Assert(!vb.has_ghost_elements(), ExcInternalError());
  Assert(v.has_ghost_elements(), ExcInternalError());

  // check local values
  {
    deallog << myid*2 << ": " << get_real_assert_zero_imag(v(myid*2)) << std::endl;
    deallog << myid*2+1 << ": " << get_real_assert_zero_imag(v(myid*2+1)) << std::endl;
  }

  Assert(get_real_assert_zero_imag(v(myid*2)) == myid*4.0, ExcInternalError());
  Assert(get_real_assert_zero_imag(v(myid*2+1)) == myid*4.0+2.0, ExcInternalError());

  ConstraintMatrix cm;
  cm.add_line(1);
  cm.add_entry(1 ,2, 3.0);
  cm.close();

  if (myid==0)
    deallog << "before: " << get_real_assert_zero_imag(vb(1)) << std::endl;
  cm.distribute(vb); // this should set x(1)= 3.0 * x(2) = 12.0
  if (myid==0)
    deallog << "after: " << get_real_assert_zero_imag(vb(1)) << std::endl;

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
