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



// test reinit()

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

  typename LA::MPI::Vector v(local_active, MPI_COMM_WORLD);
  typename LA::MPI::Vector g(local_active, local_relevant, MPI_COMM_WORLD);

  // set local values
  v(myid*2)=myid*2.0;
  v(myid*2+1)=myid*2.0+1.0;

  v.compress(VectorOperation::insert);

  g=v;

  IndexSet local_active_big(numproc*3);
  local_active_big.add_range(myid*3,myid*3+3);
  typename LA::MPI::Vector big(local_active_big, MPI_COMM_WORLD);

  typename LA::MPI::Vector x;
  Assert(!x.has_ghost_elements(), ExcInternalError());
  Assert(x.size()==0, ExcInternalError());
  x.reinit(v);
  Assert(!x.has_ghost_elements(), ExcInternalError());
  Assert(x.size()==v.size(), ExcInternalError());

  x.reinit(big);
  Assert(!x.has_ghost_elements(), ExcInternalError());
  Assert(x.size()==big.size(), ExcInternalError());

  x.reinit(v);
  Assert(!x.has_ghost_elements(), ExcInternalError());
  Assert(x.size()==v.size(), ExcInternalError());

  x.reinit(g);
  Assert(x.has_ghost_elements(), ExcInternalError());
  Assert(x.size()==g.size(), ExcInternalError());
  deallog << get_real_assert_zero_imag(x(1)) << std::endl;

  x.reinit(v);
  Assert(!x.has_ghost_elements(), ExcInternalError());
  Assert(x.size()==v.size(), ExcInternalError());

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
