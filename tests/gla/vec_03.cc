// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// assignment of ghost vectors

#include <deal.II/base/index_set.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/generic_linear_algebra.h>

#include <iostream>
#include <vector>

#include "../tests.h"

#include "gla.h"

template <class LA>
void
test()
{
  unsigned int myid    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  if (myid == 0)
    deallog << "numproc=" << numproc << std::endl;

  // each processor owns 2 indices and all
  // are ghosting Element 1 (the second)

  IndexSet local_active(numproc * 2);
  local_active.add_range(myid * 2, myid * 2 + 2);
  IndexSet local_relevant(numproc * 2);
  local_relevant.add_range(1, 2);

  typename LA::MPI::Vector vb(local_active, MPI_COMM_WORLD);
  typename LA::MPI::Vector v(local_active, local_relevant, MPI_COMM_WORLD);

  vb = 1.0;

  // set local values
  vb(myid * 2)     = myid * 2.0;
  vb(myid * 2 + 1) = myid * 2.0 + 1.0;

  vb.compress(VectorOperation::insert);
  vb *= 2.0;
  v = vb;

  Assert(vb.size() == numproc * 2, ExcInternalError());
  Assert(v.size() == numproc * 2, ExcInternalError());

  Assert(!vb.has_ghost_elements(), ExcInternalError());
  Assert(v.has_ghost_elements(), ExcInternalError());

  // check local values
  {
    deallog << myid * 2 << ": " << get_real_assert_zero_imag(v(myid * 2))
            << std::endl;
    deallog << myid * 2 + 1 << ": "
            << get_real_assert_zero_imag(v(myid * 2 + 1)) << std::endl;
  }

  Assert(get_real_assert_zero_imag(v(myid * 2)) == myid * 4.0,
         ExcInternalError());
  Assert(get_real_assert_zero_imag(v(myid * 2 + 1)) == myid * 4.0 + 2.0,
         ExcInternalError());

  using scalar_type = typename LA::MPI::BlockVector::value_type;
  AffineConstraints<scalar_type> cm(local_active,
                                    complete_index_set(local_active.size()));
  cm.add_line(1);
  cm.add_entry(1, 2, 3.0);
  cm.close();

  if (myid == 0)
    deallog << "before: " << get_real_assert_zero_imag(vb(1)) << std::endl;
  cm.distribute(vb); // this should set x(1)= 3.0 * x(2) = 12.0
  if (myid == 0)
    deallog << "after: " << get_real_assert_zero_imag(vb(1)) << std::endl;

  // done
  if (myid == 0)
    deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;
  {
    deallog.push("PETSc");
    test<LA_PETSc>();
    deallog.pop();
    deallog.push("Trilinos");
    test<LA_Trilinos>();
    deallog.pop();
  }
}
