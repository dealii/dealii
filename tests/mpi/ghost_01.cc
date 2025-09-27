// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test correct handling of ghost elements in PETScWrappers::mpi::vectors

#include <deal.II/base/index_set.h>

#include <deal.II/lac/petsc_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"


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

  PETScWrappers::MPI::Vector vb(local_active, MPI_COMM_WORLD);
  PETScWrappers::MPI::Vector v(local_active, local_relevant, MPI_COMM_WORLD);

  // set local values
  vb(myid * 2)     = myid * 2.0;
  vb(myid * 2 + 1) = myid * 2.0 + 1.0;

  vb.compress(VectorOperation::insert);
  vb *= 2.0;
  v = vb;

  Assert(!vb.has_ghost_elements(), ExcInternalError());
  Assert(v.has_ghost_elements(), ExcInternalError());

  // check local values
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      deallog << myid * 2 << ':' << get_real_assert_zero_imag(v(myid * 2))
              << std::endl;
      deallog << myid * 2 + 1 << ':'
              << get_real_assert_zero_imag(v(myid * 2 + 1)) << std::endl;
    }

  Assert(get_real_assert_zero_imag(v(myid * 2)) == myid * 4.0,
         ExcInternalError());
  Assert(get_real_assert_zero_imag(v(myid * 2 + 1)) == myid * 4.0 + 2.0,
         ExcInternalError());


  // check ghost values
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    deallog << "ghost: " << get_real_assert_zero_imag(v(1)) << std::endl;
  Assert(get_real_assert_zero_imag(v(1)) == 2.0, ExcInternalError());

  // done
  if (myid == 0)
    deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      initlog();
      deallog << std::setprecision(4);

      test();
    }
  else
    test();
}
