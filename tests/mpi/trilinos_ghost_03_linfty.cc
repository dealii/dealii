// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// while computing the l2norm on ghosted vectors is unsafe, computing
// the linfty norm should be ok. test this

#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/trilinos_vector.h>

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
  // are ghosting element 1 (the second)
  IndexSet local_active(numproc * 2);
  local_active.add_range(myid * 2, myid * 2 + 2);
  IndexSet local_relevant(numproc * 2);
  local_relevant = local_active;
  local_relevant.add_range(1, 2);

  TrilinosWrappers::MPI::Vector v(local_active, MPI_COMM_WORLD);
  TrilinosWrappers::MPI::Vector v_tmp(local_relevant, MPI_COMM_WORLD);


  // set local values
  v(myid * 2)     = myid * 2;
  v(myid * 2 + 1) = myid * 2 + 1;

  v.compress(VectorOperation::insert);

  v_tmp.reinit(v, false, true);

  bool   exc = false;
  double norm;
  try
    {
      norm = v_tmp.linfty_norm();
    }
  catch (const TrilinosWrappers::MPI::Vector::ExcTrilinosError &e)
    {
      deallog << e.get_exc_name() << std::endl;
      exc = true;
    }

  Assert(exc == false, ExcInternalError());
  Assert(norm == 2. * Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) - 1,
         ExcInternalError());
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      deallog << norm << std::endl;
      deallog << "OK" << std::endl;
    }
}



int
main(int argc, char **argv)
{
  deal_II_exceptions::disable_abort_on_exception();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

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
