// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check set of ghosts for parallel vector

#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/read_write_vector.h>

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
  IndexSet local_owned(numproc * 2);
  local_owned.add_range(myid * 2, myid * 2 + 2);
  IndexSet local_relevant(numproc * 2);
  local_relevant = local_owned;
  local_relevant.add_range(1, 2);

  LinearAlgebra::distributed::Vector<double, MemorySpace::Default> v(
    local_owned, local_relevant, MPI_COMM_WORLD);

  // set local values and check them
  LinearAlgebra::ReadWriteVector<double> rw_vector(local_owned);
  rw_vector(myid * 2)     = myid * 2.0;
  rw_vector(myid * 2 + 1) = myid * 2.0 + 1.0;
  v.import_elements(rw_vector, VectorOperation::insert);

  v *= 2.0;

  rw_vector.import_elements(v, VectorOperation::insert);
  AssertThrow(rw_vector(myid * 2) == myid * 4.0, ExcInternalError());
  AssertThrow(rw_vector(myid * 2 + 1) == myid * 4.0 + 2.0, ExcInternalError());

  // set ghost dof on all processors, compress
  // (insert mode)
  IndexSet ghost_entry(numproc * 2);
  ghost_entry.add_index(1);
  LinearAlgebra::ReadWriteVector<double> rw_ghost_entry(ghost_entry);
  rw_ghost_entry(1) = 7.;
  v.import_elements(rw_ghost_entry, VectorOperation::insert);

  rw_vector.import_elements(v, VectorOperation::insert);
  if (myid == 0)
    {
      deallog << myid * 2 << ":" << rw_vector(myid * 2) << std::endl;
      deallog << myid * 2 + 1 << ":" << rw_vector(myid * 2 + 1) << std::endl;
    }
  // import ghosts onto all procs
  v.update_ghost_values();
  rw_ghost_entry(1) = 0.;
  rw_ghost_entry.import_elements(v, VectorOperation::insert);
  AssertThrow(rw_ghost_entry(1) == 7.0, ExcInternalError());

  // check l2 norm
  const double l2_norm = v.l2_norm();
  if (myid == 0)
    deallog << "L2 norm: " << l2_norm << std::endl;

  if (myid == 0)
    deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
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
