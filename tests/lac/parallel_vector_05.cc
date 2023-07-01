// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// check that compress(add) with zero add does not change the vector entry

#include <deal.II/base/cuda.h>
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
  Assert(rw_vector(myid * 2) == myid * 4.0, ExcInternalError());
  Assert(rw_vector(myid * 2 + 1) == myid * 4.0 + 2.0, ExcInternalError());

  // set ghost dof to zero on remote processors,
  // compress
  IndexSet                               ghost_set(numproc * 2);
  LinearAlgebra::ReadWriteVector<double> ghost_vector;
  if (myid > 0)
    {
      ghost_set.add_index(1);
      ghost_vector.reinit(ghost_set);
      ghost_vector(1) = 0;
    }
  else
    ghost_vector.reinit(ghost_set);


  v.import_elements(ghost_vector, VectorOperation::add);

  // check that nothing has changed
  rw_vector.import_elements(v, VectorOperation::insert);
  Assert(rw_vector(myid * 2) == myid * 4.0, ExcInternalError());
  Assert(rw_vector(myid * 2 + 1) == myid * 4.0 + 2.0, ExcInternalError());

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
