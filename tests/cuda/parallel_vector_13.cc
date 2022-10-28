// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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


// check for ghosts on parallel vector: similar to parallel_vector_03, but
// setting where one ghost is zero and should not have an effect on vector
// entries

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

  LinearAlgebra::distributed::Vector<double, MemorySpace::CUDA> v(
    local_owned, local_relevant, MPI_COMM_WORLD);

  // set local values and check them
  LinearAlgebra::ReadWriteVector<double> rw_vector(local_owned);
  rw_vector(myid * 2)     = myid * 2.0;
  rw_vector(myid * 2 + 1) = myid * 2.0 + 1.0;

  v.import(rw_vector, VectorOperation::insert);
  v *= 2.0;
  v.add(1.0);

  rw_vector.import(v, VectorOperation::insert);
  AssertThrow(rw_vector(myid * 2) == myid * 4.0 + 1, ExcInternalError());
  AssertThrow(rw_vector(myid * 2 + 1) == myid * 4.0 + 3.0, ExcInternalError());

  // set ghost dof on all processors, compress
  // (insert mode)
  IndexSet index(numproc * 2);
  index.add_index(1);
  LinearAlgebra::ReadWriteVector<double> local_rw_vector(index);
  local_rw_vector(1) = 7;
  v.import(local_rw_vector, VectorOperation::insert);

  {
    rw_vector.import(v, VectorOperation::insert);
    deallog << myid * 2 << ":" << rw_vector(myid * 2) << std::endl;
    deallog << myid * 2 + 1 << ":" << rw_vector(myid * 2 + 1) << std::endl;
  }

  local_rw_vector(1) = -7;
  v.import(local_rw_vector, VectorOperation::insert);

  {
    rw_vector.import(v, VectorOperation::insert);
    deallog << myid * 2 << ":" << rw_vector(myid * 2) << std::endl;
    deallog << myid * 2 + 1 << ":" << rw_vector(myid * 2 + 1) << std::endl;
  }

  // import ghosts onto all procs
  v.update_ghost_values();
  local_rw_vector.import(v, VectorOperation::insert);
  AssertThrow(local_rw_vector(1) == -7.0, ExcInternalError());

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

  init_cuda(true);

  if (myid == 0)
    {
      initlog();
      deallog << std::setprecision(4);

      test();
    }
  else
    test();
}
