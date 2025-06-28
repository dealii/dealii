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


// check parallel_vector assignment and ghost values. Same vector layout as in
// parallel_vector_07.cc

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

  const unsigned int set = 200;
  AssertIndexRange(numproc, set - 2);
  const unsigned int local_size  = set - myid;
  unsigned int       global_size = 0;
  unsigned int       my_start    = 0;
  for (unsigned int i = 0; i < numproc; ++i)
    {
      global_size += set - i;
      if (i < myid)
        my_start += set - i;
    }
  // each processor owns some indices and all
  // are ghosting elements from three
  // processors (the second). some entries
  // are right around the border between two
  // processors
  IndexSet local_owned(global_size);
  local_owned.add_range(my_start, my_start + local_size);
  IndexSet local_relevant(global_size);
  local_relevant                            = local_owned;
  const std::size_t n_ghosts                = 10;
  unsigned int      ghost_indices[n_ghosts] = {
    1, 2, 13, set - 3, set - 2, set - 1, set, set + 1, set + 2, set + 3};
  local_relevant.add_indices(std::begin(ghost_indices),
                             std::end(ghost_indices));

  // v has ghosts, w has none. set some entries
  // on w, copy into v and check if they are
  // there
  LinearAlgebra::distributed::Vector<double, MemorySpace::Default> v(
    local_owned, local_relevant, MPI_COMM_WORLD);
  LinearAlgebra::distributed::Vector<double, MemorySpace::Default> w(
    local_owned, local_owned, MPI_COMM_WORLD);

  // set a few of the local elements
  LinearAlgebra::ReadWriteVector<double> rw_vector(local_owned);
  for (unsigned i = 0; i < local_size; ++i)
    rw_vector.local_element(i) = 2.0 * (i + my_start);
  w.import_elements(rw_vector, VectorOperation::insert);

  v = w;
  v.update_ghost_values();

  // Check local and ghost values with move assignment:
  {
    decltype(v) v_move, v_copy = v;
    v_move = std::move(v_copy);

    LinearAlgebra::ReadWriteVector<double> v_rw(local_relevant);
    v_rw.import_elements(v, VectorOperation::insert);
    LinearAlgebra::ReadWriteVector<double> v_move_rw(local_relevant);
    v_move_rw.import_elements(v_move, VectorOperation::insert);

    for (unsigned int i = 0; i < v_rw.locally_owned_size(); ++i)
      AssertThrow(v_move_rw.local_element(i) == v_rw.local_element(i),
                  ExcInternalError());

    // v_copy should now be empty
    AssertThrow(v_copy.locally_owned_size() == 0, ExcInternalError());
    AssertThrow(v_copy.get_partitioner()->size() == 0, ExcInternalError());
    AssertThrow(v_copy.get_partitioner()->n_ghost_indices() == 0,
                ExcInternalError());
    AssertThrow(v_copy.get_mpi_communicator() == MPI_COMM_SELF,
                ExcInternalError());
  }

  // check local values for correctness
  rw_vector.import_elements(v, VectorOperation::insert);
  for (unsigned int i = 0; i < local_size; ++i)
    AssertThrow(rw_vector.local_element(i) == 2.0 * (i + my_start),
                ExcInternalError());

  // check non-local entries on all processors
  LinearAlgebra::ReadWriteVector<double> ghost_vector(local_relevant);
  ghost_vector.import_elements(v, VectorOperation::insert);
  for (unsigned int i = 0; i < 10; ++i)
    AssertThrow(ghost_vector(ghost_indices[i]) == 2. * ghost_indices[i],
                ExcInternalError());

  // now the same again, but import ghosts automatically because v had ghosts
  // set before calling operator =
  v.reinit(local_owned, local_relevant, MPI_COMM_WORLD);
  v.update_ghost_values();
  v = w;

  // check local values for correctness
  rw_vector.import_elements(v, VectorOperation::insert);
  for (unsigned int i = 0; i < local_size; ++i)
    AssertThrow(rw_vector.local_element(i) == 2.0 * (i + my_start),
                ExcInternalError());

  // check non-local entries on all processors
  ghost_vector.import_elements(v, VectorOperation::insert);
  for (unsigned int i = 0; i < 10; ++i)
    AssertThrow(ghost_vector(ghost_indices[i]) == 2. * ghost_indices[i],
                ExcInternalError());

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
