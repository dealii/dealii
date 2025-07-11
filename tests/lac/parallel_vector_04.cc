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


// check that operator= resets ghosts, both if they have been set and if they
// have not been set

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

  // set ghost dof on remote process, no compress called. Since we don't want to
  // call compress we cannot use import
  auto partitioner = v.get_partitioner();
  if (myid > 0)
    {
      unsigned int local_index = partitioner->global_to_local(1);
      double      *values_dev  = v.get_values();
      Kokkos::deep_copy(
        Kokkos::View<double, MemorySpace::Default::kokkos_space>(values_dev +
                                                                 local_index),
        7);
    }

  unsigned int allocated_size = local_relevant.n_elements();
  Kokkos::View<double *, MemorySpace::Default::kokkos_space> v_device(
    v.get_values(), allocated_size);
  Kokkos::View<double *, Kokkos::HostSpace> v_host("v_host", allocated_size);
  Kokkos::deep_copy(v_host, v_device);

  AssertThrow(v_host[partitioner->global_to_local(myid * 2)] == myid * 4.0,
              ExcInternalError());
  AssertThrow(v_host[partitioner->global_to_local(myid * 2 + 1)] ==
                myid * 4.0 + 2.0,
              ExcInternalError());

  if (myid > 0)
    AssertThrow(v_host[partitioner->global_to_local(1)] == 7.0,
                ExcInternalError());

  // reset to zero
  v = 0;

  Kokkos::deep_copy(v_host, v_device);
  AssertThrow(v_host[partitioner->global_to_local(myid * 2)] == 0.,
              ExcInternalError());
  AssertThrow(v_host[partitioner->global_to_local(myid * 2 + 1)] == 0.,
              ExcInternalError());

  // check that everything remains zero also
  // after compress
  v.compress(VectorOperation::add);

  Kokkos::deep_copy(v_host, v_device);
  AssertThrow(v_host[partitioner->global_to_local(myid * 2)] == 0.,
              ExcInternalError());
  AssertThrow(v_host[partitioner->global_to_local(myid * 2 + 1)] == 0.,
              ExcInternalError());

  // set element 1 on owning process to
  // something nonzero
  if (myid == 0)
    {
      unsigned int local_index = partitioner->global_to_local(1);
      double      *values_dev  = v.get_values();
      Kokkos::deep_copy(Kokkos::subview(v_device, local_index), 2);
    }
  if (myid > 0)
    {
      Kokkos::deep_copy(v_host, v_device);
      AssertThrow(v_host[partitioner->global_to_local(1)] == 0.,
                  ExcInternalError());
    }

  // check that all processors get the correct
  // value again, and that it is erased by
  // operator=
  v.update_ghost_values();

  Kokkos::deep_copy(v_host, v_device);
  AssertThrow(v_host[partitioner->global_to_local(1)] == 2.,
              ExcInternalError());

  v = 0;
  Kokkos::deep_copy(v_host, v_device);
  AssertThrow(v_host[partitioner->global_to_local(1)] == 0.,
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
