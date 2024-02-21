// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check LA::Vector::compress(VectorOperation::min/max) from ghosts

#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"


template <typename Number>
double
print_value(Number *values_dev, unsigned int index)
{
  Kokkos::View<Number, Kokkos::HostSpace> cpu_value("cpu_value");
  Kokkos::deep_copy(cpu_value, Kokkos::View<Number>(values_dev + index));
  return cpu_value();
}



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

  // create vector
  LinearAlgebra::distributed::Vector<double, MemorySpace::Default> v(
    local_owned, local_relevant, MPI_COMM_WORLD);
  Kokkos::View<double *, MemorySpace::Default::kokkos_space> v_device(
    v.get_values(), local_relevant.n_elements());
  const auto &partitioner = v.get_partitioner();

  // set local values
  Kokkos::deep_copy(Kokkos::subview(v_device,
                                    partitioner->global_to_local(myid * 2)),
                    myid * 2.0);
  Kokkos::deep_copy(Kokkos::subview(v_device,
                                    partitioner->global_to_local(myid * 2 + 1)),
                    myid * 2.0 + 1.0);
  v.compress(VectorOperation::add);
  v *= 2.0;

  // check setup of vectors
  deallog << myid << ":"
          << "first owned entry: "
          << print_value(v.get_values(), partitioner->global_to_local(myid * 2))
          << std::endl;
  deallog << myid << ":"
          << "second owned entry: "
          << print_value(v.get_values(),
                         partitioner->global_to_local(myid * 2 + 1))
          << std::endl;

  // set ghost dof on owning processor and maximize
  if (myid != 0)
    Kokkos::deep_copy(Kokkos::subview(v_device,
                                      partitioner->global_to_local(1)),
                      7. * myid);
  v.compress(VectorOperation::max);

  // import ghosts onto all procs
  v.update_ghost_values();

  // check
  deallog << myid << ":"
          << "ghost entry after max from owner: "
          << print_value(v.get_values(), partitioner->global_to_local(1))
          << std::endl;

  // ghosts are set to zero
  v.zero_out_ghost_values();

  // minimize
  v.compress(VectorOperation::min);
  v.update_ghost_values();

  // check
  deallog << myid << ":"
          << "ghost entry after min from zero: "
          << print_value(v.get_values(), partitioner->global_to_local(1))
          << std::endl;

  // set ghost dof on non-owning processors and minimize
  v.zero_out_ghost_values();
  if (myid == 0)
    Kokkos::deep_copy(Kokkos::subview(v_device,
                                      partitioner->global_to_local(1)),
                      -1);
  v.compress(VectorOperation::min);
  v.update_ghost_values();

  // check
  deallog << myid << ":"
          << "ghost entry after min from : "
          << print_value(v.get_values(), partitioner->global_to_local(1))
          << std::endl;

  // set vector to 1, zeros in ghosts except on owner where -1. is set
  v.zero_out_ghost_values();
  v = 1.0;
  if (myid == 0)
    Kokkos::deep_copy(Kokkos::subview(v_device,
                                      partitioner->global_to_local(1)),
                      -1);

  // maximize
  v.compress(VectorOperation::max);
  v.update_ghost_values();

  // even if only one value is set (-1. on owner), the other values
  // contribute a "0" and maximization receives zero and returns it
  deallog << myid << ":"
          << "ghost entry after max and partly init: "
          << print_value(v.get_values(), partitioner->global_to_local(1))
          << std::endl;

  // however, if the ghost value is set on all processors, the
  // maximum is -1:
  v.zero_out_ghost_values();
  v = 1.0;
  Kokkos::deep_copy(Kokkos::subview(v_device, partitioner->global_to_local(1)),
                    -1);
  v.compress(VectorOperation::max);
  v.update_ghost_values();
  deallog << myid << ":"
          << "ghost entry after max and full init: "
          << print_value(v.get_values(), partitioner->global_to_local(1))
          << std::endl;

  // what happens in case max is called two times and all values were smaller
  // than zero
  v.zero_out_ghost_values();
  v = -1.0;
  Kokkos::deep_copy(Kokkos::subview(v_device, partitioner->global_to_local(1)),
                    -1);
  v.compress(VectorOperation::max);
  deallog << myid << ":"
          << "ghost entry after first max: "
          << print_value(v.get_values(), partitioner->global_to_local(1))
          << std::endl;
  v.compress(VectorOperation::max);
  deallog << myid << ":"
          << "ghost entry after second max: "
          << print_value(v.get_values(), partitioner->global_to_local(1))
          << std::endl;

  if (myid == 0)
    deallog << "OK" << std::endl;
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  MPILogInitAll log;

  test();
}
