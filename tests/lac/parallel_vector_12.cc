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


// check LinearAlgebra::distributed::Vector::swap

#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/read_write_vector.h>

#include <iostream>
#include <vector>

#include "../tests.h"

DeclException2(ExcNonEqual,
               double,
               double,
               << "Left compare: " << arg1 << ", right compare: " << arg2);

void
test()
{
  unsigned int myid    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int numproc = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  if (myid == 0)
    deallog << "numproc=" << numproc << std::endl;

  // vector 0:
  // global size: 20, local_size: 3 as long as
  // less than 20
  const unsigned int local_size0  = 3;
  const unsigned int global_size0 = std::min(20U, local_size0 * numproc);
  const unsigned int my_start0    = std::min(local_size0 * myid, global_size0);
  const unsigned int my_end0 = std::min(local_size0 * (myid + 1), global_size0);
  const unsigned int actual_local_size0 = my_end0 - my_start0;

  IndexSet local_owned0(global_size0);
  if (my_end0 > my_start0)
    local_owned0.add_range(static_cast<unsigned int>(my_start0),
                           static_cast<unsigned int>(my_end0));
  IndexSet local_relevant0(global_size0);
  local_relevant0 = local_owned0;
  local_relevant0.add_index(2);
  if (numproc > 2)
    local_relevant0.add_index(8);

  LinearAlgebra::distributed::Vector<double, MemorySpace::Default> v0(
    local_owned0, local_relevant0, MPI_COMM_WORLD);

  // vector1: local size 4
  const unsigned int local_size1  = 4;
  const unsigned int global_size1 = local_size1 * numproc;
  const int          my_start1    = local_size1 * myid;
  const int          my_end1      = local_size1 * (myid + 1);

  IndexSet local_owned1(global_size1);
  local_owned1.add_range(static_cast<unsigned int>(my_start1),
                         static_cast<unsigned int>(my_end1));
  IndexSet local_relevant1(global_size1);
  local_relevant1 = local_owned1;
  local_relevant1.add_index(0);
  local_relevant1.add_index(2);
  if (numproc > 2)
    {
      local_relevant1.add_index(8);
      local_relevant1.add_index(10);
    }

  LinearAlgebra::distributed::Vector<double, MemorySpace::Default> v1(
    local_owned1, local_relevant1, MPI_COMM_WORLD);

  v0 = 1;
  v1 = 2;
  // check assignment in initial state
  LinearAlgebra::ReadWriteVector<double> v0_rw(local_owned0);
  v0_rw.import_elements(v0, VectorOperation::insert);
  LinearAlgebra::ReadWriteVector<double> v1_rw(local_owned1);
  v1_rw.import_elements(v1, VectorOperation::insert);
  for (unsigned int i = 0; i < v0.locally_owned_size(); ++i)
    AssertThrow(v0_rw.local_element(i) == 1.,
                ExcNonEqual(v0_rw.local_element(i), 1.));
  for (unsigned int i = 0; i < v1.locally_owned_size(); ++i)
    AssertThrow(v1_rw.local_element(i) == 2.,
                ExcNonEqual(v1_rw.local_element(i), 2.));

  // check ghost elements in initial state
  v0.update_ghost_values();
  v1.update_ghost_values();
  LinearAlgebra::ReadWriteVector<double> v0_ghost_rw(local_relevant0);
  v0_ghost_rw.import_elements(v0, VectorOperation::insert);
  LinearAlgebra::ReadWriteVector<double> v1_ghost_rw(local_relevant1);
  v1_ghost_rw.import_elements(v1, VectorOperation::insert);
  AssertThrow(v0_ghost_rw(2) == 1., ExcNonEqual(v0_ghost_rw(2), 1.));
  if (numproc > 2)
    AssertThrow(v0_ghost_rw(8) == 1., ExcNonEqual(v0_ghost_rw(8), 2.));
  AssertThrow(v1_ghost_rw(0) == 2., ExcNonEqual(v1_ghost_rw(0), 2.));
  AssertThrow(v1_ghost_rw(2) == 2., ExcNonEqual(v1_ghost_rw(2), 2.));
  if (numproc > 2)
    {
      AssertThrow(v1_ghost_rw(8) == 2., ExcNonEqual(v1_ghost_rw(8), 2.));
      AssertThrow(v1_ghost_rw(10) == 2., ExcNonEqual(v1_ghost_rw(10), 2.));
    }
  MPI_Barrier(MPI_COMM_WORLD);
  if (myid == 0)
    deallog << "Initial set and ghost update OK" << std::endl;

  // now swap v1 and v0
  v0.swap(v1);
  AssertDimension(v0.locally_owned_size(), local_size1);
  AssertDimension(v1.locally_owned_size(), actual_local_size0);
  AssertDimension(v0.size(), global_size1);
  AssertDimension(v1.size(), global_size0);
  v1_rw.import_elements(v0, VectorOperation::insert);
  for (unsigned int i = 0; i < local_size1; ++i)
    AssertThrow(v1_rw.local_element(i) == 2.,
                ExcNonEqual(v1_rw.local_element(i), 2.));
  v0_rw.import_elements(v1, VectorOperation::insert);
  for (unsigned int i = 0; i < actual_local_size0; ++i)
    AssertThrow(v0_rw.local_element(i) == 1.,
                ExcNonEqual(v0_rw.local_element(i), 1.));
  MPI_Barrier(MPI_COMM_WORLD);
  if (myid == 0)
    deallog << "First swap OK" << std::endl;
  v0.update_ghost_values();
  v1.update_ghost_values();
  v0_ghost_rw.import_elements(v1, VectorOperation::insert);
  AssertThrow(v0_ghost_rw(2) == 1., ExcNonEqual(v0_ghost_rw(2), 1.));
  if (numproc > 2)
    AssertThrow(v0_ghost_rw(8) == 1., ExcNonEqual(v0_ghost_rw(8), 1.));
  v1_ghost_rw.import_elements(v0, VectorOperation::insert);
  AssertThrow(v1_ghost_rw(0) == 2., ExcNonEqual(v1_ghost_rw(0), 2.));
  AssertThrow(v1_ghost_rw(2) == 2., ExcNonEqual(v1_ghost_rw(2), 2.));
  if (numproc > 2)
    {
      AssertThrow(v1_ghost_rw(8) == 2., ExcNonEqual(v1_ghost_rw(8), 2.));
      AssertThrow(v1_ghost_rw(10) == 2., ExcNonEqual(v1_ghost_rw(10), 2.));
    }
  MPI_Barrier(MPI_COMM_WORLD);
  if (myid == 0)
    deallog << "Ghost values after first swap OK" << std::endl;

  // now set the vectors to some different
  // values and check the ghost values again
  v0 = 7.;
  v1 = 42.;
  v0.update_ghost_values();
  v1.update_ghost_values();
  v0_ghost_rw.import_elements(v1, VectorOperation::insert);
  AssertThrow(v0_ghost_rw(2) == 42., ExcNonEqual(v0_ghost_rw(2), 42.));
  if (numproc > 2)
    AssertThrow(v0_ghost_rw(8) == 42., ExcNonEqual(v0_ghost_rw(8), 42.));
  v1_ghost_rw.import_elements(v0, VectorOperation::insert);
  AssertThrow(v1_ghost_rw(0) == 7., ExcNonEqual(v1_ghost_rw(0), 7.));
  AssertThrow(v1_ghost_rw(2) == 7., ExcNonEqual(v1_ghost_rw(2), 7.));
  if (numproc > 2)
    {
      AssertThrow(v1_ghost_rw(8) == 7., ExcNonEqual(v1_ghost_rw(8), 7.));
      AssertThrow(v1_ghost_rw(10) == 7., ExcNonEqual(v1_ghost_rw(10), 7.));
    }
  MPI_Barrier(MPI_COMM_WORLD);
  if (myid == 0)
    deallog << "Ghost values after re-set OK" << std::endl;

  // swap with an empty vector
  LinearAlgebra::distributed::Vector<double, MemorySpace::Default> v2;
  v2.swap(v0);
  AssertDimension(v0.size(), 0);
  AssertDimension(v2.size(), global_size1);
  AssertDimension(v2.locally_owned_size(), local_size1);
  v1_rw.import_elements(v2, VectorOperation::insert);
  for (int i = my_start1; i < my_end1; ++i)
    AssertThrow(v1_rw(i) == 7., ExcNonEqual(v1_rw(i), 7.));
  MPI_Barrier(MPI_COMM_WORLD);
  if (myid == 0)
    deallog << "Second swap OK" << std::endl;
  v2 = -1.;
  v2.update_ghost_values();
  v1_ghost_rw.import_elements(v2, VectorOperation::insert);
  AssertThrow(v1_ghost_rw(0) == -1., ExcNonEqual(v1_ghost_rw(0), -1.));
  AssertThrow(v1_ghost_rw(2) == -1., ExcNonEqual(v1_ghost_rw(2), -1.));
  if (numproc > 2)
    {
      AssertThrow(v1_ghost_rw(8) == -1., ExcNonEqual(v1_ghost_rw(8), -1.));
      AssertThrow(v1_ghost_rw(10) == -1., ExcNonEqual(v1_ghost_rw(10), -1.));
    }
  MPI_Barrier(MPI_COMM_WORLD);
  if (myid == 0)
    deallog << "Ghost values after second swap OK" << std::endl;
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
