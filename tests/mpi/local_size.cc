// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// check Vector::locally_owned_size() for all supported vector types

#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/trilinos_epetra_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_tpetra_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"


template <typename VEC>
void
check_serial()
{
  const auto dofs_per_proc = 4;
  VEC        vec(dofs_per_proc);
  deallog << "type: " << Utilities::type_to_string(vec) << std::endl;
  deallog << "local size: " << vec.locally_owned_size() << std::endl;
  deallog << "size: " << vec.size() << std::endl;
}



template <typename VEC>
void
check_unghosted_parallel()
{
  const auto n_procs = dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  const auto my_rank = dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  const auto dofs_per_proc = 4;
  const auto n_dofs        = dofs_per_proc * n_procs;

  IndexSet   local_indices(n_dofs);
  const auto my_dofs_begin = dofs_per_proc * my_rank;
  const auto my_dofs_end   = dofs_per_proc * (my_rank + 1);
  local_indices.add_range(my_dofs_begin, my_dofs_end);
  local_indices.compress();

  VEC vec(local_indices, MPI_COMM_WORLD);
  deallog << "type: " << Utilities::type_to_string(vec) << std::endl;
  deallog << "index set size: " << local_indices.n_elements() << std::endl;
  deallog << "local size: " << vec.locally_owned_size() << std::endl;
  deallog << "size: " << vec.size() << std::endl;
}



template <typename VEC>
void
check_ghosted_parallel()
{
  const auto n_procs = dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  const auto my_rank = dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  const auto dofs_per_proc = 4;
  const auto n_dofs        = dofs_per_proc * n_procs;

  IndexSet   local_indices(n_dofs);
  IndexSet   ghost_indices(n_dofs);
  const auto my_dofs_begin = dofs_per_proc * my_rank;
  const auto my_dofs_end   = dofs_per_proc * (my_rank + 1);
  local_indices.add_range(my_dofs_begin, my_dofs_end);
  local_indices.compress();
  if (my_rank == 0)
    {
      ghost_indices.add_index(my_dofs_end);
    }
  else
    {
      ghost_indices.add_index(my_dofs_begin - 1);
      ghost_indices.add_index(my_dofs_end % n_dofs);
    }

  VEC vec(local_indices, ghost_indices, MPI_COMM_WORLD);
  deallog << "type: " << Utilities::type_to_string(vec) << std::endl;
  deallog << "index set size: " << local_indices.n_elements() << std::endl;
  deallog << "local size: " << vec.locally_owned_size() << std::endl;
  deallog << "size: " << vec.size() << std::endl;
}



template <typename VEC>
void
check_ghosted_parallel_block()
{
  const auto n_procs = dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  const auto my_rank = dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  const auto dofs_per_proc = 4;
  const auto n_dofs        = dofs_per_proc * n_procs;

  IndexSet   local_indices(n_dofs);
  IndexSet   ghost_indices(n_dofs);
  const auto my_dofs_begin = dofs_per_proc * my_rank;
  const auto my_dofs_end   = dofs_per_proc * (my_rank + 1);
  local_indices.add_range(my_dofs_begin, my_dofs_end);
  local_indices.compress();
  if (my_rank == 0)
    {
      ghost_indices.add_index(my_dofs_end);
    }
  else
    {
      ghost_indices.add_index(my_dofs_begin - 1);
      ghost_indices.add_index(my_dofs_end % n_dofs);
    }

  std::vector<IndexSet> local_blocks{local_indices, local_indices};
  // for variety do not ghost the second component
  std::vector<IndexSet> ghost_blocks{ghost_indices, IndexSet()};

  VEC vec(local_blocks, ghost_blocks, MPI_COMM_WORLD);
  deallog << "type: " << Utilities::type_to_string(vec) << std::endl;
  deallog << "local size: " << vec.locally_owned_size() << std::endl;
  deallog << "size: " << vec.size() << std::endl;
}



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  // non-block vectors:
  check_serial<Vector<double>>();

  check_unghosted_parallel<LinearAlgebra::EpetraWrappers::Vector>();
  check_unghosted_parallel<LinearAlgebra::TpetraWrappers::Vector<double>>();

  check_ghosted_parallel<LinearAlgebra::distributed::Vector<double>>();
  check_ghosted_parallel<PETScWrappers::MPI::Vector>();
  check_ghosted_parallel<TrilinosWrappers::MPI::Vector>();

  // block vectors:
  check_ghosted_parallel_block<
    LinearAlgebra::distributed::BlockVector<double>>();
  check_ghosted_parallel_block<PETScWrappers::MPI::BlockVector>();
  check_ghosted_parallel_block<TrilinosWrappers::MPI::BlockVector>();
}
