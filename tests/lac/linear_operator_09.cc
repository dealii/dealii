// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test that it is possible to instantiate a LinearOperator object for all
// different kinds of PETSc matrices and vectors
// TODO: A bit more tests...

#include <deal.II/base/index_set.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/linear_operator.h>

#include "../tests.h"

// Vectors:
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>

// Block Matrix and Vectors:
#include <deal.II/lac/petsc_block_sparse_matrix.h>
#include <deal.II/lac/petsc_block_vector.h>



int
main(int argc, char *argv[])
{
  using size_type = PETScWrappers::MPI::SparseMatrix::size_type;

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  const auto rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  initlog();
  deallog << std::setprecision(10);

  {
    const unsigned int np     = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
    const unsigned int n_dofs = 4;
    if (4 % np == 0 && np <= 4)
      {
        const auto dofs_per_processor = n_dofs / np;
        IndexSet   locally_owned_dofs(n_dofs);
        locally_owned_dofs.add_range(rank * dofs_per_processor,
                                     (rank + 1) * dofs_per_processor);
        locally_owned_dofs.compress();
        DynamicSparsityPattern dsp(n_dofs, n_dofs);
        for (const auto &index : locally_owned_dofs)
          dsp.add(index, index);

        PETScWrappers::MPI::SparseMatrix a;
        a.reinit(locally_owned_dofs, locally_owned_dofs, dsp, MPI_COMM_WORLD);
        for (const auto &i : locally_owned_dofs)
          for (const auto &j : locally_owned_dofs)
            a.add(i, i, 1);
        a.compress(VectorOperation::add);
        auto op_a = linear_operator<PETScWrappers::MPI::Vector>(a);

        PETScWrappers::MPI::Vector u, v;
        op_a.reinit_domain_vector(u, true);
        op_a.reinit_range_vector(v, true);
        for (auto i : u.locally_owned_elements())
          u[i] = 1;
        for (auto i : v.locally_owned_elements())
          v[i] = 1;
        u.compress(VectorOperation::insert);
        v.compress(VectorOperation::insert);

        op_a.vmult(v, u);
      }
    deallog << "SparseMatrix MPI -> OK" << std::endl;
  }

  {
    PETScWrappers::MPI::BlockSparseMatrix a;
    auto op_a = linear_operator<PETScWrappers::MPI::BlockVector>(a);
    deallog << "BlockSparseMatrix MPI -> OK" << std::endl;
  }

  deallog << "OK" << std::endl;

  return 0;
}
