// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
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

// Test that it is possible to instantiate a LinearOperator object for all
// different kinds of PETSc matrices and vectors
// TODO: A bit more tests...

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
  typedef PETScWrappers::MPI::SparseMatrix::size_type size_type;

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();
  deallog << std::setprecision(10);

  {
    unsigned int np = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
    if (4 % np == 0 && np <= 4)
      {
        PETScWrappers::MPI::SparseMatrix a(
          MPI_COMM_WORLD, 4, 4, 4 / np, 4 / np, 1);
        for (unsigned int i = 0; i < 4; ++i)
          for (unsigned int j = 0; j < 4; ++j)
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
