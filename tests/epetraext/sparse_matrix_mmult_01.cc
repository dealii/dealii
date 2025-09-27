// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check TrilinosWrappers::SparseMatrix::mmult for non-square matrices

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>

#include "../tests.h"



template <typename MATRIX>
void
out_matrix_size(const MATRIX &M, const std::string &name)
{
  deallog << name << " = " << M.m() << " x " << M.n() << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);

  initlog();

  TrilinosWrappers::SparsityPattern P_A;
  TrilinosWrappers::SparseMatrix    A;
  P_A.reinit(4, 3);
  P_A.compress();
  A.reinit(P_A);

  TrilinosWrappers::SparsityPattern P_B;
  TrilinosWrappers::SparseMatrix    B;
  P_B.reinit(3, 5);
  P_B.compress();
  B.reinit(P_B);

  out_matrix_size(A, "A");
  out_matrix_size(B, "B");

  TrilinosWrappers::SparseMatrix AB;
  A.mmult(AB, B);

  out_matrix_size(AB, "A.B");

  Assert(AB.m() == A.m(), ExcInternalError());
  Assert(AB.n() == B.n(), ExcInternalError());

  TrilinosWrappers::SparsityPattern P_C;
  TrilinosWrappers::SparseMatrix    C;
  P_C.reinit(4, 5);
  P_C.compress();
  C.reinit(P_C);

  out_matrix_size(C, "C");

  TrilinosWrappers::SparseMatrix AtC;
  A.Tmmult(AtC, C);

  out_matrix_size(AtC, "At.C");

  Assert(AtC.m() == A.n(), ExcInternalError());
  Assert(AtC.n() == C.n(), ExcInternalError());
}
