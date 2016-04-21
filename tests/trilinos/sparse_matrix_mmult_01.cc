// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// check TrilinosWrappers::SparseMatrix::mmult for non-square matrices

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>

#include "../tests.h"

using namespace dealii;


template<typename MATRIX>
void out_matrix_size(const MATRIX &M,
                     const std::string &name)
{
  deallog << name << " = " << M.m() << " x " << M.n() << std::endl;
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv);

  mpi_initlog();

  TrilinosWrappers::SparsityPattern P_A;
  TrilinosWrappers::SparseMatrix A;
  P_A.reinit(4, 3);
  P_A.compress();
  A.reinit(P_A);

  TrilinosWrappers::SparsityPattern P_B;
  TrilinosWrappers::SparseMatrix B;
  P_B.reinit(3, 5);
  P_B.compress();
  B.reinit(P_B);

  out_matrix_size(A, "A");
  out_matrix_size(B, "B");

  TrilinosWrappers::SparseMatrix AB;
  A.mmult(AB, B);

  out_matrix_size(AB, "A.B");

  Assert (AB.m() == A.m(), ExcInternalError());
  Assert (AB.n() == B.n(), ExcInternalError());

  TrilinosWrappers::SparsityPattern P_C;
  TrilinosWrappers::SparseMatrix C;
  P_C.reinit(4, 5);
  P_C.compress();
  C.reinit(P_C);

  out_matrix_size(C, "C");

  TrilinosWrappers::SparseMatrix AtC;
  A.Tmmult(AtC, C);

  out_matrix_size(AtC, "At.C");

  Assert (AtC.m() == A.n(), ExcInternalError());
  Assert (AtC.n() == C.n(), ExcInternalError());
}
