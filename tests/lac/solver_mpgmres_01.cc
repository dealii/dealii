// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check that the preconditioners are being applied cyclically

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_matrix.h>

#include <iostream>

#include "../tests.h"


// We create a wrapper around the preconditioner that will print out a name with
// each vmult so we can check that the preconditioners are being cycled
// correctly
template <typename Preconditioner>
class MyPreconditioner
{
public:
  MyPreconditioner(const std::string    &name,
                   const Preconditioner &preconditioner)
    : name(name)
    , preconditioner(preconditioner)
    , print_counter(0)
  {}

  void
  vmult(Vector<double> &dst, const Vector<double> &src) const
  {
    if (print_counter < 3)
      deallog.get_file_stream()
        << "Applying preconditioner " << name << std::endl;
    preconditioner.vmult(dst, src);
    ++print_counter;
  }

  void
  reset_counter()
  {
    print_counter = 0;
  }

private:
  const std::string     name;
  const Preconditioner &preconditioner;
  mutable unsigned int  print_counter;
};


int
main()
{
  initlog();

  const unsigned int N = 500;
  const unsigned int M = 2 * N;

  SparsityPattern sparsity_pattern(M, M, 2);
  for (unsigned int i = 0; i < N; ++i)
    sparsity_pattern.add(i, i + 1);
  sparsity_pattern.add(N, N - 1);
  sparsity_pattern.compress();
  SparseMatrix<double> matrix(sparsity_pattern);

  for (unsigned int i = 0; i < N; ++i)
    {
      matrix.diag_element(i) = 1.0 + 1 * i;
      matrix.set(i, i + 1, 1.);
    }
  matrix.set(N, N - 1, 1.);
  for (unsigned int i = N; i < M; ++i)
    matrix.diag_element(i) = 1.;

  PreconditionChebyshev<SparseMatrix<double>> prec_a;
  prec_a.initialize(matrix);
  MyPreconditioner wrapper_a("A", prec_a);

  PreconditionJacobi<SparseMatrix<double>> prec_b;
  prec_b.initialize(matrix);
  MyPreconditioner wrapper_b("B", prec_b);

  PreconditionSOR<SparseMatrix<double>> prec_c;
  prec_c.initialize(matrix);
  MyPreconditioner wrapper_c("C", prec_c);

  Vector<double> rhs(M);
  Vector<double> sol(M);
  rhs = 1.;

  SolverControl control(M, 1e-8);

  deallog.get_file_stream() << "FGMRES:\n";
  SolverFGMRES<Vector<double>> solver_fgmres(control);
  check_solver_within_range(
    solver_fgmres.solve(matrix, sol, rhs, wrapper_a, wrapper_b, wrapper_c),
    control.last_step(),
    7 - 6,
    7 + 6);

  deallog.get_file_stream() << "\nMPGMRES: (truncated)\n";
  sol = 0.;
  wrapper_a.reset_counter();
  wrapper_b.reset_counter();
  wrapper_c.reset_counter();

  SolverMPGMRES<Vector<double>>::AdditionalData data;
  data.use_truncated_mpgmres_strategy = true;
  SolverMPGMRES<Vector<double>> solver_trmpgmres(control, data);
  check_solver_within_range(
    solver_trmpgmres.solve(matrix, sol, rhs, wrapper_a, wrapper_b, wrapper_c),
    control.last_step(),
    19 - 6,
    19 + 6);


  deallog.get_file_stream() << "\nMPGMRES: (full)\n";
  sol = 0.;
  wrapper_a.reset_counter();
  wrapper_b.reset_counter();
  wrapper_c.reset_counter();
  data.use_truncated_mpgmres_strategy = false;
  SolverMPGMRES<Vector<double>> solver_mpgmres(control, data);

  check_solver_within_range(
    solver_mpgmres.solve(matrix, sol, rhs, wrapper_a, wrapper_b, wrapper_c),
    control.last_step(),
    20 - 12,
    20 + 12);
}
