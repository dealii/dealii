// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// tests that GMRES builds an orthonormal basis properly for a few difficult
// test matrices for the delayed classical Gram-Schmidt method. As opposed to
// the gmres_reorthogonalize_xx tests, this test ensures that no
// reorthogonalization is performed for that variant (because it is done
// unconditionally) and that the Hessenberg matrix is highly accurate.

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"



template <typename number>
void
test(unsigned int variant, unsigned int min_convergence_steps)
{
  const unsigned int n = 64;
  Vector<number>     rhs(n), sol(n);
  rhs = 1.;

  FullMatrix<number> matrix(n, n);
  for (unsigned int i = 0; i < n; ++i)
    for (unsigned int j = 0; j < n; ++j)
      matrix(i, j) = random_value<double>(-.1, .1);

  // put diagonal entries of different strengths. these are very challenging
  // for GMRES and will usually take a lot of iterations until the Krylov
  // subspace is complete enough
  if (variant == 0)
    for (unsigned int i = 0; i < n; ++i)
      matrix(i, i) = (i + 1);
  else if (variant == 1)
    for (unsigned int i = 0; i < n; ++i)
      matrix(i, i) = (i + 1) * (i + 1) * (i + 1) * (i + 1);
  else if (variant == 2)
    for (unsigned int i = 0; i < n; ++i)
      matrix(i, i) = 1e10 * (i + 1);
  else if (variant == 3)
    for (unsigned int i = 0; i < n; ++i)
      matrix(i, i) = 1e10 * (i + 1) * (i + 1) * (i + 1) * (i + 1);
  else if (variant == 4)
    for (unsigned int i = 0; i < n; ++i)
      matrix(i, i) = 1e30 * (i + 1);
  else if (variant == 5)
    for (unsigned int i = 0; i < n; ++i)
      matrix(i, i) = 1e30 * (i + 1) * (i + 1) * (i + 1) * (i + 1);
  else
    Assert(false, ExcMessage("Invalid variant"));
  if (std::is_same_v<number, float> == true)
    Assert(variant < 4, ExcMessage("Invalid_variant"));

  deallog.push(Utilities::int_to_string(variant, 1));

  SolverControl control(1000, 5e2 * std::numeric_limits<number>::epsilon());
  typename SolverGMRES<Vector<number>>::AdditionalData data;
  data.max_basis_size = min_convergence_steps + 2;
  data.orthogonalization_strategy =
    LinearAlgebra::OrthogonalizationStrategy::delayed_classical_gram_schmidt;

  SolverGMRES<Vector<number>> solver(control, data);
  auto print_re_orthogonalization = [](int accumulated_iterations) {
    deallog.get_file_stream() << "Re-orthogonalization enabled at step "
                              << accumulated_iterations << std::endl;
  };

  auto print_hessenberg_matrix = [](const FullMatrix<double> &H) {
    deallog.get_file_stream() << "Scaled Hessenberg matrix" << std::endl;
    // Print the first 30 entries in the Hessenberg matrix
    for (unsigned int i = 0; i < 30; ++i)
      {
        for (unsigned int j = 0; j < 30; ++j)
          {
            if (std::abs(H(i, j) / H(0, 0)) >
                (std::is_same_v<number, float> ? 1e-6 : 1e-12))
              deallog.get_file_stream() << H(i, j) / H(0, 0) << " ";
            else
              deallog.get_file_stream() << "0 ";
          }
        deallog.get_file_stream() << std::endl;
      }
  };

  auto print_orthogonality =
    [](const internal::SolverGMRESImplementation::TmpVectors<Vector<number>>
         &tmp_vector) {
      deallog.get_file_stream() << "Orthogonality quality" << std::endl;
      for (unsigned int i = 0; i < tmp_vector.size() - 3; ++i)
        {
          if (const_cast<internal::SolverGMRESImplementation::TmpVectors<
                Vector<number>> &>(tmp_vector)(i, {})
                .size() > 0)
            {
              for (unsigned int j = 0; j <= i; ++j)
                {
                  const double product = tmp_vector[i] * tmp_vector[j];
                  if (std::abs(product) >
                      (std::is_same_v<number, float> ? 5e-7 : 1e-12))
                    deallog.get_file_stream() << product << " ";
                  else
                    deallog.get_file_stream() << "0 ";
                }
              deallog.get_file_stream() << std::endl;
            }
        }
    };

  if (std::is_same_v<number, double> && variant < 4)
    solver.connect_hessenberg_slot(print_hessenberg_matrix);
  solver.connect_krylov_space_slot(print_orthogonality);

  check_solver_within_range(
    solver.solve(matrix, sol, rhs, PreconditionIdentity()),
    control.last_step(),
    min_convergence_steps,
    min_convergence_steps + 2);

  deallog.pop();
}

int
main()
{
  initlog();
  deallog << std::setprecision(10);

  deallog.push("double");
  test<double>(0, 56);
  test<double>(1, 64);
  test<double>(2, 56);
  test<double>(3, 64);
  test<double>(4, 56);
  test<double>(5, 64);
  deallog.pop();
  deallog.push("float");
  test<float>(0, 34);
  test<float>(1, 64);
  test<float>(2, 34);
  test<float>(3, 64);
  deallog.pop();
}
