// ---------------------------------------------------------------------
//
// Copyright (C) 2013 - 2018 by the deal.II authors
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


// test eigenvalue approximation by GMRES algorithm

#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"

template <class NUMBER>
void
output_eigenvalues(const std::vector<NUMBER> &eigenvalues,
                   const std::string &        text)
{
  deallog << text;
  for (unsigned int j = 0; j < eigenvalues.size(); ++j)
    {
      deallog << ' ' << eigenvalues.at(j);
    }
  deallog << std::endl;
}

template <typename number>
void
test(unsigned int variant)
{
  const unsigned int n = variant < 3 ? 64 : 20;
  Vector<number>     rhs(n), sol(n);
  rhs = 1.;

  LAPACKFullMatrix<number> matrix(n, n);

  if (variant == 0)
    for (unsigned int i = 0; i < n; ++i)
      matrix(i, i) = std::sqrt(i + 1);
  else if (variant == 1)
    for (unsigned int i = 0; i < n; ++i)
      matrix(i, i) = 1.001 * (3 + i % 5);
  else if (variant == 2)
    for (unsigned int i = 0; i < n; ++i)
      matrix(i, i) = (i % 2 ? 1. : -1.) * std::sqrt(i + 1);
  else if (variant == 3)
    for (unsigned int i = 0; i < n; ++i)
      {
        matrix(i, i) = (i % 2 ? 1. : -1.) * std::sqrt(i + 1);
        if (i < n - 1)
          matrix(i, i + 1) = 1.5 + i;
        if (i < n - 2)
          matrix(i, i + 2) = -1.65;
        matrix(i, n - 1) = 2.;
        matrix(n - 1, i) = -2.;
      }
  else
    Assert(false, ExcMessage("Invalid variant"));
  if (std::is_same<number, float>::value == true)
    Assert(variant < 4, ExcMessage("Invalid_variant"));

  deallog.push(Utilities::int_to_string(variant, 1));

  // always do exactly 5 (variant 1 which should have 5 eigenvalues) or 20
  // iterations -> set a tolerance that is too tight otherwise
  IterationNumberControl control(variant == 1 ? 5 : 20, 1e-40);

  SolverGMRES<Vector<number>> solver(control);
  solver.connect_eigenvalues_slot(
    std::bind(output_eigenvalues<std::complex<double>>,
              std::placeholders::_1,
              "Eigenvalue estimate: "));
  solver.solve(matrix, sol, rhs, PreconditionIdentity());

  if (variant < 2)
    {
      IterationNumberControl   cg_control(variant == 1 ? 6 : 21, 1e-40);
      SolverCG<Vector<number>> solver_cg(cg_control);
      solver_cg.connect_eigenvalues_slot(std::bind(output_eigenvalues<double>,
                                                   std::placeholders::_1,
                                                   "Eigenvalue estimate: "));
      sol = 0;
      solver_cg.solve(matrix, sol, rhs, PreconditionIdentity());
    }

  if (variant == 3)
    {
      matrix.compute_eigenvalues();
      std::vector<std::complex<double>> eigenvalues(n);
      for (unsigned int i = 0; i < n; ++i)
        eigenvalues[i] = matrix.eigenvalue(i);

      std::sort(eigenvalues.begin(),
                eigenvalues.end(),
                internal::SolverGMRESImplementation::complex_less_pred);

      deallog << "Actual eigenvalues:        ";
      for (unsigned int i = 0; i < n; ++i)
        deallog << ' ' << eigenvalues[i];
      deallog << std::endl;
    }
  deallog.pop();
}

int
main()
{
  initlog();
  deallog << std::setprecision(10);

  deallog.push("double");
  test<double>(0);
  test<double>(1);
  test<double>(2);
  test<double>(3);
  deallog.pop();
}
