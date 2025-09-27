// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/utilities.h>

#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_fire.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>

#include "../tests.h"



// Test to verify correctness of SolverFIRE::solve()
// The objective function is the extended Rosenbrock function.
// The Rosenbrock function is a non-convex function used as a test problem
// for optimization algorithms introduced by Howard H. Rosenbrock.
//
// f(X) = f(x_0, x_1, ..., x_{N-1})
//
//      = \sum_{i=0}^{\frac{N}{2} -1}
//
//        \left[
//                  a ( x_{2i}^2 - x_{2i+1} )^2
//                  +
//                  b ( x_{2i}   - 1        )^2
//        \right],
//
//   where N is even and a = 100 and b = 1.
//
// DOI: 10.1007/BF02196600


using vector_t = typename dealii::Vector<double>;


double
compute(vector_t &G, const vector_t &X)
{
  AssertThrow(X.size() % 2 == 0, ExcInternalError());

  double value = 0.;

  // Value of the objective function.
  for (unsigned int i = 0; i < X.size() / 2; ++i)
    value += 100 * dealii::Utilities::fixed_power<2>(X(2 * i) * X(2 * i) -
                                                     X(2 * i + 1)) +
             dealii::Utilities::fixed_power<2>(X(2 * i) - 1);

  // Gradient of the objective function.
  for (unsigned int i = 0; i < X.size() / 2; ++i)
    {
      G(2 * i) = (X(2 * i) * X(2 * i) - X(2 * i + 1)) * X(2 * i) * 400 +
                 (X(2 * i) - 1) * 2;

      G(2 * i + 1) = (X(2 * i) * X(2 * i) - X(2 * i + 1)) * -200;
    }

  return value;
}



void
check_value(const unsigned int N, const double tol)
{
  AssertThrow(N % 2 == 0, ExcInternalError());

  vector_t X(N);

  // Use this to initialize DiagonalMatrix
  X = 1.;

  // Create inverse diagonal matrix.
  DiagonalMatrix<vector_t> inv_mass;
  inv_mass.reinit(X);

  // Set initial guess.
  for (unsigned int i = 0; i < N / 2; ++i)
    {
      X(2 * i)     = -1.2;
      X(2 * i + 1) = 1.0;
    }

  auto additional_data = SolverFIRE<vector_t>::AdditionalData(0.1, 1, 1);

  SolverControl solver_control(1e5, tol);

  SolverFIRE<vector_t> fire(solver_control, additional_data);

  fire.solve(compute, X, inv_mass);

  deallog << "FIRE::Solution vector: ";

  X.print(deallog.get_file_stream());
}



int
main()
{
  std::ofstream logfile("output");
  //  logfile.setf(std::ios::fixed);
  deallog << std::setprecision(4);
  deallog.attach(logfile);

  check_value(2, 1e-14);
  check_value(10, 1e-14);
  check_value(20, 1e-14);
}
