// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// test limited memory BFGS with quadratic function
// f(x) = 0.5 x*Lx - x*f
// f'(x) = Lx - f
// where L is 1d FD Laplacian
//
// We compare results to the companion Octave file which uses BFGS
// from optim package. The output of this test is made similar to the
// Octave output.


#include <deal.II/base/logstream.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/optimization/solver_bfgs.h>

#include <fstream>
#include <iostream>

#include "../tests.h"

template <typename number>
void
test()
{
  using VectorType = Vector<number>;

  // size of the problem
  const unsigned int N = 10;

  // parameters:
  const unsigned int itmax = 100;
  const double       tol   = 1e-8;

  // 1D Laplace with zero Dirichlet BC on both sides
  FullMatrix<number> L(N);
  L = 0.;

  for (unsigned int i = 0; i < N; ++i)
    {
      if (i > 0)
        L(i, i - 1) = -1.;
      L(i, i) = 2.;
      if (i < N - 1)
        L(i, i + 1) = -1.;
    }

  // L.print_formatted(deallog.get_file_stream(), 6, false, 10);

  // RHS
  VectorType b(N);
  for (unsigned int i = 0; i < N; ++i)
    b(i) = 1.;

  // solution
  VectorType x(N);
  x = 1.;

  // safety measure to not modify L or b within Lambda.
  const FullMatrix<number> &L_const = L;
  const VectorType         &b_const = b;
  const auto                func    = [&](const VectorType &x, VectorType &g) {
    L_const.vmult(g, x);
    number res = 0.5 * (g * x) - x * b_const;
    g.add(-1, b);
    return res;
  };

  // exact line minimization for quadratic function
  /*
  f(x) := a*x**2 + b*x + c;
  g(x) := ''(diff(f(x),x));
  sol : solve([f(0)=f0, f(1)=f1, g(0)=g0],[a,b,c]);
  sol2 : solve(diff(f(u),u)=0,u);
  subst(sol,sol2);

  u=-g0/(2*(-g0+f1-f0))

  */

  bool first_step = true;

  int        iteration = 0;
  VectorType dx(x);

  const auto line_min =
    [&](number &f, VectorType &x, VectorType &g, const VectorType &p) {
      deallog << "-------------------" << std::endl
              << "Line search " << iteration++ << ':' << std::endl;

      const number g_norm_sqr = g.norm_sqr();

      deallog << "Gradient:" << std::endl;
      g.print(deallog.get_file_stream(), 5, false);

      // directional derivative
      const number df = g * p;
      Assert(df < 0, ExcInternalError());
      // do the full step
      x.add(1., p);
      // save old value
      const number f_old = f;
      // calculate new value
      f = func(x, g);
      // get the step size
      const number denom = -df + f - f_old;
      Assert(denom != 0., ExcDivideByZero());
      Assert(denom > 0, ExcInternalError());
      const number step = -df * 0.5 / denom;
      // do the step
      x.add(step - 1., p);
      f = func(x, g);

      first_step = false;
      deallog << "Function value: " << f_old << " stepsize: " << step
              << " new value: " << f << std::endl;
      deallog << "Change:" << std::endl;
      dx = p;
      dx *= step;
      dx.print(deallog.get_file_stream(), 5, false);

      // finally return the step size
      return step;
    };

  SolverControl solver_control(itmax, tol, true);
  typename SolverBFGS<VectorType>::AdditionalData data(10, false);
  SolverBFGS<VectorType>                          solver(solver_control, data);
  solver.connect_line_search_slot(line_min);
  solver.solve(func, x);

  deallog << "Limited memory BFGS solution:" << std::endl;
  x.print(deallog.get_file_stream());
}

int
main()
{
  initlog();
  deallog << std::setprecision(5);

  test<double>();
}
