//-----------------------------------------------------------
//
//    Copyright (C) 2018 by the deal.II authors
//
//    This file is part of the deal.II library.
//
//    The deal.II library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE.md at
//    the top level directory of deal.II.
//
//---------------------------------------------------------------

// test limited memory BFGS with Rosenbrock function.
// same as bfgs_05 but tests with default line search function which is the
// same asin bfgs_05.cc and therefore the number of iterations until convergence
// is exactly the same.

#include <deal.II/base/logstream.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/optimization/line_minimization.h>
#include <deal.II/optimization/solver_bfgs.h>

#include <fstream>
#include <iostream>

#include "../tests.h"

using namespace LineMinimization;

template <typename number>
void
test()
{
  auto &out = deallog.get_file_stream();
  out << std::setprecision(5) << std::fixed << std::right;

  typedef Vector<number> VectorType;

  // size of the problem
  const unsigned int N = 21;

  // parameters:
  const unsigned int itmax = 150;
  const double       gtol  = 1e-5; // gradient tolerance
  const unsigned int m_max = 3;

  // solution
  VectorType x(N), x_shifted(x), x0(x);
  x = 0.;

  // shift minimizer to this point
  VectorType location(x);
  for (unsigned int i = 0; i < N; ++i)
    location(i) = double(i) / (N - 1);

  // see
  // https://sourceforge.net/p/octave/optim/ci/default/tree/inst/rosenbrock.m#l26
  const auto rosenbrok = [&](VectorType &x, VectorType &g) {
    const unsigned int N   = x.size();
    double             res = 0.;
    g                      = 0;
    for (unsigned int i = 0; i < N; ++i)
      {
        const double xi2 = x(i) * x(i);

        if (i < N - 1)
          {
            res += 100. * dealii::Utilities::fixed_power<2>(x(i + 1) - xi2) +
                   dealii::Utilities::fixed_power<2>(1. - x(i));

            g(i) += -400. * x(i) * (x(i + 1) - xi2) - 2. * (1. - x(i));
          }

        if (i > 0)
          g(i) += 200. * (x(i) - x(i - 1) * x(i - 1));
      }
    return res;
  };



  unsigned int tot_fun_calls = 0;
  const auto   func          = [&](const VectorType &x, VectorType &g) {
    tot_fun_calls++;
    for (unsigned int i = 0; i < x.size(); ++i)
      x_shifted(i) = x(i) - location(i) + 1.;

    return rosenbrok(x_shifted, g);
  };

  const auto preconditioner = [](VectorType &                         g,
                                 const FiniteSizeHistory<VectorType> &s,
                                 const FiniteSizeHistory<VectorType> &y) {
    if (s.size() > 0)
      {
        // default preconditioning using the oldest {s,y} pair, see
        // lbfgs_recursion() in __bfgsmin.cc of "optim" Octave package.
        const unsigned int i  = s.size() - 1;
        const double       yy = y[i] * y[i];
        const double       sy = s[i] * y[i];
        Assert(yy > 0 && sy > 0, ExcInternalError());
        g *= sy / yy;
      }
  };

  SolverControl solver_control(itmax, gtol, false);
  typename SolverBFGS<VectorType>::AdditionalData data(m_max, false);
  SolverBFGS<VectorType>                          solver(solver_control, data);
  solver.connect_preconditioner_slot(preconditioner);
  solver.solve(func, x);

  deallog << "Limited memory BFGS solution:" << std::endl;
  x.print(deallog.get_file_stream());

  deallog << "Function value: " << func(x, x0) << std::endl;

  x.add(-1, location);
  deallog << "Linf error in solution: " << x.linfty_norm() << std::endl;

  deallog << "function calls: "
          << (tot_fun_calls - 1) /*one evaluation above*/ << std::endl;
}

int
main()
{
  initlog();
  deallog << std::setprecision(5);

  test<double>();
}
