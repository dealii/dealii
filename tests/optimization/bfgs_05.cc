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

// test limited memory BFGS with Rosenbrock function
// Octave with optim@1.5.2 converges with:
//
// 114 iterations
// function value: 6.74546e-13
// linf_norm =     1.4103e-06
//
// whereas we converge with
//
// 127 iterations
// function value: 1.3096e-13
// linf_norm =     1.6564e-07
// function calls: 133
//
// Scipy@1.1.0 converges with
//
// 130 iterations
// function value: 9.10208322706e-13
// linf_norm =     1.1038165888e-06
// Gradient norm:  1.72580704033e-05
// function calls: 143
//
// Note that both Octave and Scipy can not be matched exactly as
// there is no simple way to control all parameters and algorithms (i.e. line
// searches)


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
  VectorType x(N), x_shifted(x);
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

  VectorType x0(x), dx(x), g_old(x);

  bool   first_step = true;
  double f_prev     = 0.;

  const unsigned int print_n_iterations     = 5;
  unsigned int       iteration              = 0;
  unsigned int       line_search_iterations = 0;
  const auto         line_min               = [&](number &          f,
                            VectorType &      x,
                            VectorType &      g,
                            const VectorType &p) {
    if (iteration <= print_n_iterations)
      out << "------------------------------------------------" << std::endl
          << "Line search " << iteration << std::endl
          << std::endl;

    // save current solution value and gradient
    x0    = x;
    g_old = g;

    const double f0 = f;
    const double g0 = g * p;
    Assert(g0 < 0, ExcInternalError());

    // see scipy implementation
    // https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.line_search.html#scipy.optimize.line_search
    // and Eq. 2.6.8 in Fletcher 2013, Practical methods of optimization
    double df = f_prev - f;
    Assert(first_step || df >= 0., ExcInternalError());
    df = std::max(df, 100. * std::numeric_limits<double>::epsilon());
    const double a1 = (first_step ? 1. : std::min(1., -1.01 * 2. * df / g0));
    Assert(a1 > 0., ExcInternalError());
    f_prev = f;

    // 1D line-search function
    auto line_func = [&](const double &x_line) -> std::pair<double, double> {
      x = x0;
      x.add(x_line, p);
      const double f_line = func(x, g);
      const double g_line = g * p;
      f                   = f_line;

      return std::make_pair(f_line, g_line);
    };

    const auto res =
      line_search<double>(line_func, f0, g0, poly_fit<double>, a1, 0.9, 0.001);

    line_search_iterations += res.second;

    if (iteration <= print_n_iterations)
      {
        out << "function value: " << f0 << "  stepsize: " << res.first
            << std::endl
            << std::endl;

        // change:
        dx = p;
        dx *= res.first;

        const std::string s = "        ";
        for (unsigned int i = 0; i < N; ++i)
          out << s << std::setw(9) << x0(i) << s << std::setw(9) << g_old(i)
              << s << std::setw(9) << dx(i) << std::endl;
      }

    if (first_step)
      first_step = false;

    iteration++;
    return res.first;
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
  solver.connect_line_search_slot(line_min);
  solver.connect_preconditioner_slot(preconditioner);
  solver.solve(func, x);

  Assert(tot_fun_calls == line_search_iterations + 1, ExcInternalError());

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
