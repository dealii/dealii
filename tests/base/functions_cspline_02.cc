// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test CSpline value, gradient and hessian
// take example from
// https://www.gnu.org/software/gsl/manual/html_node/Interpolation-Example-programs.html


#include <deal.II/base/function_cspline.h>

#include "../tests.h"

template <int dim>
void
check()
{
  const unsigned int  n_points = 10;
  std::vector<double> x(n_points), y(n_points);
  for (unsigned int i = 0; i < n_points; ++i)
    {
      x[i] = i + 0.5 * std::sin(i);
      y[i] = i + std::cos(i * i);
    }


  // native:
  gsl_interp_accel *acc    = gsl_interp_accel_alloc();
  gsl_spline       *spline = gsl_spline_alloc(gsl_interp_cspline, n_points);

  gsl_spline_init(spline, &x[0], &y[0], n_points);

  // dealii:
  Functions::CSpline<dim> cspline(x, y);

  for (double xi = x[0]; xi <= x.back(); xi += 0.01)
    {
      const double f   = gsl_spline_eval(spline, xi, acc);
      const double df  = gsl_spline_eval_deriv(spline, xi, acc);
      const double ddf = gsl_spline_eval_deriv2(spline, xi, acc);

      const double         y   = cspline.value(Point<dim>(xi));
      const Tensor<1, dim> dy  = cspline.gradient(Point<dim>(xi));
      const double         ddy = cspline.laplacian(Point<dim>(xi));

      AssertThrow(std::fabs(f - y) <= std::fabs(f) * 1e-10, ExcInternalError());

      AssertThrow(std::fabs(df - dy[0]) <= std::fabs(df) * 1e-10,
                  ExcInternalError());

      AssertThrow(std::fabs(ddf - ddy) <= std::fabs(ddf) * 1e-10,
                  ExcInternalError());
    }
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);

  deallog << "OK" << std::endl;
}

int
main()
{
  initlog();

  check<1>();
}
