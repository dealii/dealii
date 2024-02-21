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


// CSpline
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


  std::vector<double> y_native;
  // native version
  {
    gsl_interp_accel *acc    = gsl_interp_accel_alloc();
    gsl_spline       *spline = gsl_spline_alloc(gsl_interp_cspline, n_points);

    gsl_spline_init(spline, &x[0], &y[0], n_points);

    for (double xi = x[0]; xi <= x.back(); xi += 0.01)
      {
        const double yi = gsl_spline_eval(spline, xi, acc);
        // deallog << xi << ' ' << yi << std::endl;
        y_native.push_back(yi);
      }
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
  }

  std::vector<double> y_dealii;
  {
    Functions::CSpline<dim> cspline(x, y);
    for (double xi = x[0]; xi <= x.back(); xi += 0.01)
      {
        const double yi = cspline.value(Point<dim>(xi));
        // deallog << xi << ' ' << yi << std::endl;
        y_dealii.push_back(yi);
      }
  }

  AssertThrow(
    std::equal(y_native.begin(), y_native.end(), y_dealii.begin()),
    ExcMessage("deal.II implementation of CSpline does not match native GSL."));

  deallog << "Ok" << std::endl;
}

int
main()
{
  initlog();

  check<1>();
}
