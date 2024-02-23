// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test GSL by running CSpline. Copy-paste of /base/functions_cspline
// take example from
// https://www.gnu.org/software/gsl/manual/html_node/Interpolation-Example-programs.html

#include <deal.II/base/function_cspline.h>
#include <deal.II/base/logstream.h>

#include <fstream>
#include <string>
#include <vector>

using namespace dealii;

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

  std::vector<double>     y_dealii;
  Functions::CSpline<dim> cspline(x, y);
  for (double xi = x[0]; xi <= x.back(); xi += 0.01)
    {
      const double yi = cspline.value(Point<dim>(xi));
      y_dealii.push_back(yi);
    }

  deallog << "Ok" << std::endl;
}

int
main()
{
  check<1>();
}
