// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2020 by the deal.II authors
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
  for (unsigned int i = 0; i < n_points; i++)
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
