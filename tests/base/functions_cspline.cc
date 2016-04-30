// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// CSpline
// take example from https://www.gnu.org/software/gsl/manual/html_node/Interpolation-Example-programs.html

#include "../tests.h"
#include <cmath>
#include <deal.II/base/function_cspline.h>

template <int dim>
void check()
{
  const unsigned int n_points = 10;
  std::vector<double> x(n_points), y(n_points);
  for (unsigned int i = 0; i < n_points; i++)
    {
      x[i] = i + 0.5 * std::sin (i);
      y[i] = i + std::cos (i * i);
    }


  std::vector<double> y_native;
  // native version
  {
    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, n_points);

    gsl_spline_init (spline, &x[0], &y[0], n_points);

    for (double xi = x[0]; xi <= x.back(); xi += 0.01)
      {
        const double yi = gsl_spline_eval (spline, xi, acc);
        //deallog << xi << " " << yi << std::endl;
        y_native.push_back(yi);
      }
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
  }

  std::vector<double> y_dealii;
  {
    Functions::CSpline<dim> cspline(x,y);
    for (double xi = x[0]; xi <= x.back(); xi += 0.01)
      {
        const double yi = cspline.value(Point<dim>(xi));
        //deallog << xi << " " << yi << std::endl;
        y_dealii.push_back(yi);
      }
  }

  AssertThrow(std::equal(y_native.begin(), y_native.end(), y_dealii.begin()),
              ExcMessage("deal.II implementation of CSpline does not match native GSL."));

  deallog << "Ok"<< std::endl;
}

int main()
{
  std::string logname = "output";
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  check<1>();
}

