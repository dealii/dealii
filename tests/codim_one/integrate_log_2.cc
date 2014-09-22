// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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



// integrates the function x^M ln(|x-x0|/alpha) on different
// situations, to test the QGaussLogR quadrature formula

#include "../tests.h"
#include <fstream>
#include <deal.II/base/logstream.h>

// all include files needed for the program

#include <deal.II/base/exceptions.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>

#include <fstream>
#include <string>

#include <math.h>

using namespace std;

std::ofstream logfile("output");


// Returns the following integral: /int_0^1 x^N * ln(|x-point|/alpha) dx
double log_integral(const unsigned int N, const double point, const double alpha);

int main()
{
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog<<std::fixed;
  deallog<<std::setprecision(10);


  vector<Point<1> > origins;
  origins.push_back(Point<1>());
  origins.push_back(Point<1>(1.));
  origins.push_back(Point<1>(1./3.));
  origins.push_back(Point<1>(1./2.));
  origins.push_back(Point<1>(1./5.));
  origins.push_back(Point<1>(3./5.));

  vector<double> alphas;
  alphas.push_back(1.);
  alphas.push_back(1./3.);
  alphas.push_back(3./2.);


  double exact_integral;
  double eps = 1e-10;

  for (unsigned int nos = 0; nos< origins.size(); ++nos)
    {
      for (unsigned int nas = 0; nas < alphas.size(); ++nas)
        {
          for (unsigned int power = 0; power<13; ++power)
            {
              exact_integral = log_integral(power, origins[nos][0], alphas[nas]);
              deallog << "==================================================" << endl;
              deallog << "int_0^1 x^" << power
                      << " ln(|x-" << origins[nos] << "|/"
                      << alphas[nas] << ") = "
                      << exact_integral << endl;
              for (unsigned int nq = 1; nq < 13; ++nq)
                {
                  QGaussLogR<1> quad(nq, origins[nos], alphas[nas]);
                  QGaussLogR<1> factored_quad(nq, origins[nos], alphas[nas], true);

                  double approx_integral = 0;
                  double approx_integral_factored = 0;
                  for (unsigned int q=0; q<quad.size(); ++q)
                    {
                      double qpow = pow(quad.point(q)[0], (double) power);
                      approx_integral += qpow  * quad.weight(q);
                      qpow *= log(abs( (factored_quad.point(q)-origins[nos])[0] )/alphas[nas] );
                      approx_integral_factored += qpow * factored_quad.weight(q);
                    }

                  deallog << "   Error(n=" << quad.size()
                          << ") = " << (approx_integral - exact_integral);
                  if ( abs(approx_integral - approx_integral_factored) > eps )
                    deallog << ", difference between factored and unfactored: "
                            << abs(approx_integral - approx_integral_factored);
                  deallog << endl;
                }
            }
          deallog << "==================================================" << endl;
        }
    }
}



double log_integral(const unsigned int N, const double point, const double alpha)
{
  double result = 0;
  if (point == 0)
    {
      result = - (1.+log(alpha)*N +log(alpha)) / (2 * N + N*N + 1);
    }
  else if (point == 1)
    {
      switch (N)
        {
        case 0:
          result = -0.1e1 - log(alpha);
          break;
        case 1:
          result = -log(alpha) / 0.2e1 - 0.3e1 / 0.4e1;
          break;
        case 2:
          result = -log(alpha) / 0.3e1 - 0.11e2 / 0.18e2;
          break;
        case 3:
          result = -log(alpha) / 0.4e1 - 0.25e2 / 0.48e2;
          break;
        case 4:
          result = -0.137e3 / 0.300e3 - log(alpha) / 0.5e1;
          break;
        case 5:
          result = -log(alpha) / 0.6e1 - 0.49e2 / 0.120e3;
          break;
        case 6:
          result = -log(alpha) / 0.7e1 - 0.363e3 / 0.980e3;
          break;
        case 7:
          result = -log(alpha) / 0.8e1 - 0.761e3 / 0.2240e4;
          break;
        case 8:
          result = -0.7129e4 / 0.22680e5 - log(alpha) / 0.9e1;
          break;
        case 9:
          result = -log(alpha) / 0.10e2 - 0.7381e4 / 0.25200e5;
          break;
        case 10:
          result = -log(alpha) / 0.11e2 - 0.83711e5 / 0.304920e6;
          break;
        case 11:
          result = -0.86021e5 / 0.332640e6 - log(alpha) / 0.12e2;
          break;
        case 12:
          result = -log(alpha) / 0.13e2 - 0.1145993e7 / 0.4684680e7;
          break;
        }
    }
  else
    {
      switch (N)
        {
        case 0:
          result = point * log(point) + log(0.1e1 - point) - log(alpha) - 0.1e1 - point * log(0.1e1 - point);
          break;
        case 1:
          result = pow(point, 0.2e1) * log(point) / 0.2e1 - pow(point, 0.2e1) * log(0.1e1 - point) / 0.2e1 + log(0.1e1 - point) / 0.2e1 - log(alpha) / 0.2e1 - 0.1e1 / 0.4e1 - point / 0.2e1;
          break;
        case 2:
          result = pow(point, 0.3e1) * log(point) / 0.3e1 - point / 0.6e1 - pow(point, 0.3e1) * log(0.1e1 - point) / 0.3e1 - pow(point, 0.2e1) / 0.3e1 + log(0.1e1 - point) / 0.3e1 - log(alpha) / 0.3e1 - 0.1e1 / 0.9e1;
          break;
        case 3:
          result = pow(point, 0.4e1) * log(point) / 0.4e1 + log(0.1e1 - point) / 0.4e1 - log(alpha) / 0.4e1 - pow(point, 0.4e1) * log(0.1e1 - point) / 0.4e1 - point / 0.12e2 - pow(point, 0.3e1) / 0.4e1 - pow(point, 0.2e1) / 0.8e1 - 0.1e1 / 0.16e2;
          break;
        case 4:
          result = pow(point, 0.5e1) * log(point) / 0.5e1 - 0.1e1 / 0.25e2 + log(0.1e1 - point) / 0.5e1 - log(alpha) / 0.5e1 - pow(point, 0.4e1) / 0.5e1 - point / 0.20e2 - pow(point, 0.5e1) * log(0.1e1 - point) / 0.5e1 - pow(point, 0.2e1) / 0.15e2 - pow(point, 0.3e1) / 0.10e2;
          break;
        case 5:
          result = pow(point, 0.6e1) * log(point) / 0.6e1 - pow(point, 0.3e1) / 0.18e2 - pow(point, 0.5e1) / 0.6e1 - point / 0.30e2 - pow(point, 0.2e1) / 0.24e2 - pow(point, 0.6e1) * log(0.1e1 - point) / 0.6e1 + log(0.1e1 - point) / 0.6e1 - log(alpha) / 0.6e1 - 0.1e1 / 0.36e2 - pow(point, 0.4e1) / 0.12e2;
          break;
        case 6:
          result = pow(point, 0.7e1) * log(point) / 0.7e1 - pow(point, 0.4e1) / 0.21e2 - pow(point, 0.5e1) / 0.14e2 - point / 0.42e2 - pow(point, 0.2e1) / 0.35e2 - pow(point, 0.7e1) * log(0.1e1 - point) / 0.7e1 - pow(point, 0.3e1) / 0.28e2 - 0.1e1 / 0.49e2 - pow(point, 0.6e1) / 0.7e1 + log(0.1e1 - point) / 0.7e1 - log(alpha) / 0.7e1;
          break;
        case 7:
          result = pow(point, 0.8e1) * log(point) / 0.8e1 - 0.1e1 / 0.64e2 + log(0.1e1 - point) / 0.8e1 - log(alpha) / 0.8e1 - pow(point, 0.6e1) / 0.16e2 - pow(point, 0.8e1) * log(0.1e1 - point) / 0.8e1 - pow(point, 0.5e1) / 0.24e2 - pow(point, 0.3e1) / 0.40e2 - pow(point, 0.7e1) / 0.8e1 - pow(point, 0.2e1) / 0.48e2 - point / 0.56e2 - pow(point, 0.4e1) / 0.32e2;
          break;
        case 8:
          result = pow(point, 0.9e1) * log(point) / 0.9e1 - pow(point, 0.4e1) / 0.45e2 - pow(point, 0.9e1) * log(0.1e1 - point) / 0.9e1 - pow(point, 0.6e1) / 0.27e2 - pow(point, 0.2e1) / 0.63e2 + log(0.1e1 - point) / 0.9e1 - log(alpha) / 0.9e1 - pow(point, 0.8e1) / 0.9e1 - 0.1e1 / 0.81e2 - pow(point, 0.5e1) / 0.36e2 - point / 0.72e2 - pow(point, 0.3e1) / 0.54e2 - pow(point, 0.7e1) / 0.18e2;
          break;
        case 9:
          result = pow(point, 0.10e2) * log(point) / 0.10e2 + log(0.1e1 - point) / 0.10e2 - log(alpha) / 0.10e2 - pow(point, 0.2e1) / 0.80e2 - 0.1e1 / 0.100e3 - point / 0.90e2 - pow(point, 0.4e1) / 0.60e2 - pow(point, 0.9e1) / 0.10e2 - pow(point, 0.8e1) / 0.20e2 - pow(point, 0.6e1) / 0.40e2 - pow(point, 0.5e1) / 0.50e2 - pow(point, 0.3e1) / 0.70e2 - pow(point, 0.10e2) * log(0.1e1 - point) / 0.10e2 - pow(point, 0.7e1) / 0.30e2;
          break;
        case 10:
          result = pow(point, 0.11e2) * log(point) / 0.11e2 + log(0.1e1 - point) / 0.11e2 - log(alpha) / 0.11e2 - pow(point, 0.2e1) / 0.99e2 - pow(point, 0.11e2) * log(0.1e1 - point) / 0.11e2 - 0.1e1 / 0.121e3 - pow(point, 0.3e1) / 0.88e2 - pow(point, 0.4e1) / 0.77e2 - pow(point, 0.8e1) / 0.33e2 - pow(point, 0.5e1) / 0.66e2 - pow(point, 0.10e2) / 0.11e2 - pow(point, 0.7e1) / 0.44e2 - pow(point, 0.9e1) / 0.22e2 - point / 0.110e3 - pow(point, 0.6e1) / 0.55e2;
          break;
        case 11:
          result = pow(point, 0.12e2) * log(point) / 0.12e2 - pow(point, 0.8e1) / 0.48e2 - pow(point, 0.2e1) / 0.120e3 - pow(point, 0.11e2) / 0.12e2 - pow(point, 0.4e1) / 0.96e2 - pow(point, 0.10e2) / 0.24e2 - pow(point, 0.5e1) / 0.84e2 - 0.1e1 / 0.144e3 - pow(point, 0.3e1) / 0.108e3 - point / 0.132e3 - pow(point, 0.9e1) / 0.36e2 + log(0.1e1 - point) / 0.12e2 - log(alpha) / 0.12e2 - pow(point, 0.7e1) / 0.60e2 - pow(point, 0.6e1) / 0.72e2 - pow(point, 0.12e2) * log(0.1e1 - point) / 0.12e2;
          break;
        case 12:
          result = pow(point, 0.13e2) * log(point) / 0.13e2 - pow(point, 0.4e1) / 0.117e3 - pow(point, 0.7e1) / 0.78e2 - pow(point, 0.5e1) / 0.104e3 - pow(point, 0.11e2) / 0.26e2 - pow(point, 0.13e2) * log(0.1e1 - point) / 0.13e2 - pow(point, 0.3e1) / 0.130e3 - pow(point, 0.6e1) / 0.91e2 - point / 0.156e3 - pow(point, 0.2e1) / 0.143e3 - pow(point, 0.10e2) / 0.39e2 + log(0.1e1 - point) / 0.13e2 - log(alpha) / 0.13e2 - 0.1e1 / 0.169e3 - pow(point, 0.9e1) / 0.52e2 - pow(point, 0.12e2) / 0.13e2 - pow(point, 0.8e1) / 0.65e2;
          break;
        }
    }
  return result;
}
