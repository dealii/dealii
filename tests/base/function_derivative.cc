// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2007 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check automatic differentiation

#include <deal.II/base/auto_derivative_function.h>
#include <deal.II/base/function_derivative.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/lac/vector.h>

#include <string>
#include <vector>

#include "../tests.h"


template <int dim>
void
check_derivative_order(const std::vector<Tensor<1, dim>> &gradients,
                       FunctionDerivative<dim>           &df,
                       const Quadrature<dim>             &quadrature,
                       const unsigned int                 direction,
                       const double                       order)
{
  std::vector<double> derivatives(quadrature.size());
  std::vector<double> differences(quadrature.size());

  // Compute derivatives with one
  // step size and store errors
  df.set_h(1.e-2);
  df.value_list(quadrature.get_points(), derivatives);
  for (unsigned int i = 0; i < gradients.size(); ++i)
    differences[i] = gradients[i][direction] - derivatives[i];

  df.set_h(5.e-3);
  df.value_list(quadrature.get_points(), derivatives);

  // Expected reduction
  const double expected = std::pow(.5, order);

  for (unsigned int i = 0; i < gradients.size(); ++i)
    {
      const double reduction =
        std::fabs(gradients[i][direction] - derivatives[i]) /
        std::fabs(differences[i]);
      if (reduction > 1.2 * expected || reduction < .8 * expected)
        deallog << "Derivative error " << direction << ' ' << order << ' ' << i
                << "   " << reduction << std::endl;
    }
}


template <int dim>
void
check_hessian_order(const std::vector<double> &values,
                    FunctionDerivative<dim>   &df,
                    const Quadrature<dim>     &quadrature,
                    const Point<dim>          &k,
                    const unsigned int         direction,
                    const double               order)
{
  std::vector<Tensor<1, dim>> derivatives(quadrature.size());
  std::vector<Tensor<1, dim>> differences(quadrature.size());

  const double h = (order < 3) ? 1.e-2 : 1.e-1;
  // Compute derivatives with one
  // step size and store errors
  df.set_h(h);
  df.gradient_list(quadrature.get_points(), derivatives);
  for (unsigned int i = 0; i < values.size(); ++i)
    for (unsigned int d = 0; d < dim; ++d)
      // Use what we know about
      // derivatives of the sine wave
      differences[i][d] = -values[i] * k[direction] * k[d] - derivatives[i][d];

  df.set_h(h / 2.);
  df.gradient_list(quadrature.get_points(), derivatives);

  // Expected reduction
  const double expected = std::pow(.5, order);

  for (unsigned int i = 0; i < values.size(); ++i)
    for (unsigned int d = 0; d < dim; ++d)
      {
        const double reduction =
          std::fabs(-values[i] * k[direction] * k[d] - derivatives[i][d]) /
          std::fabs(differences[i][d]);
        if (reduction > 1.2 * expected || reduction < .8 * expected)
          deallog << "Hessian error " << direction << ' ' << d << ' ' << order
                  << ' ' << i << "   " << reduction << "   " << expected
                  << "   " << differences[i][d] << std::endl;
      }
}


template <int dim>
void
check_sine(unsigned int nquad)
{
  QGauss<dim> quadrature(nquad);

  Point<dim> wave_vector;
  for (unsigned int d = 0; d < dim; ++d)
    wave_vector[d] = d + 2.;

  Functions::FourierSineFunction<dim> f(wave_vector);

  std::vector<double>         values(quadrature.size());
  std::vector<Tensor<1, dim>> gradients(quadrature.size());

  f.value_list(quadrature.get_points(), values);
  f.gradient_list(quadrature.get_points(), gradients);

  // Check derivatives in all directions
  for (unsigned int d = 0; d < dim; ++d)
    {
      deallog << "Direction " << d << std::endl;
      Point<dim> dir;
      dir[d] = 1.;
      deallog.push("Euler");
      FunctionDerivative<dim> df(f, dir, 1.e-4);
      check_derivative_order(gradients, df, quadrature, d, 2);
      check_hessian_order(values, df, quadrature, wave_vector, d, 2);
      deallog.pop();
      deallog.push("UpwindEuler");
      df.set_formula(FunctionDerivative<dim>::UpwindEuler);
      check_derivative_order(gradients, df, quadrature, d, 1);
      check_hessian_order(values, df, quadrature, wave_vector, d, 1);
      deallog.pop();
      deallog.push("FourthOrder");
      df.set_formula(FunctionDerivative<dim>::FourthOrder);
      check_derivative_order(gradients, df, quadrature, d, 4);
      check_hessian_order(values, df, quadrature, wave_vector, d, 4);
      deallog.pop();
    }
}


int
main()
{
  initlog();

  check_sine<2>(3);
  check_sine<3>(3);
}
