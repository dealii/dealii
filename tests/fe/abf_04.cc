// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// ABF elements must reproduce vector polynomials in their anisotropic space,
// including their divergence, beyond the lowest order in three dimensions.

#include <deal.II/fe/fe_abf.h>

#include <deal.II/lac/vector.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include "../tests.h"

void
test(const unsigned int degree)
{
  FE_ABF<3>  fe(degree);
  const auto value = [degree](const Point<3> &p, const unsigned int component) {
    double result = component + 1.;
    for (unsigned int d = 0; d < 3; ++d)
      result *= std::pow(1. + p[d], degree + (d == component ? 2 : 0));
    return result;
  };
  const auto divergence = [degree](const Point<3> &p) {
    double result = 0.;
    for (unsigned int component = 0; component < 3; ++component)
      {
        double term = (component + 1.) * (degree + 2.);
        for (unsigned int d = 0; d < 3; ++d)
          term *= std::pow(1. + p[d], degree + (d == component ? 1 : 0));
        result += term;
      }
    return result;
  };

  const auto                 &points = fe.get_generalized_support_points();
  std::vector<Vector<double>> samples(points.size(), Vector<double>(3));
  for (unsigned int q = 0; q < points.size(); ++q)
    for (unsigned int c = 0; c < 3; ++c)
      samples[q][c] = value(points[q], c);
  std::vector<double> coefficients(fe.n_dofs_per_cell());
  fe.convert_generalized_support_point_values_to_dof_values(samples,
                                                            coefficients);

  double value_error      = 0.;
  double divergence_error = 0.;
  for (const Point<3> &p : {Point<3>(.13, .31, .72),
                            Point<3>(.89, .27, .43),
                            Point<3>(.42, .63, .19),
                            Point<3>(0., 0., 0.),
                            Point<3>(1., 1., 1.)})
    {
      double div = 0.;
      for (unsigned int c = 0; c < 3; ++c)
        {
          double interpolated = 0.;
          for (unsigned int i = 0; i < coefficients.size(); ++i)
            {
              interpolated +=
                coefficients[i] * fe.shape_value_component(i, p, c);
              div += coefficients[i] * fe.shape_grad_component(i, p, c)[c];
            }
          const double relative_error =
            std::abs(interpolated - value(p, c)) / value(p, c);
          AssertThrow(std::isfinite(relative_error), ExcInternalError());
          value_error = std::max(value_error, relative_error);
        }
      const double relative_error =
        std::abs(div - divergence(p)) / divergence(p);
      AssertThrow(std::isfinite(relative_error), ExcInternalError());
      divergence_error = std::max(divergence_error, relative_error);
    }
  deallog << "degree " << degree << ": relative value error " << value_error
          << ", relative divergence error " << divergence_error << std::endl;
  AssertThrow(value_error < 1.e-8, ExcInternalError());
  AssertThrow(divergence_error < 1.e-7, ExcInternalError());
}

int
main()
{
  initlog();
  test(0);
  test(1);
  test(2);
}
