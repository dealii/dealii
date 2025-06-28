// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/function_tools.h>

#include <cmath>


DEAL_II_NAMESPACE_OPEN
namespace FunctionTools
{
  template <int dim>
  void
  taylor_estimate_function_bounds(
    const Function<dim>                        &function,
    const BoundingBox<dim>                     &box,
    std::pair<double, double>                  &value_bounds,
    std::array<std::pair<double, double>, dim> &gradient_bounds,
    const unsigned int                          component)
  {
    const Point<dim>     center   = box.center();
    const double         value    = function.value(center, component);
    const Tensor<1, dim> gradient = function.gradient(center, component);
    const SymmetricTensor<2, dim> hessian = function.hessian(center, component);

    // Deviation from function value at the center, based on the
    // Taylor-expansion: |f'| * dx + 1/2 * |f''| * dx^2,   (in 1d). dx is half
    // the side-length of the box.
    double taylor_bound_f = 0;

    for (unsigned int i = 0; i < dim; ++i)
      {
        const double dx_i = .5 * box.side_length(i);

        taylor_bound_f += std::abs(gradient[i]) * dx_i;

        // Deviation from value of df/dx_i at the center,
        // |f''| * dx,  (in 1d).
        double taylor_bound_dfdxi = 0;

        for (unsigned int j = 0; j < dim; ++j)
          {
            const double dx_j = .5 * box.side_length(j);

            taylor_bound_dfdxi += std::abs(hessian[i][j]) * dx_j;
            taylor_bound_f += .5 * std::abs(hessian[i][j]) * dx_i * dx_j;
          }

        gradient_bounds[i].first  = gradient[i] - taylor_bound_dfdxi;
        gradient_bounds[i].second = gradient[i] + taylor_bound_dfdxi;
      }

    value_bounds.first  = value - taylor_bound_f;
    value_bounds.second = value + taylor_bound_f;
  }

} // namespace FunctionTools

#include "base/function_tools.inst"

DEAL_II_NAMESPACE_CLOSE
