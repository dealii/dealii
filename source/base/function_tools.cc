// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2019 by the deal.II authors
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

#include <deal.II/base/function_tools.h>

#include <cmath>


DEAL_II_NAMESPACE_OPEN
namespace FunctionTools
{
  template <int dim>
  void
  taylor_estimate_function_bounds(
    const Function<dim> &                       function,
    const BoundingBox<dim> &                    box,
    std::pair<double, double> &                 value_bounds,
    std::array<std::pair<double, double>, dim> &gradient_bounds,
    const unsigned int                          component)
  {
    const Point<dim>     center   = box.center();
    const double         value    = function.value(center, component);
    const Tensor<1, dim> gradient = function.gradient(center, component);
    const SymmetricTensor<2, dim> hessian = function.hessian(center, component);

    // Deviation from function value at the center, based on the
    // Taylor-expansion: |f'| * dx + 1/2 * |f''| * dx^2,   (in 1D). dx is half
    // the side-length of the box.
    double taylor_bound_f = 0;

    for (unsigned int i = 0; i < dim; ++i)
      {
        const double dx_i = .5 * box.side_length(i);

        taylor_bound_f += std::abs(gradient[i]) * dx_i;

        // Deviation from value of df/dx_i at the center,
        // |f''| * dx,  (in 1D).
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

#include "function_tools.inst"

DEAL_II_NAMESPACE_CLOSE
