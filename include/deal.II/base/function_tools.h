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

#ifndef dealii_function_tools_h
#define dealii_function_tools_h


#include <deal.II/base/config.h>

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/function.h>

DEAL_II_NAMESPACE_OPEN

namespace FunctionTools
{
  /**
   * Estimate bounds on the value and bounds on each gradient component of a
   * Function, $f$, over a BoundingBox, by approximating it by a 2nd order
   * Taylor polynomial starting from the box center.
   *
   * Each lower and upper bound is returned as a
   * <code>std::pair<double, double></code>, such that the first entry is the
   * lower bound, $L$, and the second is the upper bound, $U$, i.e.
   * $f(x) \in [L, U]$.
   *
   * The function value, gradient, and Hessian are computed at the box center.
   * The bounds on the value of the function are then estimated as
   *
   * $f(x) \in [f(x_c) - F, f(x_c) + F]$,
   * where
   * $F = \sum_i |\partial_i f(x_c)| h_i
   *    + 1/2 \sum_i \sum_j |\partial_i \partial_j f(x_c)| h_i h_j$.
   *
   * Here, $h_i$ is half the side length of the box in the $i$th coordinate
   * direction, which is the distance we extrapolate. The bounds on the gradient
   * components are estimated similarly as
   *
   * $\partial_i f \in [\partial_i f(x_c) - G_i, \partial_i f(x_c) + G_i]$,
   * where
   * $G_i = \sum_j |\partial_i \partial_j f(x_c)| h_j$.
   *
   * If the function has more than 1 component the @p component parameter can
   * be used to specify which function component the bounds should be computed
   * for.
   */
  template <int dim>
  void
  taylor_estimate_function_bounds(
    const Function<dim> &                       function,
    const BoundingBox<dim> &                    box,
    std::pair<double, double> &                 value_bounds,
    std::array<std::pair<double, double>, dim> &gradient_bounds,
    const unsigned int                          component = 0);

} // namespace FunctionTools
DEAL_II_NAMESPACE_CLOSE

#endif /* dealii_function_tools_h */
