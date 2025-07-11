// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/function_cspline.h>
#include <deal.II/base/point.h>

#ifdef DEAL_II_WITH_GSL
#  include <gsl/gsl_spline.h>

#  include <algorithm>
#  include <cmath>


DEAL_II_NAMESPACE_OPEN
namespace Functions
{
  template <int dim>
  CSpline<dim>::CSpline(const std::vector<double> &x_,
                        const std::vector<double> &y_)
    : interpolation_points(x_)
    , interpolation_values(y_)
    , acc(gsl_interp_accel_alloc(),
          [](gsl_interp_accel *p) { gsl_interp_accel_free(p); })
    , cspline(gsl_spline_alloc(gsl_interp_cspline, interpolation_points.size()),
              [](gsl_spline *p) { gsl_spline_free(p); })
  {
    Assert(interpolation_points.size() > 0,
           ExcCSplineEmpty(interpolation_points.size()));

    Assert(interpolation_points.size() == interpolation_values.size(),
           ExcCSplineSizeMismatch(interpolation_points.size(),
                                  interpolation_values.size()));

    // check that input vector @p interpolation_points is provided in ascending order:
    for (unsigned int i = 0; i < interpolation_points.size() - 1; ++i)
      AssertThrow(interpolation_points[i] < interpolation_points[i + 1],
                  ExcCSplineOrder(i,
                                  interpolation_points[i],
                                  interpolation_points[i + 1]));

    const unsigned int n = interpolation_points.size();
    // gsl_spline_init returns something but it seems nobody knows what
    gsl_spline_init(cspline.get(),
                    interpolation_points.data(),
                    interpolation_values.data(),
                    n);
  }



  template <int dim>
  double
  CSpline<dim>::value(const Point<dim> &p, const unsigned int) const
  {
    // GSL functions may modify gsl_interp_accel *acc object (last argument).
    // This can only work in multithreaded applications if we lock the data
    // structures via a mutex.
    std::lock_guard<std::mutex> lock(acc_mutex);

    const double x = p[0];
    Assert(x >= interpolation_points.front() &&
             x <= interpolation_points.back(),
           ExcCSplineRange(x,
                           interpolation_points.front(),
                           interpolation_points.back()));

    return gsl_spline_eval(cspline.get(), x, acc.get());
  }



  template <int dim>
  Tensor<1, dim>
  CSpline<dim>::gradient(const Point<dim> &p, const unsigned int) const
  {
    // GSL functions may modify gsl_interp_accel *acc object (last argument).
    // This can only work in multithreaded applications if we lock the data
    // structures via a mutex.
    std::lock_guard<std::mutex> lock(acc_mutex);

    const double x = p[0];
    Assert(x >= interpolation_points.front() &&
             x <= interpolation_points.back(),
           ExcCSplineRange(x,
                           interpolation_points.front(),
                           interpolation_points.back()));

    const double   deriv = gsl_spline_eval_deriv(cspline.get(), x, acc.get());
    Tensor<1, dim> res;
    res[0] = deriv;
    return res;
  }



  template <int dim>
  double
  CSpline<dim>::laplacian(const Point<dim> &p, const unsigned int) const
  {
    // GSL functions may modify gsl_interp_accel *acc object (last argument).
    // This can only work in multithreaded applications if we lock the data
    // structures via a mutex.
    std::lock_guard<std::mutex> lock(acc_mutex);

    const double x = p[0];
    Assert(x >= interpolation_points.front() &&
             x <= interpolation_points.back(),
           ExcCSplineRange(x,
                           interpolation_points.front(),
                           interpolation_points.back()));

    return gsl_spline_eval_deriv2(cspline.get(), x, acc.get());
  }



  template <int dim>
  SymmetricTensor<2, dim>
  CSpline<dim>::hessian(const Point<dim> &p, const unsigned int) const
  {
    SymmetricTensor<2, dim> res;
    res[0][0] = laplacian(p);
    return res;
  }



  template <int dim>
  std::size_t
  CSpline<dim>::memory_consumption() const
  {
    // only simple data elements, so
    // use sizeof operator
    return sizeof(*this) + 2 * sizeof(double) * interpolation_values.size();
  }


  // explicit instantiations
  template class CSpline<1>;

} // namespace Functions

DEAL_II_NAMESPACE_CLOSE

#endif
