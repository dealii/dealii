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

#include <deal.II/base/point.h>
#include <deal.II/base/function_cspline.h>

#ifdef DEAL_II_WITH_GSL
#include <cmath>
#include <algorithm>

DEAL_II_NAMESPACE_OPEN
namespace Functions
{

  template <int dim>
  CSpline<dim>::CSpline(const std::vector<double> &x_,
                        const std::vector<double> &y_)
    :
    interpolation_points(x_),
    interpolation_values(y_)
  {
    Assert (interpolation_points.size() > 0,
            ExcCSplineEmpty(interpolation_points.size()));

    Assert (interpolation_points.size() == interpolation_values.size(),
            ExcCSplineSizeMismatch(interpolation_points.size(),interpolation_values.size()));

    // check that input vector @p interpolation_points is provided in ascending order:
    for (unsigned int i = 0; i < interpolation_points.size()-1; i++)
      AssertThrow(interpolation_points[i]<interpolation_points[i+1],
                  ExcCSplineOrder(i,interpolation_points[i],interpolation_points[i+1]));

    acc = gsl_interp_accel_alloc ();
    const unsigned int n = interpolation_points.size();
    cspline = gsl_spline_alloc (gsl_interp_cspline, n);
    // gsl_spline_init returns something but it seems nobody knows what
    gsl_spline_init (cspline, &interpolation_points[0], &interpolation_values[0], n);
  }

  template <int dim>
  CSpline<dim>::~CSpline()
  {
    gsl_interp_accel_free (acc);
    gsl_spline_free (cspline);
    acc = NULL;
    cspline = NULL;
  }


  template <int dim>
  double
  CSpline<dim>::value (const Point<dim>   &p,
                       const unsigned int) const
  {
    const double &x = p[0];
    Assert (x >= interpolation_points.front() && x <= interpolation_points.back(),
            ExcCSplineRange(x,interpolation_points.front(),interpolation_points.back()));

    return gsl_spline_eval (cspline, x, acc);
  }


  template <int dim>
  Tensor<1,dim>
  CSpline<dim>::gradient (const Point<dim>   &p,
                          const unsigned int) const
  {
    const double &x = p[0];
    Assert (x >= interpolation_points.front() && x <= interpolation_points.back(),
            ExcCSplineRange(x,interpolation_points.front(),interpolation_points.back()));

    const double deriv = gsl_spline_eval_deriv (cspline, x, acc);
    Tensor<1,dim> res;
    res[0] = deriv;
    return res;
  }


  template <int dim>
  std::size_t
  CSpline<dim>::memory_consumption () const
  {
    // only simple data elements, so
    // use sizeof operator
    return sizeof (*this) + 2*sizeof(double)*interpolation_values.size();
  }


// explicit instantiations
  template class CSpline<1>;

}

DEAL_II_NAMESPACE_CLOSE

#endif
