// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_function_cspline_h
#define dealii_function_cspline_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_GSL
#  include <deal.II/base/function.h>
#  include <deal.II/base/mutex.h>
#  include <deal.II/base/point.h>

#  include <gsl/gsl_spline.h>

#  include <memory>

DEAL_II_NAMESPACE_OPEN

namespace Functions
{
  DeclException1(ExcCSplineEmpty,
                 int,
                 << "Interpolation points vector size can not be <" << arg1
                 << ">.");

  DeclException2(ExcCSplineSizeMismatch,
                 int,
                 int,
                 << "The size of interpolation points <" << arg1
                 << "> is different from the size of interpolation values <"
                 << arg2 << ">.");


  DeclException3(ExcCSplineOrder,
                 int,
                 double,
                 double,
                 << "The input interpolation points are not strictly ordered : "
                 << std::endl
                 << "x[" << arg1 << "] = " << arg2 << " >= x[" << (arg1 + 1)
                 << "] = " << arg3 << '.');

  DeclException3(
    ExcCSplineRange,
    double,
    double,
    double,
    << "Spline function can not be evaluated outside of the interpolation range: "
    << std::endl
    << arg1 << " is not in [" << arg2 << ';' << arg3 << "].");

  /**
   * The cubic spline function using GNU Scientific Library.
   * The resulting curve is piecewise cubic on each interval, with matching
   * first and second derivatives at the supplied data-points. The second
   * derivative is chosen to be zero at the first point and last point.
   *
   * @note This function is only implemented for dim==1 .
   *
   * @ingroup functions
   */
  template <int dim>
  class CSpline : public Function<dim>
  {
  public:
    /**
     * Constructor which should be provided with a set of points at which
     * interpolation is to be done @p interpolation_points and a set of function
     * values @p interpolation_values .
     */
    CSpline(const std::vector<double> &interpolation_points,
            const std::vector<double> &interpolation_values);

    virtual double
    value(const Point<dim>  &point,
          const unsigned int component = 0) const override;

    virtual Tensor<1, dim>
    gradient(const Point<dim>  &p,
             const unsigned int component = 0) const override;

    virtual SymmetricTensor<2, dim>
    hessian(const Point<dim>  &p,
            const unsigned int component = 0) const override;

    virtual double
    laplacian(const Point<dim>  &p,
              const unsigned int component = 0) const override;

    /**
     * Return an estimate for the memory consumption, in bytes, of this object.
     */
    virtual std::size_t
    memory_consumption() const override;

  private:
    /**
     * Points at which interpolation is provided
     */
    const std::vector<double> interpolation_points;

    /**
     * Values of the function at interpolation points
     */
    const std::vector<double> interpolation_values;

    /**
     * GSL accelerator for spline interpolation
     */
    std::unique_ptr<gsl_interp_accel, void (*)(gsl_interp_accel *)> acc;

    /**
     * GSL cubic spline interpolator
     */
    std::unique_ptr<gsl_spline, void (*)(gsl_spline *)> cspline;

    /**
     * A mutex for accelerator object.
     */
    mutable Threads::Mutex acc_mutex;
  };
} // namespace Functions

DEAL_II_NAMESPACE_CLOSE

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif

#endif
