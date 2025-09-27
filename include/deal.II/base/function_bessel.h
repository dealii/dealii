// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_function_bessel_h
#define dealii_function_bessel_h


#include <deal.II/base/config.h>

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>

#include <vector>


DEAL_II_NAMESPACE_OPEN

namespace Functions
{
  /**
   * The Bessel functions of first kind or positive integer order.
   *
   * @ingroup functions
   */
  template <int dim>
  class Bessel1 : public Function<dim>
  {
  public:
    /**
     * Constructor. @p wave_number must be nonnegative.
     */
    Bessel1(const unsigned int order,
            const double       wave_number,
            const Point<dim>   center = Point<dim>());

    virtual double
    value(const Point<dim>  &points,
          const unsigned int component = 0) const override;

    virtual void
    value_list(const std::vector<Point<dim>> &points,
               std::vector<double>           &values,
               const unsigned int             component = 0) const override;

    virtual Tensor<1, dim>
    gradient(const Point<dim>  &p,
             const unsigned int component = 0) const override;

    virtual void
    gradient_list(const std::vector<Point<dim>> &points,
                  std::vector<Tensor<1, dim>>   &gradients,
                  const unsigned int             component = 0) const override;

  private:
    unsigned int order;
    double       wave_number;
    Point<dim>   center;
  };
} // namespace Functions

DEAL_II_NAMESPACE_CLOSE

#endif
