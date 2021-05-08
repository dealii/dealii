// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2020 by the deal.II authors
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

#ifndef dealii_function_bessel_h
#define dealii_function_bessel_h


#include <deal.II/base/config.h>

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>

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
    value(const Point<dim> & points,
          const unsigned int component = 0) const override;

    virtual void
    value_list(const std::vector<Point<dim>> &points,
               std::vector<double> &          values,
               const unsigned int             component = 0) const override;

    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;

    virtual void
    gradient_list(const std::vector<Point<dim>> &points,
                  std::vector<Tensor<1, dim>> &  gradients,
                  const unsigned int             component = 0) const override;

  private:
    unsigned int order;
    double       wave_number;
    Point<dim>   center;
  };
} // namespace Functions

DEAL_II_NAMESPACE_CLOSE

#endif
