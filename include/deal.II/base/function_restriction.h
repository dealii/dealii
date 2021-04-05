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

#ifndef dealii_function_restriction_h
#define dealii_function_restriction_h

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/tensor.h>

DEAL_II_NAMESPACE_OPEN

namespace Functions
{
  /**
   * This class takes a function in `dim + 1` dimensions and creates a new
   * function in one dimension lower by restricting one of the coordinates to a
   * given value. Mathematically this corresponds to taking a function
   * $f = f(x, y, z)$,
   * a fixed value, $Z$, and defining a new function (the restriction)
   * $g = g(x, y) = f(x, y, Z)$.
   * Using this class, this translates to
   * @code
   *   Function<3> &            function             = ...
   *   double                   z                    = ...
   *   unsigned int             restricted_direction = 2;
   *   CoordinateRestriction<2> restriction(function, restricted_direction, z);
   * @endcode
   *
   * The `dim`-dimensional coordinates on the restriction are ordered starting
   * from the restricted (`dim + 1`)-coordinate. In particular, this means that
   * if the $y$-coordinate is locked to $Y$ in 3D, the coordinates are ordered
   * as $(z, x)$ on the restriction:
   * $g = g(z, x) = f(x, Y, z)$.
   * This is the same convention as in BoundingBox::cross_section.
   */
  template <int dim>
  class CoordinateRestriction : public Function<dim>
  {
  public:
    /**
     * Constructor, takes the (`dim + 1`)-coordinate direction and the value
     * that the incoming function should be restricted to.
     *
     * A pointer to the incoming function is stored internally, so the function
     * must have a longer lifetime than the created restriction.
     */
    CoordinateRestriction(const Function<dim + 1> &function,
                          const unsigned int       direction,
                          const double             coordinate_value);

    double
    value(const Point<dim> &point, const unsigned int component) const override;

    Tensor<1, dim>
    gradient(const Point<dim> & point,
             const unsigned int component) const override;

    SymmetricTensor<2, dim>
    hessian(const Point<dim> & point,
            const unsigned int component) const override;

  private:
    // The higher-dimensional function that has been restricted.
    const SmartPointer<const Function<dim + 1>> function;

    // The (`dim + 1`)-coordinate direction that has been restricted.
    const unsigned int restricted_direction;

    // Value of the restriced coordinate.
    const double coordinate_value;
  };



  /**
   * This class creates a 1-dimensional function from a `dim + 1` dimensional
   * function by restricting `dim` of the coordinate values to a given point.
   * Mathematically this corresponds to taking a function, $f = f(x, y, z)$, and
   * a point $(Y, Z)$, and defining a new function $g = g(x) = f(x, Y, Z)$.
   * Using this class, this translates to
   * @code
   *   Function<3> &       function = ...
   *   Point<2>            point(y, z);
   *   unsigned int        open_direction = 0;
   *   PointRestriction<2> restriction(function, open_direction, point);
   * @endcode
   *
   * The coordinates of the point will be expanded in the higher-dimensional
   * functions coordinates starting from the open-direction (and wrapping
   * around). In particular, if we restrict to a point $(Z, X)$ and choose to
   * keep the y-direction open, the restriction that is created is the function
   * $g(y) = f(X, y, Z)$.
   * This is consistent with the convention in BoundingBox::cross_section.
   */
  template <int dim>
  class PointRestriction : public Function<1>
  {
  public:
    /**
     * Constructor, takes the point that the incoming function should be
     * restricted to and which (`dim + 1`)-dimensional coordinate direction
     * should be kept "open".
     *
     * A pointer to the incoming function is stored internally, so the function
     * must have a longer lifetime than the created restriction.
     */
    PointRestriction(const Function<dim + 1> &function,
                     const unsigned int       open_direction,
                     const Point<dim> &       point);

    double
    value(const Point<1> &point, const unsigned int component) const override;

    Tensor<1, 1>
    gradient(const Point<1> &   point,
             const unsigned int component) const override;

    SymmetricTensor<2, 1>
    hessian(const Point<1> &point, const unsigned int component) const override;

  private:
    // The higher-dimensional function that has been restricted.
    const SmartPointer<const Function<dim + 1>> function;

    // The (`dim + 1`)-coordinate direction that is kept "open"
    const unsigned int open_direction;

    // The point that we have restricted the above function to.
    const Point<dim> point;
  };

} // namespace Functions


namespace internal
{
  /**
   * Creates a (`dim + 1`)-dimensional point by copying over the coordinates of
   * the incoming `dim`-dimensional point and setting the "missing"
   * (`dim + 1`)-dimensional component to the incoming coordinate value.
   *
   * For example, given the input
   * $\{(x, y), 2, z \}$ this function creates the point $(x, y, z)$.
   *
   * The coordinates of the `dim`-dimensional point are written to the
   * coordinates of the (`dim + 1`)-dimensional point in the order of the
   * convention given by the function coordinate_to_one_dim_higher. Thus, the
   * order of coordinates on the lower-dimensional point are not preserved:
   * $\{(z, x), 1, y \}$ creates the point $(x, y, z)$.
   */
  template <int dim>
  Point<dim + 1>
  create_higher_dim_point(const Point<dim> & point,
                          const unsigned int component_in_dim_plus_1,
                          const double       coordinate_value);
} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif /* dealii_function_restriction_h */
