// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2017 by the deal.II authors
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

#ifndef dealii_function_spherical_h
#define dealii_function_spherical_h

#include <deal.II/base/config.h>

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>

#include <array>

DEAL_II_NAMESPACE_OPEN

namespace Functions
{
  /**
   * An abstract base class for a scalar-valued function $f=f(r,\theta,\phi)$
   * defined in spherical coordinates. This class wraps transformation of values,
   * gradients and hessians from spherical coordinates to the Cartesian coordinate
   * system used by the Function base class.
   * Therefore derived classes only need to implement those functions in
   * spherical coordinates (specifically svalue(), sgradient() and
   * shessian() ). The convention for angles is the same as in
   * GeometricUtilities::Coordinates.
   *
   * @note This function is currently only implemented for dim==3 .
   *
   * @author Denis Davydov, 2017
   */
  template <int dim>
  class Spherical : public Function<dim>
  {
  public:
    /**
     * Constructor which should be provided with @p center defining the origin
     * of the coordinate system.
     *
     * Note that components of this function are treated as entirely separate
     * quantities -- not as the components of a vector that will be
     * re-interpreted in a different coordinate system.
     */
    Spherical(const Point<dim>&  center       = Point<dim>(),
              const unsigned int n_components = 1);

    /**
     * Return the value of the function at the given point.
     *
     * This function converts the given point to spherical coordinates,
     * calls svalue() with it, and returns the result.
     */
    virtual double
    value(const Point<dim>&  point,
          const unsigned int component = 0) const override;

    /**
     * Return the gradient with respect to the Cartesian coordinates at point @p p.
     *
     * This function converts the given point to spherical coordinates,
     * calls sgradient() with it, and converts the result into Cartesian
     * coordinates.
     */
    virtual Tensor<1, dim>
    gradient(const Point<dim>&  p,
             const unsigned int component = 0) const override;

    /**
     * Return the Hessian with respect to the Cartesian coordinates at point @p p.
     *
     * This function converts the given point to spherical coordinates,
     * calls sgradient and shessian() with it, and converts the result into
     * Cartesian coordinates.
     */
    virtual SymmetricTensor<2, dim>
    hessian(const Point<dim>&  p,
            const unsigned int component = 0) const override;

    std::size_t
    memory_consumption() const;

  private:
    /**
     * Return the value at point @p sp. Here, @p sp is provided in spherical
     * coordinates.
     */
    virtual double
    svalue(const std::array<double, dim>& sp,
           const unsigned int             component) const;

    /**
     * Return the gradient in spherical coordinates.
     *
     * The returned object should contain derivatives in the following order:
     * $\{ f_{,r},\, f_{,\theta},\, f_{,\phi}\}$.
     */
    virtual std::array<double, dim>
    sgradient(const std::array<double, dim>& sp,
              const unsigned int             component) const;

    /**
     * Return the Hessian in spherical coordinates.
     *
     * The returned object should contain derivatives in the following order:
     * $\{ f_{,rr},\, f_{,\theta\theta},\, f_{,\phi\phi},\, f_{,r\theta},\, f_{,r\phi},\, f_{,\theta\phi}\}$.
     */
    virtual std::array<double, 6>
    shessian(const std::array<double, dim>& sp,
             const unsigned int             component) const;

    /**
     * A vector from the origin to the center of spherical coordinate system.
     */
    const Tensor<1, dim> coordinate_system_offset;
  };
} // namespace Functions

DEAL_II_NAMESPACE_CLOSE

#endif
