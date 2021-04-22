// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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

#ifndef dealii_function_level_set_h
#define dealii_function_level_set_h

#include <deal.II/base/config.h>

#include <deal.II/base/function.h>

DEAL_II_NAMESPACE_OPEN

namespace Functions
{
  namespace LevelSet
  {
    /**
     * Signed-distance level set function of a sphere:
     * $\psi(x) = \| x - x^c \| - R$.
     * Here, $x^c$ is the center of the sphere and $R$ is its radius. This
     * function is thus zero on the sphere, negative "inside" the ball having
     * the sphere as its boundary, and positive in the rest of
     * $\mathbb{R}^{dim}$.
     *
     * This function has gradient and Hessian equal to
     * $\partial_i \psi(x) = (x - x^c)/\| x - x^c \|$,
     * $\partial_i \partial_j \psi =
     * \delta_{ij}/\| x - x^c \| - (x_i - x_i^c)(x_j - x_j^c)/\| x - x^c \|^3$,
     * where $\delta_{ij}$ is the Kronecker delta function.
     *
     * @ingroup functions
     */
    template <int dim>
    class Sphere : public Function<dim>
    {
    public:
      /**
       * Constructor, takes the center and radius of the sphere.
       */
      Sphere(const Point<dim> &center = Point<dim>(), const double radius = 1);

      double
      value(const Point<dim> & point,
            const unsigned int component = 0) const override;

      /**
       * @copydoc Function::gradient()
       *
       * @note The gradient is singular at the center of the sphere. If the
       * incoming @p point is too close to the center, a floating-point
       * exception may be thrown or entries in the gradient may be +inf/-inf
       * or +nan/-nan, depending on how the point is located relative to the
       * singularity.
       */
      Tensor<1, dim>
      gradient(const Point<dim> & point,
               const unsigned int component = 0) const override;

      /**
       * @copydoc Function::hessian()
       *
       * @note The Hessian is singular at the center of the sphere. If the
       * incoming @p point is too close to the center, a floating-point
       * exception may be thrown or entries in the Hessian may be +inf/-inf
       * or +nan/-nan, depending on how the point is located relative to the
       * singularity.
       */
      SymmetricTensor<2, dim>
      hessian(const Point<dim> & point,
              const unsigned int component = 0) const override;

    private:
      const Point<dim> center;
      const double     radius;
    };


    /**
     * Signed level set function of a plane in $\mathbb{R}^{dim}$:
     * $\psi(x) = n \cdot (x - x_p)$.
     * Here, $n$ is the plane normal and $x_p$ is a point in the plane.
     * Thus, with respect to the direction of the normal, this function is
     * positive above the plane, zero in the plane, and negative below the
     * plane. If the normal is normalized, $\psi$ will be the signed distance to
     * the closest point in the plane.
     *
     * @ingroup functions
     */
    template <int dim>
    class Plane : public Function<dim>
    {
    public:
      /**
       * Constructor, takes a point in the plane and the plane normal.
       */
      Plane(const Point<dim> &point, const Tensor<1, dim> &normal);

      double
      value(const Point<dim> & point,
            const unsigned int component = 0) const override;

      Tensor<1, dim>
      gradient(const Point<dim> &,
               const unsigned int component = 0) const override;

      SymmetricTensor<2, dim>
      hessian(const Point<dim> &,
              const unsigned int component = 0) const override;

    private:
      const Point<dim>     point_in_plane;
      const Tensor<1, dim> normal;
    };

  } // namespace LevelSet
} // namespace Functions

DEAL_II_NAMESPACE_CLOSE

#endif /* CE8BBB3E_B726_40A7_B963_561AE7B84973 */
