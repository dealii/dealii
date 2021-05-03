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

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function_level_set.h>

DEAL_II_NAMESPACE_OPEN

namespace Functions
{
  namespace LevelSet
  {
    template <int dim>
    Sphere<dim>::Sphere(const Point<dim> &center, const double radius)
      : center(center)
      , radius(radius)
    {
      Assert(radius > 0, ExcMessage("Radius must be positive."))
    }



    template <int dim>
    double
    Sphere<dim>::value(const Point<dim> & point,
                       const unsigned int component) const
    {
      AssertIndexRange(component, this->n_components);
      (void)component;

      return point.distance(center) - radius;
    }



    template <int dim>
    Tensor<1, dim>
    Sphere<dim>::gradient(const Point<dim> & point,
                          const unsigned int component) const
    {
      AssertIndexRange(component, this->n_components);
      (void)component;

      const Tensor<1, dim> center_to_point = point - center;
      const Tensor<1, dim> grad = center_to_point / center_to_point.norm();
      return grad;
    }



    template <int dim>
    SymmetricTensor<2, dim>
    Sphere<dim>::hessian(const Point<dim> & point,
                         const unsigned int component) const
    {
      AssertIndexRange(component, this->n_components);
      (void)component;

      const Tensor<1, dim> center_to_point = point - center;
      const double         distance        = center_to_point.norm();

      const SymmetricTensor<2, dim> hess =
        unit_symmetric_tensor<dim>() / distance -
        symmetrize(outer_product(center_to_point, center_to_point)) /
          std::pow(distance, 3);

      return hess;
    }



    template <int dim>
    Plane<dim>::Plane(const Point<dim> &point, const Tensor<1, dim> &normal)
      : point_in_plane(point)
      , normal(normal)
    {
      Assert(normal.norm() > 0, ExcMessage("Plane normal must not be 0."))
    }



    template <int dim>
    double
    Plane<dim>::value(const Point<dim> & point,
                      const unsigned int component) const
    {
      AssertIndexRange(component, this->n_components);
      (void)component;

      return normal * (point - point_in_plane);
    }



    template <int dim>
    Tensor<1, dim>
    Plane<dim>::gradient(const Point<dim> &, const unsigned int component) const
    {
      AssertIndexRange(component, this->n_components);
      (void)component;

      return normal;
    }



    template <int dim>
    SymmetricTensor<2, dim>
    Plane<dim>::hessian(const Point<dim> &, const unsigned int component) const
    {
      AssertIndexRange(component, this->n_components);
      (void)component;

      return SymmetricTensor<2, dim>();
    }

  } // namespace LevelSet
} // namespace Functions

#include "function_level_set.inst"

DEAL_II_NAMESPACE_CLOSE
