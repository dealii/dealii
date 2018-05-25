// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2017 by the deal.II authors
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

#include <deal.II/non_matching/immersed_surface_quadrature.h>

DEAL_II_NAMESPACE_OPEN
namespace NonMatching
{
  template <int dim>
  ImmersedSurfaceQuadrature<dim>::ImmersedSurfaceQuadrature(
    const std::vector<Point<dim>> &    points,
    const std::vector<double> &        weights,
    const std::vector<Tensor<1, dim>> &normals) :
    Quadrature<dim>(points, weights),
    normals(normals)
  {
    AssertDimension(weights.size(), points.size());
    AssertDimension(normals.size(), points.size());
    for (auto normal : normals)
      {
        (void)normal;
        Assert(std::abs(normal.norm() - 1.0) < 1e-9,
               ExcMessage("Normal is not normalized."));
      }
  }



  template <int dim>
  void
  ImmersedSurfaceQuadrature<dim>::push_back(const Point<dim> &    point,
                                            const double          weight,
                                            const Tensor<1, dim> &normal)
  {
    this->quadrature_points.push_back(point);
    this->weights.push_back(weight);
    this->normals.push_back(normal);
    Assert(std::abs(normal.norm() - 1.0) < 1e-9,
           ExcMessage("Normal is not normalized."));
  }



  template <int dim>
  const Tensor<1, dim> &
  ImmersedSurfaceQuadrature<dim>::normal_vector(const unsigned int i) const
  {
    AssertIndexRange(i, this->size());
    return normals[i];
  }



  template <int dim>
  const std::vector<Tensor<1, dim>> &
  ImmersedSurfaceQuadrature<dim>::get_normal_vectors() const
  {
    return normals;
  }



  template class ImmersedSurfaceQuadrature<1>;
  template class ImmersedSurfaceQuadrature<2>;
  template class ImmersedSurfaceQuadrature<3>;

} // namespace NonMatching
DEAL_II_NAMESPACE_CLOSE
