// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/non_matching/immersed_surface_quadrature.h>


DEAL_II_NAMESPACE_OPEN
namespace NonMatching
{
  template <int dim, int spacedim>
  ImmersedSurfaceQuadrature<dim, spacedim>::ImmersedSurfaceQuadrature(
    const std::vector<Point<dim>>          &points,
    const std::vector<double>              &weights,
    const std::vector<Tensor<1, spacedim>> &normals)
    : Quadrature<dim>(points, weights)
    , normals(normals)
  {
    AssertDimension(weights.size(), points.size());
    AssertDimension(normals.size(), points.size());
    for (const auto &normal : normals)
      {
        (void)normal;
        Assert(std::abs(normal.norm() - 1.0) < 1e-9,
               ExcMessage("Normal is not normalized."));
      }
  }



  template <int dim, int spacedim>
  inline void
  ImmersedSurfaceQuadrature<dim, spacedim>::clear()
  {
    this->quadrature_points.clear();
    this->weights.clear();
    this->normals.clear();
  }



  template <int dim, int spacedim>
  void
  ImmersedSurfaceQuadrature<dim, spacedim>::push_back(
    const Point<dim>          &point,
    const double               weight,
    const Tensor<1, spacedim> &normal)
  {
    this->quadrature_points.push_back(point);
    this->weights.push_back(weight);
    this->normals.push_back(normal);
    Assert(std::abs(normal.norm() - 1.0) < 1e-9,
           ExcMessage("Normal is not normalized."));
  }



  template <int dim, int spacedim>
  const Tensor<1, spacedim> &
  ImmersedSurfaceQuadrature<dim, spacedim>::normal_vector(
    const unsigned int i) const
  {
    AssertIndexRange(i, this->size());
    return normals[i];
  }



  template <int dim, int spacedim>
  const std::vector<Tensor<1, spacedim>> &
  ImmersedSurfaceQuadrature<dim, spacedim>::get_normal_vectors() const
  {
    return normals;
  }



  template class ImmersedSurfaceQuadrature<1, 1>;
  template class ImmersedSurfaceQuadrature<2, 2>;
  template class ImmersedSurfaceQuadrature<3, 3>;
  template class ImmersedSurfaceQuadrature<0, 1>;
  template class ImmersedSurfaceQuadrature<1, 2>;
  template class ImmersedSurfaceQuadrature<2, 3>;

} // namespace NonMatching
DEAL_II_NAMESPACE_CLOSE
