// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/signaling_nan.h>

#include <deal.II/fe/mapping_related_data.h>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace FEValuesImplementation
  {
    template <int dim, int spacedim>
    void
    MappingRelatedData<dim, spacedim>::initialize(
      const unsigned int n_quadrature_points,
      const UpdateFlags  flags)
    {
      if (flags & update_quadrature_points)
        this->quadrature_points.resize(
          n_quadrature_points,
          Point<spacedim>(numbers::signaling_nan<Tensor<1, spacedim>>()));

      if (flags & update_JxW_values)
        this->JxW_values.resize(n_quadrature_points,
                                numbers::signaling_nan<double>());

      if (flags & update_jacobians)
        this->jacobians.resize(
          n_quadrature_points,
          numbers::signaling_nan<DerivativeForm<1, dim, spacedim>>());

      if (flags & update_jacobian_grads)
        this->jacobian_grads.resize(
          n_quadrature_points,
          numbers::signaling_nan<DerivativeForm<2, dim, spacedim>>());

      if (flags & update_jacobian_pushed_forward_grads)
        this->jacobian_pushed_forward_grads.resize(
          n_quadrature_points, numbers::signaling_nan<Tensor<3, spacedim>>());

      if (flags & update_jacobian_2nd_derivatives)
        this->jacobian_2nd_derivatives.resize(
          n_quadrature_points,
          numbers::signaling_nan<DerivativeForm<3, dim, spacedim>>());

      if (flags & update_jacobian_pushed_forward_2nd_derivatives)
        this->jacobian_pushed_forward_2nd_derivatives.resize(
          n_quadrature_points, numbers::signaling_nan<Tensor<4, spacedim>>());

      if (flags & update_jacobian_3rd_derivatives)
        this->jacobian_3rd_derivatives.resize(n_quadrature_points);

      if (flags & update_jacobian_pushed_forward_3rd_derivatives)
        this->jacobian_pushed_forward_3rd_derivatives.resize(
          n_quadrature_points, numbers::signaling_nan<Tensor<5, spacedim>>());

      if (flags & update_inverse_jacobians)
        this->inverse_jacobians.resize(
          n_quadrature_points,
          numbers::signaling_nan<DerivativeForm<1, spacedim, dim>>());

      if (flags & update_boundary_forms)
        this->boundary_forms.resize(
          n_quadrature_points, numbers::signaling_nan<Tensor<1, spacedim>>());

      if (flags & update_normal_vectors)
        this->normal_vectors.resize(
          n_quadrature_points, numbers::signaling_nan<Tensor<1, spacedim>>());
    }



    template <int dim, int spacedim>
    std::size_t
    MappingRelatedData<dim, spacedim>::memory_consumption() const
    {
      return (
        MemoryConsumption::memory_consumption(JxW_values) +
        MemoryConsumption::memory_consumption(jacobians) +
        MemoryConsumption::memory_consumption(jacobian_grads) +
        MemoryConsumption::memory_consumption(jacobian_pushed_forward_grads) +
        MemoryConsumption::memory_consumption(jacobian_2nd_derivatives) +
        MemoryConsumption::memory_consumption(
          jacobian_pushed_forward_2nd_derivatives) +
        MemoryConsumption::memory_consumption(jacobian_3rd_derivatives) +
        MemoryConsumption::memory_consumption(
          jacobian_pushed_forward_3rd_derivatives) +
        MemoryConsumption::memory_consumption(inverse_jacobians) +
        MemoryConsumption::memory_consumption(quadrature_points) +
        MemoryConsumption::memory_consumption(normal_vectors) +
        MemoryConsumption::memory_consumption(boundary_forms));
    }
  } // namespace FEValuesImplementation
} // namespace internal

#include "fe/mapping_related_data.inst"

DEAL_II_NAMESPACE_CLOSE
