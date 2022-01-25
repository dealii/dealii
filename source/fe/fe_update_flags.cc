// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2021 by the deal.II authors
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


#include <deal.II/fe/fe_update_flags.h>

DEAL_II_NAMESPACE_OPEN

constexpr UpdateFlags UpdateFlags::update_default =
  internal::make_update_flags(0);
constexpr UpdateFlags UpdateFlags::update_values =
  internal::make_update_flags(0x0001);
constexpr UpdateFlags UpdateFlags::update_gradients =
  internal::make_update_flags(0x0002);
constexpr UpdateFlags UpdateFlags::update_hessians =
  internal::make_update_flags(0x0004);
constexpr UpdateFlags UpdateFlags::update_3rd_derivatives =
  internal::make_update_flags(0x0008);
constexpr UpdateFlags UpdateFlags::update_boundary_forms =
  internal::make_update_flags(0x0010);
constexpr UpdateFlags UpdateFlags::update_quadrature_points =
  internal::make_update_flags(0x0020);
constexpr UpdateFlags UpdateFlags::update_JxW_values =
  internal::make_update_flags(0x0040);
constexpr UpdateFlags UpdateFlags::update_normal_vectors =
  internal::make_update_flags(0x0080);
constexpr UpdateFlags UpdateFlags::update_jacobians =
  internal::make_update_flags(0x0100);
constexpr UpdateFlags UpdateFlags::update_jacobian_grads =
  internal::make_update_flags(0x0200);
constexpr UpdateFlags UpdateFlags::update_inverse_jacobians =
  internal::make_update_flags(0x0400);
constexpr UpdateFlags UpdateFlags::update_covariant_transformation =
  internal::make_update_flags(0x0800);
constexpr UpdateFlags UpdateFlags::update_contravariant_transformation =
  internal::make_update_flags(0x1000);
constexpr UpdateFlags UpdateFlags::update_transformation_values =
  internal::make_update_flags(0x2000);
constexpr UpdateFlags UpdateFlags::update_transformation_gradients =
  internal::make_update_flags(0x4000);
constexpr UpdateFlags UpdateFlags::update_volume_elements =
  internal::make_update_flags(0x10000);
constexpr UpdateFlags UpdateFlags::update_jacobian_pushed_forward_grads =
  internal::make_update_flags(0x100000);
constexpr UpdateFlags UpdateFlags::update_jacobian_2nd_derivatives =
  internal::make_update_flags(0x200000);
constexpr UpdateFlags
  UpdateFlags::update_jacobian_pushed_forward_2nd_derivatives =
    internal::make_update_flags(0x400000);
constexpr UpdateFlags UpdateFlags::update_jacobian_3rd_derivatives =
  internal::make_update_flags(0x800000);
constexpr UpdateFlags update_jacobian_pushed_forward_3rd_derivatives =
  internal::make_update_flags(0x1000000);
constexpr UpdateFlags update_piola =
  update_volume_elements | update_contravariant_transformation;
constexpr UpdateFlags update_mapping =
  update_quadrature_points | update_JxW_values | update_jacobians |
  update_jacobian_grads | update_jacobian_pushed_forward_grads |
  update_jacobian_2nd_derivatives |
  update_jacobian_pushed_forward_2nd_derivatives |
  update_jacobian_3rd_derivatives |
  update_jacobian_pushed_forward_3rd_derivatives | update_inverse_jacobians |
  update_boundary_forms | update_normal_vectors |
  update_covariant_transformation | update_contravariant_transformation |
  update_transformation_values | update_transformation_gradients |
  update_volume_elements;

constexpr UpdateFlags update_default   = internal::make_update_flags(0);
constexpr UpdateFlags update_values    = internal::make_update_flags(0x0001);
constexpr UpdateFlags update_gradients = internal::make_update_flags(0x0002);
constexpr UpdateFlags update_hessians  = internal::make_update_flags(0x0004);
constexpr UpdateFlags update_3rd_derivatives =
  internal::make_update_flags(0x0008);
constexpr UpdateFlags update_boundary_forms =
  internal::make_update_flags(0x0010);
constexpr UpdateFlags update_quadrature_points =
  internal::make_update_flags(0x0020);
constexpr UpdateFlags update_JxW_values = internal::make_update_flags(0x0040);
constexpr UpdateFlags update_normal_vectors =
  internal::make_update_flags(0x0080);
constexpr UpdateFlags update_jacobians = internal::make_update_flags(0x0100);
constexpr UpdateFlags update_jacobian_grads =
  internal::make_update_flags(0x0200);
constexpr UpdateFlags update_inverse_jacobians =
  internal::make_update_flags(0x0400);
constexpr UpdateFlags update_covariant_transformation =
  internal::make_update_flags(0x0800);
constexpr UpdateFlags update_contravariant_transformation =
  internal::make_update_flags(0x1000);
constexpr UpdateFlags update_transformation_values =
  internal::make_update_flags(0x2000);
constexpr UpdateFlags update_transformation_gradients =
  internal::make_update_flags(0x4000);
constexpr UpdateFlags update_volume_elements =
  internal::make_update_flags(0x10000);
constexpr UpdateFlags update_jacobian_pushed_forward_grads =
  internal::make_update_flags(0x100000);
constexpr UpdateFlags update_jacobian_2nd_derivatives =
  internal::make_update_flags(0x200000);
constexpr UpdateFlags
  UpdateFlags::update_jacobian_pushed_forward_2nd_derivatives =
    internal::make_update_flags(0x400000);
constexpr UpdateFlags update_jacobian_3rd_derivatives =
  internal::make_update_flags(0x800000);
constexpr UpdateFlags
  UpdateFlags::update_jacobian_pushed_forward_3rd_derivatives =
    internal::make_update_flags(0x1000000);
constexpr UpdateFlags update_piola =
  update_volume_elements | update_contravariant_transformation;
constexpr UpdateFlags update_mapping =
  update_quadrature_points | update_JxW_values | update_jacobians |
  update_jacobian_grads | update_jacobian_pushed_forward_grads |
  update_jacobian_2nd_derivatives |
  update_jacobian_pushed_forward_2nd_derivatives |
  update_jacobian_3rd_derivatives |
  update_jacobian_pushed_forward_3rd_derivatives | update_inverse_jacobians |
  update_boundary_forms | update_normal_vectors |
  update_covariant_transformation | update_contravariant_transformation |
  update_transformation_values | update_transformation_gradients |
  update_volume_elements;

DEAL_II_NAMESPACE_CLOSE
