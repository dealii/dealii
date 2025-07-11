// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_mapping_internal_h
#define dealii_mapping_internal_h

#include <deal.II/base/config.h>

#include <deal.II/base/derivative_form.h>
#include <deal.II/base/tensor.h>

// Implementations of transformations used by several Mapping classes (such as
// MappingFE, MappingQ, and MappingFEField, and MappingManifold)

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  /**
   * Map the gradient of a covariant vector field. For more information see the
   * overload of Mapping::transform() which maps 2-differential forms from the
   * reference cell to the physical cell.
   */
  template <int dim, int spacedim, typename Number>
  Tensor<3, spacedim, Number>
  apply_covariant_gradient(
    const DerivativeForm<1, dim, spacedim, Number> &covariant,
    const DerivativeForm<2, dim, spacedim, Number> &input);

  /**
   * Map the Hessian of a contravariant vector field. For more information see
   * the overload of Mapping::transform() which maps 3-differential forms from
   * the reference cell to the physical cell.
   */
  template <int dim, int spacedim, typename Number>
  Tensor<3, spacedim, Number>
  apply_contravariant_hessian(
    const DerivativeForm<1, dim, spacedim, Number> &covariant,
    const DerivativeForm<1, dim, spacedim, Number> &contravariant,
    const Tensor<3, dim, Number>                   &input);

  /**
   * Map the Hessian of a covariant vector field. For more information see the
   * overload of Mapping::transform() which maps 3-differential forms from the
   * reference cell to the physical cell.
   */
  template <int dim, int spacedim, typename Number>
  Tensor<3, spacedim, Number>
  apply_covariant_hessian(
    const DerivativeForm<1, dim, spacedim, Number> &covariant,
    const Tensor<3, dim, Number>                   &input);

  /**
   * Map the Hessian of a Piola vector field. For more information see the
   * overload of Mapping::transform() which maps 3-differential forms from the
   * reference cell to the physical cell.
   */
  template <int dim, int spacedim, typename Number>
  Tensor<3, spacedim, Number>
  apply_piola_hessian(
    const DerivativeForm<1, dim, spacedim, Number> &covariant,
    const DerivativeForm<1, dim, spacedim, Number> &contravariant,
    const Number                                   &volume_element,
    const Tensor<3, dim, Number>                   &input);
} // namespace internal

namespace internal
{
  template <int dim, int spacedim, typename Number>
  Tensor<3, spacedim, Number>
  apply_covariant_gradient(
    const DerivativeForm<1, dim, spacedim, Number> &covariant,
    const DerivativeForm<2, dim, spacedim, Number> &input)
  {
    Tensor<3, spacedim, Number> output;
    for (unsigned int i = 0; i < spacedim; ++i)
      for (unsigned int j = 0; j < spacedim; ++j)
        {
          double tmp[dim];
          for (unsigned int K = 0; K < dim; ++K)
            {
              tmp[K] = covariant[j][0] * input[i][0][K];
              for (unsigned int J = 1; J < dim; ++J)
                tmp[K] += covariant[j][J] * input[i][J][K];
            }
          for (unsigned int k = 0; k < spacedim; ++k)
            {
              output[i][j][k] = covariant[k][0] * tmp[0];
              for (unsigned int K = 1; K < dim; ++K)
                output[i][j][k] += covariant[k][K] * tmp[K];
            }
        }

    return output;
  }



  template <int dim, int spacedim, typename Number>
  inline Tensor<3, spacedim, Number>
  apply_contravariant_hessian(
    const DerivativeForm<1, dim, spacedim, Number> &covariant,
    const DerivativeForm<1, dim, spacedim, Number> &contravariant,
    const Tensor<3, dim, Number>                   &input)
  {
    Tensor<3, spacedim, Number> output;
    for (unsigned int i = 0; i < spacedim; ++i)
      {
        Number tmp1[dim][dim];
        for (unsigned int J = 0; J < dim; ++J)
          for (unsigned int K = 0; K < dim; ++K)
            {
              tmp1[J][K] = contravariant[i][0] * input[0][J][K];
              for (unsigned int I = 1; I < dim; ++I)
                tmp1[J][K] += contravariant[i][I] * input[I][J][K];
            }
        for (unsigned int j = 0; j < spacedim; ++j)
          {
            Number tmp2[dim];
            for (unsigned int K = 0; K < dim; ++K)
              {
                tmp2[K] = covariant[j][0] * tmp1[0][K];
                for (unsigned int J = 1; J < dim; ++J)
                  tmp2[K] += covariant[j][J] * tmp1[J][K];
              }
            for (unsigned int k = 0; k < spacedim; ++k)
              {
                output[i][j][k] = covariant[k][0] * tmp2[0];
                for (unsigned int K = 1; K < dim; ++K)
                  output[i][j][k] += covariant[k][K] * tmp2[K];
              }
          }
      }

    return output;
  }



  template <int dim, int spacedim, typename Number>
  inline Tensor<3, spacedim, Number>
  apply_covariant_hessian(
    const DerivativeForm<1, dim, spacedim, Number> &covariant,
    const Tensor<3, dim, Number>                   &input)
  {
    Tensor<3, spacedim, Number> output;
    for (unsigned int i = 0; i < spacedim; ++i)
      {
        Number tmp1[dim][dim];
        for (unsigned int J = 0; J < dim; ++J)
          for (unsigned int K = 0; K < dim; ++K)
            {
              tmp1[J][K] = covariant[i][0] * input[0][J][K];
              for (unsigned int I = 1; I < dim; ++I)
                tmp1[J][K] += covariant[i][I] * input[I][J][K];
            }
        for (unsigned int j = 0; j < spacedim; ++j)
          {
            Number tmp2[dim];
            for (unsigned int K = 0; K < dim; ++K)
              {
                tmp2[K] = covariant[j][0] * tmp1[0][K];
                for (unsigned int J = 1; J < dim; ++J)
                  tmp2[K] += covariant[j][J] * tmp1[J][K];
              }
            for (unsigned int k = 0; k < spacedim; ++k)
              {
                output[i][j][k] = covariant[k][0] * tmp2[0];
                for (unsigned int K = 1; K < dim; ++K)
                  output[i][j][k] += covariant[k][K] * tmp2[K];
              }
          }
      }

    return output;
  }



  template <int dim, int spacedim, typename Number>
  inline Tensor<3, spacedim, Number>
  apply_piola_hessian(
    const DerivativeForm<1, dim, spacedim, Number> &covariant,
    const DerivativeForm<1, dim, spacedim, Number> &contravariant,
    const Number                                   &volume_element,
    const Tensor<3, dim, Number>                   &input)
  {
    Tensor<3, spacedim, Number> output;
    for (unsigned int i = 0; i < spacedim; ++i)
      {
        Number factor[dim];
        for (unsigned int I = 0; I < dim; ++I)
          factor[I] = contravariant[i][I] * (1. / volume_element);
        Number tmp1[dim][dim];
        for (unsigned int J = 0; J < dim; ++J)
          for (unsigned int K = 0; K < dim; ++K)
            {
              tmp1[J][K] = factor[0] * input[0][J][K];
              for (unsigned int I = 1; I < dim; ++I)
                tmp1[J][K] += factor[I] * input[I][J][K];
            }
        for (unsigned int j = 0; j < spacedim; ++j)
          {
            Number tmp2[dim];
            for (unsigned int K = 0; K < dim; ++K)
              {
                tmp2[K] = covariant[j][0] * tmp1[0][K];
                for (unsigned int J = 1; J < dim; ++J)
                  tmp2[K] += covariant[j][J] * tmp1[J][K];
              }
            for (unsigned int k = 0; k < spacedim; ++k)
              {
                output[i][j][k] = covariant[k][0] * tmp2[0];
                for (unsigned int K = 1; K < dim; ++K)
                  output[i][j][k] += covariant[k][K] * tmp2[K];
              }
          }
      }

    return output;
  }
} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
